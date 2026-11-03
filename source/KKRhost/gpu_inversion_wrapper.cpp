#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuComplex.h>
#include <stdlib.h>
#include <stdio.h>

// On-device dense LU solve for the Chebyshev radial solver (rll_local_solutions).
// Factor A (n x n, device) once and solve two RHS B1,B2 (each n x nrhs, device) in place.
// Replaces host zgetrf/zgetrs so slv/yrll/zrll stay device-resident (no managed-memory
// page migration). Called thousands of times/iteration -> persistent handle + workspace.
static cusolverDnHandle_t g_rll_handle = NULL;
static cuDoubleComplex*   g_rll_work   = NULL;
static int                g_rll_lwork  = 0;
static int*               g_rll_ipiv   = NULL;
static int                g_rll_ipiv_n = 0;
static int*               g_rll_info   = NULL;

extern "C" void gpu_rll_solve_c(void* A, void* B1, void* B2, int n, int nrhs)
{
    if (!g_rll_handle) {
        cusolverDnCreate(&g_rll_handle);
        cudaMalloc((void**)&g_rll_info, sizeof(int));
    }
    int lwork = 0;
    cusolverDnZgetrf_bufferSize(g_rll_handle, n, n, (cuDoubleComplex*)A, n, &lwork);
    if (lwork > g_rll_lwork) {
        if (g_rll_work) cudaFree(g_rll_work);
        cudaMalloc((void**)&g_rll_work, sizeof(cuDoubleComplex)*lwork);
        g_rll_lwork = lwork;
    }
    if (n > g_rll_ipiv_n) {
        if (g_rll_ipiv) cudaFree(g_rll_ipiv);
        cudaMalloc((void**)&g_rll_ipiv, sizeof(int)*n);
        g_rll_ipiv_n = n;
    }
    cusolverStatus_t st = cusolverDnZgetrf(g_rll_handle, n, n, (cuDoubleComplex*)A, n, g_rll_work, g_rll_ipiv, g_rll_info);
    if (st != CUSOLVER_STATUS_SUCCESS) fprintf(stderr, "cusolverDnZgetrf failed (%d)\n", (int)st);
    cusolverDnZgetrs(g_rll_handle, CUBLAS_OP_N, n, nrhs, (cuDoubleComplex*)A, n, g_rll_ipiv, (cuDoubleComplex*)B1, n, g_rll_info);
    cusolverDnZgetrs(g_rll_handle, CUBLAS_OP_N, n, nrhs, (cuDoubleComplex*)A, n, g_rll_ipiv, (cuDoubleComplex*)B2, n, g_rll_info);
    cudaDeviceSynchronize();
}

extern "C" void gpu_inversion_batched_c(
    void* gllke_batch,
    void* gtemp_batch,
    int* ipvt_batch,
    int* info_batch,
    int alm,
    int batchSize) 
{
    cublasHandle_t handle;
    cublasStatus_t status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "CUBLAS initialization failed!\n");
        return;
    }

    // Allocate host arrays of pointers to matrices
    cuDoubleComplex** A_host_ptrs = (cuDoubleComplex**)malloc(batchSize * sizeof(cuDoubleComplex*));
    cuDoubleComplex** B_host_ptrs = (cuDoubleComplex**)malloc(batchSize * sizeof(cuDoubleComplex*));

    cuDoubleComplex* A_base = (cuDoubleComplex*)gllke_batch;
    cuDoubleComplex* B_base = (cuDoubleComplex*)gtemp_batch;

    for (int i = 0; i < batchSize; i++) {
        A_host_ptrs[i] = A_base + i * alm * alm;
        B_host_ptrs[i] = B_base + i * alm * alm;
    }

    // Allocate device arrays of pointers to matrices
    cuDoubleComplex** A_dev_ptrs = NULL;
    cuDoubleComplex** B_dev_ptrs = NULL;
    
    cudaError_t err1 = cudaMalloc((void**)&A_dev_ptrs, batchSize * sizeof(cuDoubleComplex*));
    cudaError_t err2 = cudaMalloc((void**)&B_dev_ptrs, batchSize * sizeof(cuDoubleComplex*));
    if (err1 != cudaSuccess || err2 != cudaSuccess) {
        fprintf(stderr, "CUDA malloc for pointer arrays failed!\n");
        free(A_host_ptrs);
        free(B_host_ptrs);
        cublasDestroy(handle);
        return;
    }

    // Copy pointer arrays from host to device
    cudaMemcpy(A_dev_ptrs, A_host_ptrs, batchSize * sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);
    cudaMemcpy(B_dev_ptrs, B_host_ptrs, batchSize * sizeof(cuDoubleComplex*), cudaMemcpyHostToDevice);

    // Call cublasZgetrfBatched (batched LU factorization).
    // getrf's infoArray IS a DEVICE array of size batchSize (per-matrix singularity flag),
    // so passing the device pointer info_batch is correct here.
    status = cublasZgetrfBatched(handle, alm, A_dev_ptrs, alm, ipvt_batch, info_batch, batchSize);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "cublasZgetrfBatched failed (status=%d)!\n", (int)status);
    }

    // Call cublasZgetrsBatched (batched LU solve: A * X = B, B is replaced by X).
    // IMPORTANT: getrsBatched's `info` argument is a HOST pointer to a SINGLE int
    // (it reports an invalid-argument error, NOT per-matrix singularity). Passing a
    // device pointer here makes cuBLAS dereference device memory on the host -> SIGSEGV.
    // Use a host int, not info_batch.
    int info_getrs = 0;
    status = cublasZgetrsBatched(handle, CUBLAS_OP_N, alm, alm, A_dev_ptrs, alm, ipvt_batch, B_dev_ptrs, alm, &info_getrs, batchSize);
    if (status != CUBLAS_STATUS_SUCCESS || info_getrs != 0) {
        fprintf(stderr, "cublasZgetrsBatched failed (status=%d, info=%d)!\n", (int)status, info_getrs);
    }

    // Ensure the batched factor/solve have completed before we read the results,
    // and surface any asynchronous launch/execution error instead of crashing later.
    cudaError_t syncErr = cudaDeviceSynchronize();
    if (syncErr != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize after batched inversion failed: %s\n", cudaGetErrorString(syncErr));
    }

    // Check the per-matrix getrf info (device array) for singular factorizations,
    // which would otherwise silently propagate garbage into the Green function.
    int* info_host = (int*)malloc(batchSize * sizeof(int));
    if (info_host) {
        cudaMemcpy(info_host, info_batch, batchSize * sizeof(int), cudaMemcpyDeviceToHost);
        for (int i = 0; i < batchSize; i++) {
            if (info_host[i] != 0) {
                fprintf(stderr, "WARNING: batched LU factorization %d reported info=%d (singular U?)\n", i, info_host[i]);
            }
        }
        free(info_host);
    }

    // Copy back solved B (= A^{-1}, in gtemp_batch) to gllke_batch, where Fortran reads it.
    // Cast the element count to size_t so the byte size does not overflow int for large cells.
    cudaMemcpy(gllke_batch, gtemp_batch, (size_t)batchSize * alm * alm * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);

    // Free device and host memory
    cudaFree(A_dev_ptrs);
    cudaFree(B_dev_ptrs);
    free(A_host_ptrs);
    free(B_host_ptrs);

    cublasDestroy(handle);
}
