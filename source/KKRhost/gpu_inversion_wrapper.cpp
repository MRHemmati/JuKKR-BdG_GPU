#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#include <stdlib.h>
#include <stdio.h>

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

    // Call cublasZgetrfBatched (batched LU factorization)
    status = cublasZgetrfBatched(handle, alm, A_dev_ptrs, alm, ipvt_batch, info_batch, batchSize);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "cublasZgetrfBatched failed!\n");
    }

    // Call cublasZgetrsBatched (batched LU solve: A * X = B, B is replaced by X)
    status = cublasZgetrsBatched(handle, CUBLAS_OP_N, alm, alm, A_dev_ptrs, alm, ipvt_batch, B_dev_ptrs, alm, info_batch, batchSize);
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "cublasZgetrsBatched failed!\n");
    }

    // Copy back solved B (which is gtemp_batch) to gllke_batch
    // Since Fortran expects the result in gllke, we copy gtemp_batch back to gllke_batch on the device
    cudaMemcpy(gllke_batch, gtemp_batch, batchSize * alm * alm * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice);

    // Free device and host memory
    cudaFree(A_dev_ptrs);
    cudaFree(B_dev_ptrs);
    free(A_host_ptrs);
    free(B_host_ptrs);

    cublasDestroy(handle);
}
