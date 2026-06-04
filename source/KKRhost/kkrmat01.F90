!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_kkrmat01

  private
  public :: kkrmat01

contains

  !-------------------------------------------------------------------------------
  !> Summary: Performs k-space integration to determine scattering path operator
  !> \[\tau = \left(g\left(\mathbf{k},e\right)-t^{-1}\right)^{-1}\]
  !>  and Greens function of the real system -> \[GS\]
  !> Author: 
  !> Category: KKRhost, k-points, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Modifications according to H. Hoehler ( July 2002)
  !> define Fourier transformation as
  !>
  !> $$ G\left(\mu,\mu'\right)_{L,L'}=\frac{1}{2}\left[\sum_n G^{n,0}\left(\mu,\mu'\right)_{L,L'}\exp\left(-iKR^n\right)+ \sum_n G^{n,0}\left(\mu,\mu'\right)_{L,L'}\exp\left(-iK\left(-R\right)^n\right)\right]  $$
  !>
  !> this operation has to be done to satisfy the point symmetry;
  !> the application of the Fourier transformation is just an
  !> approximation for the tb system, since the translational invariance
  !> is not satisfied --> force it by R, -R

  !> @note
  !> - New version 10.99: up -> left , down -> right, for decimation
  !> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine kkrmat01(bzkp, nofks, gs, volcub, tinvll, rrot, nshell, nsdia, alat, nsymat, naez, cls, nacls, naclsmax, rr, ezoa, atom, nsh1, nsh2, ginp, rbasis, rcls, tinvbup, &
    tinvbdown, vacflag, nlbasis, nrbasis, factl, icheck, invmod, ideci, srrel, irrel, nrrel, dtrefll, dtmatll, dginp, refpot, lly_grtr, tracet, cfctor, lly) ! LLY
#ifdef CPP_MPI
    use :: mpi

    use :: mod_types, only: t_mpi_c_grid
    use :: mod_mympi, only: myrank, nranks, master, distribute_linear_on_tasks
#else
    use :: mod_mympi, only: myrank, nranks, master
#endif
#ifdef CPP_HYBRID
    use :: omp_lib
#endif
#ifdef CPP_TIMING
    use :: mod_timing, only: timing_start, timing_pause, timing_stop
#endif
    use :: mod_types, only: t_inc
    use :: mod_runoptions, only: print_program_flow, use_Chebychev_solver, use_qdos, use_virtual_atoms, &
      write_green_imp, write_rhoq_input, decouple_spin_cheby, calc_complex_bandstructure, symmetrize_gmat, use_deci_onebulk
    use :: mod_rhoqtools, only: rhoq_find_kmask, rhoq_saveg, rhoq_write_tau0, rhoq_read_mu0_scoef
    use :: global_variables, only: nembd1, nembd2, nsheld, nclsd, naclsd, lmmaxd, nprincd, nrd, nrefd, lmgf0d, krel, ndim_slabinv, alm, almgf0 
    use :: mod_constants, only: czero, cone, nsymaxd, ci,pi
    use :: mod_datatypes, only: dp
    use :: mod_decimate, only: decimate
    use :: mod_dlke0, only: dlke0
    use :: mod_inversion, only: inversion, inversion_gpu_batched
    use :: mod_cinit, only: cinit

    implicit none
    ! .. Input variables
    integer, intent (in) :: lly    !! LLY <> 0 --> use Lloyds formula
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: nofks  !! number of k-points
    integer, intent (in) :: nsdia   
    integer, intent (in) :: ideci   
    integer, intent (in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, intent (in) :: nsymat
    integer, intent (in) :: invmod !! Inversion scheme
    integer, intent (in) :: nlbasis !! Number of basis layers of left host (repeated units)
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    integer, intent (in) :: naclsmax !! maximal number of atoms in screening cluster
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    integer, dimension (nembd2), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (nsheld), intent (in) :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nsheld), intent (in) :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (nembd2), intent (in) :: refpot !! Ref. pot. card  at position ! REFPOT(NAEZD+NEMBD)
    integer, dimension (naclsd, nembd2), intent (in) :: atom !! Atom at site in cluster
    integer, dimension (naclsd, nembd2), intent (in) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (2, lmmaxd), intent (in) :: nrrel
    integer, dimension (naez/nprincd, naez/nprincd), intent (in) :: icheck
    integer, dimension (2, 2, lmmaxd), intent (in) :: irrel
    real (kind=dp), dimension (*), intent (in) :: volcub 
    real (kind=dp), dimension (3, 0:nrd), intent (in) :: rr !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension (3, *), intent (in) :: bzkp !! array of k-points in irreducable part of BZ
    real (kind=dp), dimension (3, nembd2), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension (48, 3, nsheld), intent (in) :: rrot
    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: factl
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: tinvll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1), intent (in) :: tinvbup
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nrefd), intent (in) :: dtrefll ! LLY dtref/dE
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: dtmatll ! LLY  dt/dE (should be av.-tmatrix in CPA)
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1), intent (in) :: tinvbdown
    complex (kind=dp), dimension (lmgf0d*naclsmax, lmgf0d, *), intent (in) :: ginp ! Gref
    complex (kind=dp), dimension (lmgf0d*naclsmax, lmgf0d, *), intent (in) :: dginp ! LLY dGref/dE
    complex (kind=dp), dimension (2, 2, lmmaxd), intent (in) :: srrel
    logical, dimension (2), intent (in) :: vacflag
    ! .. In/Out variables
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd, *), intent (inout) :: gs
    ! .. Output variables
    complex (kind=dp), intent (out) :: lly_grtr ! Trace Eq.5.38 PhD Thiess  (integrated) ! LLY Lloyd
    ! .. Local variables
    integer :: i_stat, i_all
    integer :: ikm1, ikm2, is, n1, n2, j1, j2, i2
    integer :: iq1, iq2, ioff1, ioff2, joff1, joff2
    integer :: i, i1, ilm, isym, iu, j, jlm, il1, kpt, lm, lm1, lm2, ns, il2, jl1, jl2
    real (kind=dp) :: zktr
    complex (kind=dp) :: lly_grtr_k ! Trace Eq.5.38 PhD Thiess  (k-dependent) ! LLY Lloyd
    complex (kind=dp) :: csum1, csum2, trace, tracet ! LLY Lloyd
    complex (kind=dp) :: carg, citpi, cfctor
    real (kind=dp), dimension (3) :: kp
    real (kind=dp), dimension (6) :: bzkpk
    real (kind=dp), dimension (3, 0:nrd) :: rrm
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: g
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gaux1 ! LLY
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gaux2 ! LLY
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gaux3 ! LLY
    complex (kind=dp), dimension (nsymaxd, nsheld) :: etaikr
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez) :: t_aux ! LLY auxiliary array for t-matrix manipulation
    ! .. Local allocatable arrays
    complex (kind=dp), dimension (:, :), allocatable :: gllke, gllkem, gllken
    complex (kind=dp), dimension (:, :), allocatable :: dgllke, dgllkem, dgllken, grefllke ! LLY
    complex (kind=dp), dimension (:, :), allocatable :: gllke0v, gllke0v2, gllketv ! for VIRTUAL ATOMS
    complex (kind=dp), dimension (:, :), allocatable :: gllketv_new ! for VIRTUAL ATOMS
    complex (kind=dp), dimension (:, :), allocatable :: gllke0, gllke0m

#ifdef CPP_GPU
    complex (kind=dp), dimension (:, :, :), allocatable :: gllken_batch, gllkem_batch, gllke_batch, gtemp_batch
    complex (kind=dp), dimension (:, :, :), allocatable :: dgllken_batch, dgllkem_batch, dgllke_batch
    complex (kind=dp), dimension (:, :, :), allocatable :: gllken1_batch, gllkem1_batch
    complex (kind=dp), dimension (:, :, :), allocatable :: gllke0_batch, gllke0m_batch
    complex (kind=dp), dimension (:, :, :), allocatable :: gllke01_batch, gllke0m1_batch
    integer, dimension (:, :), allocatable :: ipvt_batch
    integer, dimension (:), allocatable :: info_batch
    complex (kind=dp), dimension (:, :), allocatable :: etaikr_batch
    complex (kind=dp), dimension (:, :, :), allocatable :: grefllke_batch
    integer :: batchSize
    real (kind=dp) :: tpi_gpu, convpu_gpu
    complex (kind=dp) :: arg1, arg2, arg3, eikr_gpu, tt_gpu, eta_val
    complex (kind=dp), dimension (6) :: kp_gpu
    integer :: idx, k
#endif
    ! .. Parameters
    complex (kind=dp), parameter :: cmi = -ci !! negative imaginary part \[-i\]

#ifdef CPP_MPI
    integer :: ntot1
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt
#endif
    integer :: k_start, k_end
    ! !#ifdef CPP_MPI
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd) :: work
    integer :: ierr, iwork
    ! !#endif
    integer :: mu, nscoef, imin
    integer, allocatable :: iatomimp(:)

    real (kind=dp), allocatable :: rhoq_kmask(:, :) ! only in reduced number of kpts
    integer, allocatable :: kmask(:) ! logical array over all kpts (determine if kpt=1,nofks is in reduced set)
    integer :: mythread


    ! NDIM=LMGF0D*NAEZ
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Array sizes definitions
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (print_program_flow .and. (t_inc%i_write>0)) write (1337, *) '>>> kkrmat1: loop over k-points'

    citpi = cmi*2.0_dp*pi    ! = -i*2*PI

    do ns = 1, nshell
      do iu = 1, nsymaxd
        call cinit(lmmaxd*lmmaxd, gs(1,1,iu,ns))
      end do
    end do

    lly_grtr = czero               ! LLY Lloyd
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Array allocations
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (gllke(alm,alm), stat=i_stat)
    ! LLY Lloyd
    if (lly/=0) then
      allocate (dgllke(alm,alm), stat=i_stat)
      allocate (grefllke(alm,alm), stat=i_stat)
    end if

    if (use_virtual_atoms) then
      allocate (gllke0v(alm,alm), stat=i_stat)
      allocate (gllke0v2(alm,alm), stat=i_stat)
      allocate (gllketv(alm,lmmaxd), stat=i_stat)
      allocate (gllketv_new(lmmaxd,alm), stat=i_stat)
    end if                         ! ( use_virtual_atoms ) THEN
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of array allocations
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! K-points loop
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (write_rhoq_input) then
      call rhoq_read_mu0_scoef(iatomimp, mu, nscoef, imin)
      call rhoq_find_kmask(nofks, k_end, bzkp(1:3,1:nofks), kmask, rhoq_kmask)
    end if                         ! write_rhoq_input

#ifdef CPP_MPI
    ! MPI:
    if (.not. use_qdos) then

      if (write_rhoq_input) then
        ntot1 = k_end
      else
        ntot1 = nofks
      end if

      if (myrank==master) write (1337, *) 'kkrmat k loop:', nofks, t_mpi_c_grid%nranks_ie
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, ntot1, ntot_pt, ioff_pt, .true.)

      k_start = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
      k_end = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1 = ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot_pt1 = ntot_pt
      t_mpi_c_grid%ioff_pt1 = ioff_pt

    else                           ! .not.use_qdos

      k_start = 1
      k_end = nofks

    end if                         ! .not.use_qdos
#else
    k_start = 1
    if (.not. write_rhoq_input) k_end = nofks
#endif

    ! k-loop not needed for GREENIMP-case
    if (write_green_imp) then
      if (myrank==master) write (*, *) 'Skipping kloop in kkrmat'
      k_start = 1
      k_end = 0
    end if

    ! Print header of statusbar for k-loop
    if (t_inc%i_write>0) then
      write (1337, '("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
      write (1337, fmt=100, advance='no') ! Beginning of statusbar
    end if


    ! ----------------------------------------------------------------------
    ! allocattions of work arrays
    if (krel==0) then
      allocate (gllken(almgf0,almgf0), stat=i_stat)
      gllken(:, :) = czero
      allocate (gllkem(almgf0,almgf0), stat=i_stat)
      gllkem(:, :) = czero
      if (lly/=0) then
        allocate (dgllken(almgf0,almgf0), stat=i_stat)
        dgllken(:, :) = czero
        allocate (dgllkem(almgf0,almgf0), stat=i_stat)
        dgllkem(:, :) = czero
      end if
    else
      allocate (gllke0(almgf0,almgf0), stat=i_stat)
      allocate (gllke0m(almgf0,almgf0), stat=i_stat)
    end if                         ! (KREL.EQ.0)
    ! ----------------------------------------------------------------------

#ifdef CPP_GPU
    batchSize = max(0, k_end - k_start + 1)
    if (batchSize > 0) then
      allocate (gllke_batch(alm, alm, batchSize), stat=i_stat)
      allocate (gtemp_batch(alm, alm, batchSize), stat=i_stat)
      allocate (ipvt_batch(alm, batchSize), stat=i_stat)
      allocate (info_batch(batchSize), stat=i_stat)
      allocate (etaikr_batch(nsymaxd * nsheld, batchSize), stat=i_stat)
      if (krel==0) then
        allocate (gllken_batch(almgf0, almgf0, batchSize), stat=i_stat)
        allocate (gllkem_batch(almgf0, almgf0, batchSize), stat=i_stat)
        if (lly/=0) then
          allocate (dgllken_batch(almgf0, almgf0, batchSize), stat=i_stat)
          allocate (dgllkem_batch(almgf0, almgf0, batchSize), stat=i_stat)
          allocate (dgllke_batch(alm, alm, batchSize), stat=i_stat)
          allocate (grefllke_batch(alm, alm, batchSize), stat=i_stat)
        end if
        if (symmetrize_gmat) then
          allocate (gllken1_batch(almgf0, almgf0, batchSize), stat=i_stat)
          allocate (gllkem1_batch(almgf0, almgf0, batchSize), stat=i_stat)
        end if
      else
        allocate (gllke0_batch(almgf0, almgf0, batchSize), stat=i_stat)
        allocate (gllke0m_batch(almgf0, almgf0, batchSize), stat=i_stat)
        if (symmetrize_gmat) then
          allocate (gllke01_batch(almgf0, almgf0, batchSize), stat=i_stat)
          allocate (gllke0m1_batch(almgf0, almgf0, batchSize), stat=i_stat)
        end if
      end if

      !$acc enter data create(gllke_batch, gtemp_batch, ipvt_batch, info_batch, etaikr_batch)
      if (krel==0) then
        !$acc enter data create(gllken_batch, gllkem_batch)
        if (lly/=0) then
          !$acc enter data create(dgllken_batch, dgllkem_batch, dgllke_batch, grefllke_batch)
        end if
        if (symmetrize_gmat) then
          !$acc enter data create(gllken1_batch, gllkem1_batch)
        end if
      else
        !$acc enter data create(gllke0_batch, gllke0m_batch)
        if (symmetrize_gmat) then
          !$acc enter data create(gllke01_batch, gllke0m1_batch)
        end if
      end if
    end if
#endif
    ! ----------------------------------------------------------------------


#ifdef CPP_HYBRID
    ! $omp parallel default(shared) &
    ! $omp private(kpt, ns, i, j, isym, carg, i1, zktr, i2, iq1, iq2, ioff1) &
    ! $omp private(ioff2, joff1, joff2, ikm1, ikm2, csum1, is, n1, n2) &
    ! $omp private(j1, csum2, j2, il1, il2, lm1, lm2, gaux1, gaux2) &
    ! $omp private(jl1, jl2, gaux3, ilm, jlm, mythread) &
    ! $omp reduction(+:trace)
    mythread = omp_get_thread_num()
#else
    mythread = 0
#endif

    !$acc data copyin(bzkp(1:3, 1:nofks), &
    !$acc&            ginp(1:lmgf0d*naclsmax, 1:lmgf0d, 1:nclsd), &
    !$acc&            dginp(1:lmgf0d*naclsmax, 1:lmgf0d, 1:nclsd)) &
    !$acc& present(rr(1:3, 0:nrd), ezoa(1:naclsd, 1:nembd2), &
    !$acc&         atom(1:naclsd, 1:nembd2), cls(1:nembd2), nacls(1:nclsd), &
    !$acc&         rcls(1:3, 1:naclsd, 1:nclsd)) &
    !$acc& create(rrm(1:3, 0:nrd))

    ! Pre-compute rrm = -rr on the GPU (rr already resident from main1b data region)
    !$acc parallel loop collapse(2) present(rr(1:3, 0:nrd), rrm(1:3, 0:nrd))
    do i = 1, 3
      do j = 1, nrd
        rrm(i, j) = -rr(i, j)
      end do
    end do
    !$acc parallel loop present(rrm(1:3, 0:nrd))
    do i = 1, 3
      rrm(i, 0) = 0.0_dp
    end do

    ! kpts loop
#ifdef CPP_GPU
    if (.not. use_virtual_atoms .and. ideci == 0 .and. krel == 0) then
      tpi_gpu = 8.0_dp*atan(1.0_dp)
      convpu_gpu = alat/tpi_gpu

      if (batchSize > 0) then
        !$acc parallel loop collapse(3) present(gllken_batch, gllkem_batch, gllke_batch)
        do k = 1, batchSize
          do j = 1, almgf0
            do i = 1, almgf0
              gllken_batch(i, j, k) = czero
              gllkem_batch(i, j, k) = czero
            end do
          end do
        end do
        !$acc parallel loop collapse(3) present(gllke_batch)
        do k = 1, batchSize
          do j = 1, alm
            do i = 1, alm
              gllke_batch(i, j, k) = czero
            end do
          end do
        end do
        if (symmetrize_gmat) then
          !$acc parallel loop collapse(3) present(gllken1_batch, gllkem1_batch)
          do k = 1, batchSize
            do j = 1, almgf0
              do i = 1, almgf0
                gllken1_batch(i, j, k) = czero
                gllkem1_batch(i, j, k) = czero
              end do
            end do
          end do
        end if

        !$acc parallel loop collapse(1) &
        !$acc& present(bzkp, ginp, rr, ezoa, atom, cls, nacls, rcls) &
        !$acc& present(gllken_batch, gllkem_batch, gllke_batch, etaikr_batch) &
        !$acc& present(rrm, tinvll) &
        !$acc& private(kp_gpu, carg, zktr, tt_gpu, eikr_gpu, arg1, arg2, arg3) &
        !$acc& private(i, m, lm2, lm1, ic, im, am, isym, ns, j, i1, idx)
        do kpt = k_start, k_end
          idx = kpt - k_start + 1
          kp_gpu(1:3) = bzkp(1:3, kpt)
          if (calc_complex_bandstructure) then
            kp_gpu(4:6) = bzkp(4:6, kpt)
          else
            kp_gpu(4:6) = czero
          end if

          ! Phase factors etaikr
          do ns = 1, nshell
            i = nsh1(ns)
            j = nsh2(ns)
            do isym = 1, nsymat
              if (ns <= nsdia) then
                etaikr_batch((isym-1)*nshell + ns, idx) = volcub(kpt)
              else
                carg = czero
                do i1 = 1, 3
                  zktr = rrot(isym, i1, ns) - rbasis(i1, j) + rbasis(i1, i)
                  zktr = kp_gpu(i1)*zktr
                  carg = carg + zktr
                end do
                etaikr_batch((isym-1)*nshell + ns, idx) = volcub(kpt)*exp(carg*citpi)
              end if
            end do
          end do

          ! Fourier transform for gllken_batch (using rr)
          do i = 1, naez
            do m = 1, naclsmax
              do lm2 = 1, lmgf0d
                if (m <= nacls(cls(i))) then
                  if (atom(m, i) >= 0) then
                    ic = cls(i)
                    if (use_deci_onebulk) then
                      arg1 = -ci*tpi_gpu*rcls(1, m, ic)
                      arg2 = -ci*tpi_gpu*rcls(2, m, ic)
                      arg3 = -ci*tpi_gpu*rcls(3, m, ic)
                    else
                      arg1 = -ci*tpi_gpu*rr(1, ezoa(m, i))
                      arg2 = -ci*tpi_gpu*rr(2, ezoa(m, i))
                      arg3 = -ci*tpi_gpu*rr(3, ezoa(m, i))
                    end if

                    tt_gpu = kp_gpu(1)*arg1 + kp_gpu(2)*arg2 + kp_gpu(3)*arg3
                    if (calc_complex_bandstructure) then
                      tt_gpu = tt_gpu + ci*(kp_gpu(4)*arg1+kp_gpu(5)*arg2+kp_gpu(6)*arg3)
                    end if

                    eikr_gpu = exp(tt_gpu)*convpu_gpu

                    im = 1 + (m-1)*lmgf0d
                    am = 1 + (atom(m, i)-1)*lmgf0d
                    do lm1 = 1, lmgf0d
                      gllken_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) = &
                        gllken_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) + eikr_gpu * ginp(im + lm1 - 1, lm2, ic)
                    end do
                  end if
                end if
              end do
            end do
          end do

          if (symmetrize_gmat) then
            do i = 1, naez
              do m = 1, naclsmax
                do lm2 = 1, lmgf0d
                  if (m <= nacls(cls(i))) then
                    if (atom(m, i) >= 0) then
                      ic = cls(i)
                      if (use_deci_onebulk) then
                        arg1 = -ci*tpi_gpu*rcls(1, m, ic)
                        arg2 = -ci*tpi_gpu*rcls(2, m, ic)
                        arg3 = -ci*tpi_gpu*rcls(3, m, ic)
                      else
                        arg1 = -ci*tpi_gpu*rr(1, ezoa(m, i))
                        arg2 = -ci*tpi_gpu*rr(2, ezoa(m, i))
                        arg3 = -ci*tpi_gpu*rr(3, ezoa(m, i))
                      end if

                      tt_gpu = -kp_gpu(1)*arg1 - kp_gpu(2)*arg2 - kp_gpu(3)*arg3
                      if (calc_complex_bandstructure) then
                        tt_gpu = tt_gpu - ci*(kp_gpu(4)*arg1+kp_gpu(5)*arg2+kp_gpu(6)*arg3)
                      end if

                      eikr_gpu = exp(tt_gpu)*convpu_gpu

                      im = 1 + (m-1)*lmgf0d
                      am = 1 + (atom(m, i)-1)*lmgf0d
                      do lm1 = 1, lmgf0d
                        gllken1_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) = &
                          gllken1_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) + eikr_gpu * ginp(im + lm1 - 1, lm2, ic)
                      end do
                    end if
                  end if
                end do
              end do
            end do

            do j = 1, naez
              do i = 1, naez
                do lm2 = 1, lmgf0d
                  do lm1 = 1, lmgf0d
                    im = (i-1)*lmgf0d + lm1
                    jn = (j-1)*lmgf0d + lm2
                    gllken_batch(im, jn, idx) = (gllken_batch(im, jn, idx) + gllken1_batch(jn, im, idx)) * 0.5_dp
                  end do
                end do
              end do
            end do
          end if

          ! Fourier transform for gllkem_batch (using rrm)
          do i = 1, naez
            do m = 1, naclsmax
              do lm2 = 1, lmgf0d
                if (m <= nacls(cls(i))) then
                  if (atom(m, i) >= 0) then
                    ic = cls(i)
                    if (use_deci_onebulk) then
                      arg1 = -ci*tpi_gpu*rcls(1, m, ic)
                      arg2 = -ci*tpi_gpu*rcls(2, m, ic)
                      arg3 = -ci*tpi_gpu*rcls(3, m, ic)
                    else
                      arg1 = -ci*tpi_gpu*rrm(1, ezoa(m, i))
                      arg2 = -ci*tpi_gpu*rrm(2, ezoa(m, i))
                      arg3 = -ci*tpi_gpu*rrm(3, ezoa(m, i))
                    end if

                    tt_gpu = kp_gpu(1)*arg1 + kp_gpu(2)*arg2 + kp_gpu(3)*arg3
                    if (calc_complex_bandstructure) then
                      tt_gpu = tt_gpu + ci*(kp_gpu(4)*arg1+kp_gpu(5)*arg2+kp_gpu(6)*arg3)
                    end if

                    eikr_gpu = exp(tt_gpu)*convpu_gpu

                    im = 1 + (m-1)*lmgf0d
                    am = 1 + (atom(m, i)-1)*lmgf0d
                    do lm1 = 1, lmgf0d
                      gllkem_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) = &
                        gllkem_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) + eikr_gpu * ginp(im + lm1 - 1, lm2, ic)
                    end do
                  end if
                end if
              end do
            end do
          end do

          if (symmetrize_gmat) then
            do i = 1, naez
              do m = 1, naclsmax
                do lm2 = 1, lmgf0d
                  if (m <= nacls(cls(i))) then
                    if (atom(m, i) >= 0) then
                      ic = cls(i)
                      if (use_deci_onebulk) then
                        arg1 = -ci*tpi_gpu*rcls(1, m, ic)
                        arg2 = -ci*tpi_gpu*rcls(2, m, ic)
                        arg3 = -ci*tpi_gpu*rcls(3, m, ic)
                      else
                        arg1 = -ci*tpi_gpu*rrm(1, ezoa(m, i))
                        arg2 = -ci*tpi_gpu*rrm(2, ezoa(m, i))
                        arg3 = -ci*tpi_gpu*rrm(3, ezoa(m, i))
                      end if

                      tt_gpu = -kp_gpu(1)*arg1 - kp_gpu(2)*arg2 - kp_gpu(3)*arg3
                      if (calc_complex_bandstructure) then
                        tt_gpu = tt_gpu - ci*(kp_gpu(4)*arg1+kp_gpu(5)*arg2+kp_gpu(6)*arg3)
                      end if

                      eikr_gpu = exp(tt_gpu)*convpu_gpu

                      im = 1 + (m-1)*lmgf0d
                      am = 1 + (atom(m, i)-1)*lmgf0d
                      do lm1 = 1, lmgf0d
                        gllkem1_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) = &
                          gllkem1_batch(am + lm1 - 1, (i-1)*lmgf0d + lm2, idx) + eikr_gpu * ginp(im + lm1 - 1, lm2, ic)
                      end do
                    end if
                  end if
                end do
              end do
            end do

            do j = 1, naez
              do i = 1, naez
                do lm2 = 1, lmgf0d
                  do lm1 = 1, lmgf0d
                    im = (i-1)*lmgf0d + lm1
                    jn = (j-1)*lmgf0d + lm2
                    gllkem_batch(im, jn, idx) = (gllkem_batch(im, jn, idx) + gllkem1_batch(jn, im, idx)) * 0.5_dp
                  end do
                end do
              end do
            end do
          end if

          ! Combine gllken_batch and gllkem_batch into gllke_batch
          do i2 = 1, alm
            do i1 = 1, alm
              gllke_batch(i1, i2, idx) = (gllken_batch(i1, i2, idx) + gllkem_batch(i2, i1, idx)) * 0.5_dp
            end do
          end do

          ! Subtract tinvll
          do i1 = 1, naez
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                il1 = lmmaxd*(i1-1) + lm1
                il2 = lmmaxd*(i1-1) + lm2
                gllke_batch(il1, il2, idx) = gllke_batch(il1, il2, idx) - tinvll(lm1, lm2, i1)
              end do
            end do
          end do

        end do

        ! Call batched inversion on GPU
        call inversion_gpu_batched(gllke_batch, gtemp_batch, ipvt_batch, info_batch, batchSize)

        ! Parallel accumulation into gs
        !$acc parallel loop collapse(4) present(gs, gllke_batch, etaikr_batch) &
        !$acc& private(eta_val, idx, i, j, ilm, jlm)
        do ns = 1, nshell
          do isym = 1, nsymat
            do lm2 = 1, lmmaxd
              do lm1 = 1, lmmaxd
                do kpt = k_start, k_end
                  idx = kpt - k_start + 1
                  i = nsh1(ns)
                  j = nsh2(ns)
                  ilm = lmmaxd*(i-1) + lm1
                  jlm = lmmaxd*(j-1) + lm2
                  eta_val = etaikr_batch((isym-1)*nshell + ns, idx)
                  gs(lm1, lm2, isym, ns) = gs(lm1, lm2, isym, ns) + eta_val * gllke_batch(ilm, jlm, idx)
                end do
              end do
            end do
          end do
        end do
      end if

    else
#endif
    ! kpts loop
    do kpt = k_start, k_end
      gllke(:, :) = czero
      if (lly/=0) dgllke(:, :) = czero
      kp(1:3) = bzkp(1:3, kpt)

      ! overwrite kpt in case of rhoqtest (take only reduced set of kpts)
      if (write_rhoq_input) kp(1:3) = rhoq_kmask(1:3, kpt)

      etaikr(1:nsymat, 1:nshell) = volcub(kpt)
      ! -------------------------------------------------------------------------
      ! First NAEZ/NATYP elements of GS() are site-diagonal
      ! -------------------------------------------------------------------------
#ifdef CPP_HYBRID
      ! $omp do
#endif
      do ns = nsdia + 1, nshell
        i = nsh1(ns)
        j = nsh2(ns)
        do isym = 1, nsymat
          carg = czero
          do i1 = 1, 3
            zktr = rrot(isym, i1, ns) - rbasis(i1, j) + rbasis(i1, i)
            zktr = kp(i1)*zktr
            carg = carg + zktr
          end do
          etaikr(isym, ns) = etaikr(isym, ns)*exp(carg*citpi)
        end do
      end do
#ifdef CPP_HYBRID
      ! $omp end do
#endif

      bzkpk(1:3) = kp(1:3)
      bzkpk(4:6) = 0.0_dp

      ! -------------------------------------------------------------------------
      ! Fourier transformation
      ! -------------------------------------------------------------------------
#ifdef CPP_TIMING
      if (mythread==0 .and. t_inc%i_time>0) call timing_start('main1b - fourier')
#endif

      ! (rrm is now pre-computed on GPU before the k-loop, no CPU assignment needed)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! KREL .EQ. 0/1
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (krel==0) then

        ! $omp single
        gllken(:, :) = czero
        call dlke0(gllken, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, bzkpk, rcls, ginp)

        gllkem(:, :) = czero
        call dlke0(gllkem, alat, naez, cls, nacls, naclsmax, rrm, ezoa, atom, bzkpk, rcls, ginp)
        ! ----------------------------------------------------------------------
        ! LLY Lloyd
        ! Fourier for dGref/dE for Lloyds formula, repeat the above allocation
        ! and Fourier transform for the derivatives.
        ! ----------------------------------------------------------------------
        if (lly/=0) then
          dgllken(:, :) = czero
          call dlke0(dgllken, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, bzkpk, rcls, dginp)

          dgllkem(:, :) = czero
          call dlke0(dgllkem, alat, naez, cls, nacls, naclsmax, rrm, ezoa, atom, bzkpk, rcls, dginp)
        end if
        ! $omp end single
        ! ----------------------------------------------------------------------
        ! LLY Lloyd
        ! --------------------------------------------------------------------
        if (.not. use_Chebychev_solver .or. decouple_spin_cheby) then
          ! $omp single
          do i2 = 1, alm
            do i1 = 1, alm
              gllke(i1, i2) = (gllken(i1,i2)+gllkem(i2,i1))*0.5_dp
              if (lly/=0) then
                dgllke(i1, i2) = (dgllken(i1,i2)+dgllkem(i2,i1))*0.5_dp ! LLY Lloyd
              end if
            end do
          end do
          ! $omp end single
        else                       ! (.NOT.use_Chebychev_solver)
          ! $omp single
          do i2 = 1, almgf0
            do i1 = 1, almgf0
              gllken(i1, i2) = (gllken(i1,i2)+gllkem(i2,i1))*0.5_dp
              if (lly/=0) then
                dgllken(i1, i2) = (dgllken(i1,i2)+dgllkem(i2,i1))*0.5_dp ! LLY Lloyd
              end if
            end do
          end do
          ! $omp end single

          ! bigger GLLKE matrix and rearrange with atom block
#ifdef CPP_HYBRID
          ! $omp do
#endif
          do iq1 = 1, naez
            do iq2 = 1, naez
              ioff1 = lmmaxd*(iq1-1)
              joff1 = lmgf0d*(iq1-1)
              ioff2 = lmmaxd*(iq2-1)
              joff2 = lmgf0d*(iq2-1)
              do lm1 = 1, lmgf0d
                do lm2 = 1, lmgf0d
                  gllke(ioff1+lm1, ioff2+lm2) = gllken(joff1+lm1, joff2+lm2)
                  gllke(ioff1+lm1+lmgf0d, ioff2+lm2+lmgf0d) = gllken(joff1+lm1, joff2+lm2)
                  if (lly/=0) then ! LLY Lloyd
                    dgllke(ioff1+lm1, ioff2+lm2) = & ! LLY Lloyd
                      dgllken(joff1+lm1, joff2+lm2) ! LLY Lloyd
                    dgllke(ioff1+lm1+lmgf0d, ioff2+lm2+lmgf0d) = & ! LLY Lloyd
                      dgllken(joff1+lm1, joff2+lm2) ! LLY Lloyd
                  end if           ! (LLY.NE.0)                                       ! LLY Lloyd
                end do
              end do
            end do
          end do
#ifdef CPP_HYBRID
          ! $omp end do
#endif
        end if                     ! (.NOT.use_Chebychev_solver)
        ! ----------------------------------------------------------------------
        ! LLY Lloyd At this point DGLLKE contains the Fourier transform of the dGref/dE
        ! ----------------------------------------------------------------------
      else                         ! (KREL.EQ.0)
        ! ----------------------------------------------------------------------
        ! LLY Lloyd Not implementing Lloyds formula for KREL=1 (Dirac ASA)
        ! ----------------------------------------------------------------------
        ! $omp single
        gllke0(:, :) = czero
        gllke0m(:, :) = czero
        call dlke0(gllke0, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, bzkpk, rcls, ginp)
        call dlke0(gllke0m, alat, naez, cls, nacls, naclsmax, rrm, ezoa, atom, bzkpk, rcls, ginp)

        do i2 = 1, almgf0
          do i1 = 1, almgf0
            gllke0(i1, i2) = (gllke0(i1,i2)+gllke0m(i2,i1))*0.5_dp
          end do
        end do
        ! $omp end single
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Double the GLLKE0 matrix and transform to the REL representation
        ! ==> GLLKE
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CPP_HYBRID
        ! $omp do
#endif
        do iq1 = 1, naez
          do iq2 = 1, naez
            ioff1 = lmmaxd*(iq1-1)
            joff1 = lmgf0d*(iq1-1)
            ioff2 = lmmaxd*(iq2-1)
            joff2 = lmgf0d*(iq2-1)
            ! ---------------------------------------------------------------
            do ikm2 = 1, lmmaxd
              do ikm1 = 1, lmmaxd
                csum1 = czero
                do is = 1, 2
                  n1 = nrrel(is, ikm1)
                  n2 = nrrel(is, ikm2)
                  do i1 = 1, n1
                    j1 = irrel(i1, is, ikm1) + joff1
                    csum2 = czero
                    do i2 = 1, n2
                      j2 = irrel(i2, is, ikm2) + joff2
                      csum2 = csum2 + gllke0(j1, j2)*srrel(i2, is, ikm2)
                    end do
                    csum1 = csum1 + conjg(srrel(i1,is,ikm1))*csum2
                  end do
                end do
                gllke(ioff1+ikm1, ioff2+ikm2) = csum1
              end do
            end do
            ! ----------------------------------------------------------------
          end do
        end do
#ifdef CPP_HYBRID
        ! $omp end do
#endif
      end if                       ! (KREL.EQ.0)
#ifdef CPP_TIMING
      if (mythread==0 .and. t_inc%i_time>0) call timing_pause('main1b - fourier')
#endif

      if (lly/=0) grefllke(1:alm, 1:alm) = gllke(1:alm, 1:alm) ! LLY Save k-dependent Gref

      if (ideci==1) then
        call decimate(gllke, naez, tinvbup, tinvbdown, vacflag, factl, nlbasis, nrbasis)
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Construct the matrix M=[-(t)^-1 + G^r] and store it
      ! in the same matrix GLLKE where G^r was stored.
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not. use_virtual_atoms) then
        ! $omp single
        ! Run on CPU since gllke assembly and subsequent inversion are on CPU
        ! (avoids host-device data roundtrips for trivial size-81 matrix subtraction)
        do i1 = 1, naez
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              il1 = lmmaxd*(i1-1) + lm1
              il2 = lmmaxd*(i1-1) + lm2
              gllke(il1, il2) = gllke(il1, il2) - tinvll(lm1, lm2, i1)
            end do
          end do
        end do
        ! $omp end single

#ifdef CPP_TIMING
        if (mythread==0 .and. t_inc%i_time>0) call timing_start('main1b - inversion')
#endif
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Perform the inversion of matrix M
        ! the output is the scattering path operator TAU stored in GLLKE
        ! Actually -TAU, because TAU = (Deltat^-1 - Gref)^-1
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! $omp single
        if (lly/=0 .and. invmod/=3) then           ! If LLY, full inversion is needed
          call inversion(gllke, 0, icheck) ! LLY
        else ! MdSD: godfrin options are set to full inversion if LLY is true
          call inversion(gllke, invmod, icheck)
        end if
        ! $omp end single
#ifdef CPP_TIMING
        if (mythread==0 .and. t_inc%i_time>0) call timing_pause('main1b - inversion')
#endif
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LLY Lloyd
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (lly/=0) then
          ! -------------------------------------------------------------------
          ! LLY  Prepare quantities for Lloyds formula.
          ! LLY  Needed is Trace[ (1-Gref * Deltat)^-1 * d(1-Gref * Deltat)/dE ] (PhD Thiess Eq.5.38)
          ! LLY  where Deltat = t-tref. This is re-written as:
          ! LLY  -Trace[ Tau * ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 ) ]
          ! LLY  where Tau is the scattering path operator Tau = (Deltat^-1 - Gref)^-1
          ! LLY  (negative of array GLLKE) and (t-tref)^-1 is in array TINVLL.
          ! LLY  The quantities Gref, dGref/dE, dt/dE have been prepared by main1a.
          ! LLY  Quantity dtref/dE is in array DTREFLL
          ! -------------------------------------------------------------------

          ! First set up (dt/dE - dtref/dE) Deltat^-1, store in array t_aux
          ! $omp single
          do i1 = 1, naez
            ! GAUX1 = dt/dE-dtref/dE
            gaux1(1:lmmaxd, 1:lmmaxd) = (1.0_dp/cfctor)*(dtmatll(1:lmmaxd,1:lmmaxd,i1)-dtrefll(1:lmmaxd,1:lmmaxd,refpot(i1)))
            gaux2(1:lmmaxd, 1:lmmaxd) = tinvll(1:lmmaxd, 1:lmmaxd, i1)
            ! T_AUX = (dt/dE-dtref/dE)* Deltat^-1
            call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, gaux1, lmmaxd, gaux2, lmmaxd, czero, gaux3, lmmaxd)
            t_aux(1:lmmaxd, 1:lmmaxd, i1) = gaux3(1:lmmaxd, 1:lmmaxd)
          end do
          ! -------------------------------------------------------------------
          ! Now perform dGref/dE + Gref * t_aux
          ! (Gref is ALM*ALM ; t_aux site-diagonal LMMAXD*LMMAXD)
          ! -------------------------------------------------------------------
          do j1 = 1, naez          ! Loop over columns of Gref
            jl1 = lmmaxd*(j1-1) + 1
            jl2 = lmmaxd*(j1-1) + lmmaxd
            gaux3(1:lmmaxd, 1:lmmaxd) = t_aux(1:lmmaxd, 1:lmmaxd, j1)
            do i1 = 1, naez        ! Loop over rows of Gref
              il1 = lmmaxd*(i1-1) + 1
              il2 = lmmaxd*(i1-1) + lmmaxd
              ! Copy to small matrices
              gaux1(1:lmmaxd, 1:lmmaxd) = grefllke(il1:il2, jl1:jl2)
              gaux2(1:lmmaxd, 1:lmmaxd) = dgllke(il1:il2, jl1:jl2)
              ! GAUX2 = GAUX2 + GAUX1 * T_AUX
              call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, gaux1, lmmaxd, gaux3, lmmaxd, cone, gaux2, lmmaxd)
              ! Copy back to large matrix, use again array DGLLKE
              ! (the I1-J1 block is not needed any more)
              dgllke(il1:il2, jl1:jl2) = gaux2(1:lmmaxd, 1:lmmaxd)
            end do
          end do
          ! $omp end single

          ! full matrix multiple
          ! ALLOCATE(GLLKE0(ALM,ALM))
          ! GLLKE0=CZERO
          ! DO I1=1,NAEZ
          ! DO LM1 = 1,LMMAXD
          ! DO LM2 = 1,LMMAXD
          ! IL1 = LMMAXD*(I1-1)+LM1
          ! IL2 = LMMAXD*(I1-1)+LM2
          ! GLLKE0(IL1,IL2)= T_AUX(LM1,LM2,I1)
          ! ENDDO
          ! ENDDO
          ! ENDDO
          ! CALL ZGEMM('N','N',ALM,ALM,ALM,CONE,GREFLLKE,
          ! &                 ALM,GLLKE0,ALM,CONE,DGLLKE,ALM)
          ! DEALLOCATE(GLLKE0)

          ! Now array DGLLKE contains
          ! ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 )
          ! Build trace of tau * DGLLKE, -tau is conained in GLLKE.
          trace = czero
          do il1 = 1, alm
            do il2 = 1, alm
              trace = trace + gllke(il1, il2)*dgllke(il2, il1)
            end do
          end do
          lly_grtr_k = trace
        end if                     ! (LLY.NE.0)
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LLY Lloyd
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else                         ! .not. use_virtual_atoms
        ! LLY Lloyd formula not built in yet for viratoms
        gllke0v(1:alm, 1:alm) = gllke(1:alm, 1:alm)

        do i1 = 1, naez
          il1 = (i1-1)*lmmaxd + 1
          ! GLLKETV = -GLLKE0V * TINVLL,
          ! where TINVLL contains (t-tref) and not 1/(t-tref) in case of opt VIRATOMS
          ! tref=0 for each vir. atom.
          call zgemm('N', 'N', ndim_slabinv, lmmaxd, lmmaxd, -cone, gllke0v(1,il1), alm, tinvll(1,1,i1), lmgf0d, czero, gllketv(1,1), alm)
          call zcopy(alm*lmmaxd, gllketv(1,1), 1, gllke0v2(1,il1), 1)
        end do
        ! ----------------------------------------------------------------------
        ! Solve (1-gt)G=g instead of [Gref - t^-1]^-1 for viratoms
        ! because for a virtual atom t=0, t^-1 undefined.
        ! ----------------------------------------------------------------------
        call gtdyson(gllke0v2, gllke, ndim_slabinv, alm, alm)
      end if                       ! .not. use_virtual_atoms
      ! -------------------------------------------------------------------------
      ! Global sum on array gs
      ! -------------------------------------------------------------------------
      ! no omp at this loop because of rhoq output
      ! $omp single
      do ns = 1, nshell
        i = nsh1(ns)
        j = nsh2(ns)
        ilm = lmmaxd*(i-1) + 1
        jlm = lmmaxd*(j-1)

        do lm = 1, lmmaxd
          call zcopy(lmmaxd, gllke(ilm,jlm+lm), 1, g(1,lm), 1)
        end do
        ! ----------------------------------------------------------------------
        do isym = 1, nsymat
          do lm2 = 1, lmmaxd
            do lm1 = 1, lmmaxd
              gs(lm1, lm2, isym, ns) = gs(lm1, lm2, isym, ns) + etaikr(isym, ns)*g(lm1, lm2)
            end do
          end do
          if (write_rhoq_input) then
            call rhoq_saveg(nscoef, rhoq_kmask, kpt, k_end, kp, i, j, mu, imin, iatomimp, lmmaxd, g)
          end if
        end do                     ! isym
        ! ----------------------------------------------------------------------
      end do                       ! ns
      ! $omp end single
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LLY Lloyd Integration
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (lly/=0) then
        lly_grtr = lly_grtr + lly_grtr_k*volcub(kpt)*nsymat
      end if

      ! Update statusbar
      if (mythread==0) then
#ifdef CPP_MPI
        ! $omp critical
        if ((((k_end-k_start)/50)==0 .or. mod(kpt-k_start,(k_end-k_start)/50)==0) .and. t_inc%i_write>0) write (1337, fmt=110, advance='no')
        ! $omp end critical
#else
        !if (((nofks/50)==0 .or. mod(kpt,nofks/50)==0) .and. t_inc%i_write>0) write (1337, fmt=110, advance='no')
        if (t_inc%i_write>0) then
          if (nofks/50==0) then
            write (1337, fmt=110, advance='no')
          else
            if (mod(kpt,nofks/50)==0) then
              write (1337, fmt=110, advance='no')
            end if
          end if
        end if
#endif
      end if                       ! mythread==0

    end do                         ! KPT = 1,NOFKS   end K-points loop
#ifdef CPP_GPU
    end if
#endif
    !$acc end data
100 format ('                 |')  ! status bar
110 format ('|')                   ! status bar
    ! $omp critical
    if (t_inc%i_write>0) write (1337, *) ! finalize status bar
    ! $omp end critical

    ! ----------------------------------------------------------------------
    ! deallocattions of work arrays
    if (krel==0) then
      if (mythread==0) then
        deallocate (gllken, stat=i_stat)
        deallocate (gllkem, stat=i_stat)
      end if
      if (lly/=0) then
        if (mythread==0) then
          deallocate (dgllken, stat=i_stat)
          deallocate (dgllkem, stat=i_stat)
        end if
      end if
    else
      if (mythread==0) then
        deallocate (gllke0, stat=i_stat)
        deallocate (gllke0m, stat=i_stat)
      end if
    end if                         ! (KREL.EQ.0)
    ! ----------------------------------------------------------------------

#ifdef CPP_HYBRID
    ! $omp end parallel
#endif

    if (write_rhoq_input) then
#ifdef CPP_TIMING
      call timing_start('main1b - kkrmat01 - writeout_rhoq')
#endif
      call rhoq_write_tau0(nofks, nshell, nsh1, nsh2, nsymat, nscoef, mu, iatomimp, kmask, lmmaxd, bzkp, imin)
#ifdef CPP_TIMING
      call timing_stop('main1b - kkrmat01 - writeout_rhoq')
#endif
    end if                         ! write_rhoq_input

    tracet = czero
    if (lly==2) then
      ! Add trace of (t-tref)^-1 * d(t-tref)/dE. Remember that in this case
      ! Tr(alpha^-1 d alpha/dE) should be subtracted and
      ! Tr(alpha_ref^-1 d alpha_ref/dE) should be added.
      do i1 = 1, naez
        gaux1(1:lmmaxd, 1:lmmaxd) = cfctor* & ! (1.D0/CFCTOR) *
          (dtmatll(1:lmmaxd,1:lmmaxd,i1)-dtrefll(1:lmmaxd,1:lmmaxd,refpot(i1)))
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            tracet = tracet + gaux1(lm1, lm2)*tinvll(lm2, lm1, i1)
          end do
        end do
      end do
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! deallocate arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i_all = -product(shape(gllke))*kind(gllke)
#ifdef CPP_GPU
    if (batchSize > 0) then
      !$acc exit data delete(gllke_batch, gtemp_batch, ipvt_batch, info_batch, etaikr_batch)
      deallocate(gllke_batch, gtemp_batch, ipvt_batch, info_batch, etaikr_batch, stat=i_stat)
      if (allocated(gllken_batch)) then
        !$acc exit data delete(gllken_batch, gllkem_batch)
        deallocate(gllken_batch, gllkem_batch, stat=i_stat)
        if (allocated(dgllken_batch)) then
          !$acc exit data delete(dgllken_batch, dgllkem_batch, dgllke_batch, grefllke_batch)
          deallocate(dgllken_batch, dgllkem_batch, dgllke_batch, grefllke_batch, stat=i_stat)
        end if
        if (allocated(gllken1_batch)) then
          !$acc exit data delete(gllken1_batch, gllkem1_batch)
          deallocate(gllken1_batch, gllkem1_batch, stat=i_stat)
        end if
      else if (allocated(gllke0_batch)) then
        !$acc exit data delete(gllke0_batch, gllke0m_batch)
        deallocate(gllke0_batch, gllke0m_batch, stat=i_stat)
        if (allocated(gllke01_batch)) then
          !$acc exit data delete(gllke01_batch, gllke0m1_batch)
          deallocate(gllke01_batch, gllke0m1_batch, stat=i_stat)
        end if
      end if
    end if
#endif

    deallocate (gllke, stat=i_stat)
    ! LLY Lloyd
    if (lly/=0) then
      i_all = -product(shape(dgllke))*kind(dgllke)
      deallocate (dgllke, stat=i_stat)
      i_all = -product(shape(grefllke))*kind(grefllke)
      deallocate (grefllke, stat=i_stat)
    end if
    ! ----------------------------------------------------------------------------

#ifdef CPP_MPI
    if (.not. use_qdos) then
      do ns = 1, nshell
        iwork = lmmaxd*lmmaxd*nsymaxd
        work = czero
        call mpi_allreduce(gs(1,1,1,ns), work, iwork, mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_ie, ierr)
        call zcopy(iwork, work, 1, gs(1,1,1,ns), 1)
      end do

      if (lly/=0) then
        iwork = 1
        work = czero
        call mpi_allreduce(lly_grtr, work, iwork, mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_ie, ierr)
        call zcopy(iwork, work(1,1,1), 1, lly_grtr, 1)
      end if

      if (lly==2) then
        iwork = 1
        work = czero
        call mpi_allreduce(tracet, work, iwork, mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_ie, ierr)
        call zcopy(iwork, work, 1, tracet, 1)
      end if
    end if                         ! .not.use_qdos
#endif

    if (print_program_flow .and. (t_inc%i_write>0)) write (1337, *) '<<< KKRMAT1'

  end subroutine kkrmat01

  !-------------------------------------------------------------------------------
  !> Summary: Solve the Dyson equation \((1-g t)  G = g\)
  !> Category: KKRhost, k-points, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine gtdyson(gtmat, gmat, ndim, lmgf0d, ngd)

    use :: mod_constants
    use :: mod_datatypes, only: dp

    implicit none
    ! .. Input variables
    integer, intent (in) :: ngd
    integer, intent (in) :: ndim
    integer, intent (in) :: lmgf0d !! (LMAX+1)**2
    ! .. In/Out variables
    complex (kind=dp), dimension (ngd, lmgf0d), intent (inout) :: gmat
    complex (kind=dp), dimension (ngd, ngd), intent (inout) :: gtmat
    ! .. Local variables
    integer :: i, info
    integer, dimension (ngd) :: ipvt
    ! ..

    do i = 1, ndim
      gtmat(i, i) = cone + gtmat(i, i) ! GTMAT= 1 - G * T
    end do
    ! ----------------------------------------------------------------------------
    ! SOLVE THE SYSTEM OF LINEAR EQUATIONS
    ! ----------------------------------------------------------------------------
    call zgetrf(ndim, ndim, gtmat, ngd, ipvt, info)
    call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, gmat, ngd, info)
    return

  end subroutine gtdyson

end module mod_kkrmat01
