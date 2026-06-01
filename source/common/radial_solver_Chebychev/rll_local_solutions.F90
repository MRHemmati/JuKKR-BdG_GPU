!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the local regular solutions
!> Author: 
!> Calculation of the local regular solutions
!>
!> 1. Prepare `VJLR`, `VNL`, `VHLR`, which appear in the integrands `TAU(K,IPAN)` 
!> is used instead of `TAU(K,IPAN)**2`, which directly gives `RLL(r)` and `SLL(r)` 
!> multiplied with `r`. `TAU` is the radial mesh.
!> 2. Prepare the source terms `YR`, `ZR`, `YI`, `ZI` because of the conventions used
!> by Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
!> a factor \(\sqrt(E)\) is included in the source terms this factor is removed by 
!> the definition of `ZSLC1SUM` given below
!> \begin{equation}
!> vjlr = \kappa J V = \kappa r j V
!> \end{equation}
!> \begin{equation}
!> vhlr = \kappa H V = \kappa r h V
!> \end{equation}
!> i.e. prepare terms \(\kappa JDV\), \(\kappa HDV\) appearing in 5.11, 5.12.
!>
!> Then determine the local solutions solve the equations 
!> \begin{equation}
!> SLV \times YRLL=S
!> \end{equation}
!> and 
!> \begin{equation}
!> SLV \times ZRLL=C
!> \end{equation}
!> and
!> \begin{equation}
!> SRV \times YILL=C
!> \end{equation}
!> and
!> \begin{equation}
!> SRV \times ZILL=S
!> \end{equation}
!> i.e., solve system \(A U=J\), see eq. 5.68.
!------------------------------------------------------------------------------------
module mod_rll_local_solutions

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the local regular solutions
  !> Author: 
  !> Category: single-site, KKRhost
  !> Deprecated: False 
  !> Calculation of the local regular solutions
  !>
  !> 1. Prepare `VJLR`, `VNL`, `VHLR`, which appear in the integrands `TAU(K,IPAN)` 
  !> is used instead of `TAU(K,IPAN)**2`, which directly gives `RLL(r)` and `SLL(r)` 
  !> multiplied with `r`. `TAU` is the radial mesh.
  !> 2. Prepare the source terms `YR`, `ZR`, `YI`, `ZI` because of the conventions used
  !> by Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
  !> a factor \(\sqrt(E)\) is included in the source terms this factor is removed by 
  !> the definition of `ZSLC1SUM` given below
  !> \begin{equation}
  !> vjlr = \kappa J V = \kappa r j V
  !> \end{equation}
  !> \begin{equation}
  !> vhlr = \kappa H V = \kappa r h V
  !> \end{equation}
  !> i.e. prepare terms \(\kappa JDV\), \(\kappa HDV\) appearing in 5.11, 5.12.
  !>
  !> Then determine the local solutions solve the equations 
  !> \begin{equation}
  !> SLV \times YRLL=S
  !> \end{equation}
  !> and 
  !> \begin{equation}
  !> SLV \times ZRLL=C
  !> \end{equation}
  !> and
  !> \begin{equation}
  !> SRV \times YILL=C
  !> \end{equation}
  !> and
  !> \begin{equation}
  !> SRV \times ZILL=S
  !> \end{equation}
  !> i.e., solve system \(A U=J\), see eq. 5.68.
  !-------------------------------------------------------------------------------
  subroutine rll_local_solutions(vll,tau,drpan2,cslc1,slc1sum,mrnvy,mrnvz,mrjvy,    &
    mrjvz,yrf,zrf,ncheb,ipan,lmsize,lmsize2,nrmax,nvec,jlk_index,hlk,jlk,hlk2,jlk2, &
    gmatprefactor,cmoderll,lbessel,use_sratrick1)

    use :: mod_datatypes, only: dp
    use :: mod_sll_local_solutions, only: svpart
    use :: mod_constants, only: cone,czero

    implicit none
    integer :: ncheb               ! number of chebyshev nodes
    integer :: lmsize              ! lm-components * nspin
    integer :: lmsize2             ! lmsize * nvec
    integer :: nvec                ! spinor integer
    ! nvec=1 non-rel, nvec=2 for sra and dirac
    integer :: nrmax               ! total number of rad. mesh points

    integer :: lbessel, use_sratrick1 ! dimensions etc., needed only for host code interface

    ! running indices
    integer :: ivec, ivec2
    integer :: l1, l2, lm1, lm2, lm3
    integer :: info, icheb2, icheb, ipan, mn, nplm

    ! source terms
    complex (kind=dp) :: gmatprefactor ! prefactor of green function
    ! non-rel: = kappa = sqrt e

    complex (kind=dp) :: hlk(lbessel, nrmax), jlk(lbessel, nrmax), hlk2(lbessel, nrmax), jlk2(lbessel, nrmax)


    integer :: jlk_index(2*lmsize)


    character (len=1) :: cmoderll
    ! cmoderll ="1" : op( )=identity       for reg. solution
    ! cmoderll ="T" : op( )=transpose in L for reg. solution

    complex (kind=dp) :: vll(lmsize*nvec, lmsize*nvec, nrmax) ! potential term in 5.7

    complex (kind=dp) :: mrnvy(lmsize, lmsize), mrnvz(lmsize, lmsize), mrjvy(lmsize, lmsize), mrjvz(lmsize, lmsize), yrf(lmsize2, lmsize, 0:ncheb), zrf(lmsize2, lmsize, 0:ncheb)
    complex (kind=dp) :: slv(0:ncheb, lmsize2, 0:ncheb, lmsize2), slv1(0:ncheb, lmsize, 0:ncheb, lmsize), yrll1(0:ncheb, lmsize, lmsize), zrll1(0:ncheb, lmsize, lmsize), &
      yrll2(0:ncheb, lmsize, lmsize), zrll2(0:ncheb, lmsize, lmsize), yrll(0:ncheb, lmsize2, lmsize), zrll(0:ncheb, lmsize2, lmsize), vjlr(lmsize, lmsize2, 0:ncheb), &
      vhlr(lmsize, lmsize2, 0:ncheb), vjlr_yrll1(lmsize, lmsize), vhlr_yrll1(lmsize, lmsize), vjlr_zrll1(lmsize, lmsize), vhlr_zrll1(lmsize, lmsize), yrll1temp(lmsize, lmsize), &
      zrll1temp(lmsize, lmsize)

    complex (kind=dp) :: jlmkmn(0:ncheb, lmsize2, 0:ncheb), hlmkmn(0:ncheb, lmsize2, 0:ncheb)

    ! chebyshev arrays
    complex (kind=dp) :: zslc1sum(0:ncheb)
    real (kind=dp) :: drpan2
    real (kind=dp) :: cslc1(0:ncheb, 0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
      tau(0:ncheb), &              ! Radial mesh point
      slc1sum(0:ncheb), taucslcr, tau_icheb
    complex (kind=dp) :: gf_tau_icheb
    complex (kind=dp) :: tmp_nvy, tmp_jvy, tmp_nvz, tmp_jvz

    integer :: ipiv(0:ncheb, lmsize2)
    integer :: use_sratrick

    if (cmoderll/='1' .and. cmoderll/='T') then
      stop '[rll-loc] mode not known'
    end if

    ! GPU: enter OpenACC data region
    !$acc data copyin(tau, jlk_index, jlk2, vll, hlk2, jlk, hlk, cslc1, slc1sum) &
    !$acc copy(yrf, zrf, mrnvy, mrnvz, mrjvy, mrjvz) &
    !$acc create(vjlr, vhlr, yrll, zrll, yrll1, zrll1, yrll2, zrll2, slv, slv1, jlmkmn, hlmkmn, yrll1temp, zrll1temp, vhlr_yrll1, vhlr_zrll1, vjlr_yrll1, vjlr_zrll1, zslc1sum)

    if (lmsize==1) then
      use_sratrick = 0
    else
      use_sratrick = use_sratrick1
    end if

    ! initialization

    ! GPU: array initialization on GPU
    !$acc kernels
    vhlr = czero
    vjlr = czero

    if (use_sratrick==0) then
      yrll = czero
      zrll = czero
    else
      yrll1 = czero
      zrll1 = czero
      yrll2 = czero
      zrll2 = czero
    end if
    !$acc end kernels

    ! ---------------------------------------------------------------------
    ! 1. prepare VJLR, VNL, VHLR, which appear in the integrands
    ! TAU(K,IPAN) is used instead of TAU(K,IPAN)**2, which directly gives
    ! RLL(r) and SLL(r) multiplied with r. TAU is the radial mesh.

    ! 2. prepare the source terms YR, ZR, YI, ZI
    ! because of the conventions used by
    ! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
    ! a factor sqrt(E) is included in the source terms
    ! this factor is removed by the definition of ZSLC1SUM given below

    ! vjlr = \kappa * J * V = \kappa * r * j *V
    ! vhlr = \kappa * H * V = \kappa * r * h *V

    ! i.e. prepare terms kappa*J*DV, kappa*H*DV appearing in 5.11, 5.12.

    ! GPU: Setup vjlr and vhlr loops
    !$acc parallel loop private(mn, tau_icheb, gf_tau_icheb, l1, l2)
    do icheb = 0, ncheb
      mn = ipan*ncheb + ipan - icheb
      tau_icheb = tau(icheb)
      gf_tau_icheb = gmatprefactor*tau_icheb
      if (cmoderll=='1') then
        do ivec2 = 1, nvec
          do lm2 = 1, lmsize
            do ivec = 1, nvec
              do lm1 = 1, lmsize
                l1 = jlk_index(lm1+lmsize*(ivec-1))
                vjlr(lm1, lm2+lmsize*(ivec2-1), icheb) = vjlr(lm1, lm2+lmsize*(ivec2-1), icheb) + gf_tau_icheb*jlk2(l1, mn)*vll(lm1+lmsize*(ivec-1), lm2+lmsize*(ivec2-1), mn)
                vhlr(lm1, lm2+lmsize*(ivec2-1), icheb) = vhlr(lm1, lm2+lmsize*(ivec2-1), icheb) + gf_tau_icheb*hlk2(l1, mn)*vll(lm1+lmsize*(ivec-1), lm2+lmsize*(ivec2-1), mn)
              end do
            end do
          end do
        end do
      else ! cmoderll=='T'
        do ivec2 = 1, nvec
          do lm2 = 1, lmsize
            do ivec = 1, nvec
              do lm1 = 1, lmsize
                l1 = jlk_index(lm1+lmsize*(ivec-1))
                vjlr(lm1, lm2+lmsize*(ivec2-1), icheb) = vjlr(lm1, lm2+lmsize*(ivec2-1), icheb) + gf_tau_icheb*jlk2(l1, mn)*vll(lm2+lmsize*(ivec-1), lm1+lmsize*(ivec2-1), mn)
                vhlr(lm1, lm2+lmsize*(ivec2-1), icheb) = vhlr(lm1, lm2+lmsize*(ivec2-1), icheb) + gf_tau_icheb*hlk2(l1, mn)*vll(lm2+lmsize*(ivec-1), lm1+lmsize*(ivec2-1), mn)
              end do
            end do
          end do
        end do                     ! nvec
      end if

      ! calculation of the J (and H) matrix according to equation 5.69 (2nd eq.)
      if (use_sratrick==0) then
        do ivec = 1, nvec          ! index for large/small component
          do lm1 = 1, lmsize
            l1 = jlk_index(lm1+lmsize*(ivec-1))
            yrll(icheb, lm1+lmsize*(ivec-1), lm1) = tau_icheb*jlk(l1, mn)
            zrll(icheb, lm1+lmsize*(ivec-1), lm1) = tau_icheb*hlk(l1, mn)
          end do
        end do                     ! ivec=1,nvec
      else if (use_sratrick==1) then
        do lm1 = 1, lmsize
          l1 = jlk_index(lm1)
          l2 = jlk_index(lm1+lmsize)
          yrll1(icheb, lm1, lm1) = tau_icheb*jlk(l1, mn)
          zrll1(icheb, lm1, lm1) = tau_icheb*hlk(l1, mn)
          yrll2(icheb, lm1, lm1) = tau_icheb*jlk(l2, mn)
          zrll2(icheb, lm1, lm1) = tau_icheb*hlk(l2, mn)
        end do
      end if
    end do                         ! icheb

    ! calculation of A in 5.68
    if (use_sratrick==0) then
      ! GPU: setup slv kernel matrix
      !$acc parallel loop collapse(2) private(taucslcr, mn, lm1, l1)
      do icheb2 = 0, ncheb
        do icheb = 0, ncheb
          taucslcr = tau(icheb)*cslc1(icheb, icheb2)*drpan2
          mn = ipan*ncheb + ipan - icheb
          do lm2 = 1, lmsize2
            do ivec = 1, nvec
              do lm3 = 1, lmsize
                lm1 = lm3 + (ivec-1)*lmsize
                l1 = jlk_index(lm1)
                slv(icheb, lm1, icheb2, lm2) = taucslcr*(jlk(l1,mn)*vhlr(lm3,lm2,icheb2)-hlk(l1,mn)*vjlr(lm3,lm2,icheb2))
              end do
            end do
          end do
        end do
      end do
      ! GPU: add diagonal identity
      !$acc parallel loop collapse(2)
      do lm1 = 1, lmsize2
        do icheb = 0, ncheb
          slv(icheb, lm1, icheb, lm1) = slv(icheb, lm1, icheb, lm1) + 1.0_dp
        end do
      end do
    else if (use_sratrick==1) then
      ! GPU: setup jlmkmn and hlmkmn
      !$acc parallel loop collapse(2) private(taucslcr, mn, l1)
      do icheb2 = 0, ncheb
        do icheb = 0, ncheb
          taucslcr = tau(icheb)*cslc1(icheb, icheb2)*drpan2
          mn = ipan*ncheb + ipan - icheb
          do lm1 = 1, lmsize
            l1 = jlk_index(lm1)
            jlmkmn(icheb, lm1, icheb2) = -taucslcr*jlk(l1, mn)
            hlmkmn(icheb, lm1, icheb2) = -taucslcr*hlk(l1, mn)
          end do
        end do
      end do

      ! GPU: call svpart in parallel
      !$acc parallel loop collapse(2)
      do lm2 = 1, lmsize
        do icheb2 = 0, ncheb
          call svpart(slv1(0,1,icheb2,lm2), jlmkmn(0,1,icheb2), hlmkmn(0,1,icheb2), vhlr(1,lm2,icheb2), vjlr(1,lm2,icheb2), ncheb, lmsize)
        end do
      end do
      ! GPU: add diagonal identity
      !$acc parallel loop collapse(2)
      do lm1 = 1, lmsize
        do icheb = 0, ncheb
          slv1(icheb, lm1, icheb, lm1) = slv1(icheb, lm1, icheb, lm1) + 1.0_dp
        end do
      end do

    else
      stop '[rll-loc] error in inversion'
    end if

    ! -------------------------------------------------------
    ! determine the local solutions
    ! solve the equations SLV*YRLL=S and SLV*ZRLL=C
    ! and SRV*YILL=C and SRV*ZILL=S
    ! i.e., solve system A*U=J, see eq. 5.68.

    if (use_sratrick==0) then

      nplm = (ncheb+1)*lmsize2

      if (cmoderll/='0') then
        call zgetrf(nplm, nplm, slv, nplm, ipiv, info)
        if (info/=0) stop 'rll-loc: zgetrf'
        call zgetrs('n', nplm, lmsize, slv, nplm, ipiv, yrll, nplm, info)
        call zgetrs('n', nplm, lmsize, slv, nplm, ipiv, zrll, nplm, info)
      end if

    else if (use_sratrick==1) then

      nplm = (ncheb+1)*lmsize

      call zgetrf(nplm, nplm, slv1, nplm, ipiv, info)
      call zgetrs('n', nplm, lmsize, slv1, nplm, ipiv, yrll1, nplm, info)
      call zgetrs('n', nplm, lmsize, slv1, nplm, ipiv, zrll1, nplm, info)

      do icheb2 = 0, ncheb

        do lm2 = 1, lmsize
          do lm1 = 1, lmsize
            yrll1temp(lm1, lm2) = yrll1(icheb2, lm1, lm2)
            zrll1temp(lm1, lm2) = zrll1(icheb2, lm1, lm2)
          end do
        end do

        call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, vhlr(1,1,icheb2), lmsize, yrll1temp, lmsize, czero, vhlr_yrll1, lmsize)
        call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, vhlr(1,1,icheb2), lmsize, zrll1temp, lmsize, czero, vhlr_zrll1, lmsize)
        call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, vjlr(1,1,icheb2), lmsize, yrll1temp, lmsize, czero, vjlr_yrll1, lmsize)
        call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, vjlr(1,1,icheb2), lmsize, zrll1temp, lmsize, czero, vjlr_zrll1, lmsize)

        do icheb = 0, ncheb
          taucslcr = -tau(icheb)*cslc1(icheb, icheb2)*drpan2
          mn = ipan*ncheb + ipan - icheb
          do lm2 = 1, lmsize
            do lm3 = 1, lmsize
              lm1 = lm3 + lmsize
              l1 = jlk_index(lm1)

              yrll2(icheb, lm3, lm2) = yrll2(icheb, lm3, lm2) + taucslcr*(jlk(l1,mn)*vhlr_yrll1(lm3,lm2)-hlk(l1,mn)*vjlr_yrll1(lm3,lm2))

              zrll2(icheb, lm3, lm2) = zrll2(icheb, lm3, lm2) + taucslcr*(jlk(l1,mn)*vhlr_zrll1(lm3,lm2)-hlk(l1,mn)*vjlr_zrll1(lm3,lm2))
            end do
          end do
        end do                     ! icheb

      end do                       ! icheb2

    else
      stop '[rll-loc] error in inversion'
    end if

    ! GPU: Reorient indices for later use on GPU
    if (use_sratrick==0) then
      !$acc parallel loop collapse(3)
      do icheb = 0, ncheb
        do lm2 = 1, lmsize
          do lm1 = 1, lmsize2
            yrf(lm1, lm2, icheb) = yrll(icheb, lm1, lm2)
            zrf(lm1, lm2, icheb) = zrll(icheb, lm1, lm2)
          end do
        end do
      end do

    else if (use_sratrick==1) then
      !$acc parallel loop collapse(3)
      do icheb = 0, ncheb
        do lm2 = 1, lmsize
          do lm1 = 1, lmsize
            yrf(lm1, lm2, icheb) = yrll1(icheb, lm1, lm2)
            zrf(lm1, lm2, icheb) = zrll1(icheb, lm1, lm2)
            yrf(lm1+lmsize, lm2, icheb) = yrll2(icheb, lm1, lm2)
            zrf(lm1+lmsize, lm2, icheb) = zrll2(icheb, lm1, lm2)
          end do
        end do
      end do

    end if

    ! GPU: Setup zslc1sum and perform local matrix multiplication reductions on GPU
    !$acc parallel loop
    do icheb = 0, ncheb
      zslc1sum(icheb) = slc1sum(icheb)*drpan2
    end do

    !$acc parallel loop collapse(2) private(tmp_nvy, tmp_jvy, tmp_nvz, tmp_jvz, lm3, icheb)
    do lm2 = 1, lmsize
      do lm1 = 1, lmsize
        tmp_nvy = czero
        tmp_jvy = czero
        tmp_nvz = czero
        tmp_jvz = czero
        do icheb = 0, ncheb
          do lm3 = 1, lmsize2
            tmp_nvy = tmp_nvy + zslc1sum(icheb) * vhlr(lm1, lm3, icheb) * yrf(lm3, lm2, icheb)
            tmp_jvy = tmp_jvy + zslc1sum(icheb) * vjlr(lm1, lm3, icheb) * yrf(lm3, lm2, icheb)
            tmp_nvz = tmp_nvz + zslc1sum(icheb) * vhlr(lm1, lm3, icheb) * zrf(lm3, lm2, icheb)
            tmp_jvz = tmp_jvz + zslc1sum(icheb) * vjlr(lm1, lm3, icheb) * zrf(lm3, lm2, icheb)
          end do
        end do
        mrnvy(lm1, lm2) = tmp_nvy
        mrjvy(lm1, lm2) = tmp_jvy
        mrnvz(lm1, lm2) = tmp_nvz
        mrjvz(lm1, lm2) = tmp_jvz
      end do
    end do

    !$acc end data
  end subroutine rll_local_solutions


end module mod_rll_local_solutions
