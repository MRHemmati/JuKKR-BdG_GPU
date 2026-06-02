!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dlke0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver for lattice fourier transform
  !> Author: 
  !> Category: KKRhost, k-points, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Function, set up in the spin-independent non-relativstic
  !> (l,m_l)-representation
  !-------------------------------------------------------------------------------
  subroutine dlke0(gllke, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, bzkp, rcls, ginp)

    use :: mod_runoptions, only: calc_complex_bandstructure, symmetrize_gmat, use_deci_onebulk
    use :: global_variables, only: lmgf0d, almgf0, naclsd, nrd
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: ci
    implicit none

    real (kind=dp) :: alat
    integer :: naez, naclsmax
    complex (kind=dp) :: ginp(lmgf0d*naclsmax, lmgf0d, *), gllke(almgf0, *)
    real (kind=dp) :: bzkp(*), rcls(3, naclsd, *), rr(3, 0:nrd)
    integer :: atom(naclsd, *), cls(*), ezoa(naclsd, *), nacls(*)
    ! ..
    integer :: i, ic, m, lm1, lm2, im, am, jn, j
    real (kind=dp) :: tpi, convpu
    complex (kind=dp) :: tt, eikr
    complex (kind=dp) :: arg(3)
    complex (kind=dp) :: kp(6)
    complex (kind=dp), dimension(:,:), allocatable :: gllke1

    tpi = 8.0_dp*atan(1.0_dp)
    convpu = alat/tpi

    gllke(1:almgf0, 1:almgf0) = (0.0_dp, 0.0_dp)

    kp(1) = bzkp(1)
    kp(2) = bzkp(2)
    kp(3) = bzkp(3)
    if (calc_complex_bandstructure) then
      kp(4) = bzkp(4)
      kp(5) = bzkp(5)
      kp(6) = bzkp(6)
    else
      kp(4) = (0.0_dp, 0.0_dp)
      kp(5) = (0.0_dp, 0.0_dp)
      kp(6) = (0.0_dp, 0.0_dp)
    end if

    !$acc parallel loop collapse(3) copyin(kp, rr, ezoa, atom, cls, nacls, rcls, ginp) copy(gllke) private(ic, arg, tt, eikr, im, am)
    do i = 1, naez
      do m = 1, nacls(cls(i))
        do lm2 = 1, lmgf0d
          if (atom(m, i) >= 0) then
            ic = cls(i)
            if (use_deci_onebulk) then
              arg(1) = -ci*tpi*rcls(1, m, ic)
              arg(2) = -ci*tpi*rcls(2, m, ic)
              arg(3) = -ci*tpi*rcls(3, m, ic)
            else
              arg(1) = -ci*tpi*rr(1, ezoa(m, i))
              arg(2) = -ci*tpi*rr(2, ezoa(m, i))
              arg(3) = -ci*tpi*rr(3, ezoa(m, i))
            end if

            tt = kp(1)*arg(1) + kp(2)*arg(2) + kp(3)*arg(3)
            if (calc_complex_bandstructure) then
              tt = tt + ci*(kp(4)*arg(1)+kp(5)*arg(2)+kp(6)*arg(3))
            end if

            eikr = exp(tt)*convpu

            im = 1 + (m-1)*lmgf0d
            am = 1 + (atom(m, i)-1)*lmgf0d
            do lm1 = 1, lmgf0d
              gllke(am + lm1 - 1, (i-1)*lmgf0d + lm2) = gllke(am + lm1 - 1, (i-1)*lmgf0d + lm2) + eikr * ginp(im + lm1 - 1, lm2, ic)
            end do
          end if
        end do
      end do
    end do

    if (symmetrize_gmat) then
      allocate(gllke1(almgf0, almgf0))
      gllke1(:, :) = (0.0_dp, 0.0_dp)

      kp(1) = -bzkp(1)
      kp(2) = -bzkp(2)
      kp(3) = -bzkp(3)
      if (calc_complex_bandstructure) then
        kp(4) = -bzkp(4)
        kp(5) = -bzkp(5)
        kp(6) = -bzkp(6)
      end if

      !$acc parallel loop collapse(3) copyin(kp, rr, ezoa, atom, cls, nacls, rcls, ginp) copy(gllke1) private(ic, arg, tt, eikr, im, am)
      do i = 1, naez
        do m = 1, nacls(cls(i))
          do lm2 = 1, lmgf0d
            if (atom(m, i) >= 0) then
              ic = cls(i)
              if (use_deci_onebulk) then
                arg(1) = -ci*tpi*rcls(1, m, ic)
                arg(2) = -ci*tpi*rcls(2, m, ic)
                arg(3) = -ci*tpi*rcls(3, m, ic)
              else
                arg(1) = -ci*tpi*rr(1, ezoa(m, i))
                arg(2) = -ci*tpi*rr(2, ezoa(m, i))
                arg(3) = -ci*tpi*rr(3, ezoa(m, i))
              end if

              tt = kp(1)*arg(1) + kp(2)*arg(2) + kp(3)*arg(3)
              if (calc_complex_bandstructure) then
                tt = tt + ci*(kp(4)*arg(1)+kp(5)*arg(2)+kp(6)*arg(3))
              end if

              eikr = exp(tt)*convpu

              im = 1 + (m-1)*lmgf0d
              am = 1 + (atom(m, i)-1)*lmgf0d
              do lm1 = 1, lmgf0d
                gllke1(am + lm1 - 1, (i-1)*lmgf0d + lm2) = gllke1(am + lm1 - 1, (i-1)*lmgf0d + lm2) + eikr * ginp(im + lm1 - 1, lm2, ic)
              end do
            end if
          end do
        end do
      end do

      !$acc parallel loop collapse(2) copy(gllke) copyin(gllke1)
      do j = 1, naez
        do i = 1, naez
          do lm2 = 1, lmgf0d
            do lm1 = 1, lmgf0d
              im = (i-1)*lmgf0d + lm1
              jn = (j-1)*lmgf0d + lm2
              gllke(im, jn) = (gllke(im, jn) + gllke1(jn, im)) * 0.5_dp
            end do
          end do
        end do
      end do

      deallocate(gllke1)
    end if

    return

  end subroutine dlke0

end module mod_dlke0
