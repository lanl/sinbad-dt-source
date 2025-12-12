! © 2025. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos Nationa!l Laboratory (LANL), which is operated by Triad National
! Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC,
! and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce, prepare.
! derivative works, distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.

module srcdx_mod
    ! Use subroutines/funcitons
    use mcnp_interfaces_mod, only: expirx
    use source_mod, only: xinterp2d

    ! Use Parameters
    use mcnp_params, only: dknd, zero, half, one, two
    use source_mod, only: npts1, npts2, ma, mbv, m1, m2v, ed, emin

    ! Use Varibales
    use pblcom, only: pbl
    use tskcom, only: var_reduce_thread
    use source_mod, only: coslab, sinlab, pe1, dedx, ean, uione, vione, wione

    ! MCNP_DEBUG variables defined on idum and rdum cards
    use mcnp_debug, only: idum, rdum

    implicit none
    private

    public :: srcdx

contains

! *********************************************************
    subroutine srcdx

        integer :: i, iflag
        real(dknd) :: eb, cth, a, smd, sth, mb, m2, pana, t1n, t2n, t3n, &
                      x, x1, x2, x3, x4, x5, x6, x7, x8, x9
        real(dknd), dimension(npts1) :: pl

        ! Save variabels; they are computed only once.
        integer, save :: isrcdx_flag = 0
        real(dknd), save :: sml
        real(dknd), save :: thlab(npts2, npts1)
!$omp threadprivate(isrcdx_flag, sml, thlab)

        integer, parameter :: nterm = 3

        intrinsic :: acos, abs, sqrt

        iflag = idum(1)          ! =1 for DT reaction, =2 for DD reaction
        eb = rdum(1)*1e6_dknd    ! deuteron beam energy in eV

        if (isrcdx_flag == 0) then
            isrcdx_flag = 1
            call compute_sml_thlab(nterm, eb, sml, thlab) ! compute sml, thlab
        end if

        cth = uione*pbl%r%u + vione*pbl%r%v + wione*pbl%r%w
        a = acos(cth)
        do i = 1, npts1  ! energy bins
            call xinterp2d(ed, thlab(:, i), pe1, npts1, npts2, ed(i), a, pl(i))
            if (dedx(i) /= zero) then
                pl(i) = pl(i)/dedx(i)
            end if
        end do
        call xint1d(ed, pl, npts1, nterm, emin, eb, smd)
        var_reduce_thread%psc = smd*sml

        ! Energy calculation for psuedo-particle travelling directly to
        ! next event estimators
        mb = mbv(iflag)
        m2 = m2v(iflag)
        x1 = ma*ma
        x2 = cth*cth
        x3 = mb*mb
        x4 = m1*m1
        x5 = m2*m2
        x = ean*ean - x1
        if (x > zero) then
            pana = sqrt(x)
            x6 = x
        else
            pana = zero
            x6 = zero
        end if
        x = one - x2
        if (x > zero) then
            sth = sqrt(x)
            x7 = x
        else
            sth = zero
            x7 = zero
        end if
        x8 = mb*ean
        x9 = ean + mb

        x = x9*x9 - x6*x7
        if (x /= zero) then
            t1n = one/x
        else
            t1n = zero
        end if
        t2n = x9*(x8 + half*(x1 + x3 + x4 - x5))
        x = (x8 + half*(x1 + x3 - x4 - x5))**2 - x4*(x5 - x6*x7)
        if (x > zero) then
            t3n = pana*cth*sqrt(x)
        else
            t3n = zero
        end if

        x = t1n*(t2n + t3n) - m1
        if (x > zero) then
            pbl%r%erg = x
        end if

    end subroutine srcdx

! *********************************************************
! Private Subroutines/Functions
! *********************************************************
    subroutine compute_sml_thlab(nterm, eb, sml, thlab)
        ! Compute sml and thalab

        integer, intent(in)     :: nterm
        real(dknd), intent(in)  :: eb
        real(dknd), intent(out) :: sml
        real(dknd), intent(out) :: thlab(npts2, npts1)

        integer :: i, j
        real(dknd) :: pe2(npts2, npts1)
        real(dknd) :: thlab1(npts2), pl(npts2)
        real(dknd) :: pyi(npts1)

        intrinsic :: acos

        do j = 1, npts1 ! energy bins
            do i = 1, npts2  ! angle bins
                thlab(i, j) = acos(coslab(j, i))
                if (dedx(j) /= zero) then
                    pe2(i, j) = pe1(j, i)/dedx(j)*sinlab(j, i)
                else
                    pe2(i, j) = zero
                end if
            end do
        end do

        do j = 1, npts1  ! energy bins
            thlab1(:) = thlab(:, j)
            pl(:) = pe2(:, j)
            call xint1d(thlab1, pl, npts2, nterm, thlab1(1), thlab1(npts2), pyi(j))
        end do
        call xint1d(ed, pyi, npts1, nterm, emin, eb, sml)
        if (sml /= zero) then
            sml = one/sml
        end if

    end subroutine compute_sml_thlab

! *********************************************************
    subroutine xint1d(e, f, npts, nterms, e1, e2, tot)
        ! 1D numberical integration

        integer, intent(in)     :: npts, nterms
        real(dknd), intent(in)  :: e1, e2
        real(dknd), intent(in) :: e(npts), f(npts)
        real(dknd), intent(out) :: tot

        integer :: i, il, ih, k, ila, iha
        real(dknd) :: x(npts), y(npts)

        tot = zero

        do i = 1, npts
            if (e1 < e(i)) exit
        end do
        il = i

        do i = 1, npts
            if (e2 < e(i)) exit
        end do
        ih = i - 1

        if (ih >= il) then
            ila = il - nterms/2
            if (ila < 0) ila = 1
            if (ila + nterms - 1 > npts) ila = npts - nterms + 1
            iha = ih - nterms/2
            if (iha < 1) iha = 1
            if (iha + nterms - 1 > npts) iha = npts - nterms + 1
            call iteg(e, f, nterms, ila, e1, e(il), tot)
            call iteg(e, f, nterms, iha, e(ih), e2, tot)

            k = 0
            do i = il, ih
                k = k + 1
                x(k) = e(i)
                y(k) = f(i)
            end do
            if (k > 0) then
                tot = tot + area(x, y, k, nterms)
            end if

        else
            ila = il - nterms/2
            if (ila <= 0) ila = 1
            if (ila + nterms - 1 > npts) ila = npts - nterms + 1
            call iteg(e, f, nterms, ila, e1, e2, tot)
        end if

    end subroutine xint1d

! *********************************************************
    function area(x, y, npts, nterms)
        ! Compute area

        integer, intent(in) :: npts, nterms
        real(dknd), intent(in) :: x(:), y(:)

        integer :: i, j, i1, idelta, imax, neven
        real(dknd) :: area, tot, x1, x2

        tot = zero
        if ((npts - nterms) <= 0) then
            x1 = x(1)
            x2 = x(npts)
            call iteg(x, y, npts, 1, x1, x2, tot)
            area = tot
            return
        end if

        neven = 2*(nterms/2)
        idelta = nterms/2 - 1

        j = nterms - idelta
        x1 = x(1)
        x2 = x(j)
        call iteg(x, y, nterms, 1, x1, x2, tot)

        i1 = npts - nterms + 1
        j = i1 + idelta
        imax = i1 - 1

        if ((nterms - neven) <= 0) then
            x1 = x(j)
            x2 = x(npts)
            call iteg(x, y, nterms, i1, x1, x2, tot)
            do i = 2, imax
                j = i + idelta
                x1 = x(j)
                x2 = x(j + 1)
                call iteg(x, y, nterms, i, x1, x2, tot)
            end do

        else
            x1 = (x(j) + x(j + 1))/two
            x2 = x(npts)
            call iteg(x, y, nterms, i1, x1, x2, tot)
            do i = 2, imax
                j = i + idelta
                x1 = (x(j + 1) + x(j))/two
                x2 = (x(j + 2) + x(j + 1))/two
                call iteg(x, y, nterms, i, x1, x2, tot)
            end do
        end if

        area = tot

    end function area

! *********************************************************
    subroutine iteg(x, y, nterms, i1, x1, x2, tot)
        ! Integration helper

        integer, intent(in) :: nterms, i1
        real(dknd), intent(in) :: x1, x2
        real(dknd), intent(in) :: x(:), y(:)
        real(dknd), intent(inout) :: tot

        integer :: j, i, k, imid
        real(dknd) :: det, deltax, dx1, dx2, a
        real(dknd) :: array(nterms, nterms)

        array(1:nterms, 1) = one
        do k = 2, nterms
            array(1, k) = one
            do j = 2, nterms
                i = i1 + j - 1
                deltax = x(i) - x(i1)
                array(j, k) = deltax**(k - 1)
            end do
        end do
        call dmatinv(array, nterms, det)

        if (det == zero) then
            imid = i1 + nterms/2
            tot = tot + y(imid)*(x2 - x1)
        else
            dx1 = x1 - x(i1)
            dx2 = x2 - x(i1)
            do k = 1, nterms
                i = k + i1 - 1
                do j = 1, nterms
                    a = y(i)*array(j, k)
                    tot = tot + (a/j)*(dx2**j - dx1**j)
                end do
            end do
        end if

    end subroutine iteg

! *********************************************************
    subroutine dmatinv(matrix, norder, det)
        ! Compute matrix inverse

        integer, intent(in) :: norder
        real(dknd), intent(in)  :: matrix(norder, norder)
        real(dknd), intent(out) :: det

        integer :: i, j, k, ik, jk
        real(dknd) :: amax, tmp
        real(dknd) :: array(norder, norder)

        intrinsic :: abs

        det = one
        array(:, :) = matrix(:, :)
        do k = 1, norder
            amax = zero
            ! Find the largest array(i,j)
            do j = k, norder
                do i = k, norder
                    if (abs(amax) <= abs(array(i, j))) then
                        amax = array(i, j)
                        ik = i
                        jk = j
                    end if
                end do
            end do
            if (amax == zero) then
                det = zero
                return
            end if

            ! At this point amax /= zero, ik >= k and jk >= k, array(ik, jk) = amax

            if (ik > k) then
                ! Set array(k, :) to array(ik, ;) and array(ik, :) to negative of old array(k, :)
                do j = 1, norder
                    tmp = array(k, j)
                    array(k, j) = array(ik, j)
                    array(ik, j) = -tmp
                end do
            end if

            if (jk > k) then
                ! Set array(:, k) to array(:, jk) and array(:, jk) to negative of old array(:, k)
                do i = 1, norder
                    tmp = array(i, k)
                    array(i, k) = array(i, jk)
                    array(i, jk) = -tmp
                end do
            end if

            ! Set array(i, k) = -array(i, k)/amax for each i where i /= k
            do i = 1, norder
                if (i == k) cycle
                array(i, k) = -array(i, k)/amax
            end do

            ! Compute array(i, j) for each i and j where i /= k or j /= k
            do j = 1, norder
                if (j == k) cycle
                do i = 1, norder
                    if (i == k) cycle
                    array(i, j) = array(i, j) + array(i, k)*array(k, j)
                end do
            end do

            ! Set array(k, j) = array(k, j)/amax for each j where j /= k
            do j = 1, norder
                if (j == k) cycle
                array(k, j) = array(k, j)/amax
            end do

            array(k, k) = one/amax
            det = det*amax
        end do

    end subroutine dmatinv

end module srcdx_mod
