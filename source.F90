! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos Nationa!l Laboratory (LANL), which is operated by Triad National
! Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC,
! and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce, prepare.
! derivative works, distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.

module source_mod
    ! Use subroutines/functions
    use mcnp_interfaces_mod, only: expirx

    ! Use Parameters
    use mcnp_params, only: dknd, very_small, hundredth, half, pie, &
                           zero, one, two, four, five, ten

    ! Use Variables
    use mcnp_global, only: ncl
    use mcnp_random, only: rng_history
    use fixcom, only: mxa
    use pblcom, only: pbl

    ! MCNP_DEBUG variables defined on idum and rdum cards
    use mcnp_debug, only: idum, rdum

    implicit none
    private

    public :: source, xinterp2d

    integer, public, parameter :: npts1 = 38   ! number of energy bins
    integer, public, parameter :: npts2 = 37   ! number of angle bins

    ! Masses from NIST Reference on Constants, Units and Uncertainties Web site:
    ! https://physics.nist.gov/cuu/Constants/index.html
    !    ma = Deuteron mass, MeV/c^2
    !    mbv(1) = Triton mass, MeV/c^2; mbv(2) = Deuteron mass, MeV/c^2
    !    m1 = Neutron mass, MeV/c^2
    !    m2v(1) = Alpha mass, MeV/c^2, m2v(2) = Helion (He-3) mass, MeV/c^2
    real(dknd), public, parameter :: ma = 1875.612793_dknd
    real(dknd), public, parameter :: mbv(1:2) = [2808.920906_dknd, 1875.612793_dknd]
    real(dknd), public, parameter :: m1 = 939.565346_dknd
    real(dknd), public, parameter :: m2v(1:2) = [3727.379109_dknd, 2808.391383_dknd]
    ! Min energy, MeV
    real(dknd), public, parameter :: emin = hundredth  ! 10 keV
    ! Max energy, MeV
    real(dknd), public, parameter :: emax = ten        ! 10 MeV
    ! Energy bins, MeV
    real(dknd), public, parameter, dimension(npts1) :: &
        ed = [0.01_dknd, 0.02_dknd, 0.03_dknd, 0.04_dknd, 0.05_dknd, &
              0.06_dknd, 0.07_dknd, 0.08_dknd, 0.09_dknd, 0.10_dknd, &
              0.108_dknd, 0.12_dknd, 0.14_dknd, 0.16_dknd, 0.18_dknd, &
              0.20_dknd, 0.22_dknd, 0.24_dknd, 0.26_dknd, 0.28_dknd, &
              0.30_dknd, 0.32_dknd, 0.34_dknd, 0.36_dknd, 0.38_dknd, &
              0.40_dknd, 0.50_dknd, 0.75_dknd, 1.0_dknd, 2.0_dknd, &
              3.0_dknd, 4.0_dknd, 5.0_dknd, 6.0_dknd, 7.0_dknd, &
              8.0_dknd, 9.0_dknd, 10.0_dknd]

! ----------------------------------------------------------------
! The following variables must be threadprivate since these
! variables are defined and changed in the OMP parallel region.
! They are computed only once when the source subroutine is called.
! -----------------------------------------------------------------

    ! Scalars that are used in srcdx
    real(dknd), public :: ean, uione, vione, wione
!$omp threadprivate(ean, uione, vione, wione)

    ! Arrays that are used in srcdx
    real(dknd), public, dimension(npts1) :: dedx
    real(dknd), public, dimension(npts1, npts2) :: coslab, sinlab, pe1
!$omp threadprivate(dedx, coslab, sinlab, pe1)

! -----------
! Parameters
! -----------
    integer, parameter :: mit = 50  ! max number of calling rot before reseting

contains

! *********************************************************
    subroutine source

!      The original SOURCE subroutine in SINBAD was created by A. Milocco and
!      is dicussed in the following reference:
!      [1] Alberto Milocco. MCNPX/MCNP5 Routine for Simulating D-T Neutron Source
!      in Ti-T Targets. Tech. rep. IJS-DP-9988. Institute Josef Stefan, Reactor Physics
!      Department, Ljubljana, 2008.
!
!      The method for computing the low energy deuteron scattering was created by
!      Biersack and Haggmark, and is described in the following reference:
!      [2] J.P. Biersack and L.G. Haggmark. “A Monte Carlo computer program for the
!      transport of energetic ions in amorphous targets”. In: Nuclear Instruments and
!      Methods 174.1 (1980), pp. 257–269. issn: 0029-554X.
!      doi: https://doi.org/10.1016/0029-554X(80)90440-1.
!      url: https://www.sciencedirect.com/science/article/pii/0029554X80904401.
!
!      This subroutine is for D-T/D-D source anisotropy interaction includes:
!      * ion scattering
!      * 37 angles
!      * 38 energy points
!      * minimum beam energy  10 keV
!      * maximun energy 10 MeV (NOT BREAK-UP)
!
!      The source routine is called to set following source particle (pbl) values:
!      pbl%r%x, pbl%r%y, pbl%r%z, pbl%i%icl, pbl%i%jsu, pbl%i%erg, pbl%r%wgt,
!      pbl%r%tme, pbl%i%ipt, pbl%r%u, pbl%r%v, pbl%r%w.
!      The subroutine srcdx may also be needed.
!
!      INPUT PARAMETERS ON RDUM CARD:
!      * rdum(1)=deuteron beam energy, MeV
!      * rdum(2)=Tritium/Titanium atomic ratio
!      * rdum(3)= initial position in x coordinate, cm
!      * rdum(4)= position in y coordinate for neutron, cm
!      * rdum(5)= initial position in z coordinate, cm
!      * rdum(6)=deuteron beam width, cm
!
!      INPUT PARAMETERS ON IDUM CARD:
!      * idum(1)=flag, 1 for D-T reaction, 2 for D-D reaction
!      * idum(2)=source cell defined in an MCNP input file
!      * idum(3)=maximum iterations for rejection; will set to 1000 if idum(3) < 1000
!
!      comments:
!      * rdum(4) should be positive
!      * D-D not tested.

! ****************
! Local Variables
! ****************
        integer :: i, j, k, l, n, irtt, icycle, rej
        integer :: iflag, cell_number, cell_index, max_iter

        real(dknd) :: eb, xt, ry, mb, m2, ebb, rtitol, yylab, &
                      ma2, mb2, m12, m22, mambm, mambp, mambm2, mambp2, &
                      ztm, atm, rhod, a, c, thl, ea, pa, s, s2, s3, s4, s5, s6, &
                      gammacm, rlambda, rlambda1, betacm, beta1, rg1, theta1, &
                      p1t1, p1t2, p1t3, p1t, e1t1, e1t2, e1t3, e1t, &
                      eold, eps, eeg, pmax1, ffpath, ein, exout, dee, stdee, &
                      p, b, sqe, c2, r, rr, ex1, ex2, ex3, ex4, v, v1, fr, fr1, qq, &
                      roc, cc, aa, ff, delta, ct, st, cu, psi, deny, uzero, &
                      ep, t2, a1, pn, ctn, stn, ean, pana, ph, sinlab1, &
                      t1n, t2n, t3n, e1n, ry2, x1, z1

        real(dknd), dimension(2)  :: stoich
        real(dknd), dimension(3) :: am, z
        real(dknd), dimension(11) :: pl

        real(dknd), allocatable :: coscm(:), sincm(:)
        real(dknd), allocatable :: xscman(:, :), xslaban(:, :)

!*****************
! Saved Variables
!*****************
        ! The following variables are computed only once and fixed
        integer, save :: isource_flag = 0
        real(dknd), save :: epeak, atmrho, ffpa, ffpf, epsdg, stbohr, pemax
!$omp threadprivate(isource_flag, epeak, atmrho, ffpa, ffpf, epsdg, stbohr, pemax)

        ! The following variables will be varied.
        real(dknd), save :: emev, e, distalb, uione1, vione1, wione1
!$omp threadprivate(emev, e, distalb, uione1, vione1, wione1)

        ! The following variables are computed only once and fixed.
        real(dknd), save, dimension(2:3) :: m2m1, ec, a2, f
        real(dknd), save, dimension(npts2) :: thcm
        real(dknd), save, dimension(npts1, npts2) :: pe
!$omp threadprivate(m2m1, ec, a2, f, thcm, pe)

! ************
! Parameters
! ************
        real(dknd), parameter :: fuzz = 1.e-7_dknd
        real(dknd), parameter :: calb = 1000._dknd/(two*pie)

        ! Legendre coefficients from ENDF/B-VIII.0 (MF=6, MT=50)
        real(dknd), dimension(10, npts1, 2), parameter :: &
            coeff = reshape([ &
                            ! D-T
                            7.327700e-5_dknd, 1.608200e-4_dknd, -7.473100e-7_dknd, -4.203800e-8_dknd, 5.11080e-10_dknd, &
                            1.45780e-11_dknd, 1.15120e-13_dknd, -7.95610e-15_dknd, 1.08850e-14_dknd, -8.57510e-15_dknd, &
                            2.984100e-4_dknd, 1.653600e-4_dknd, -8.807900e-7_dknd, -5.931300e-8_dknd, 1.057200e-9_dknd, &
                            3.96900e-11_dknd, 3.45020e-13_dknd, 0.000000e+0_dknd, 8.45920e-15_dknd, -6.54030e-15_dknd, &
                            4.680500e-4_dknd, 1.519000e-4_dknd, -8.589000e-7_dknd, -6.668500e-8_dknd, 2.001500e-9_dknd, &
                            7.50340e-11_dknd, 7.03140e-13_dknd, 1.47510e-14_dknd, 4.41030e-15_dknd, -3.96070e-15_dknd, &
                            6.115700e-4_dknd, 1.258900e-4_dknd, -7.121700e-7_dknd, -6.472400e-8_dknd, 3.359900e-9_dknd, &
                            1.18000e-10_dknd, 1.13570e-12_dknd, 3.80350e-14_dknd, 6.91890e-16_dknd, -4.17330e-16_dknd, &
                            7.418600e-4_dknd, 8.968800e-5_dknd, -4.540100e-7_dknd, -5.319400e-8_dknd, 5.159600e-9_dknd, &
                            1.67550e-10_dknd, 1.60020e-12_dknd, 7.21140e-14_dknd, -3.53700e-15_dknd, 2.98190e-15_dknd, &
                            8.659800e-4_dknd, 4.467600e-5_dknd, -9.224000e-8_dknd, -3.139600e-8_dknd, 7.478300e-9_dknd, &
                            2.26230e-10_dknd, 2.14090e-12_dknd, 1.20150e-13_dknd, -8.93660e-15_dknd, 7.51230e-15_dknd, &
                            9.881700e-4_dknd, -8.173400e-6_dknd, 3.681600e-7_dknd, 1.734200e-9_dknd, 1.046000e-8_dknd, &
                            3.01800e-10_dknd, 2.94820e-12_dknd, 1.86840e-13_dknd, -1.40100e-14_dknd, 1.22060e-14_dknd, &
                            1.111100e-3_dknd, -6.811900e-5_dknd, 9.240100e-7_dknd, 4.757600e-8_dknd, 1.432400e-8_dknd, &
                            4.08440e-10_dknd, 4.44840e-12_dknd, 2.86040e-13_dknd, -2.01520e-14_dknd, 1.68560e-14_dknd, &
                            1.236500e-3_dknd, -1.345600e-4_dknd, 1.573500e-6_dknd, 1.077800e-7_dknd, 1.936600e-8_dknd, &
                            5.67800e-10_dknd, 7.35900e-12_dknd, 4.34100e-13_dknd, -2.47950e-14_dknd, 2.32730e-14_dknd, &
                            1.365600e-3_dknd, -2.070000e-4_dknd, 2.315700e-6_dknd, 1.842600e-7_dknd, 2.596100e-8_dknd, &
                            8.09830e-10_dknd, 1.27710e-11_dknd, 6.63720e-13_dknd, -3.01930e-14_dknd, 2.94700e-14_dknd, &
                            1.472000e-3_dknd, -2.689600e-4_dknd, 2.976200e-6_dknd, 2.585700e-7_dknd, 3.265600e-8_dknd, &
                            1.088800e-9_dknd, 1.99090e-11_dknd, 9.34890e-13_dknd, -3.38570e-14_dknd, 3.41840e-14_dknd, &
                            1.637300e-3_dknd, -3.681400e-4_dknd, 4.078800e-6_dknd, 3.946700e-7_dknd, 4.569400e-8_dknd, &
                            1.707100e-9_dknd, 3.76940e-11_dknd, 1.55960e-12_dknd, -3.67380e-14_dknd, 4.21480e-14_dknd, &
                            1.929600e-3_dknd, -5.486300e-4_dknd, 6.219800e-6_dknd, 6.979200e-7_dknd, 7.805200e-8_dknd, &
                            3.528500e-9_dknd, 9.77880e-11_dknd, 3.55240e-12_dknd, -3.23910e-14_dknd, 5.85190e-14_dknd, &
                            2.243500e-3_dknd, -7.460400e-4_dknd, 8.751700e-6_dknd, 1.115600e-6_dknd, 1.287300e-7_dknd, &
                            6.868400e-9_dknd, 2.21540e-10_dknd, 7.61890e-12_dknd, -3.62360e-16_dknd, 7.62440e-14_dknd, &
                            2.578800e-3_dknd, -9.583600e-4_dknd, 1.169200e-5_dknd, 1.671400e-6_dknd, 2.045400e-7_dknd, &
                            1.249700e-8_dknd, 4.49160e-10_dknd, 1.52360e-11_dknd, 8.08060e-14_dknd, 9.97730e-14_dknd, &
                            2.934900e-3_dknd, -1.183800e-3_dknd, 1.506100e-5_dknd, 2.390400e-6_dknd, 3.133300e-7_dknd, &
                            2.136800e-8_dknd, 8.34240e-10_dknd, 2.85190e-11_dknd, 2.46860e-13_dknd, 1.30540e-13_dknd, &
                            3.310900e-3_dknd, -1.420900e-3_dknd, 1.888200e-5_dknd, 3.299200e-6_dknd, 4.639200e-7_dknd, &
                            3.462300e-8_dknd, 1.445300e-9_dknd, 5.03420e-11_dknd, 5.49000e-13_dknd, 1.71650e-13_dknd, &
                            3.705600e-3_dknd, -1.668200e-3_dknd, 2.317800e-5_dknd, 4.425500e-6_dknd, 6.660400e-7_dknd, &
                            5.359500e-8_dknd, 2.367100e-9_dknd, 8.44660e-11_dknd, 1.05620e-12_dknd, 2.29510e-13_dknd, &
                            4.117900e-3_dknd, -1.924400e-3_dknd, 2.797500e-5_dknd, 5.797600e-6_dknd, 9.302900e-7_dknd, &
                            7.981200e-8_dknd, 3.702200e-9_dknd, 1.35680e-10_dknd, 1.86310e-12_dknd, 3.12680e-13_dknd, &
                            4.546600e-3_dknd, -2.188400e-3_dknd, 3.329800e-5_dknd, 7.445000e-6_dknd, 1.268100e-6_dknd, &
                            1.150000e-7_dknd, 5.571900e-9_dknd, 2.09910e-10_dknd, 3.09040e-12_dknd, 4.30890e-13_dknd, &
                            4.990700e-3_dknd, -2.459100e-3_dknd, 3.917500e-5_dknd, 9.397300e-6_dknd, 1.691500e-6_dknd, &
                            1.610700e-7_dknd, 8.117700e-9_dknd, 3.14380e-10_dknd, 4.89030e-12_dknd, 6.00750e-13_dknd, &
                            5.449100e-3_dknd, -2.735600e-3_dknd, 4.563300e-5_dknd, 1.168500e-5_dknd, 2.213400e-6_dknd, &
                            2.201300e-7_dknd, 1.150200e-8_dknd, 4.57720e-10_dknd, 7.45520e-12_dknd, 8.39580e-13_dknd, &
                            5.920700e-3_dknd, -3.016900e-3_dknd, 5.270000e-5_dknd, 1.433800e-5_dknd, 2.847200e-6_dknd, &
                            2.945000e-7_dknd, 1.590900e-8_dknd, 6.50130e-10_dknd, 1.10170e-11_dknd, 1.17270e-12_dknd, &
                            6.404600e-3_dknd, -3.302100e-3_dknd, 6.040300e-5_dknd, 1.738900e-5_dknd, 3.607100e-6_dknd, &
                            3.866800e-7_dknd, 2.154700e-8_dknd, 9.03450e-10_dknd, 1.58550e-11_dknd, 1.62950e-12_dknd, &
                            6.899800e-3_dknd, -3.590400e-3_dknd, 6.877300e-5_dknd, 2.086800e-5_dknd, 4.507600e-6_dknd, &
                            4.993400e-7_dknd, 2.864600e-8_dknd, 1.231300e-9_dknd, 2.23040e-11_dknd, 2.24820e-12_dknd, &
                            7.405600e-3_dknd, -3.881100e-3_dknd, 7.783800e-5_dknd, 2.480600e-5_dknd, 5.563900e-6_dknd, &
                            6.353700e-7_dknd, 3.746100e-8_dknd, 1.649400e-9_dknd, 3.07550e-11_dknd, 3.07300e-12_dknd, &
                            1.006600e-2_dknd, -5.344500e-3_dknd, 1.346600e-4_dknd, 5.248800e-5_dknd, 1.374100e-5_dknd, &
                            1.779300e-6_dknd, 1.180300e-7_dknd, 5.788200e-9_dknd, 1.21570e-10_dknd, 1.25660e-11_dknd, &
                            1.728600e-2_dknd, -8.669600e-3_dknd, 3.825000e-4_dknd, 2.007700e-4_dknd, 6.702200e-5_dknd, &
                            1.074200e-5_dknd, 8.822900e-7_dknd, 5.303000e-8_dknd, 1.368200e-9_dknd, 1.70200e-10_dknd, &
                            2.474200e-2_dknd, -1.076300e-2_dknd, 8.355400e-4_dknd, 5.062900e-4_dknd, 1.954500e-4_dknd, &
                            3.630700e-5_dknd, 3.479400e-6_dknd, 2.421600e-7_dknd, 7.217900e-9_dknd, 1.053100e-9_dknd, &
                            5.173400e-2_dknd, 4.230700e-3_dknd, 5.734300e-3_dknd, 3.878700e-3_dknd, 2.020600e-3_dknd, &
                            5.367600e-4_dknd, 7.651200e-5_dknd, 7.532800e-6_dknd, 3.165600e-7_dknd, 7.121200e-8_dknd, &
                            7.356900e-2_dknd, 6.055900e-2_dknd, 1.518900e-2_dknd, 9.433400e-3_dknd, 6.518100e-3_dknd, &
                            2.105300e-3_dknd, 3.916400e-4_dknd, 4.665300e-5_dknd, 2.383600e-6_dknd, 7.229100e-7_dknd, &
                            9.546600e-2_dknd, 1.265700e-1_dknd, 2.491100e-2_dknd, 1.347000e-2_dknd, 1.373100e-2_dknd, &
                            4.868700e-3_dknd, 1.102600e-3_dknd, 1.549100e-4_dknd, 9.306900e-6_dknd, 3.533800e-6_dknd, &
                            1.101100e-1_dknd, 1.545900e-1_dknd, 3.149000e-2_dknd, 1.782300e-2_dknd, 2.421400e-2_dknd, &
                            8.843700e-3_dknd, 2.249500e-3_dknd, 3.911600e-4_dknd, 2.835700e-5_dknd, 1.254200e-5_dknd, &
                            1.196600e-1_dknd, 1.451700e-1_dknd, 3.155400e-2_dknd, 2.927100e-2_dknd, 3.851500e-2_dknd, &
                            1.451900e-2_dknd, 3.925400e-3_dknd, 8.666000e-4_dknd, 7.659600e-5_dknd, 3.814500e-5_dknd, &
                            1.167600e-1_dknd, 1.309000e-1_dknd, 2.144400e-2_dknd, 4.427300e-2_dknd, 5.337200e-2_dknd, &
                            2.192300e-2_dknd, 6.188700e-3_dknd, 1.724300e-3_dknd, 1.839500e-4_dknd, 1.031500e-4_dknd, &
                            1.015800e-1_dknd, 1.313500e-1_dknd, 8.771600e-3_dknd, 5.615700e-2_dknd, 6.396600e-2_dknd, &
                            3.087800e-2_dknd, 8.775300e-3_dknd, 3.079100e-3_dknd, 6.463300e-4_dknd, 3.051800e-4_dknd, &
                            9.051100e-2_dknd, 1.326300e-1_dknd, -1.724800e-3_dknd, 6.077300e-2_dknd, 7.008100e-2_dknd, &
                            4.026600e-2_dknd, 1.063400e-2_dknd, 4.861800e-3_dknd, 1.104600e-3_dknd, 6.656100e-4_dknd, &
                            8.831100e-2_dknd, 1.308000e-1_dknd, -1.163300e-2_dknd, 6.132700e-2_dknd, 7.390400e-2_dknd, &
                            5.016900e-2_dknd, 1.212100e-2_dknd, 6.803000e-3_dknd, 1.505000e-3_dknd, 1.061600e-3_dknd, &
                            ! D-D
                            0.0_dknd, 4.5426e-2_dknd, 0.0_dknd, 7.0971e-5_dknd, 0.0_dknd, &
                            3.5316e-9_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 5.2119e-2_dknd, 0.0_dknd, 1.2521e-4_dknd, 0.0_dknd, &
                            9.5527e-9_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 5.8303e-2_dknd, 0.0_dknd, 1.8575e-4_dknd, 0.0_dknd, &
                            1.9418e-8_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 6.4068e-2_dknd, 0.0_dknd, 2.5274e-4_dknd, 0.0_dknd, &
                            3.3838e-8_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 6.9460e-2_dknd, 0.0_dknd, 3.2584e-4_dknd, 0.0_dknd, &
                            5.3474e-8_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 7.4508e-2_dknd, 0.0_dknd, 4.0466e-4_dknd, 0.0_dknd, &
                            7.8944e-8_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 7.9237e-2_dknd, 0.0_dknd, 4.8880e-4_dknd, 0.0_dknd, &
                            1.1082e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 8.3670e-2_dknd, 0.0_dknd, 5.7785e-4_dknd, 0.0_dknd, &
                            1.4965e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 8.7829e-2_dknd, 0.0_dknd, 6.7147e-4_dknd, 0.0_dknd, &
                            1.9592e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 9.1731e-2_dknd, 0.0_dknd, 7.6933e-4_dknd, 0.0_dknd, &
                            2.5013e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 9.4337e-2_dknd, 0.0_dknd, 8.5637e-4_dknd, 0.0_dknd, &
                            3.1499e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 9.8247e-2_dknd, 0.0_dknd, 9.8692e-4_dknd, 0.0_dknd, &
                            4.1229e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.0476e-1_dknd, 0.0_dknd, 1.2045e-3_dknd, 0.0_dknd, &
                            5.7444e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.1046e-1_dknd, 0.0_dknd, 1.4371e-3_dknd, 0.0_dknd, &
                            7.8882e-7_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.1534e-1_dknd, 0.0_dknd, 1.6846e-3_dknd, 0.0_dknd, &
                            1.0554e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.2022e-1_dknd, 0.0_dknd, 1.9322e-3_dknd, 0.0_dknd, &
                            1.3220e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.2397e-1_dknd, 0.0_dknd, 2.2043e-3_dknd, 0.0_dknd, &
                            1.7122e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.2772e-1_dknd, 0.0_dknd, 2.4764e-3_dknd, 0.0_dknd, &
                            2.1025e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.3107e-1_dknd, 0.0_dknd, 2.7595e-3_dknd, 0.0_dknd, &
                            2.5642e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.3402e-1_dknd, 0.0_dknd, 3.0534e-3_dknd, 0.0_dknd, &
                            3.0973e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.3697e-1_dknd, 0.0_dknd, 3.3473e-3_dknd, 0.0_dknd, &
                            3.6305e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.3934e-1_dknd, 0.0_dknd, 3.6616e-3_dknd, 0.0_dknd, &
                            4.3265e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.4171e-1_dknd, 0.0_dknd, 3.9759e-3_dknd, 0.0_dknd, &
                            5.0224e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.4386e-1_dknd, 0.0_dknd, 4.3001e-3_dknd, 0.0_dknd, &
                            5.8102e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.4580e-1_dknd, 0.0_dknd, 4.6342e-3_dknd, 0.0_dknd, &
                            6.6899e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.4774e-1_dknd, 0.0_dknd, 4.9682e-3_dknd, 0.0_dknd, &
                            7.5696e-6_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.5513e-1_dknd, 0.0_dknd, 6.7831e-3_dknd, 0.0_dknd, &
                            1.3570e-5_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.6587e-1_dknd, 0.0_dknd, 1.2154e-2_dknd, 0.0_dknd, &
                            4.0835e-5_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.7109e-1_dknd, 0.0_dknd, 1.8687e-2_dknd, 0.0_dknd, &
                            9.2837e-5_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.7641e-1_dknd, 0.0_dknd, 5.4104e-2_dknd, 0.0_dknd, &
                            7.7734e-4_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.8281e-1_dknd, 0.0_dknd, 9.6315e-2_dknd, 0.0_dknd, &
                            2.8960e-3_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 1.9856e-1_dknd, 0.0_dknd, 1.3677e-1_dknd, 0.0_dknd, &
                            7.4280e-3_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.2182e-1_dknd, 0.0_dknd, 1.7078e-1_dknd, 0.0_dknd, &
                            1.5224e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.4748e-1_dknd, 0.0_dknd, 1.9501e-1_dknd, 0.0_dknd, &
                            2.6667e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.6829e-1_dknd, 0.0_dknd, 2.0589e-1_dknd, 0.0_dknd, &
                            4.1448e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.7499e-1_dknd, 0.0_dknd, 1.9917e-1_dknd, 0.0_dknd, &
                            5.8470e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.5869e-1_dknd, 0.0_dknd, 1.7305e-1_dknd, 0.0_dknd, &
                            7.5678e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, &
                            0.0_dknd, 2.1701e-1_dknd, 0.0_dknd, 1.3419e-1_dknd, 0.0_dknd, &
                            9.0279e-2_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd, 0.0_dknd], &
                            [10, npts1, 2])

        ! Stopping force from SRIM2006 (kev/mg/cm)
        real(dknd), dimension(npts1, 2), parameter :: &
            cdedx = reshape([ &
                            ! D-Ti
                            143.4_dknd, 195.8_dknd, 234.2_dknd, 264.5_dknd, 289.1_dknd, &
                            309.4_dknd, 326.2_dknd, 340.1_dknd, 351.7_dknd, 361.1_dknd, &
                            368.8_dknd, 375.0_dknd, 383.5_dknd, 387.9_dknd, 389.3_dknd, &
                            388.4_dknd, 385.5_dknd, 381.5_dknd, 376.6_dknd, 371.1_dknd, &
                            365.2_dknd, 358.7_dknd, 352.2_dknd, 345.7_dknd, 339.1_dknd, &
                            332.5_dknd, 302.0_dknd, 244.8_dknd, 207.6_dknd, 139.9_dknd, &
                            109.2_dknd, 91.0_dknd, 78.7_dknd, 69.7_dknd, 62.8_dknd, &
                            57.3_dknd, 52.7_dknd, 49.0_dknd, &
                            ! D-Ti also for Ti-D target
                            143.4_dknd, 195.8_dknd, 234.2_dknd, 264.5_dknd, 289.1_dknd, &
                            309.4_dknd, 326.2_dknd, 340.1_dknd, 351.7_dknd, 361.1_dknd, &
                            368.8_dknd, 375.0_dknd, 383.5_dknd, 387.9_dknd, 389.3_dknd, &
                            388.4_dknd, 385.5_dknd, 381.5_dknd, 376.6_dknd, 371.1_dknd, &
                            365.2_dknd, 358.7_dknd, 352.2_dknd, 345.7_dknd, 339.1_dknd, &
                            332.5_dknd, 302.0_dknd, 244.8_dknd, 207.6_dknd, 139.9_dknd, &
                            109.2_dknd, 91.0_dknd, 78.7_dknd, 69.7_dknd, 62.8_dknd, &
                            57.3_dknd, 52.7_dknd, 49.0_dknd], &
                            [npts1, 2])

        !   D-CU; These data were from the orginal SINBAD. To use, need to add `_dknd`
        !   &86.403,115.32,135.72,150.89, 162.8,172.28,180.11,186.58,192.07,196.67,
        !   &199.81,204.03,209.32,213.13,215.86, 217.6,218.59,219.05,219.09,218.78,
        !   &218.11,217.21,216.13,214.91,213.57,212.11,203.94, 182.9,164.69,118.91,
        !   &94.497, 79.76, 69.56,61.992,56.097,51.353, 47.45,44.167/

        real(dknd), dimension(npts1, 2), parameter :: &
            tdedx = reshape([ &
                            ! D-T
                            516.7_dknd, 707.4_dknd, 841.0_dknd, 940.1_dknd, 1014.2_dknd, &
                            1069.9_dknd, 1109.6_dknd, 1137.5_dknd, 1155.3_dknd, 1165.2_dknd, &
                            1167.1_dknd, 1165.0_dknd, 1147.9_dknd, 1118.8_dknd, 1084.7_dknd, &
                            1046.7_dknd, 1007.6_dknd, 968.8_dknd, 930.9_dknd, 894.3_dknd, &
                            859.3_dknd, 826.6_dknd, 795.7_dknd, 766.5_dknd, 739.1_dknd, &
                            713.4_dknd, 607.5_dknd, 447.8_dknd, 357.8_dknd, 212.9_dknd, &
                            170.1_dknd, 134.8_dknd, 111.6_dknd, 95.5_dknd, 83.7_dknd, &
                            74.7_dknd, 67.6_dknd, 61.8_dknd, &
                            ! D-D
                            779.1_dknd, 1062.0_dknd, 1261.2_dknd, 1409.7_dknd, 1520.7_dknd, &
                            1603.1_dknd, 1663.6_dknd, 1705.2_dknd, 1730.9_dknd, 1745.6_dknd, &
                            1748.7_dknd, 1745.3_dknd, 1719.0_dknd, 1676.8_dknd, 1624.6_dknd, &
                            1567.5_dknd, 1508.9_dknd, 1451.0_dknd, 1394.6_dknd, 1339.9_dknd, &
                            1287.0_dknd, 1238.2_dknd, 1191.7_dknd, 1148.1_dknd, 1107.2_dknd, &
                            1068.8_dknd, 909.86_dknd, 670.66_dknd, 535.86_dknd, 318.89_dknd, &
                            254.63_dknd, 201.9_dknd, 167.18_dknd, 143.07_dknd, 125.36_dknd, &
                            111.85_dknd, 101.15_dknd, 92.495_dknd], &
                            [npts1, 2])

        ! Integral XS from ENDF/B-VIII.0 (MF=3, MT=50)
        real(dknd), dimension(npts1, 2), parameter :: &
            ed2c = reshape([ &
                           ! D-T
                           1.732600e-3_dknd, 5.967800e-2_dknd, 2.809600e-1_dknd, 7.163600e-1_dknd, 1.370200e+0_dknd, &
                           2.202100e+0_dknd, 3.113400e+0_dknd, 3.957100e+0_dknd, 4.590200e+0_dknd, 4.937000e+0_dknd, &
                           5.014200e+0_dknd, 4.876900e+0_dknd, 4.302700e+0_dknd, 3.647500e+0_dknd, 3.072300e+0_dknd, &
                           2.605100e+0_dknd, 2.233200e+0_dknd, 1.936900e+0_dknd, 1.698500e+0_dknd, 1.504600e+0_dknd, &
                           1.344800e+0_dknd, 1.211600e+0_dknd, 1.099400e+0_dknd, 1.003900e+0_dknd, 9.218200e-1_dknd, &
                           8.507400e-1_dknd, 6.048500e-1_dknd, 3.377800e-1_dknd, 2.315800e-1_dknd, 1.132000e-1_dknd, &
                           8.999200e-2_dknd, 8.446100e-2_dknd, 8.058200e-2_dknd, 7.434800e-2_dknd, 6.797900e-2_dknd, &
                           6.123400e-2_dknd, 5.568000e-2_dknd, 5.137200e-2_dknd, &
                           ! D-D
                           8.8216e-6_dknd, 2.7734e-4_dknd, 1.1747e-3_dknd, 2.6785e-3_dknd, 4.6105e-3_dknd, &
                           6.8050e-3_dknd, 9.1418e-3_dknd, 1.1541e-2_dknd, 1.3950e-2_dknd, 1.6339e-2_dknd, &
                           1.8124e-2_dknd, 2.0801e-2_dknd, 2.5264e-2_dknd, 2.9419e-2_dknd, 3.3266e-2_dknd, &
                           3.7113e-2_dknd, 4.0411e-2_dknd, 4.3709e-2_dknd, 4.6779e-2_dknd, 4.9621e-2_dknd, &
                           5.2463e-2_dknd, 5.4925e-2_dknd, 5.7387e-2_dknd, 5.9688e-2_dknd, 6.1829e-2_dknd, &
                           6.3970e-2_dknd, 7.2709e-2_dknd, 8.6567e-2_dknd, 9.3566e-2_dknd, 1.0425e-1_dknd, &
                           1.0388e-1_dknd, 1.0046e-1_dknd, 9.6621e-2_dknd, 9.2854e-2_dknd, 8.9473e-2_dknd, &
                           8.5752e-2_dknd, 8.3193e-2_dknd, 8.2559e-2_dknd], &
                           [npts1, 2])

        intrinsic :: abs, acos, atan, cos, exp, log, max, sin, sqrt

        ! The following calcultions are performed when each source is called,
        ! as we do not want to add `save` to these variables.
        ! These calculations are not expensive.

        ! Set values from idum card
        iflag = idum(1)           ! 1 for DT reaction, 2 for DD
        cell_number = idum(2)     ! source cell
        max_iter = idum(3)        ! Maximum number of iterations for rejection

        ! Set values from rdum card
        eb = rdum(1)              ! Deuteron  beam energy, MeV
        xt = rdum(2)              ! T/Ti atom ratio
        ry = rdum(6)              ! Deuteron beam width, cm

        ! Check iflag value
        if (iflag == 1 .or. iflag == 2) then
            mb = mbv(iflag)   ! Target mass, triton or deuteron for D-T or D-D, respectively
            m2 = m2v(iflag)   ! Outgoing particle mass, alpha or triton for D-T or D-D, respectively
        else
            call expirx(1, 'source', 'bad idum(1); it must be 1 or 2')
        end if

        ! Change the problem cell number (i.e, cell number defiend in the MCNP input file) to
        ! the program cell number (i.e, the cell number calculated by MCNP when processing the input file).
        ! For a problem with very large number of cells, put the source cell to be the 1st cell number
        ! in the MCNP input file to speed up the calculations.
        cell_index = 0
        do i = 1, mxa ! mxa is the number of cells in the MCNP input file
            if (ncl(i) == cell_number) then
                cell_index = i
                exit
            end if
        end do
        if (cell_index == 0) then
            call expirx(1, 'source', 'bad idum(2); it must be a cell number')
        end if

        ! Set the max iterations to 1000 or the user value if it is larger
        max_iter = max(max_iter, 1000)

        ! Check deuteron beam energy
        if (eb <= zero) then
            call expirx(1, 'source', 'bad rdum(1); it must be positive')
        end if
        ! Check T/Ti atom ratio
        if (xt <= zero) then
            call expirx(1, 'source', 'bad rdum(2); it must be positive')
        end if
        ! Check position in y direction for source
        if (rdum(4) <= zero) then
            call expirx(1, 'source', 'bad rdum(4); it must be positive')
        end if
        ! Check deuteron beam with
        if (ry <= zero) then
            call expirx(1, 'source', 'bad rdum(6); it must be positive')
        end if

        ! ebb is beam energy converted to eV; eb = rdum(1)
        ebb = eb*1.e6_dknd

        ! xt is T/Ti atom ratio; xt = rdum(2)
        rtitol = one/(one + xt)

        ! set the cutoff distance
        yylab = rdum(4)*1.e8_dknd

        ma2 = ma*ma  ! ma = Deuteron mass
        mb2 = mb*mb  ! mb = Target mass
        m12 = m1*m1  ! m1 = Neutroon mass, triton or deuteron for D-T or D-D, respectively
        m22 = m2*m2  ! m2 = Outgoing mass, alpha for DT reaction or triton for D-D reaction
        mambm = ma - mb
        mambp = ma + mb
        mambm2 = mambm*mambm
        mambp2 = mambp*mambp

        if (isource_flag == 0) then
            isource_flag = 1

            ! Specify z num and atomic mass for D-T reaction, 1: D, 2: Ti, 3: T
            z(1) = one
            z(2) = 22._dknd
            z(3) = one
            am(1) = 2.0141010_dknd
            am(2) = 47.867_dknd
            select case (iflag)
            case (1) ! T target
                am(3) = 3.016049_dknd
            case (2) ! D target
                am(3) = 2.0141010_dknd
            end select

            ! Set energy
            emev = zero  ! varied, MeV
            epeak = 1.0e5_dknd*z(1)**0.67_dknd*am(1) ! fixed, eV

            ! Compute stoich(1), ratio of Ti atm to total in molecule and same for T
            stoich(1) = rtitol       ! ratio of Ti atm to total
            stoich(2) = one - rtitol ! ratio of T atm to total

            ! Compute average z num for the molecule
            ztm = z(2)*stoich(1) + z(3)*stoich(2)
            ! Compute avg atomic mass for the molecule
            atm = am(2)*stoich(1) + am(3)*stoich(2)

            rhod = 4.51_dknd ! target density

            ! Compute the atm density for the molecule
            atmrho = rhod*0.6022_dknd/atm

            ! screeening length of ZBL potential
            a = am(1)/atm
            c = 0.5292_dknd*0.8853_dknd  ! 0.5292 is bohr radius
            ffpa = c/(z(1)**0.23_dknd + ztm**0.23_dknd)          ! Eq 2 in the D-T paper,in ref [1]
            ffpf = ffpa*atm/(z(1)*ztm*14.4_dknd*(am(1) + atm))   ! almost equal to reduced CM energy, without mult by energy
            epsdg = five*ffpf*(one + a)**2/(four*a)              ! the right part of this equation is just 1/gamma in [2]
            ! bohr straggling of de/dx. number "12" normalizes dist.
            ! bohr straggling eq divided by dx and rooted. Units of MeV/cm
            stbohr = 117._dknd*z(1)*sqrt(ztm*atmrho)

            ! atom-atom scattering parameter
            do j = 2, 3
                m2m1(j) = am(1)/am(j)
                ec(j) = four*m2m1(j)/(one + m2m1(j))**2
                a2(j) = c/(z(1)**0.23_dknd + z(j)**0.23_dknd)
                f(j) = a2(j)/(z(1)*z(j)*14.41_dknd*(one + m2m1(j)))
            end do

            ! dE/dx in target, stopping power
            ! Calculate dE/dx as the sum of the individual stopping powers weighted
            do i = 1, npts1  ! npts1 = number of energy bins
                dedx(i) = cdedx(i, iflag)*am(2)/(am(2) + am(3)*xt) + tdedx(i, iflag)*am(3)*xt/(am(2) + am(3)*xt)
            end do
            dedx = 0.01_dknd*rhod*dedx

            ! ****************
            ! This is the relativistic kinematics from the appendix of [1]
            ! In the kinematics equations, 'a' represents the incident particle, the deuteron,
            ! 'b' represents the target, the triton before the reaction. Following the reaction,
            ! the outgoing neutron is represented by '1' and the alhpa by '2'.

            ! In the SOURCE routine, this nomenclature is used for the masses and energies
            ! of the respective particles, e.g., ma for the deuteron mass, mb for the
            ! triton mass, ea for the deuteron energy, and e1 for the neutron energy.

            ! Calculate the cross section for the range of angles and energies of interest
            c = pie/180._dknd
            thl = -five
            pl(1) = one  ! 1st term of Legendre polynomial

            allocate (coscm(npts2))
            allocate (sincm(npts2))
            allocate (xscman(npts1, npts2))
            allocate (xslaban(npts1, npts2))

            do n = 1, npts2  ! npts2 = number of angle bins
                thl = thl + five

                ! Populate thcm with angles in [0,180], spaced by 5 degree increments
                ! The angles are converted to radians before being stored
                ! thcm was defined with fixed dimansion and save.
                thcm(n) = c*thl

                ! The cosine and sine of each angle are saved in coscm and sincm, respectively
                coscm(n) = cos(thcm(n))
                sincm(n) = sqrt(one - coscm(n)*coscm(n))

                ! Legendre polynomials
                pl(2) = coscm(n)

                ! Recursively generate the legendre polynomials
                DO l = 2, 10
                    pl(l + 1) = (pl(l)*coscm(n)*(2*l - 1) - pl(l - 1)*(l - 1))/real(l, dknd)
                end do

                do j = 1, npts1  ! npts1 = number of energy bins
                    ! Center of mass xs calculation
                    ! Uses cross section expansion, ((2L + 1)/2) * PL(mu) * sigl(E)
                    ! Calculates the cross section for each angle at all energies
                    ! This is for the CM frame, the lab result is desired and
                    ! calculated next
                    xscman(j, n) = half*pl(1)
                    do l = 2, 11
                        xscman(j, n) = xscman(j, n) + ((two*(l - one) + one)*half*pl(l)*coeff(l - 1, j, iflag))
                    end do

                    ! The Deuteron energy and momentum is calculated for each energy bin
                    ea = ma + ed(j)           ! ma is Deuteron mass in MeV/c**2
                    pa = sqrt(ea*ea - ma2)    ! ma2 = ma*ma

                    s = ma2 + mb2 + two*mb*ea ! mb2 = mb*mb; mb = Taget mass

                    ! gammacm = (s - ma*ma + mb*mb)/(two*mb*sqrt(s))
                    !         = (two*mb*mb + two*mb*ea)/(two*mb*sqrt(s))
                    !         = (mb + ea)/sqrt(s)
                    gammacm = (mb + ea)/sqrt(s)

                    rg1 = zero
                    ! rlambda = (s - (ma + mb)*(ma + mb))*(s - (ma - mb)*(ma - mb))
                    s2 = s - mambp2   ! mambp2 = (ma+mb)**2
                    s3 = s - mambm2   ! mambm2 = (ma-mb)**2
                    rlambda = s2*s3
                    if (rlambda > zero) then
                        ! s4 = s - ma*ma + mb*mb = 2*mb*mb + 2*mb*ea - 2*mb*(mb+ea)
                        s4 = two*mb*(mb + ea)
                        betacm = sqrt(rlambda)/s4

                        ! m1 = Neutron mass, m2 = Outgoing particle mass
                        s5 = s - m12 - m22  ! s5 = s - m1*m1 - m2*m2
                        s6 = s + m12 - m22  ! s6 = s + m1*m1 - m2*m2
                        rlambda1 = s5*s5 - four*m12*m22
                        if (rlambda1 > zero .and. s6 /= zero) then
                            beta1 = sqrt(rlambda1)/s6
                            rg1 = betacm/beta1
                        end if
                    end if  ! rlambda > zero

                    s = gammacm*(coscm(n) + rg1)
                    if (s /= zero) then
                        theta1 = atan(sincm(n)/s)
                    else
                        theta1 = zero
                    end if
                    if (n > 1) then
                        if (theta1 <= zero) theta1 = theta1 + pie
                    end if

                    ! coslab and sinlab were definded with fixed dimension and save
                    s = cos(theta1)
                    coslab(j, n) = s
                    sinlab(j, n) = sqrt(one - s*s)

                    s = (ea + mb)*(ea + mb) - pa*pa*coslab(j, n)*coslab(j, n)
                    if (s /= zero) then
                        p1t1 = one/s
                        p1t2 = pa*coslab(j, n)*(mb*ea + 0.5_dknd*(ma2 + mb2 + m12 - m22))

                        s3 = (mb*ea + 0.5_dknd*(ma2 + mb2 - m12 - m22))**2
                        s4 = s3 - m12*m22 - m12*pa*pa*sinlab(j, n)*sinlab(j, n)
                        if (s4 > zero) then
                            p1t3 = (ea + mb)*sqrt(s4)
                            p1t = p1t1*(p1t2 + p1t3)
                        else
                            p1t = p1t1*p1t2
                        end if
                    else
                        p1t1 = zero
                        p1t = zero
                    end if

                    e1t1 = p1t1
                    if (e1t1 /= zero) then
                        e1t2 = (ea + mb)*(mb*ea + 0.5_dknd*(ma2 + mb2 + m12 - m22))
                        s = (mb*ea + 0.5_dknd*(ma2 + mb2 - m12 - m22))**2 - m12*m22 - m12*pa*pa*sinlab(j, n)*sinlab(j, n)
                        if (s > zero) then
                            e1t3 = pa*coslab(j, n)*sqrt(s)
                            e1t = e1t1*(e1t2 + e1t3)
                        else
                            e1t = e1t1*e1t2
                        end if
                    else
                        e1t = zero
                    end if

                    if (p1t /= zero .and. e1t /= zero) then
                        s = mb*pa*(p1t*(ea + mb) - e1t*pa*coslab(j, n))
                        if (s /= zero) then
                            xslaban(j, n) = xscman(j, n)*p1t*p1t/s
                        else
                            xslaban(j, n) = zero
                        end if
                    else
                        xslaban(j, n) = zero
                    end if

                end do  ! j loop
            end do      ! n loop

            pemax = zero
            do j = 1, npts2
                do i = 1, npts1
                    ! The probability dist for D-T reaction to occur
                    pe(i, j) = calb*xscman(i, j)*ed2c(i, iflag)/dedx(i)
                    if (pe(i, j) > pemax) pemax = pe(i, j)
                end do
            end do

            ! pe1 is used in srcdx
            do j = 1, npts2
                do i = 1, npts1
                    ! The probability dist for D-T reaction to occur
                    pe1(i, j) = calb*xslaban(i, j)*ed2c(i, iflag)
                end do
            end do

            deallocate (coscm)
            deallocate (sincm)
            deallocate (xscman)
            deallocate (xslaban)
        end if  ! isource_flag == 0

        ! Do rejection for a set maximum number of iterations
        do rej = 1, max_iter
            ! emev, e, distalb, uione1, vione1, wione1 were defined with save.
            ! These variables are varied and computed via Monte Carlo method
            ! (i.e., using physics + random numbers).
            if (emev <= emin) then
                ! Reset values

                ! irtt is used in the rot subroutine
                ! it is decreased each call and renormalizes params at 0
                irtt = mit

                ! ebb is the beam energy in eV
                e = ebb

                distalb = zero  ! distance
                uione1 = zero   ! x-direction
                vione1 = one    ! y-direction
                wione1 = zero   ! z-direction
            end if  ! emev <= emin

            ! Save the current energy
            eold = e  ! in eV

            ! Save current direction
            uione = uione1
            vione = vione1
            wione = wione1

            ! free flight path
            ! eps is set as beam energy times ffpf
            ! ffpf is the reduced CM energy divided by the beam energy
            ! thus eps is the reduced CM energy
            ! ffpf was computed and fixed when the subroutine was called at the 1st time
            eps = e*ffpf

            ! eeg is the same as xi in the paper, following is right side of eq 13 in [1]
            ! epsdg was computed and saved when this subroutine was called at the 1st time
            eeg = sqrt(eps*epsdg)

            ! The following is the left side of eq 13 in [1]
            ! ffpa was computed and saved when this subroutine was called at the 1st time
            pmax1 = ffpa/(eeg + sqrt(eeg) + 0.125_dknd*eeg**0.1_dknd)

            ! Calculate the free flight path using pmax1
            ! atmrho was computed and saved when this subroutine was called at the 1st time
            ffpath = one/(pie*pmax1*pmax1*atmrho)

            ! create ein, takes the value of e which is reduced or reset to ebb
            ! as the particle is tracked. ebb is beam energy defined by rdum(1)
            ein = e/1.0e6_dknd  ! MeV

            ! interpolate a stopping power, exout, based on the current particle energy, ein
            ! dedx was computed and saved (global variables are saved by default) when this subroutine
            ! was called at the 1st time
            call interp(ed, dedx, npts1, ein, exout)
            ! calculate the energy loss from continuous slowing with the stopping power
            dee = exout*ffpath

            ! bohr straggling
            ! uses the bohr straggling val and sqrt(ffpath) as the delta x
            ! atmrho was computed and saved when this subroutine was called at the 1st time
            stdee = (rng_history%rn() - 0.5_dknd)*stbohr*sqrt(ffpath)

            ! If the energy e is less than epeak, reduce stdee by the fraction
            ! epeak was computed when this subroutine was called at the 1st time and
            ! it was defined with save.
            if (e <= epeak) stdee = stdee*e/epeak

            ! dee + stdee  = the energy loss from slowing down and bohr straggling
            ! if it is greater than 0, reduce the particle energy by that amount
            s = dee + stdee
            if (s > zero) e = e - s
            if (e <= zero) then
                emev = zero
                cycle
            end if

            ! use pmax1 and a random num to get a random p
            p = pmax1*sqrt(rng_history%rn())

            ! collision atom
            ! Sample collision atom
            ! Set j, 2 for titanium, 3 for tritium. Used immediately as well as later
            ! rtitol is ratio of Ti atm to total
            if (rng_history%rn() < rtitol) then
                j = 2  ! Titanium sampled
            else
                j = 3  ! Tritium sampled
            end if

            ! calculate the reduced impact parameter for the corresponding atom
            ! a2 was was computed and saved when this subroutine was called at the 1st time
            b = p/a2(j)

            ! set the eps, the reduced CM energy, to that of the corresponding atom
            ! f was computed and saved when this subroutine was called at the 1st time
            eps = e*f(j) ! eps > zero

            ! Rutherford scattering
            ! The reduced CM energy is greater than 10 MeV, do rutherford scattering
            if (eps > emax) then
                ! s2, sin^2(theta/2) is found using eq 10 in [1]
                s2 = one/(one + (one + b*(one + b))*(two*eps*b)**2)
                c2 = one - s2

            else
                ! Magic formula, described in detail in [2]
                ! r is used as the reduced atomic radius below
                ! Set r to reduced impact param; r /= zero
                r = b
                rr = -2.7_dknd*log(abs(eps*b) + fuzz)
                if (rr >= b) then
                    rr = -2.7_dknd*log(abs(eps*rr) + fuzz)
                    if (rr >= b) then
                        r = rr
                    end if
                end if

                ! below is zbl universal int. potential
                do icycle = 1, 10
                    if (r > 200._dknd) r = 200._dknd

                    ! universal potential function, eq 3 in [1]
                    ex1 = 0.18175_dknd*exp(-3.1998_dknd*r)
                    ex2 = 0.50986_dknd*exp(-0.94229_dknd*r)
                    ex3 = 0.28022_dknd*exp(-0.40290_dknd*r)
                    ex4 = 0.02817_dknd*exp(-0.20162_dknd*r)

                    ! Interatomic potential function ( phi (ex1 to 4 summed) / r)
                    ! The elementary charge is set to 1 (natural units)
                    v = (ex1 + ex2 + ex3 + ex4)/r

                    ! This is the correct derivative of v, product rule.
                    v1 = -(v + 3.1998_dknd*ex1 + 0.94229_dknd*ex2 + 0.4029_dknd*ex3 + 0.20162_dknd*ex4)/r

                    ! The following represent terms that appear in the scattering integral, eq 5 in [1]
                    fr = b*b/r + v*r/eps - r
                    fr1 = -b*b/(r*r) + (v + v1*r)/eps - one
                    if (fr /= zero .and. fr1 /= zero) then
                        qq = fr/fr1
                        r = r - qq
                        if (r == zero) exit
                        if (abs(qq/r) > 0.0001_dknd) then
                            cycle
                        else
                            exit
                        end if
                    else
                        exit
                    end if
                end do  ! icycle

                ! roc is Rc from eq 8 in [1]
                roc = -two*(eps - v)/v1

                ! the reduced CM energy is rooted
                sqe = sqrt(eps)  !  eps > zero

                ! The following calculations are for delta, which is described in [1]
                ! Represents beta in [2]
                cc = (0.011615_dknd + sqe)/(0.0071222_dknd + sqe)

                ! Represents A in [2] and cc was used
                ! It may be in a bug in an old code in SINBAD since c was used
                aa = two*eps*(one + (0.99229_dknd/sqe))*b**cc

                ! whole term is G inverse in [2]
                ff = (sqrt(aa*aa + one) - aa)*((9.3066_dknd + eps)/(14.813_dknd + eps))

                ! calculate delta for the magic formula
                delta = (r - b)*aa*ff/(ff + one)  ! ff > zero
                if (r /= -roc) then
                    c = (b + delta + roc)/(r + roc)
                    c2 = c*c      ! cos(theta/2)**2
                    s2 = one - c2 ! sin(theta/2)**2
                else
                    c2 = zero
                    s2 = one
                end if

            end if ! eps > emax

            ! convert the cm scattering angle from the previous calculations to lab
            ! Use double angle formula to get ct = cos(theta)
            ct = two*c2 - one

            ! Use trig identity to get st = sin(theta)
            st = sqrt(1.0001_dknd - ct*ct)

            ! Add a term for the conversion to lab frame
            cu = ct + m2m1(j)

            ! psi is the scattering angle in the lab frame
            if (cu /= zero) then
                psi = atan(st/cu)
                if (psi < zero) psi = psi + pie
            else
                psi = zero
            end if

            ! use the cm angle to calculate the energy transfer to the target (eq 6 from [1])
            ! ec is exactly equivalent to gamma in eq 6
            ! ec was computed and saved when this routine was called at the 1st time.
            deny = ec(j)*s2*e

            ! remove the energy transferred from the D energy
            e = e - deny
            if (e <= zero) then
                emev = zero
                cycle
            end if

            ! calculate the cosine of the lab scattering angle
            uzero = cos(psi)

            ! call rot to sample a new set of directions using the scat cos. sample azimuthal uniformly
            call rot(uzero, uione, vione, wione, uione1, vione1, wione1, irtt)

            ! Turn e into MeV
            emev = e/1.e6_dknd
            if (emev < emin) cycle

            ! calculates the total change in energy for the flight
            ! Adds the random sampled energy to emev to get ep.
            ep = emev + rng_history%rn()*(eold/1.e6_dknd - emev)
            if (ep <= zero) then
                emev = zero
                cycle
            end if

            ! thickness condition
            ! increment distance travelled in y dir
            distalb = distalb + ffpath*vione
            if (distalb < zero .or. distalb > yylab) then
                emev = zero
                cycle
            end if

            ! reaction
            ! sample a random cos from -1 to 1
            t2 = two*rng_history%rn() - one

            ! get the angle from the random cosine, outgoing angle
            a1 = acos(t2)

            ! call the 2D interp subroutine
            ! uses pe, the probability distribution for a D-T reaction to occur
            ! takes in the value ep for the energy, calculated angle and gets pn
            call xinterp2D(ed, thcm, pe, npts1, npts2, ep, a1, pn)

            ! pemax was computed and saved when this routine was called at the 1st time.
            if (rng_history%rn()*pemax <= pn) exit
        end do  ! rej loop

        if (rej > max_iter) then
            call expirx(1, 'source', 'bad DT reaction rejection sampling; increase idum(3)')
        end if

        ! Set the val of the random cosine that was used.
        ! -one <= t2 <= one
        ctn = t2                ! cosine
        stn = sqrt(one - t2*t2) ! sine

        ! Relativistic kinematics
        ! sum of Deuteron rest mass energy (ma) and KE (ep); ep > zero
        ean = ma + ep

        ! momentum calc for the D using the energy momentum eq and natural units
        pana = sqrt(ean*ean - ma2) ! pana > 0 since ean > ma

        if (ctn == -one) then
            uzero = ctn
        else if (ctn == one) then
            uzero = ctn
        else
            ! s mandelstam var, eq 1 in appendix
            ! ma2 = ma*ma, ma = Deuteron mass
            ! mb2 = mb*mb, mb = Target mass
            s = ma2 + mb2 + two*mb*ean

            ! gamma is the lorentz factor
            ! gammacm = (s - ma2 + mb2)/(two*mb*sqrt(s))
            !         = (two*mb*mb + two*mb*ean)/(two*mb*sqt(s))
            !         =  (mb + ean)/sqrt(s)
            gammacm = (mb + ean)/sqrt(s)

            ! kinematical function lambda, eq 12 in appendix
            ! mambp2 = (ma + mb)**2, mambm2 = (ma - mb)**2
            rg1 = zero
            rlambda = (s - mambp2)*(s - mambm2)
            if (rlambda > zero) then
                ! s2 = s - ma*ma + mb*mb = 2*mb*mb + 2*mb*ean = 2*mb*(mb+ean)
                s2 = two*mb*(mb + ean)
                betacm = sqrt(rlambda)/s2

                ! m12 = m1*m1, m1 = Neutron mass
                ! m22 = m2*m2, m2 = Outgoing particle mass
                s3 = s - m12 - m22
                s4 = s + m12 - m22
                rlambda1 = s3*s3 - four*m12*m22
                if (rlambda1 > zero .and. s4 /= zero) then
                    beta1 = sqrt(rlambda1)/s4
                    rg1 = betacm/beta1
                end if
            end if

            ! calculate lab angle corresponding to randomly sampled cm angle
            s = gammacm*(ctn + rg1)
            if (s /= zero) then
                ph = atan(stn/s)
                if (ph < zero) ph = ph + pie
                uzero = cos(ph)
            else
                uzero = zero
            end if

        end if  ! ctn == -one

        ! Energy calc from eq 22 in the appendix of [1], uses the angle calculated from kinematics
        s = (ean + mb)*(ean + mb) - (pana*pana)*(uzero*uzero)
        if (s /= zero) then
            ! calculate the sin of the lab angle that was found
            sinlab1 = sqrt(one - uzero*uzero)

            t1n = one/s
            t2n = (ean + mb)*(mb*ean + half*(ma2 + mb2 + m12 - m22))

            ! calculate the sin of the lab angle that was found
            s2 = (mb*ean + half*(ma2 + mb2 - m12 - m22))**2 - m12*(m22 - pana*pana*sinlab1*sinlab1)
            if (s2 > zero) then
                t3n = pana*uzero*sqrt(s2)
                e1n = t1n*(t2n + t3n) ! energy calc, eq 22 in appendix of [1]
            else
                e1n = t1n*t2n
            end if
        else
            e1n = zero
        end if

        ! Subtracts the neutron mass energy from the resulting energy to get KE
        ! m1 = Neutron mass, MeV/c^2
        if (e1n > m1) then
            pbl%r%erg = e1n - m1
        else
            pbl%r%erg = very_small
        end if

        ! Call to rot to compute the neutron direction
        ! Uses lab angle corresponding to randomly sampled CM angle
        ! Take in direction of the incident D and adds the scattering angle to get
        ! direction for the outgoing neutron.
        call rot(uzero, uione, vione, wione, pbl%r%u, pbl%r%v, pbl%r%w, irtt)

        ry2 = ry*ry  ! Deuteron beam width squared
        do rej = 1, max_iter
            ! Randomly sample an x1 and z1 from -ry to ry
            x1 = (ry - two*ry*rng_history%rn())
            z1 = (ry - two*ry*rng_history%rn())
            rr = x1*x1 + z1*z1
            if (rr <= ry2) exit
        end do 
        if (rej > max_iter) then
            call expirx(1, 'source', 'bad in position rejection sampling; incrase idum(3)')
        end if

        ! Set particle starting position
        pbl%r%x = rdum(3) + x1
        pbl%r%y = rdum(4)
        pbl%r%z = rdum(5) + z1

        ! Set other required particle parameters
        pbl%i%icl = cell_index
        pbl%i%jsu = 0
        pbl%i%ipt = 1  ! neutron
        pbl%r%tme = zero
        pbl%r%wgt = one

    end subroutine source

! *********************************************************
    subroutine xinterp2d(x1, x2, y, npts1, npts2, xin1, xin2, yout)

        integer, intent(in) :: npts1, npts2
        real(dknd), intent(in) :: xin1, xin2
        real(dknd), intent(in), dimension(npts1) :: x1
        real(dknd), intent(in), dimension(npts2) :: x2
        real(dknd), intent(in), dimension(npts1, npts2) :: y
        real(dknd), intent(out) :: yout

        integer :: i, j, ii, jj
        real(dknd) :: v1, v2, v3, v4, denom

        real(dknd), parameter :: eps = 0.00001_dknd

        intrinsic :: abs, max

        ii = 0
        do i = 1, npts1 - 1
            if (x1(i + 1) <= x1(i)) then
                call expirx(1, 'source', 'xinterp2D: x1 is not monotonically increasing')
            end if
            if ((x1(i) <= xin1) .and. (xin1 <= x1(i + 1))) then
                ii = i
                exit
            end if
        end do
        if (ii == 0) then
            call expirx(1, 'source', 'xinterp2D: ii = 0')
        end if

        jj = 0
        do j = 1, npts2 - 1
            if (x2(j + 1) <= x2(j)) then
                call expirx(1, 'source', 'xinterp2D: x2 is not monotonically increasing')
                return
            end if
            if ((x2(j) <= xin2) .and. (xin2 <= x2(j + 1))) then
                jj = j
                exit
            end if
        end do
        if (jj == 0) then
            call expirx(1, 'source', 'xinterp2D: jj = 0')
        end if

        denom = ((x1(ii + 1) - x1(ii))*(x2(jj + 1) - x2(jj))) ! denom /= zero
        denom = one/denom
        v1 = x1(ii + 1) - xin1
        v2 = x2(jj + 1) - xin2
        v3 = xin2 - x2(jj)
        v4 = xin1 - x1(ii)
        yout = v1*v2*denom*y(ii, jj) + &
               v1*v3*denom*y(ii, jj + 1) + &
               v4*v2*denom*y(ii + 1, jj) + &
               v4*v3*denom*y(ii + 1, jj + 1)

    end subroutine xinterp2d

! *********************************************************
! Private subroutines/functions
! *********************************************************
    subroutine interp(x, y, npts, xin, yout)

        integer, intent(in) :: npts
        real(dknd), intent(in), dimension(npts) :: x
        real(dknd), intent(in), dimension(npts) :: y
        real(dknd), intent(in) :: xin
        real(dknd), intent(out) :: yout

        real(dknd), parameter :: eps = 0.00001_dknd

        integer :: i
        real(dknd) :: ax, axin

        intrinsic :: abs, max

        do i = 2, npts
            if (x(i) <= x(i - 1)) then
                call expirx(1, 'interp', 'x is not monotonically increasing')
            end if
        end do

        axin = abs(xin)
        do i = 1, npts
            ax = abs(x(i))
            if (abs(x(i) - xin) < eps*max(ax, axin)) then
                yout = y(i)
                return
            end if
            if (i > 1 .and. x(i) > xin) then
                yout = y(i - 1) + (y(i) - y(i - 1))/(x(i) - x(i - 1))*(xin - x(i - 1))
                return
            end if
        end do
        call expirx(1, 'source', 'interp: 1d interpolation error')

    end subroutine interp

! *********************************************************
    subroutine rot(c, uold, vold, wold, uuu, vvv, www, irtt)

        integer, intent(inout) :: irtt
        real(dknd), intent(in) :: c
        real(dknd), intent(in) :: uold
        real(dknd), intent(in) :: vold
        real(dknd), intent(in) :: wold
        real(dknd), intent(inout) :: uuu
        real(dknd), intent(inout) :: vvv
        real(dknd), intent(inout) :: www

        real(dknd), parameter :: aone = 1.0001_dknd
        real(dknd), parameter :: mval = 0.9_dknd

        real(dknd) :: ac, t1, t2, r, s, t

        intrinsic :: abs, sqrt

        ac = abs(c)
        if (ac > aone) then
            call expirx(1, 'rot', 'invalid cosine scattering')
        else if (ac >= one) then  !   1.0 <= ac <= aone
            uuu = c*uold
            vvv = c*vold
            www = c*wold
        else  !   zero =< ac < one
            do
                t1 = two*rng_history%rn() - one
                t2 = two*rng_history%rn() - one
                r = t1*t1 + t2*t2
                if (abs(r) <= very_small) cycle
                if (r <= one) exit
            end do
            r = sqrt((one - c*c)/r)
            t1 = t1*r
            t2 = t2*r
            if (abs(wold) > mval) then
                s = sqrt(uold*uold + wold*wold)
                t = one/s
                uuu = uold*c + (t1*uold*vold + t2*wold)*t
                vvv = vold*c - t1*s
                www = wold*c + (t1*wold*vold - t2*uold)*t
            else
                s = sqrt(uold*uold + vold*vold)
                if (s <= very_small) then
                    uuu = uold*c + t1*uold*wold - t2*vold
                    vvv = vold*c + t1*vold*wold + t2*uold
                else
                    t = one/s
                    uuu = uold*c + (t1*uold*wold - t2*vold)*t
                    vvv = vold*c + (t1*vold*wold + t2*uold)*t
                end if
                www = wold*c - t1*s
            end if
        end if

        ! renormalize to prevent error buildup.
        irtt = irtt - 1
        if (irtt /= 0) return
        irtt = mit
        s = sqrt(uuu*uuu + vvv*vvv + www*www)
        if (s > very_small) then
            t = one/s
            uuu = uuu*t
            vvv = vvv*t
            www = www*t
        end if

    end subroutine rot

end module source_mod
