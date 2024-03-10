
!=======================================
!
! H4 hydrogen bond correction for semiempirical methods - Version 1.2D0
!
! Includes H-H repulsion correction from D3H4 scheme
!
! Reference: J. Rezac, P. Hobza J. Chem. Theory Comput. 8, 141-151 (2012)
!            http:!dx%doi%org/10.1021D0/ct200751e
!
! Code copyright (c) 2011 Jan Rezac - see the license at the end of the file
! E-mail rezac@uochb%cas%cz
!
! The code contains only the parameters for PM6 method. Other sets of parameters
! are listed in the paper.
!
! Compile using any C comiler, linking math library. Example for gcc:
! gcc -lm -o h_bonds4 h_bonds4.c
!
! Use: The program reads coordinates in .xyz format from standard input,
! prints the correction energy (in kcal/mol) and gradient (kcal/mol/A)
! to standard output:
! ./h_bonds4 .lt. geometry%xyz
!
!=======================================

!=======================================
! Changelog
!
! 1.1D0 (Jan 18, 2012) - Code for H-H repulsion added
! 1.2D0 (Jan 24, 2012) - Analytical derivatives
!=======================================

!=======================================
! Define type used for coordinates
!=======================================


module my_types

    implicit none
    ! Cartesian coordinates
    type :: coord_t
        real*8 ::       x
        real*8 ::       y
        real*8 ::       z
    end type coord_t

    ! Atom: coordinates + proton number
    type :: atom_t
        real*8 ::       x
        real*8 ::       y
        real*8 ::       z
        integer ::      e
    end type atom_t

    !=======================================
    ! Tabular data
    !=======================================
    ! Element names
    character(LEN=5), dimension(0:118) :: element_names = (/"X", "H", "He", "Li", "Be", "B", "C", "N",&
                                                "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",&
                                                "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",&
                                                "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",&
                                                "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",&
                                                "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",&
                                                "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",&
                                                "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",&
                                                "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",&
                                                "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub", "Uut", "Uuq", "Uup",&
                                                "Uuh", "Uus", "Uuo"/)
    ! Covalent radii (indexed by proton number), zero -.gt. not available
    real*8, dimension(0:118) :: covalent_radii = (/0.0D0, 0.37D0, 0.32D0, 1.34D0, 0.9D0, 0.82D0, 0.77D0,&
                                                0.75D0, 0.73D0, 0.71D0, 0.69D0, 1.54D0, 1.3D0, 1.18D0, 1.11D0, 1.06D0, 1.02D0, 0.99D0, 0.97D0,&
                                                1.96D0, 1.74D0, 1.44D0, 1.36D0, 1.25D0, 1.27D0, 1.39D0, 1.25D0, 1.26D0, 1.21D0, 1.38D0, 1.31D0,&
                                                1.26D0, 1.22D0, 1.19D0, 1.16D0, 1.14D0, 1.1D0, 2.11D0, 1.92D0, 1.62D0, 1.48D0, 1.37D0, 1.45D0,&
                                                1.56D0, 1.26D0, 1.35D0, 1.31D0, 1.53D0, 1.48D0, 1.44D0, 1.41D0, 1.38D0, 1.35D0, 1.33D0, 1.3D0,&
                                                2.25D0, 1.98D0, 1.69D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,&
                                                0.0D0, 0.0D0, 1.6D0, 1.5D0, 1.38D0, 1.46D0, 1.59D0, 1.28D0, 1.37D0, 1.28D0, 1.44D0, 1.49D0, 0.0D0,&
                                                0.0D0, 1.46D0, 0.0D0, 0.0D0, 1.45D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,&
                                                0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,&
                                                0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)

    !=======================================
    ! Constants
    !=======================================
    integer, parameter :: HYDROGEN=1
    integer, parameter :: CARBON=6
    integer, parameter :: NITROGEN=7
    integer, parameter :: OXYGEN=8

    real*8, parameter :: M_PI=3.14159265359D0

    !=======================================
    ! Cutoffs
    !=======================================

    ! Cutoff for the correction, more distant donor-acceptor pairs do not contribute
    real*8, parameter :: HB_R_CUTOFF=5.5D0
    ! Short-range cutoff, closer donor-acceptor pairs do not contribute
    real*8, parameter :: HB_R_0=1.5D0
    ! max. X-H covalent bond distance
    real*8, parameter :: MAX_XH_BOND=1.15D0

    !=======================================
    ! Parameters (For PM6-D3)
    !=======================================

    ! H4 correction
    real*8, parameter :: para_oh_o = 2.32D0
    real*8, parameter :: para_oh_n = 3.10D0
    real*8, parameter :: para_nh_o = 1.07D0
    real*8, parameter :: para_nh_n = 2.01D0
    real*8, parameter :: multiplier_wh_o = 0.42D0
    real*8, parameter :: multiplier_nh4 = 3.61D0
    real*8, parameter :: multiplier_coo = 1.41D0

    contains
    !=======================================
    ! Auxiliary functions
    !=======================================

    !------------------------------------------------------------------------------
    ! Power function
    real*8 function pow (a,b)
            real*8 :: a
            integer :: b
            pow = a**b
            return
    end function

    !------------------------------------------------------------------------------
    ! Distance between two atoms
    real*8 function distance (a,b)
            type(atom_t) :: a,b
            distance = dsqrt((a%x - b%x)*(a%x - b%x) + (a%y - b%y)*(a%y - b%y) + (a%z - b%z)*(a%z - b%z))
            return
    end function

    !------------------------------------------------------------------------------
    ! Angle between three atoms A-B-C
    real*8 function atomangle (a,b,c)
            ! Two vectors ...
            type(atom_t) :: a,b,c
            real*8 :: ux, uy, uz
            real*8 :: vx, vy, vz
            real*8 :: abs1, abs2, dot, cosi
            ux = a%x - b%x
            uy = a%y - b%y
            uz = a%z - b%z
            vx = c%x - b%x
            vy = c%y - b%y
            vz = c%z - b%z
            abs1 = dsqrt(ux*ux + uy*uy + uz*uz)
            abs2 = dsqrt(vx*vx + vy*vy + vz*vz)
            if (abs1 .eq. 0.0D0 .or. abs2 .eq. 0.0D0) then
                write (*,*) "Coordinate for angle calculation can not be zero"
                call exit(1)
            end if
            dot = ux*vx + uy*vy +uz*vz
            cosi = dot/abs1/abs2
            ! Numerical issue: can be close to 1 but out of interval:
            if (cosi .lt. -1.0D0) cosi = -1.0D0
            if (cosi .gt. 1.0D0) cosi = 1.0D0
            atomangle=dacos(cosi)
            return
    end function

    !------------------------------------------------------------------------------
    ! Continuous valence contribution of a pair of atoms
    real*8 function cvalence_contribution(a, b)
            type(atom_t) :: a,b
            real*8 :: r
            real*8 :: ri, rj
            real*8 :: r0, r1
            real*8 :: x
            ri = covalent_radii(a%e)
            rj = covalent_radii(b%e)
            r0 = ri + rj
            r1 = r0 * 1.6D0
            r = distance(a,b)
            if (r .eq. 0.0D0) then
                cvalence_contribution = 0.0D0
                return
            end if
            if (r .ge. r1) then
                cvalence_contribution = 0.0D0
                return
            end if
            if (r .le. r0) then
                cvalence_contribution = 1.0D0
                return
            end if
            x = (r - r0) / (r1 - r0)
            cvalence_contribution = 1.0D0 - (-20.0D0*pow(x,7) + 70.0D0*pow(x,6) -84.0D0*pow(x,5) + 35.0D0*pow(x,4))
            return
    end function

    !------------------------------------------------------------------------------
    ! Continuous valence contribution of a pair of atoms - derivative in the internal coordinate
    real*8 function cvalence_contribution_d(a, b)
            type(atom_t) :: a,b
            real*8 :: r
            real*8 :: ri, rj
            real*8 :: r0, r1
            real*8 :: x
            ri = covalent_radii(a%e)
            rj = covalent_radii(b%e)
            r0 = ri + rj
            r1 = r0 * 1.6D0
            r = distance(a,b)
            if (r .eq. 0.0D0) then
                cvalence_contribution_d = 0.0D0
                return
            end if
            if (r .ge. r1) then
                cvalence_contribution_d = 0.0D0
                return
            end if
            if (r .le. r0) then
                cvalence_contribution_d = 0.0D0
                return
            end if
            x = (r - r0) / (r1 - r0)
            cvalence_contribution_d = -(-140.0D0*pow(x,6) + 420.0D0*pow(x,5) - 420.0D0*pow(x,4) + 140.0D0*pow(x,3)) / (r1 - r0)
            return
    end function

    !------------------------------------------------------------------------------
    ! Write array of coordinates (gradient) to standard output
    subroutine gradient_write(natom, crd)
            type(coord_t),dimension(:) :: crd
            integer :: i, natom
            do i = 1, natom
                write(*,'(F14.6,1X,F14.6,1X,F14.6)') crd(i)%x, crd(i)%y, crd(i)%z
            end do
            return
    end subroutine

    !------------------------------------------------------------------------------
    ! Multiplication of coordinate vector by a number
    type(coord_t) function coord_scale(coord, factor)
            type(coord_t) :: coord,resul
            real*8 :: factor
            resul%x = coord%x * factor
            resul%y = coord%y * factor
            resul%z = coord%z * factor
            coord_scale = resul
            return
    end function

    !------------------------------------------------------------------------------
    ! Coordinate vector addition to array of coordinates used to construct the total gradient from atomic contributions
    subroutine coord_add(coord, i, add)
            type(coord_t),dimension(:), intent(INOUT) :: coord
            type(coord_t) :: add
            integer :: i
            coord(i)%x = coord(i)%x + add%x
            coord(i)%y = coord(i)%y + add%y
            coord(i)%z = coord(i)%z + add%z
            return
    end subroutine

    !=======================================
    ! H4 correction calculation
    !=======================================
    real*8 function energy_corr_h4(natom, geo, grad, do_grad)
        real*8 :: e_corr_sum = 0.0D0 ! correction energy
        logical :: do_grad
        integer :: natom
        type(atom_t), dimension(:) :: geo
        type(coord_t), dimension(:), intent(INOUT) :: grad

        ! H-bond description
        integer :: d_i, a_i, h_i ! donor, acceptor, hydrogen indices
        real*8 :: rda ! donor-acceptor distance
        real*8 :: rdh, rah
        real*8 :: angle


        ! Energy terms:
        real*8 :: e_para
        real*8 :: e_bond_switch
        real*8 :: e_radial
        real*8 :: e_angular
        real*8 :: e_scale_w
        real*8 :: e_scale_chd
        real*8 :: e_scale_cha
        real*8 :: e_corr

        ! Derivatives
        real*8 :: d_radial
        type(coord_t) :: d_radial_d, d_radial_a

        real*8 :: d_angular
        type(coord_t) :: d_angular_d, d_angular_h, d_angular_a

        real*8 :: d_bs
        type(coord_t) :: d_bs_d, d_bs_a, d_bs_h

        type(coord_t) :: g

        ! Scaling derivatives
        real*8 :: sign_wat

        integer :: o1
        integer :: o2
        integer :: cc

        real*8 :: cv_o1
        real*8 :: cv_o2
        real*8 :: cv_cc

        real*8 :: f_o1
        real*8 :: f_o2
        real*8 :: f_cc

        ! Auxiliary variables
        integer :: i, j, k ! iteration counters
        real*8 :: rih, rjh ! auxiliary distances
        real*8 :: x, xd, xd2, a, d
        real*8 :: slope, v, fv, fv2
        real*8 :: rdhs, ravgs
        real*8 :: sign, hydrogens, others
        real*8 :: cdist, odist

        e_corr_sum = 0.0D0
        ! Iterate over donor/acceptor pairs
        do i = 1, natom
            if (geo(i)%e .eq. NITROGEN .or. geo(i)%e .eq. OXYGEN) then
                do j = 1, i-1
                    if (geo(j)%e .eq. NITROGEN .or. geo(j)%e .eq. OXYGEN) then
                        ! Calculate donor-acceptor distance
                        rda = distance(geo(i), geo(j))
                        ! Continue only when in range where correction acts
                        if (rda .gt. HB_R_0 .and. rda .lt. HB_R_CUTOFF) then
                            ! Iterate over hydrogens
                            do h_i = 1, natom
                                if (geo(h_i)%e .eq. HYDROGEN) then
                                    ! Distances to hydrogen
                                    rih = distance(geo(i), geo(h_i))
                                    rjh = distance(geo(j), geo(h_i))

                                    angle = M_PI - atomangle(geo(i), geo(h_i), geo(j))
                                    ! if (rih*rih + rjh*rjh .lt. rda*rda) then
                                    if (angle .lt. M_PI/2.0D0) then
                                        ! Here, we have filterd out everything but corrected H-bonds
                                        ! Determine donor and acceptor - donor is the closer one
                                        if (rih .le. rjh) then
                                            d_i = i
                                            a_i = j
                                            rdh = rih
                                            rah = rjh
                                        else
                                            d_i = j
                                            a_i = i
                                            rdh = rjh
                                            rah = rih
                                        end if

                                        ! Radial term
                                        e_radial = -0.00303407407407313510D0 * pow(rda,7) + &
                                                    0.07357629629627092382D0 * pow(rda,6) + &
                                                    -0.70087111111082800452D0 * pow(rda,5) + &
                                                    3.25309629629461749545D0 * pow(rda,4) + &
                                                    -7.20687407406838786983D0 * pow(rda,3) + &
                                                    5.31754666665572184314D0 * pow(rda,2) + &
                                                    3.40736000001102778967D0 * rda + &
                                                    -4.68512000000450434811D0

                                        ! Radial gradient
                                        if (do_grad) then
                                            ! In rDA coordinate
                                            d_radial = -0.02123851851851194655D0 * pow(rda,6) + &
                                                        0.44145777777762551519D0 * pow(rda,5) + &
                                                        -3.50435555555413991158D0 * pow(rda,4) + &
                                                        13.01238518517846998179D0 * pow(rda,3) + &
                                                        -21.62062222220516360949D0 * pow(rda,2) + &
                                                        10.63509333331144368628D0 * rda + &
                                                        3.40736000001102778967D0

                                            ! Cartesian gradients on D and A atoms
                                            d_radial_d%x = (geo(d_i)%x - geo(a_i)%x)/rda * d_radial
                                            d_radial_d%y = (geo(d_i)%y - geo(a_i)%y)/rda * d_radial
                                            d_radial_d%z = (geo(d_i)%z - geo(a_i)%z)/rda * d_radial

                                            d_radial_a%x = -d_radial_d%x
                                            d_radial_a%y = -d_radial_d%y
                                            d_radial_a%z = -d_radial_d%z
                                        end if

                                        ! Angular term
                                        a = angle/(M_PI/2.0D0)
                                        x = -20.0D0*pow(a,7) + 70.0D0*pow(a,6) - 84.0D0*pow(a,5) + 35.0D0*pow(a,4)
                                        e_angular = 1.0D0 - x*x

                                        ! Angular gradient
                                        if (do_grad) then
                                            xd = (-140.0D0*pow(a,6) + 420.0D0*pow(a,5) - 420.0D0*pow(a,4) + 140.0D0*pow(a,3)) / (M_PI/2.0D0)
                                            d_angular = -xd * 2.0D0 * x

                                            ! Dot product of bond vectors
                                            d = (geo(d_i)%x - geo(h_i)%x)*(geo(a_i)%x - geo(h_i)%x) + &
                                                (geo(d_i)%y - geo(h_i)%y)*(geo(a_i)%y - geo(h_i)%y) + &
                                                (geo(d_i)%z - geo(h_i)%z)*(geo(a_i)%z - geo(h_i)%z)

                                            x = -d_angular / dsqrt(1.0D0 - (d*d) / (rdh*rdh) / (rah*rah))

                                            ! Donor atom
                                            d_angular_d.x = x * -((geo(a_i)%x - geo(h_i)%x)/rdh/rah - (geo(d_i)%x - geo(h_i)%x)*d/pow(rdh,3)/rah)
                                            d_angular_d.y = x * -((geo(a_i)%y - geo(h_i)%y)/rdh/rah - (geo(d_i)%y - geo(h_i)%y)*d/pow(rdh,3)/rah)
                                            d_angular_d.z = x * -((geo(a_i)%z - geo(h_i)%z)/rdh/rah - (geo(d_i)%z - geo(h_i)%z)*d/pow(rdh,3)/rah)
                                            ! Acceptor atom
                                            d_angular_a.x = x * -((geo(d_i)%x - geo(h_i)%x)/rdh/rah - (geo(a_i)%x - geo(h_i)%x)*d/rdh/pow(rah,3))
                                            d_angular_a.y = x * -((geo(d_i)%y - geo(h_i)%y)/rdh/rah - (geo(a_i)%y - geo(h_i)%y)*d/rdh/pow(rah,3))
                                            d_angular_a.z = x * -((geo(d_i)%z - geo(h_i)%z)/rdh/rah - (geo(a_i)%z - geo(h_i)%z)*d/rdh/pow(rah,3))
                                            ! Hydrogen
                                            d_angular_h%x = -d_angular_d%x - d_angular_a%x
                                            d_angular_h%y = -d_angular_d%y - d_angular_a%y
                                            d_angular_h%z = -d_angular_d%z - d_angular_a%z
                                        end if

                                        ! Energy coefficient
                                        if (geo(d_i)%e .eq. OXYGEN .and. geo(a_i)%e .eq. OXYGEN)     e_para = para_oh_o
                                        if (geo(d_i)%e .eq. OXYGEN .and. geo(a_i)%e .eq. NITROGEN)   e_para = para_oh_n
                                        if (geo(d_i)%e .eq. NITROGEN .and. geo(a_i)%e .eq. OXYGEN)   e_para = para_nh_o
                                        if (geo(d_i)%e .eq. NITROGEN .and. geo(a_i)%e .eq. NITROGEN) e_para = para_nh_n

                                        ! Bond switching
                                        if (rdh .gt. 1.15D0) then
                                            rdhs = rdh - 1.15D0
                                            ravgs = 0.5D0*rdh + 0.5D0*rah - 1.15D0
                                            x = rdhs/ravgs
                                            e_bond_switch = 1.0D0-(-20.0D0*pow(x,7) + 70.0D0*pow(x,6) - 84.0D0*pow(x,5) + 35.0D0*pow(x,4))

                                            ! Gradient
                                            if (do_grad) then
                                                d_bs = -(-140.0D0*pow(x,6) + 420.0D0*pow(x,5) - 420.0D0*pow(x,4) + 140.0D0*pow(x,3))

                                                xd = d_bs / ravgs
                                                xd2 = 0.5D0 * d_bs * -x / ravgs

                                                d_bs_d%x = (geo(d_i)%x - geo(h_i)%x)/rdh * xd + (geo(d_i)%x - geo(h_i)%x)/rdh * xd2
                                                d_bs_d%y = (geo(d_i)%y - geo(h_i)%y)/rdh * xd + (geo(d_i)%y - geo(h_i)%y)/rdh * xd2
                                                d_bs_d%z = (geo(d_i)%z - geo(h_i)%z)/rdh * xd + (geo(d_i)%z - geo(h_i)%z)/rdh * xd2

                                                d_bs_a%x = (geo(a_i)%x - geo(h_i)%x)/rah * xd2
                                                d_bs_a%y = (geo(a_i)%y - geo(h_i)%y)/rah * xd2
                                                d_bs_a%z = (geo(a_i)%z - geo(h_i)%z)/rah * xd2

                                                d_bs_h%x = -d_bs_d%x + -d_bs_a%x
                                                d_bs_h%y = -d_bs_d%y + -d_bs_a%y
                                                d_bs_h%z = -d_bs_d%z + -d_bs_a%z
                                            end if
                                        else
                                            ! No switching, no gradient
                                            e_bond_switch = 1.0D0
                                            if (do_grad) then
                                                d_bs_d%x = 0.0D0; d_bs_d%y = 0.0D0; d_bs_d%z = 0.0D0
                                                d_bs_a%x = 0.0D0; d_bs_a%y = 0.0D0; d_bs_a%z = 0.0D0
                                                d_bs_h%x = 0.0D0; d_bs_h%y = 0.0D0; d_bs_h%z = 0.0D0
                                            end if
                                        end if

                                        ! Water scaling
                                        e_scale_w = 1.0D0
                                        if (geo(d_i)%e .eq. OXYGEN .and. geo(a_i)%e .eq. OXYGEN) then
                                            ! Count hydrogens and other atoms in vicinity
                                            hydrogens = 0.0D0
                                            others = 0.0D0
                                            do k = 1, natom
                                                if (geo(k)%e .eq. HYDROGEN) then
                                                    hydrogens = hydrogens + cvalence_contribution(geo(d_i), geo(k))
                                                else
                                                    others = others + cvalence_contribution(geo(d_i), geo(k))
                                                end if
                                            end do

                                            ! If it is water
                                            if (hydrogens .ge. 1.0D0 ) then
                                                sign_wat = 1.0D0
                                                slope = multiplier_wh_o - 1.0D0
                                                v = hydrogens
                                                fv = 0.0D0
                                                if (v .gt. 1.0D0 .and. v .le. 2.0D0) then
                                                    fv = v - 1.0D0
                                                    sign_wat = 1.0D0
                                                end if
                                                if (v .gt. 2.0D0 .and. v .lt. 3.0D0) then
                                                    fv = 3.0D0 - v
                                                    sign_wat = -1.0D0
                                                end if
                                                fv2 = 1.0D0 - others
                                                if (fv2 .lt. 0.0D0) fv2 = 0.0D0

                                                e_scale_w = 1.0D0 + slope * fv * fv2
                                            end if
                                        end if

                                        ! Charged groups
                                        e_scale_chd = 1.0D0
                                        e_scale_cha = 1.0D0

                                        ! Scaled groups: NR4+
                                        if (geo(d_i)%e .eq. NITROGEN) then
                                            slope = multiplier_nh4 - 1.0D0
                                            v = 0.0D0
                                            do k = 1, natom
                                                v = v + cvalence_contribution(geo(d_i), geo(k))
                                            end do
                                            if (v .gt. 3.0D0) then
                                                v = v - 3.0D0
                                            else
                                                v = 0.0D0
                                            end if
                                            e_scale_chd = 1.0D0 + slope * v
                                        end if

                                        ! Scaled groups: COO-
                                        f_o1 = 0.0D0
                                        f_o2 = 0.0D0
                                        f_cc = 0.0D0

                                        o1 = a_i
                                        o2 = -1
                                        cc = -1
                                        if (geo(a_i)%e .eq. OXYGEN) then
                                            slope = multiplier_coo - 1.0D0

                                            ! Search for closest C atom
                                            cdist = 9.9D9
                                            cv_o1 = 0.0D0
                                            do k = 1, natom
                                                v = cvalence_contribution(geo(o1), geo(k))
                                                cv_o1 = cv_o1 + v ! Sum O1 valence
                                                if (v .gt. 0.0D0 .and. geo(k)%e .eq. CARBON .and. distance(geo(o1), geo(k)) .lt. cdist) then
                                                    cdist =  distance(geo(o1), geo(k))
                                                    cc = k
                                                end if
                                            end do

                                            ! If C found, look for the second O
                                            if (cc .ne. -1) then
                                                odist = 9.9D9
                                                cv_cc = 0.0D0
                                                do k = 1, natom
                                                    v = cvalence_contribution(geo(cc), geo(k))
                                                    cv_cc = cv_cc + v
                                                    if (v .gt. 0.0D0 .and. k .ne. o1 .and. geo(k)%e .eq. OXYGEN .and. distance(geo(cc), geo(k)) .lt. odist) then
                                                        odist =  distance(geo(cc), geo(k))
                                                        o2 = k
                                                    end if
                                                end do
                                            end if


                                            ! O1-C-O2 triad:
                                            if (o2 .ne. -1) then
                                                ! Get O2 valence
                                                cv_o2 = 0.0D0
                                                do k = 1, natom
                                                    cv_o2 = cv_o2 + cvalence_contribution(geo(o2), geo(k))
                                                end do

                                                f_o1 = 1.0D0 - dabs(1.0D0 - cv_o1)
                                                if (f_o1 .lt. 0.0D0) f_o1 = 0.0D0

                                                f_o2 = 1.0D0 - dabs(1.0D0 - cv_o2)
                                                if (f_o2 .lt. 0.0D0) f_o2 = 0.0D0

                                                f_cc = 1.0D0 - dabs(3.0D0 - cv_cc)
                                                if (f_cc .lt. 0.0D0) f_cc = 0.0D0

                                                e_scale_cha = 1.0D0 + slope * f_o1 * f_o2 * f_cc
                                            end if
                                        end if

                                        ! Final energy
                                        e_corr = e_para * e_radial * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha
                                        e_corr_sum = e_corr_sum + e_corr


                                        ! Total gradient
                                        ! radial
                                        call coord_add(grad, d_i, coord_scale(d_radial_d, e_para * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha))
                                        call coord_add(grad, a_i, coord_scale(d_radial_a, e_para * e_angular * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha))
                                        ! angular
                                        call coord_add(grad, d_i, coord_scale(d_angular_d, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha))
                                        call coord_add(grad, a_i, coord_scale(d_angular_a, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha))
                                        call coord_add(grad, h_i, coord_scale(d_angular_h, e_para * e_radial * e_bond_switch * e_scale_w * e_scale_chd * e_scale_cha))
                                        ! bond_switch
                                        call coord_add(grad, d_i, coord_scale(d_bs_d, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha))
                                        call coord_add(grad, a_i, coord_scale(d_bs_a, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha))
                                        call coord_add(grad, h_i, coord_scale(d_bs_h, e_para * e_radial * e_angular * e_scale_w * e_scale_chd * e_scale_cha))
                                        ! water scaling
                                        if (do_grad .and. e_scale_w .ne. 1.0D0) then
                                            slope = multiplier_wh_o - 1.0D0
                                            do k = 1, natom
                                                if (k .ne. d_i) then
                                                    x = distance(geo(d_i), geo(k))
                                                    if (geo(k)%e .eq. HYDROGEN) then
                                                        xd = cvalence_contribution_d(geo(d_i), geo(k)) * sign_wat
                                                        g%x = (geo(d_i)%x - geo(k)%x) * -xd/x * slope
                                                        g%y = (geo(d_i)%y - geo(k)%y) * -xd/x * slope
                                                        g%z = (geo(d_i)%z - geo(k)%z) * -xd/x * slope
                                                        call coord_add(grad, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha))
                                                        call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha))
                                                    else
                                                        xd = cvalence_contribution_d(geo(d_i), geo(k))
                                                        g%x = (geo(d_i)%x - geo(k)%x) * xd/x * slope
                                                        g%y = (geo(d_i)%y - geo(k)%y) * xd/x * slope
                                                        g%z = (geo(d_i)%z - geo(k)%z) * xd/x * slope
                                                        call coord_add(grad, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha))
                                                        call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_cha))
                                                    end if
                                                end if
                                            end do
                                        end if
                                        ! scaled groups: NR4+
                                        if (do_grad .and. e_scale_chd .ne. 1.0D0) then
                                            slope = multiplier_nh4 - 1.0D0
                                            do k = 1, natom
                                                if (k .ne. d_i) then
                                                    x = distance(geo(d_i), geo(k))
                                                    xd = cvalence_contribution_d(geo(d_i), geo(k))
                                                    g%x = (geo(d_i)%x - geo(k)%x) * -xd/x * slope
                                                    g%y = (geo(d_i)%y - geo(k)%y) * -xd/x * slope
                                                    g%z = (geo(d_i)%z - geo(k)%z) * -xd/x * slope
                                                    call coord_add(grad, d_i, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w))
                                                    call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_cha * e_scale_w))
                                                end if
                                            end do
                                        end if
                                        ! scaled groups: COO-
                                        if (do_grad .and. f_o1 * f_o2 * f_cc .ne. 0.0D0) then
                                            slope = multiplier_coo - 1.0D0
                                            ! Atoms around O1
                                            do k = 1, natom
                                                if (k .ne. o1) then
                                                    xd = cvalence_contribution_d(geo(o1), geo(k))
                                                    if (xd .ne. 0.0D0) then
                                                        x =  distance(geo(o1), geo(k))
                                                        if (cv_o1 .gt. 1.0D0) xd = xd * -1.0D0
                                                        xd = xd * f_o2 * f_cc
                                                        g%x = (geo(o1).x - geo(k)%x) * -xd/x * slope
                                                        g%y = (geo(o1).y - geo(k)%y) * -xd/x * slope
                                                        g%z = (geo(o1).z - geo(k)%z) * -xd/x * slope
                                                        call coord_add(grad, o1, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                        call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                    end if
                                                end if
                                            end do
                                            slope = multiplier_coo - 1.0D0
                                            ! Atoms around O2
                                            do k = 1, natom
                                                if (k .ne. o2) then
                                                    xd = cvalence_contribution_d(geo(o2), geo(k))
                                                    if (xd .ne. 0.0D0) then
                                                        x =  distance(geo(o2), geo(k))
                                                        if (cv_o2 .gt. 1.0D0) xd = xd * -1.0D0
                                                        xd = xd * f_o1 * f_cc
                                                        g%x = (geo(o2).x - geo(k)%x) * -xd/x * slope
                                                        g%y = (geo(o2).y - geo(k)%y) * -xd/x * slope
                                                        g%z = (geo(o2).z - geo(k)%z) * -xd/x * slope
                                                        call coord_add(grad, o2, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                        call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                    end if
                                                end if
                                            end do
                                            slope = multiplier_coo - 1.0D0
                                            do k = 1, natom
                                                if (k .ne. cc) then
                                                    xd = cvalence_contribution_d(geo(cc), geo(k))
                                                    if (xd .ne. 0.0D0) then
                                                        x =  distance(geo(cc), geo(k))
                                                        if (cv_cc .gt. 3.0D0) xd = xd * -1.0D0
                                                        xd = xd * f_o1 * f_o2
                                                        g%x = (geo(cc)%x - geo(k)%x) * -xd/x * slope
                                                        g%y = (geo(cc)%y - geo(k)%y) * -xd/x * slope
                                                        g%z = (geo(cc)%z - geo(k)%z) * -xd/x * slope
                                                        call coord_add(grad, cc, coord_scale(g, -e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                        call coord_add(grad, k, coord_scale(g, e_para * e_radial * e_angular * e_bond_switch * e_scale_chd * e_scale_w))
                                                    end if
                                                end if
                                            end do
                                        end if
                                    end if
                                end if
                            end do
                        end if
                    end if
                end do
            end if
        end do
        energy_corr_h4 = e_corr_sum
        return
    end function energy_corr_h4

end module my_types


subroutine HCORR(do_grad)
    use my_types
    use LIMIT, only: LM1, LM1M
    implicit none

    !=======================================
    ! Main
    !=======================================

     COMMON /INOPT2/ IN2(300)
     COMMON /ATOMC / COORD(3,LM1)
     COMMON /ATOMS / NUMAT,NAT(LM1)
     COMMON /CGRAD / CG(3,LM1+LM1M)
     COMMON /CHCORR/ EHCORR

    INTEGER NUMAT,NAT,IN2,IOP,IHCORR
    REAL*8 COORD,CG,EHCORR

    real*8 :: energy_h4
    integer :: i
    type(atom_t),dimension(:),allocatable :: geometry
    type(coord_t),dimension(:),allocatable :: gradient


    LOGICAL do_grad

    IOP=IN2(2)
    IHCORR=IN2(226)

    allocate(gradient(NUMAT),geometry(NUMAT))
    do i=1, NUMAT
        gradient(i)%x = 0.0D0
        gradient(i)%y = 0.0D0
        gradient(i)%z = 0.0D0
    end do
    do i=1, NUMAT
        geometry(i)%x = COORD(1,i)
        geometry(i)%y = COORD(2,i)
        geometry(i)%z = COORD(3,i)
        geometry(i)%e = NAT(i)
    end do

    ! H4 correction
    energy_h4 = energy_corr_h4(NUMAT, geometry, gradient, do_grad)

    IF (IN2(72) .ge. -1 .and. .not. do_grad) THEN
        WRITE(*,*) ""
        WRITE(*,'(A,f11.4)') 'H-Corr        E (kcal/mol)= ',energy_h4
    END IF
    EHCORR = energy_h4


    if (do_grad) then
    ! Total gradient - add and print
        IF (IN2(41) .ge. 1) THEN
            WRITE(*,*) ""
            WRITE(*,'(A)') "     H bonding correction to gradient"
            WRITE(*,'(a5,3(a15))') "NI","X","Y","Z"
            DO i=1,NUMAT
                WRITE(*,'(i5,3(f15.8))') i, gradient(i)%x, gradient(i)%y, gradient(i)%z
            ENDDO
        END IF
        DO i=1,NUMAT
            cg(1,i)=cg(1,i) + gradient(i)%x
            cg(2,i)=cg(2,i) + gradient(i)%y
            cg(3,i)=cg(3,i) + gradient(i)%z
        ENDDO
    end if

end subroutine
!=======================================
!
! Copyright (c) 2011 Jan Rezac
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
! IN THE SOFTWARE.
!=======================================




