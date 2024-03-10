module hdlclib
  use global
  use matrixlib
  use primitive
  use constraint
  implicit none

!------------------------------------------------------------------------------
! Type definitions for HDLC objects and initialisation
! Dimensions of integer matrices: (value, entry)
! Dimensions of Cartesians matrix: (component, atom), unlike old HDLCopt
!
! Fields of an HDLC object:
!
! name                    identification as in resn()
! natom                   number of atoms
! start                   atom sequence number of the first atom
! lgmatok                 .true. if G not linear dependent
! x       (natom)         X coordinates
! y       (natom)         Y coordinates
! z       (natom)         Z coordinates
! ut      matrix(ndlc,np) U transposed
! np                      number of primitives
! nconn                   number of connections
! nbend                   number of bends, linears, impropers
! nrots                   number of torsions
! ncons                   number of constraints
! iconn   (2,nconn)       connectivity
! ibend   (4,nbend)       bends, linears, impropers
! irots   (4,nrots)       torsions
! icons   (6,ncons)       constraints
! vcons   (ncons)         values of constraints
! oldxyz  matrix(3*natom,1) old xyz values (used for cartesian fitting)
! biasv                   bias vector, currently not used
! err_cnt                 counter of HDLC breakdowns
! next    hdlc_obj        next HDLC object
! prev    hdlc_obj        previous HDLC object
!------------------------------------------------------------------------------

  type hdlc_obj
     integer name, natom, start
     logical lgmatok
     real (kind=8), pointer, dimension(:) :: x, y, z
     type(matrix), pointer :: ut
     integer np, nconn, nbend, nrots, ncons
     integer, pointer, dimension(:,:) :: iconn
     integer, pointer, dimension(:,:) :: ibend
     integer, pointer, dimension(:,:) :: irots
     integer, pointer, dimension(:,:) :: icons
     real (kind=8), pointer, dimension(:) :: vcons
     type(matrix), pointer :: oldxyz
     real (kind=8), pointer, dimension(:) :: biasv
     integer err_cnt
     type(hdlc_obj), pointer :: next, prev
  end type hdlc_obj

!------------------------------------------------------------------------------
! Module variables
! - checkpointing of hdlc and all residues: hdlc_wr_hdlc / hdlc_rd_hdlc
!
! lhdlc    .true. if there are internal (HDLC or DLC) coordinates
! internal .true. if 3N-6 DLC instead of 3N HDLC (internals only)
! ngroups  number of HDLC residues
! first    pointer to the first HDLC object
! last     pointer to the last HDLC object
!------------------------------------------------------------------------------

  type hdlc_ctrl
     logical lhdlc, internal
     integer ngroups
     type(hdlc_obj), pointer :: first, last
  end type hdlc_ctrl

  type(hdlc_ctrl) hdlc

contains

!==============================================================================
! routines formerly known as method routines of HDLC objects
!==============================================================================

!------------------------------------------------------------------------------
! subroutine hdlc_create
!
! method: create
!
! Creates a new HDLC residue from cartesian coordinates and (optionally)
! connectivity and constraints.
!------------------------------------------------------------------------------

  subroutine hdlc_create (residue, xyz, con, cns, name, start)

! args
    integer name, start
    type(hdlc_obj), pointer :: residue
    type(matrix), pointer :: xyz, cns
    type(int_matrix), pointer :: con

! local vars
    logical failed
    integer idum, n6, ni
    type(matrix), pointer, save :: bprim, bhdlc, ighdlc

! begin, do linked list stuff
    if (ctrl%printl.ge.1) then
       write (stdout,'(3x,a,i4)') 'Generating HDLC for residue ', name
    end if
    allocate (residue)
    if (hdlc%ngroups .eq. 0) then
       hdlc%first => residue
       hdlc%last => residue
       nullify (residue%prev)
    else
       hdlc%last%next => residue
       residue%prev => hdlc%last
       hdlc%last => residue
    end if
    nullify (residue%next)
    hdlc%ngroups = hdlc%ngroups + 1

! clear pointers of the new residue and other values
    nullify (residue%ibend); nullify (residue%irots); nullify (residue%icons)
    nullify (residue%vcons); nullify (residue%oldxyz); nullify (residue%biasv)
    nullify (residue%x); nullify (residue%y); nullify (residue%z)
    nullify (residue%ut)
    residue%ncons = 0; residue%nconn = 0; residue%nbend = 0; residue%nrots = 0

! store residue information
    residue%name = name
    residue%start = start

! store the coordinates as components
    residue%natom = matrix_dimension (xyz, 1)/3
    allocate (residue%x(residue%natom))
    allocate (residue%y(residue%natom))
    allocate (residue%z(residue%natom))
    call hdlc_linear_checkin (xyz, residue%x, residue%y, residue%z)

! if there are connections - fill in connections
    if (associated(con)) then
       residue%nconn = int_matrix_dimension (con, 2)
       allocate (residue%iconn(2,residue%nconn))
       call hdlc_con_checkin (con, residue%iconn)

! fill in connections required for constraints
       if (associated(cns)) then
          residue%ncons = matrix_dimension (cns, 2)
          allocate (residue%icons(6,residue%ncons))
          allocate (residue%vcons(residue%ncons))
          call ci_cons (cns, residue%ncons, residue%icons, residue%vcons, &
               residue%nconn, residue%iconn)
       end if

! generate angular primitives from connectivity
       call valcoor (residue%nconn, residue%iconn, residue%natom, residue%x, &
            residue%y, residue%z, residue%nbend, residue%ibend, &
            residue%nrots, residue%irots)

! if there are no connections yet - apply total connection scheme
    else
       nullify (residue%iconn)
       call connect_all (residue%natom, residue%nconn, residue%iconn)
    end if ! (associated(con))

! done connections, now do angular constraints
    if (residue%ncons.gt.0) then
       call ck_angle_cons (residue%ncons, residue%icons, residue%nconn, &
            residue%ibend, residue%nbend, residue%irots, residue%nrots)
    end if

! print out primitives - ugly patch for SGI machines
    if (.not.associated(residue%ibend)) allocate (residue%ibend(4,0))
    if (.not.associated(residue%irots)) allocate (residue%irots(4,0))
    call valcoor_print (residue%nconn, residue%iconn, residue%nbend, &
         residue%ibend, residue%nrots, residue%irots, residue%x, &
         residue%y, residue%z, residue%natom)

! generate the primitive B matrix
    call hdlc_make_bprim (residue%natom, residue%x, residue%y, residue%z, &
         residue%nconn, residue%iconn, residue%nbend, residue%ibend, &
         residue%nrots, residue%irots, residue%np, bprim, ni)

! generate the Ut matrix
    call hdlc_make_ut (residue%ut, bprim)

! orthogonalise against constrained space
    if (residue%ncons.gt.0) then
       call ortho_cons (residue%ut, residue%ncons, residue%icons)
    end if

! generate HDLC B matrix
    call hdlc_make_bhdlc (bprim, bhdlc, residue%ut)

! generate HDLC inverse G matrix - only as check - set residue%lgmatok
    call hdlc_make_ighdlc (bhdlc, ighdlc, failed)
    call hdlc_report_failure (residue, failed, 'create')

! here would be the initialisation of the internal biasing vector
    allocate (residue%biasv(residue%np))
    call vector_initialise (residue%biasv, residue%np, 0.0_8)

! generate matrix holding the Cartesians of the previous step
    residue%oldxyz => matrix_create (3*residue%natom, 1, 'old XYZ')
    idum = matrix_copy (xyz, residue%oldxyz)

! initialise cyclic failure counter
    residue%err_cnt = 0

! clean up
    idum = matrix_destroy (bprim)
    idum = matrix_destroy (bhdlc)
    idum = matrix_destroy (ighdlc)

! end
  end subroutine hdlc_create

!------------------------------------------------------------------------------
! subroutine hdlc_destroy
!
! method: destry
!
! Destroys a HDLC residue
!------------------------------------------------------------------------------

  subroutine hdlc_destroy (residue, first, last)

! args
    type(hdlc_obj), pointer :: residue, first, last

! local vars
    integer idum
    type(hdlc_obj), pointer :: next, prev

! begin
    if (ctrl%printl.ge.2) then
       write (stdout,'(a,i4,/)') 'Destroying HDLC residue ', residue%name
    end if
    if (associated (residue%next)) then
       if (associated (residue%prev)) then  ! not first, not last
          next => residue%next
          prev => residue%prev
          prev%next => next
          next%prev => prev
       else                                 ! first, not last
          next => residue%next
          first => next
          nullify (next%prev)
       end if
    else                                    ! not first, last
       if (associated (residue%prev)) then
          prev => residue%prev
          last => prev
          nullify (prev%next)
       else                                 ! first, last
          nullify (first)
          nullify (last)
       end if
    end if
    idum = matrix_destroy (residue%ut)
    deallocate (residue%x); deallocate (residue%y); deallocate (residue%z)
    deallocate (residue%iconn)
    if (associated (residue%ibend)) deallocate (residue%ibend)
    if (associated (residue%irots)) deallocate (residue%irots)
    if (associated (residue%icons)) then
       deallocate (residue%icons); deallocate (residue%vcons)
    end if
    deallocate (residue%oldxyz); deallocate (residue%biasv)
    deallocate (residue)
    hdlc%ngroups = hdlc%ngroups - 1

! end
  end subroutine hdlc_destroy

!------------------------------------------------------------------------------
! subroutine coord_cart_to_hdlc
!
! method: getint
!
! Convert Cartesian coordinates to HDLC
!
! Note: if this subroutine is called by coord_hdlc_to_cart,
! - set lback to .true.
! - prim is preallocated on call, prim is set with primitives on return
! - the Cartesians are not stored to res%oldxyz
!------------------------------------------------------------------------------

  subroutine coord_cart_to_hdlc (res, xyz, chdlc, prim, lback)

! args
    logical lback
    real(kind=8), dimension(:), pointer :: prim
    type(hdlc_obj), pointer :: res
    type(matrix), pointer :: xyz, chdlc

! local params
    real(kind=8) one
    parameter (one = 1.0_8)

! local vars
    integer i, idum, j, n6, nip
    real(kind=8) a1, a2, dx, dy, dz, fact, r
    real(kind=8), dimension(:), allocatable :: x, y, z, p, iut

! begin, the following exception should never occur
    if (res%natom .ne. matrix_dimension(xyz,1)/3) then
       write (stdout,'(A,I4,A,I4,A,I4)') 'Residue ', res%name, &
            ', natom: ', res%natom, '; coordinates, natom: ', &
            matrix_dimension (xyz,1)/3
       call hdlc_errflag ('Size mismatch', 'abort')
    end if

! allocate temporary space and separate Cartesian components
    allocate (x(res%natom))
    allocate (y(res%natom))
    allocate (z(res%natom))
    allocate (p(3*res%natom))
    call hdlc_linear_checkin (xyz, x, y, z)

! init memory for primitive internals and rows of UT
    if (.not. lback) allocate (prim(res%np))
    allocate (iut(res%np))
    nip = 0
    if (hdlc%internal) then
       n6 = 3*res%natom - 6
    else
       n6 = 3*res%natom
    end if

! set up scaling factor for Cartesians
    if (ctrl%cfact .eq. 0.0_8) then
       fact = 1.0_8 / real (res%natom, 8)
    else if (ctrl%cfact .lt. 0.0_8) then
       fact = -ctrl%cfact / real (res%natom, 8)
    else
       fact = ctrl%cfact
    end if

! compute all stretches
    do i = 1,res%nconn
       nip = nip+1
       dx = x(res%iconn(1,i)) - x(res%iconn(2,i))
       dy = y(res%iconn(1,i)) - y(res%iconn(2,i))
       dz = z(res%iconn(1,i)) - z(res%iconn(2,i))
       r = sqrt (dx*dx + dy*dy + dz*dz)
       prim(nip) = r
    end do

! compute all bends
    do i = 1,res%nbend
       if (res%ibend(1,i) .gt. 0) then 
          nip = nip+1

! bends: simple bend
          if (res%ibend(4,i) .eq. 0) then
             prim(nip) = vangle ( &
                  x(res%ibend(1,i)), x(res%ibend(2,i)), x(res%ibend(3,i)), &
                  y(res%ibend(1,i)), y(res%ibend(2,i)), y(res%ibend(3,i)), &
                  z(res%ibend(1,i)), z(res%ibend(2,i)), z(res%ibend(3,i)))

! bends: dihedrals
          else if (res%ibend(4,i) .gt. 0) then
             prim(nip) = vdihedral ( &
                  x(res%ibend(1,i)), y(res%ibend(1,i)), z(res%ibend(1,i)), &
                  x(res%ibend(2,i)), y(res%ibend(2,i)), z(res%ibend(2,i)), &
                  x(res%ibend(3,i)), y(res%ibend(3,i)), z(res%ibend(3,i)), &
                  x(res%ibend(4,i)), y(res%ibend(4,i)), z(res%ibend(4,i)))
             prim(nip) = prim(nip) + res%biasv(nip)

! bends: linear bends (two bends in each of two planes, see hdlc_make_bprim)
          else if (res%ibend(4,i) .le. -1) then
             a1 = vangle ( &
                  x(res%ibend(1,i)),x(res%ibend(2,i)),x(res%ibend(2,i)), &
                  y(res%ibend(1,i)),y(res%ibend(2,i)),y(res%ibend(2,i))+one, &
                  z(res%ibend(1,i)),z(res%ibend(2,i)),z(res%ibend(2,i)))
             a2 = vangle ( &
                  x(res%ibend(2,i)),    x(res%ibend(2,i)),x(res%ibend(3,i)), &
                  y(res%ibend(2,i))+one,y(res%ibend(2,i)),y(res%ibend(3,i)), &
                  z(res%ibend(2,i)),    z(res%ibend(2,i)),z(res%ibend(3,i)))
!            prim(nip) = a1 + a2
             prim(nip) = y(res%ibend(2,i))
             nip = nip + 1
             if (res%ibend(4,i) .eq. -1) then
                a1 = vangle ( &
                   x(res%ibend(1,i)),x(res%ibend(2,i)),x(res%ibend(2,i)), &
                   y(res%ibend(1,i)),y(res%ibend(2,i)),y(res%ibend(2,i)), &
                   z(res%ibend(1,i)),z(res%ibend(2,i)),z(res%ibend(2,i))+one)
                a2 = vangle ( &
                   x(res%ibend(2,i)),    x(res%ibend(2,i)),x(res%ibend(3,i)), &
                   y(res%ibend(2,i)),    y(res%ibend(2,i)),y(res%ibend(3,i)), &
                   z(res%ibend(2,i))+one,z(res%ibend(2,i)),z(res%ibend(3,i)))
             else
                a1 = vangle ( &
                   x(res%ibend(1,i)),x(res%ibend(2,i)),x(res%ibend(2,i))+one, &
                   y(res%ibend(1,i)),y(res%ibend(2,i)),y(res%ibend(2,i)), &
                   z(res%ibend(1,i)),z(res%ibend(2,i)),z(res%ibend(2,i)))
                a2 = vangle ( &
                   x(res%ibend(2,i))+one,x(res%ibend(2,i)),x(res%ibend(3,i)), &
                   y(res%ibend(2,i)),    y(res%ibend(2,i)),y(res%ibend(3,i)), &
                   z(res%ibend(2,i)),    z(res%ibend(2,i)),z(res%ibend(3,i)))
             end if
!            prim(nip) = a1 + a2
             if (res%ibend(4,i) .eq. -1) then
                prim(nip) = z(res%ibend(2,i))
             else
                prim(nip) = x(res%ibend(2,i))
             end if
          end if ! (res%ibend(4,i) .eq. 0) ... else if ...
       end if ! (res%ibend(4,i) .eq. 0)
    end do ! (i = 1,res%nbend)

! dihedrals
    do i = 1,res%nrots
       nip = nip + 1
       prim(nip) = vdihedral ( &
            x(res%irots(1,i)), y(res%irots(1,i)), z(res%irots(1,i)), &
            x(res%irots(2,i)), y(res%irots(2,i)), z(res%irots(2,i)), &
            x(res%irots(3,i)), y(res%irots(3,i)), z(res%irots(3,i)), &
            x(res%irots(4,i)), y(res%irots(4,i)), z(res%irots(4,i)))
       prim(nip) = prim(nip) + res%biasv(nip)
    end do

! Cartesians if required
    if (.not. hdlc%internal) then
       do i = 1,res%natom
          nip = nip + 1
          prim(nip) = x(i)*fact
          nip = nip + 1
          prim(nip) = y(i)*fact
          nip = nip + 1
          prim(nip) = z(i)*fact
       end do
    end if

! if nip is not res%np, something went really wrong
    if (nip .ne. res%np) then
       write (stdout,'(/,A,I4,A,I4)') 'Error, nip: ', nip, ', np: ', res%np
       call hdlc_errflag ('Error converting Cartesians to HDLC', 'stop')
    end if

!//////////////////////////////////////////////////////////////////////////////
! the primitives are available in prim(1..nip) at this point
!//////////////////////////////////////////////////////////////////////////////

! to find HDLC, multiply each row of UT by the primitive internals and sum up
    do i = 1,n6
       p(i) = 0.0_8
       idum = matrix_get_row (res%ut, iut, i)
       do j = 1,res%np
          p(i) = p(i) + iut(j)*prim(j) 
       end do
!      write (stdout,'(a,i4,a,f20.14)') 'HDLC coordinate ', i, ': ', p(i)
    end do

! if in internals, set last 6 to zero
    if (hdlc%internal) then
       do i = n6+1,n6+6
          p(i) = 0.0_8
       end do
    end if

! free memory for primitive internals and rows of UT
    deallocate (iut)
    if (.not. lback) deallocate (prim)

! store HDLC to the returned matrix and hold Cartesians in the residue
    idum = matrix_set (chdlc, p)
    if (.not. lback) idum = matrix_copy (xyz, res%oldxyz)

! clean up
    deallocate (p)
    deallocate (z)
    deallocate (y)
    deallocate (x)

  end subroutine coord_cart_to_hdlc

!------------------------------------------------------------------------------
! subroutine grad_cart_to_hdlc
!
! method: intgrd
!
! Convert Cartesian gradient to HDLC
!------------------------------------------------------------------------------

  subroutine grad_cart_to_hdlc (res, xyz, gxyz, ghdlc)

! args
    type(hdlc_obj), pointer :: res
    type(matrix), pointer :: xyz, gxyz, ghdlc

! local vars
    logical failed
    integer i, idum, n, n6, ni
    real(kind=8), dimension(:), allocatable :: x, y, z
    real(kind=8), dimension(:,:), pointer :: ghdlc_dat
    type(matrix), pointer :: bprim, bhdlc, bthdlc, ighdlc
    save bprim, bhdlc, ighdlc

! begin, the following exception should never occur
    if (res%natom .ne. matrix_dimension(gxyz,1)/3) then
       write (stdout,'(A,I4,A,I4,A,I4)') 'Residue ', res%name, &
            ', natom: ', res%natom, '; coordinates, natom: ', &
            matrix_dimension (gxyz, 1)/3
       call hdlc_errflag ('Size mismatch', 'abort')
    end if

! allocate temporary space and separate Cartesian components
    allocate (x(res%natom))
    allocate (y(res%natom))
    allocate (z(res%natom))
    call hdlc_linear_checkin (xyz, x, y, z)

! number of HDLC internals: n6, number of Cartesians: n 
    n = res%natom*3
    if (hdlc%internal) then
       n6 = n - 6
    else
       n6 = n
    end if

! generate a primitive B matrix (force generation of a new matrix)
    if (associated(bprim)) idum = matrix_destroy (bprim)
    call hdlc_make_bprim (res%natom, x, y, z, res%nconn, res%iconn, &
         res%nbend, res%ibend, res%nrots, res%irots, res%np, bprim, ni)

! generate delocalises B matrix (force again)
    if (associated(bhdlc)) idum = matrix_destroy (bhdlc)
    call hdlc_make_bhdlc (bprim, bhdlc, res%ut)

! generate HDLC inverse G matrix
    if (associated(ighdlc)) idum = matrix_destroy (ighdlc)
    call hdlc_make_ighdlc (bhdlc, ighdlc, failed)
    call hdlc_report_failure (res, failed, 'intgrd')

! make HDLC Bt**-1
    bthdlc => matrix_create (n6, n, "HDLC_BT")
    idum = matrix_multiply (1.0_8, ighdlc, bhdlc, 0.0_8, bthdlc)

! set HDLC gradient = (Bt**-1) * Cartesian gradient
    idum = matrix_multiply (1.0_8, bthdlc, gxyz, 0.0_8, ghdlc) 

! clean up
    idum = matrix_destroy (bprim)
    idum = matrix_destroy (bhdlc)
    idum = matrix_destroy (bthdlc)
    idum = matrix_destroy (ighdlc)
    deallocate (z)
    deallocate (y)
    deallocate (x)
!   call timePrint('END-INTGRAD')

  end subroutine grad_cart_to_hdlc

!------------------------------------------------------------------------------
! subroutine coord_hdlc_to_cart
!
! method: getcrt
!
! Convert HDLC to Cartesian coordinates
!
! Arrays/matrices of Cartesians:
! xyz:                 Cartesians returned by this routine
! res%oldxyz:          Cartesians stored by coords_cart_to_hdlc
! res%x, res%y, res%z: Cartesians used to generate UT (by hdlc_create)
!------------------------------------------------------------------------------

  subroutine coord_hdlc_to_cart (res, xyz, chdlc)

! args
    type(hdlc_obj), pointer :: res
    type(matrix), pointer :: xyz, chdlc

! local params
    real(kind=8) pi, poor
    parameter (pi = 3.1415926535897931_8)
    parameter (poor = 1.0E-6_8)

! local vars
    logical failed, flipped
    integer failc, i, idum, iter, j, k, n, n6, ni
    real(kind=8) absm, cexit, dx, dy, mag, trust
    real(kind=8), dimension(:), allocatable :: tmp_dif, x, y, z
    real(kind=8), dimension(:), pointer :: prim, oldprim, tmp_prim
    real(kind=8), dimension(:,:), pointer :: hdlc_data, olhdlc_data
    type(matrix), pointer :: dif, xdif, bthdlc, xyzbak, olhdlc, &
         bhdlc, ighdlc, bprim, scrdif
    save bhdlc, ighdlc, bprim

! begin, the following exception should never occur
    if (res%natom .ne. matrix_dimension(xyz,1)/3) then
       write (stdout,'(A,I4,A,I4,A,I4)') 'Residue ', res%name, &
            ', natom: ', res%natom, '; coordinates, natom: ', &
            matrix_dimension (xyz,1)/3
       call hdlc_errflag ('Size mismatch', 'abort')
    end if

! this parameter needs to be set upon every call
    cexit = 1.0E-12_8

! say what we are doing
    if (ctrl%printl.ge.5) write (stdout,'(7X,A,/)') &
         'Entering fitting algorithm'
    if (ctrl%printl.ge.4) then
       write (stdout,'(5X,A,E8.1)') &
            'Initial maximum exit error in HDLC is ', cexit
       write (stdout,'(5X,A,E8.1)') &
            'Convergence considered poor if error greater than ', poor
    end if

! n6 is the number of HDLC in the residue
    n = 3*res%natom
    if (hdlc%internal) then
       n6 = n - 6
    else
       n6 = n
    end if

! generate some work matrices
    dif    => matrix_create (n6, 1, 'HDLC diff')
    xdif   => matrix_create (n, 1, 'Cart diff')
    bthdlc => matrix_create (n6, n, 'HDLC_BT')
    xyzbak => matrix_create (n, 1, 'Cart backup')
    olhdlc  => matrix_create (n6, 1, 'old HDLC')

! allocate memory for primitives
    allocate (prim(res%np))
    allocate (oldprim(res%np))

! allocate arrays holding the Cartesian components
    allocate (x(res%natom))
    allocate (y(res%natom))
    allocate (z(res%natom))

! intial guess of Cartesians are the old Cartesians
    idum = matrix_copy (res%oldxyz, xyz)

! backup old Carts in case the iterative B matrix method fails and least square
! SD fitting is used
    idum = matrix_copy (res%oldxyz, xyzbak)

! set currect Carts to old Carts
    idum = matrix_copy (res%oldxyz, xyz)

! the following deallocations to force generation should never be necessary
    if (associated(bhdlc)) idum = matrix_destroy (bhdlc)
    if (associated(ighdlc)) idum = matrix_destroy (ighdlc)
    if (associated(bprim)) idum = matrix_destroy (bprim)

!//////////////////////////////////////////////////////////////////////////////
! Initial trust multiplier
!
! The algorithm makes non-linear fitting steps based on the B matrix.
! Often the step direction is OK, but the magnitude two large.
! The trust multiplier is used to reduce the step size upon step failure.
!//////////////////////////////////////////////////////////////////////////////

    trust = 1.0_8

! set iteration counter
    iter = 0

! set up fail counter, give up if ten consecutive steps fail to improve the fit
    failc = 0

!//////////////////////////////////////////////////////////////////////////////
! Initialise absolute error memory
!
! The absolute error is the measure by which the algorithm tests convergence,
! but the maximum element of the error vector is the exit criterion
!
! dx:   maximum of the absolute of the elements of the error vector
! absm: dx of the previous iteration
!
! dy:   absolute of the error vector (root of scalar product by itself)
! mag:  dy of the previous iteration
!//////////////////////////////////////////////////////////////////////////////

    absm = 0.0_8
    mag = 0.0_8

!//////////////////////////////////////////////////////////////////////////////
! start of main loop
!//////////////////////////////////////////////////////////////////////////////

100 continue
    call hdlc_linear_checkin (xyz, x, y, z)
    k = 0
    flipped = .false.

!//////////////////////////////////////////////////////////////////////////////
! Test if torsion flipped through 180 degrees:
!
! Torsion flip of 360 degrees due to identification of -179 degrees and 181
! degrees can cause failure of the non-linear fit
! trust radius would be > 180 degrees
! 
! The bias vector corrects for this
!//////////////////////////////////////////////////////////////////////////////

    do while (flipped .or. k.eq.0)
       k = k+1
       if (k.gt.2) call hdlc_errflag ( &
            'Code error in hdlclib$hdlc_cart_to_hdlc', 'abort')
       flipped = .false.

! get HDLC from current Cartesians xyz
       call coord_cart_to_hdlc (res, xyz, olhdlc, prim, .true.)
       if (iter.gt.0) then
          if (.not. hdlc%internal) then
             j = res%np - res%natom*3
          else
             j = res%np
          end if
          do j = res%nconn+1,j
             if (prim(j) .gt. oldprim(j)+pi) then 
                res%biasv(j) = res%biasv(j) - pi*2.0_8
                flipped = .true.
             else if (prim(j) .lt. oldprim(j)-pi) then
                res%biasv(j) = res%biasv(j) + pi*2.0_8
                flipped = .true.
             end if
          end do
       end if ! (iter.gt.0)
    end do ! (flipped .or. k.eq.0)
    tmp_prim => oldprim
    oldprim => prim
    prim => tmp_prim

! if pure internals, set last six elements of olhdlc to last six el. of chdlc
    if (hdlc%internal) then
       hdlc_data => matrix_data_array (chdlc)
       olhdlc_data => matrix_data_array (olhdlc)
       do i = n6+1,n6
          olhdlc_data(i,1) = hdlc_data(i,1)
       end do
    end if
      
! generate a primitive B matrix from current cartesians
    call hdlc_make_bprim (res%natom, x, y, z, &
         res%nconn, res%iconn, res%nbend, res%ibend, &
         res%nrots, res%irots, res%np, bprim, ni)

! generate HDLC B matrix
    call hdlc_make_bhdlc (bprim, bhdlc, res%ut)

! generate HDLC inverse G matrix
    call hdlc_make_ighdlc (bhdlc, ighdlc, failed)
    call hdlc_report_failure (res, failed, 'getcrt')

! make HDLC Bt**-1
    idum = matrix_multiply (1.0_8, ighdlc, bhdlc, 0.0_8, bthdlc)

! make [Bt**-1]t
    idum = matrix_transpose (bthdlc)

! get difference between input HDLC and HDLC of the current iteration (olhdlc)
    idum = matrix_copy (olhdlc, dif)
    idum = matrix_scale (dif, -1.0_8)
    idum = matrix_add (dif, dif, chdlc)
    dx = matrix_absmax (dif)
    dy = matrix_length (dif)

! relax limit after four and ten cycles
    if (iter.eq.4 .or. iter.eq.10) then
       cexit = cexit * 1.0e2_8
       if (ctrl%printl.ge.4) write (stdout,'(5x,a,e8.1)') &
            'Relaxing maximum exit error to ', cexit
    end if

! report convergence info if requested
    if (ctrl%printl.ge.4) write (stdout,'(5x,a,i3,a,g9.3,a,g9.3,a,f6.3)') &
         'Step ', iter, ', max error= ', dx, ', length of error= ', dy, &
         ', step scaling factor= ', trust
    call flushout

! set up on first iteration
    if (iter.eq.0) then
       absm = dx + 1.0_8
       mag = dy + 1.0_8
    end if

! test maximum absolute element of difference vector for convergence
    if (dx.lt.cexit) then
       if (mag.lt.dy) then
          idum = matrix_copy (xyzbak, xyz)
          if (ctrl%printl.ge.4) write (stdout,'(5x,a)') 'Recalling best step'
       end if
       goto 200
    end if

!//////////////////////////////////////////////////////////////////////////////
! *** NOT IMPLEMENTED BUT PLANNED ***
!
! if there are 10 unreasonably small steps 
! or a step that increases the error
! go over to a much more stable least squares fit method
!//////////////////////////////////////////////////////////////////////////////

! reject step if error increases
    if (mag.lt.dy) then
       if (ctrl%printl.ge.4) write (stdout,'(/,5x,a,/)') 'Rejecting step'
       trust = trust * 0.5_8
       idum = matrix_copy (xyzbak, xyz)
       if (iter.gt.n6) then 
          if (ctrl%printl.ge.4) write (stdout,'(/,5x,a,/)') &
               'Cannot make step - guessing'
          if (dx.gt.poor) res%lgmatok = .false.
          goto 200
       end if

! error is reduced - too slow convergence detects linear dependence
    else
       if(iter.gt.n6) then 
          if (ctrl%printl.ge.4) write (stdout,'(/,5x,a,/)') &
               'Step converged too slowly'
          if (dx.gt.poor) res%lgmatok = .false.
          goto 200
       end if

! step accepted
       absm = dx
       mag = dy
       failc = 0
       trust = min (1.0_8, trust*1.25_8)

! record best position yet
       idum = matrix_copy (xyz, xyzbak)

! set xyz = oldxyz + [Bt**-1] * trust * [HDLC-HDLC(xyz(i))]
! if internal, create a scratch copy of dif with the last six elements missing
       if (hdlc%internal) then
          scrdif => matrix_create (n6, 1, 'scrdif')
          allocate (tmp_dif(n))
          idum = matrix_get (dif, tmp_dif)
          idum = matrix_set (scrdif, tmp_dif)
          idum = matrix_multiply (1.0_8, bthdlc, scrdif, 0.0_8, xdif)
          idum = matrix_scale (xdif, trust)
          idum = matrix_add (xyz, xdif, res%oldxyz)
          idum = matrix_destroy (scrdif)
          deallocate (tmp_dif)
       else
          idum = matrix_multiply (1.0_8, bthdlc, dif, 0.0_8, xdif)
          idum = matrix_scale (xdif, trust)
          idum = matrix_add (xyz, xdif, res%oldxyz)
       end if

! now xyz contains a guess - recompute the HDLC
    end if ! (mag.lt.dy)
    idum = matrix_copy (xyz, res%oldxyz)
    idum = matrix_transpose (bthdlc)
    iter = iter + 1
    goto 100

!//////////////////////////////////////////////////////////////////////////////
! end of main loop - prepare for return 
!//////////////////////////////////////////////////////////////////////////////

200 continue

! if print level>=3, convergence has already been reported
    if (ctrl%printl.eq.3) write (stdout,'(5x,a,i3,a,/,5x,g9.3,a,/)') &
         'Converged Cartesians in ', iter, ' steps to ', dx, &
         ' maximum component of error vector'

! clean up
    deallocate (x); deallocate (y); deallocate (z)
    nullify (tmp_prim); deallocate (prim); deallocate (oldprim)
    idum = matrix_destroy (dif); idum = matrix_destroy (xdif)
    idum = matrix_destroy (ighdlc)
    idum = matrix_destroy (bthdlc); idum = matrix_destroy (bhdlc)
    idum = matrix_destroy (bprim)
    idum = matrix_destroy (xyzbak)
    idum = matrix_destroy (olhdlc)

  end subroutine coord_hdlc_to_cart

!------------------------------------------------------------------------------
! subroutine hdlc_split_cons
!
! Wrapper to split_cons in library constraint: unwraps variables from residue.
! Separates the values of the constrained degrees of freedom and stores them
! to residue%vcons if required
!
! Arguments:
! residue: residue (in)
! hdlc:    matrix of full HDLC coordinates (in)
!          matrix of active HDLC coordinates (out)
! lstore: constrained HDLC coordinates are stored to vcons if (lstore) (in)
!------------------------------------------------------------------------------

  subroutine hdlc_split_cons (residue, hdlc, lstore)
    logical lstore
    type(hdlc_obj), pointer :: residue
    type(matrix), pointer :: hdlc
    call split_cons (hdlc, lstore, residue%natom, residue%ncons, residue%vcons)
  end subroutine hdlc_split_cons

!------------------------------------------------------------------------------
! subroutine hdlc_rest_cons
!
! Wrapper to rest_cons in library constraint: unwraps variables from residue.
! Restores the values of the constrained degrees of freedom from residue%vcons
!
! Arguments:
! residue: residue (in)
! hdlc:    matrix of active HDLC coordinates (in)
!          matrix of all HDLC coordinates (out)
!------------------------------------------------------------------------------

  subroutine hdlc_rest_cons (residue, hdlc)
    type(hdlc_obj), pointer :: residue
    type(matrix), pointer :: hdlc
    call rest_cons (hdlc, residue%natom, residue%ncons, residue%vcons)
  end subroutine hdlc_rest_cons

!==============================================================================
! routines doing the calculations for the delocalisations
!==============================================================================

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_make_bprim
!
! Makes the B matrix for primitive internal coordinates
!
! On input
! ========
! n  (int)           = number of atoms
! x,y,z (dble)       = cartesian co-rodiantes
! nbend (int)        = number of bends (this is linears and wags as well)
! nrots (int)        = number of rotations
! nconn (int)        = number of connections
! bprim              = pointer to a matrix
!
! if bprim is not associated, a new matrix is made
! if bprim is associated, it is expected to point to a previous matrix
!
! iconn(nconn) (int) = connections array
!                      a,b a-b connection
!
! ibend(nbend) (int) = array from valcoor
!                      a,b,c,0  bend a-b-c
!                      a,b,c,-1 line relative to yz plane
!                      a,b,c,-2 line relative to xy plane
!                      a,b,c,d  dihedral a-b-c-d
!                      (the dihedral is no longer supported by valcoor
!                      but is here for user defined coords)
!
! irots(nrots) (int) = array from valcoor
!                      a,b,c,d  dihedral about b-c
!
! On output
! =========
! ni (int)     = number of primitive internal coordinates
! bprim        = pointer to primitive redundant B matrix
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_make_bprim (natom, x, y, z, nconn, iconn, nbend, ibend, &
       nrots, irots, np, bprim, ni)

! args
    integer natom, nconn, nbend, nrots, np, ni
    integer iconn(2,nconn), ibend(4,nbend), irots(4,nrots)
    real(kind=8), dimension(*) :: x, y, z
    type(matrix), pointer :: bprim

! local params
    integer ib(4), noint, k, idum
    real(kind=8) :: cutoff
    parameter (cutoff = 1.0E-30_8)

! local vars
    logical lbumout
    integer i, j
    real(kind=8), dimension(:,:), allocatable :: brow
    real(kind=8) :: b(3,4), d_zero, d_one, fact, cz(12)

! begin
    data lbumout / .false. /
    data d_zero, d_one / 0.0_8, 1.0_8 /
    if (ctrl%printl.ge.4) write (stdout,'(5X,A)') 'Entering B matrix generator'

! count the number of primitives
    if (hdlc%internal) then
       ni = nconn + nbend + nrots
    else
       ni = 3*natom + nconn + nbend + nrots
    end if

! count linear bends as two primitives
    do i = 1,nbend
       if (ibend(4,i).lt.0 .and. ibend(1,i).ne.0) then
          ni = ni + 1
          if (ctrl%printl.ge.4) write (stdout,'(5X,A,I5,A)') &
               'Linear bend ', i, ' counts as two primitives'
       end if
    end do

! report
    if (ctrl%printl.ge.4) then
       if (hdlc%internal) then
          write (stdout,'(5X,A,I5,/)') 'Number of internal coordinates = ', ni
       else
          write (stdout,'(5X,A,I5)') 'Number of primitives = ', ni
          write (stdout,'(5X,A,I5,A,/)')'Of which are ', 3*natom, ' Cartesians'
       end if
    end if

! allocate memory for one row and the B matrix if required
    allocate (brow(3,natom))
    if (.not. associated(bprim)) then
       if (ctrl%printl.ge.5) write (stdout,'(7X,A)') 'Allocating new B matrix'
       bprim => matrix_create (ni, 3*natom, 'B matrix')
    end if

! adjust the scaling factors for the cartesians (old code: start of ...bprim1)
    if (ctrl%cfact .eq. d_zero) then
       fact = d_one / real (natom, 8)
    elseif (ctrl%cfact.lt.d_zero) then
       fact = -ctrl%cfact / real (natom, 8)
    else
       fact = ctrl%cfact
    end if
    noint = 0

! B matrix is made one row at a time by making brow
    do i = 1,natom
       do j = 1,3
          brow(j,i) = d_zero
       end do
    end do

! loop over all connections making b matrix elements
    do i = 1,nconn
       noint = i
       cz(1) = x(iconn(1,i))
       cz(2) = y(iconn(1,i))
       cz(3) = z(iconn(1,i))
       cz(4) = x(iconn(2,i))
       cz(5) = y(iconn(2,i))
       cz(6) = z(iconn(2,i))
       call str_dlc (1, 1, 2, b, ib, cz)

! construct row of bmatrix
       do j = 1,2
          do k = 1,3
             brow(k,iconn(j,i)) = b(k,j) 
          end do
       end do
       idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
       do j = 1,2
          do k = 1,3
             brow(k,iconn(j,i)) = d_zero
          end do
       end do

! end of loop over connections
    end do

! loop over all bends making bends, wags and linears
    do i = 1,nbend
       if (ibend(1,i).gt.0) then
          noint = noint+1

! 'normal' bend
          if (ibend(4,i).eq.0) then
             cz(1) = x(ibend(1,i))
             cz(2) = y(ibend(1,i))
             cz(3) = z(ibend(1,i))
             cz(4) = x(ibend(2,i))
             cz(5) = y(ibend(2,i))
             cz(6) = z(ibend(2,i))
             cz(7) = x(ibend(3,i))
             cz(8) = y(ibend(3,i))
             cz(9) = z(ibend(3,i))
             call bend_dlc (1, 1, 2, 3, b, ib, cz)

! construct row of bmatrix
             do j = 1,3
                do k = 1,3
                   brow(k,ibend(j,i)) = b(k,j)
                end do
             end do
             idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
             do j = 1,3
                do k = 1,3
                   brow(k,ibend(j,i)) = d_zero
                end do
             end do

! improper dihedral
!
! Impropers can be constructed thus:
!
!        1         1
!        |         |
!        2    =>   2
!       / \         \  <- axis of rotation
!      4   3     4---3
!
! This is not a true wag but it should not make much
! difference in delocalised internals and requires less
! code!
          elseif (ibend(4,i).gt.0) then
             cz(1) = x(ibend(1,i))
             cz(2) = y(ibend(1,i))
             cz(3) = z(ibend(1,i))
             cz(4) = x(ibend(2,i))
             cz(5) = y(ibend(2,i))
             cz(6) = z(ibend(2,i))
             cz(7) = x(ibend(3,i))
             cz(8) = y(ibend(3,i))
             cz(9) = z(ibend(3,i))
             cz(10) = x(ibend(4,i))
             cz(11) = y(ibend(4,i))
             cz(12) = z(ibend(4,i))
             call tors_dlc (1, 1, 2, 3, 4, b, ib, cz)

! construct row of bmatrix
             do j = 1,4
                do k = 1,3
                   brow(k,ibend(j,i)) = b(k,j)
                end do
             end do
             idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
             do j = 1,4
                do k = 1,3
                   brow(k,ibend(j,i)) = d_zero
                end do
             end do

! 'linear' bend w.r.t. yz plane
!
! make a point translated through 1A from atom atom 2 in y, find both angles
! formed and contruct B matrix from contributions from real atoms, then
! repeat with y translation - two coordinates will be formed this way
!
! A   C    A Y    Y C
!  \ /  =   \| +  |/  
!   B        B    B
!
! first set up ib
          elseif (ibend(4,i).le.-1) then
             cz(1) = x(ibend(1,i))
             cz(2) = y(ibend(1,i))
             cz(3) = z(ibend(1,i))
             cz(4) = x(ibend(2,i))
             cz(5) = y(ibend(2,i))
             cz(6) = z(ibend(2,i))
             cz(7) = x(ibend(2,i))
             cz(8) = y(ibend(2,i))+1.0_8
             cz(9) = z(ibend(2,i))
             call bend_dlc (1, 1, 2, 3, b, ib, cz)

! copy atom 1 into real bmat
             do j = 1,3
                brow(j,ibend(1,i)) = b(j,1)
             end do

! perform oposite angle
             cz(1) = x(ibend(2,i))
             cz(2) = y(ibend(2,i))+1.0_8
             cz(3) = z(ibend(2,i))
             cz(7) = x(ibend(3,i))
             cz(8) = y(ibend(3,i))
             cz(9) = z(ibend(3,i))
             call bend_dlc (1, 1, 2, 3, b, ib, cz)

! copy atom 3 into bmat and add atoms 3 & 1 to make 2 
!            do j = 1,3
!               brow(j,ibend(3,i)) = b(j,3) 
!               brow(j,ibend(2,i)) = -1.0_8*(brow(j,ibend(1,i))+b(j,3))
!            end do
             brow(2,ibend(2,i))=1.0_8

! construct row of bmatrix
             idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
             do j = 1,3
                do k = 1,3
                   brow(k,ibend(j,i)) = d_zero
                end do
             end do

! increment ic counter
             noint = noint+1

! do other plane
!            in z if -1 or x if -2
             cz(1) = x(ibend(1,i))
             cz(2) = y(ibend(1,i))
             cz(3) = z(ibend(1,i))
             cz(4) = x(ibend(2,i))
             cz(5) = y(ibend(2,i))
             cz(6) = z(ibend(2,i))
             cz(8) = y(ibend(2,i))      
             if (ibend(4,i).eq.-1) then
                cz(7) = x(ibend(2,i))
                cz(9) = z(ibend(2,i))+1.0_8
             else
                cz(7) = x(ibend(2,i))+1.0_8
                cz(9) = z(ibend(2,i))
             end if
             call bend_dlc (1, 1, 2, 3, b, ib, cz)

! copy atom 1 into real bmatrix
             do j = 1,3
                brow(j,ibend(1,i)) = b(j,1)
             end do

! do opposite angle
             cz(2) = y(ibend(2,i))
             if (ibend(4,i).eq.-1) then
                cz(1) = x(ibend(2,i))
                cz(3) = z(ibend(2,i))+1.0_8
             else
                cz(1) = x(ibend(2,i))+1.0_8
                cz(3) = z(ibend(2,i))
             end if
             cz(7) = x(ibend(3,i))
             cz(8) = y(ibend(3,i))
             cz(9) = z(ibend(3,i))
             call bend_dlc (1, 1, 2, 3, b, ib, cz)

! copy atom 3 into bmatrix and add atoms 3 & 1 to make 2
!            do j = 1,3
!               brow(j,ibend(3,i))=b(j,3)
!               brow(j,ibend(2,i))=brow(j,ibend(1,i))+b(j,3)
!            end do
             if (ibend(4,i).eq.-1) then
                brow(3,ibend(2,i)) = 1.0_8
             else
                brow(1,ibend(2,i)) = 1.0_8
             end if

! construct row of bmatrix
             idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
             do j = 1,3
                do k = 1,3
                   brow(k,ibend(j,i)) = d_zero
                end do
             end do
          end if ! if (ibend(4,i).eq.0) then ... elseif ...
       end if ! if (ibend(1,i).gt.0) then

! end of loop over all bends etc.
    end do ! do i = 1,nbend

! put in dihedrals
    do i = 1,nrots
       noint = noint+1
       cz(1) = x(irots(1,i))
       cz(2) = y(irots(1,i))
       cz(3) = z(irots(1,i))
       cz(4) = x(irots(2,i))
       cz(5) = y(irots(2,i))
       cz(6) = z(irots(2,i))
       cz(7) = x(irots(3,i))
       cz(8) = y(irots(3,i))
       cz(9) = z(irots(3,i))
       cz(10) = x(irots(4,i))
       cz(11) = y(irots(4,i))
       cz(12) = z(irots(4,i))
       call tors_dlc (1, 1, 2, 3, 4, b, ib, cz)

! construct row of bmatrix
       do j = 1,4
          do k = 1,3
             brow(k,irots(j,i)) = b(k,j)
          end do
       end do
       idum = matrix_set_row (bprim, brow, noint)

! re-zero brow
       do j = 1,4
          do k = 1,3
             brow(k,irots(j,i)) = d_zero
          end do
       end do

! end of loop over dihedrals
    end do

! put in cartesians
    if (.not. hdlc%internal) then
       do i = 1,natom
          do j = 1,3
             noint = noint+1
             brow(j,i) = fact
             idum = matrix_set_row (bprim, brow, noint)
             brow(j,i) = 0.0_8
          end do
       end do
    end if

! set the number of primitives and print the matrix if requested
    np = noint
    if (ctrl%printl.ge.5) i = matrix_print (bprim)

!//////////////////////////////////////////////////////////////////////////////
! Helper routines for hdlc_make_bprim
!
! The following three subroutines are adapted from the normal coordinate
! analysis program of Schachtschneider, Shell development
!//////////////////////////////////////////////////////////////////////////////

  contains
    subroutine str_dlc (noint,i,j,b,ib,c)
      integer i, j, ib(4,*), noint
      real(kind=8) c(*), b(3,4,1)
!
      integer iaind, jaind, m
      real(kind=8) dzero, rij(3), dijsq
!
      data dzero/0.0_8/
!
      iaind=0
      jaind=3
      ib(1,noint)=i
      ib(2,noint)=j
      dijsq = dzero
      do m=1,3
         rij(m)=c(m+jaind)-c(m+iaind)
         dijsq=dijsq+rij(m)**2
      end do
      do m=1,3
         b(m,1,noint) = -rij(m)/sqrt(dijsq)
         b(m,2,noint)=-b(m,1,noint)
      end do
      return
!
    end subroutine str_dlc
!
    subroutine bend_dlc (noint,i,j,k,b,ib,c)
      integer i, j, k, ib(4,1), noint
      real(kind=8) b(3,4,1), c(*)
!
      integer iaind, jaind, kaind, m
      real(kind=8) rji(3), rjk(3), eji(3), ejk(3), djisq, djksq, dzero, done, &
           sinj, dotj, dji, djk
!
      data dzero/0.0_8/,done/1.0_8/
!
      iaind=0
      jaind=3
      kaind=6
      ib(1,noint)=i
      ib(2,noint)=j
      ib(3,noint)=k
      djisq=dzero
      djksq=dzero
      do m=1,3
         rji(m)=c(m+iaind)-c(m+jaind)
         rjk(m)=c(m+kaind)-c(m+jaind)
         djisq = djisq + rji(m)**2
         djksq = djksq + rjk(m)**2
      end do
      dji = sqrt(djisq)
      djk = sqrt(djksq)
      dotj = dzero
      do m=1,3
         eji(m)=rji(m)/dji
         ejk(m)=rjk(m)/djk
         dotj=dotj+eji(m)*ejk(m)
      end do
      sinj = sqrt(done-dotj**2)
      if(sinj.lt.1.0e-20_8)go to 145
      do m=1,3
         b(m,3,noint)=((dotj*ejk(m)-eji(m)))/(djk*sinj)
         b(m,1,noint)=((dotj*eji(m)-ejk(m)))/(dji*sinj)
         b(m,2,noint)=-b(m,1,noint)-b(m,3,noint)
      end do
!
      return
 145  continue
      do m=1,3
         b(m,3,noint)=0.0_8
         b(m,1,noint)=0.0_8
         b(m,2,noint)=0.0_8
      end do
      return
!
    end subroutine bend_dlc
!
    subroutine tors_dlc(noint,i,j,k,l,b,ib,c)
      integer i, j, k, l, ib(4,*), noint
      real(kind=8) b(3,4,1), c(*)
!
      integer iaind, jaind, kaind, laind, m
      real(kind=8) dzero, done, dij, dijsq, djk, djksq, dkl, dklsq, dotpj, &
           dotpk, sinpj, sinpk, smi, smj, sml, f1, f2
      real(kind=8) rij(3), rjk(3), rkl(3), eij(3), ejk(3), ekl(3), cr1(3), &
           cr2(3)
!
      data dzero/0.0_8/,done/1.0_8/
!
      iaind=0
      jaind=3
      kaind=6
      laind=9
      ib(1,noint)=i
      ib(2,noint)=j
      ib(3,noint)=k
      ib(4,noint)=l
      dijsq=dzero
      djksq=dzero
      dklsq=dzero
      do m=1,3
         rij(m)=c(m+jaind)-c(m+iaind)
         dijsq=dijsq+rij(m)**2
         rjk(m)=c(m+kaind)-c(m+jaind)
         djksq=djksq+rjk(m)**2
         rkl(m)=c(m+laind)-c(m+kaind)
         dklsq=dklsq+rkl(m)**2
      end do
      dij = sqrt(dijsq)
      djk = sqrt(djksq)
      dkl = sqrt(dklsq)
      do m=1,3
         eij(m)=rij(m)/dij
         ejk(m)=rjk(m)/djk
         ekl(m)=rkl(m)/dkl
      end do
      cr1(1)=eij(2)*ejk(3)-eij(3)*ejk(2)
      cr1(2)=eij(3)*ejk(1)-eij(1)*ejk(3)
      cr1(3)=eij(1)*ejk(2)-eij(2)*ejk(1)
      cr2(1)=ejk(2)*ekl(3)-ejk(3)*ekl(2)
      cr2(2)=ejk(3)*ekl(1)-ejk(1)*ekl(3)
      cr2(3)=ejk(1)*ekl(2)-ejk(2)*ekl(1)
      dotpj = -(eij(1)*ejk(1)+eij(2)*ejk(2)+eij(3)*ejk(3))
      dotpk = -(ejk(1)*ekl(1)+ejk(2)*ekl(2)+ejk(3)*ekl(3))
      sinpj = sqrt(done-dotpj**2)
      sinpk = sqrt(done-dotpk**2)
      do 164 m=1,3
      smi=-cr1(m)/(dij*sinpj*sinpj)
      b(m,1,noint)=smi
      f1=0.0_8
      f2=0.0_8
      if(sinpj.gt.1.0e-20_8) f1=(cr1(m)*(djk-dij*dotpj))/(djk*dij*sinpj*sinpj)
      if(sinpk.gt.1.0e-20_8) f2=(dotpk*cr2(m))/(djk*sinpk*sinpk)
      smj=f1-f2
      b(m,2,noint)=smj
      sml=0.0_8
      if(sinpk.gt.1.0e-20_8)sml= cr2(m)/(dkl*sinpk*sinpk)
      b(m,4,noint)=sml
      b(m,3,noint)=(-smi-smj-sml)
  164 continue
      return
    end subroutine tors_dlc

!//////////////////////////////////////////////////////////////////////////////
! end of helper routines for hdlc_make_bprim
!//////////////////////////////////////////////////////////////////////////////

  end subroutine hdlc_make_bprim

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_make_ut
!
! Constructs Ut, the transpose of the non-redundant eigenvectors of the
! redundant G matrix
!
! Alexander J Turner Dec 1997
!
!     See:
! J_Baker, A_Kessi and B_Delley
! J.Chem.Phys 105,(1),1 July 1996
!
! Input
! =====
!
! bprim:  primitive B matrix
!
! Output
! ======
!
! ut: transposed U matrix
!
! Paramters
! =========
!
! Cutoff decides if an eigenvalue of G is zero
!
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_make_ut (ut, bprim)

! args
    type(matrix), pointer :: ut, bprim

! local params
    real(kind=8) cutoff
    parameter (cutoff = 1.0E-8_8)

! local vars
    integer i, idum, j, n, n6, nii
    real(kind=8), dimension(:), allocatable :: temp
    type(matrix), pointer :: g_mat, btprim, r_mat, v_mat

! begin, nii: number of primitives
    if (ctrl%printl.ge.4) write(stdout,'(5X,A)') 'Entering UT matrix generator'
    nii = matrix_dimension (bprim, 1)

! n: number of cartesians
    n = matrix_dimension (bprim, 2)

! n6: number of delocalised internal coordinates
    if (hdlc%internal) then
       n6 = n-6
    else
       n6 = n
    end if

! test for insufficiant primitives
    if (nii.lt.n6) then
       write (stdout,'(A,I5,A,I5,A)') 'There are ', nii, ' primitives and ', &
            n6,' required!'
       call hdlc_errflag ('Insufficiant primitive coordinates', 'abort')
    end if

! allocate space
    g_mat  => matrix_create (nii, nii, 'G matrix')
    btprim => matrix_create (nii, n,  'Bt matrix')

! form G
    idum = matrix_copy (bprim, btprim)
    idum = matrix_transpose (btprim)
    idum = matrix_multiply (1.0_8, bprim, btprim, 0.0_8, g_mat)

! allocate work matrices for the diagonalisation; Bt is not used any longer
    idum = matrix_destroy (btprim)

! Next two lines change due to problems with LAPACK routine dsyevx
!   r_mat => matrix_create (nii, n6, 'R vectors')
!   v_mat => matrix_create (n6, 1,   'R roots ')
    r_mat => matrix_create (nii, nii, 'R vectors')
    v_mat => matrix_create (nii, 1,   'R roots ')

! diagonalise and form U 
!   call timeSet ()
    idum = matrix_diagonalise (g_mat, r_mat, v_mat, .false.)
    if (idum .ne. 0) then
       call hdlc_errflag ('G matrix could not be diagonalised', 'stop')
    end if
!   call timePrint ('G prim diag')

! allocate Ut matrix; the G matrix not used any longer
    idum = matrix_destroy (g_mat)
    ut  => matrix_create (n6, nii, 'Ut matrix')

! check for faulty eigenvalue structure, use temp as scratch
    allocate (temp(nii))
    idum = matrix_get (v_mat, temp)
    if (ctrl%printl.ge.4) write (stdout,'(5x,a,/)') &
         'Eigenvalue structure of primitive G matrix'
    j = 0
    do i = 1,nii
       if (ctrl%printl.ge.4) write (stdout,'(5x,a,i5,a,f13.8)') &
            'Eigenvalue ', i, ' = ', temp(i)
       if (abs(temp(i)).gt.cutoff) j = j + 1
    end do
    if (ctrl%printl.ge.3) write (stdout,'(/,5x,a,i5,a,i5,a,/)') &
         'The system has ', n6, ' degrees of freedom, and ', j, &
         ' non-zero eigenvalues'
    if (j.lt.n6) then
       if (ctrl%ctfirst.eq.4 .and. hdlc%internal) then
          call hdlc_errflag ( &
          'DLC based on total connection cannot be used for a planar system', &
          'stop')
       else
          call hdlc_errflag ( &
               'Too few delocalised coordinates with non-zero values', 'stop')
       end if
    end if

! generate Ut (the eigenvectors come out of matrix_diagonalise in columns!)
    do i = 1,n6
       j = matrix_get_column (r_mat, temp, i)
       j = matrix_set_row (ut, temp,i)
    end do

! clean up some memory
    idum = matrix_destroy (v_mat)
    idum = matrix_destroy (r_mat)
    deallocate (temp)

! end
    return
  end subroutine hdlc_make_ut

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_make_bhdlc
!
! Construct delocalised B matrix
!
! Input
! =====
! ut    : Ut matrix
! bprim : primitive B matrix
! bhdlc : unassociated => generate B matrix
!         associated   => pre-existing B matrix
!
! Output
! ======
! bhdlc : HDLC B matrix
!
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_make_bhdlc (bprim, bhdlc, ut)

! args
    type(matrix), pointer :: bhdlc, bprim, ut

! local vars
    integer idum, n, np, n6

! begin, get dimensions of B matrix
    np = matrix_dimension (bprim, 1)
    n = matrix_dimension (bprim, 2)
    if (hdlc%internal) then
       n6 = n - 6
    else
       n6 = n
    end if

! create B matrix if needed
    if (.not. associated(bhdlc)) then
       bhdlc => matrix_create (n6, n, 'B HDLC')
       if (ctrl%printl.ge.5) write (stdout,'(7X,A,/)') 'New HDLC B matrix'
    end if

! do maths
    idum = matrix_multiply (1.0_8, ut, bprim, 0.0_8, bhdlc)
    if (ctrl%printl.ge.5) idum = matrix_print (bhdlc)

  end subroutine hdlc_make_bhdlc

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_make_ighdlc
!
! Construct inverse G matrix from delocalised B matrix
!
! Input
! =====
! bhdlc  : HDLC B matrix
! ighdlc : unassociated => generate iG matrix
!          associated   => pre-existing iG matrix
!
! Output
! ======
! ighdlc : inverse HDLC G matrix
! failed : true if det|iG| < cutoff
!
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_make_ighdlc (bhdlc, ighdlc, failed)

! args
    logical failed
    type(matrix), pointer :: bhdlc, ighdlc

! local params
    real(kind=8) :: cutoff
    parameter (cutoff = 1.0E-8_8)

! local vars
    integer i, idum, n, n6
    real(kind=8) :: det
    type(matrix), pointer :: bthdlc

! begin, get number of HDLC (n6) and number of cartesians (n)
    n6 = matrix_dimension (bhdlc, 1)
    n = matrix_dimension (bhdlc, 2)

! generate iG matrix if needed
    if (.not. associated(ighdlc)) then
       ighdlc => matrix_create (n6, n6, 'HDLC iG')
       if (ctrl%printl.ge.5) write (stdout,'(7X,A,/)') 'New HDLC iG matrix'
    end if

! make B transposed
    bthdlc => matrix_create (n6, n, 'HDLC BT')
    idum = matrix_copy (bhdlc, bthdlc)
    idum = matrix_transpose (bthdlc)

! do maths
    idum = matrix_multiply (1.0_8, bhdlc, bthdlc, 0.0_8, ighdlc)
    if (ctrl%printl.ge.5) idum = matrix_print (ighdlc)
    idum = matrix_invert (ighdlc, det, .true.)

! clean up
    idum = matrix_destroy (bthdlc)

! check for linear dependency in iG, scale det to be size-independent
    i = matrix_dimension (ighdlc, 1)
    if (ctrl%printl.ge.4) then
       write (stdout,'(7X,A,E13.5,/)') 'The HDLC G matrix determinant is ', det
    end if
    det = det**(1.0_8/real(i,8))
    failed = (det.lt.cutoff)
  end subroutine hdlc_make_ighdlc

!==============================================================================
! checkpointing
!==============================================================================

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_wr_hdlc
!
! Write all HDLC objects to iunit
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_wr_hdlc (iunit, s_hdlc, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(hdlc_ctrl) :: s_hdlc

! local vars
    integer i, j, igroup
    type(hdlc_obj), pointer :: res

! begin
    lerr = .false.

! write header, return if no HDLC residues
    if (lform) then
       write (iunit,*,err=98) s_hdlc%lhdlc, s_hdlc%internal, s_hdlc%ngroups
    else
       write (iunit,err=98) s_hdlc%lhdlc, s_hdlc%internal, s_hdlc%ngroups
    end if
    if (.not. s_hdlc%lhdlc) return

! loop over all HDLC objects
    res => s_hdlc%first
    do igroup = 1,s_hdlc%ngroups

! write scalars
       if (lform) then
          write (iunit,'(8I8)',err=98) res%name, res%natom, res%start, &
               res%np, res%nconn, res%nbend, res%nrots, res%ncons
          write (iunit,*,err=98) res%err_cnt, res%lgmatok          
       else
          write (iunit,'(8I8)',err=98) res%name, res%natom, res%start, &
               res%np, res%nconn, res%nbend, res%nrots, res%ncons
          write (iunit,*,err=98) res%err_cnt, res%lgmatok
       end if

! write arrays
       if (lform) then
          write (iunit,'(5F16.10)',err=98) (res%x(i), i=1,res%natom)
          write (iunit,'(5F16.10)',err=98) (res%y(i), i=1,res%natom)
          write (iunit,'(5F16.10)',err=98) (res%z(i), i=1,res%natom)
          write (iunit,'(2I8)',err=98) ((res%iconn(i,j), i=1,2), j=1,res%nconn)
          write (iunit,'(4I8)',err=98) ((res%ibend(i,j), i=1,4), j=1,res%nbend)
          write (iunit,'(4I8)',err=98) ((res%irots(i,j), i=1,4), j=1,res%nrots)
          write (iunit,'(6I8)',err=98) ((res%icons(i,j), i=1,6), j=1,res%ncons)
          write (iunit,'(5F16.10)',err=98) (res%vcons(i), i=1,res%ncons)
       else
          write (iunit,err=98) (res%x(i), i=1,res%natom)
          write (iunit,err=98) (res%y(i), i=1,res%natom)
          write (iunit,err=98) (res%z(i), i=1,res%natom)
          write (iunit,err=98) ((res%iconn(i,j), i=1,2), j=1,res%nconn)
          write (iunit,err=98) ((res%ibend(i,j), i=1,4), j=1,res%nbend)
          write (iunit,err=98) ((res%irots(i,j), i=1,4), j=1,res%nrots)
          write (iunit,err=98) ((res%icons(i,j), i=1,6), j=1,res%ncons)
          write (iunit,err=98) (res%vcons(i), i=1,res%ncons)
       end if

! write UT matrix - see matrixlib
       call hdlc_wr_matrix (iunit, res%ut, lform, lerr)

! end loop over all HDLC residues
       res => res%next
    end do

! I/O failure jump point
    return
98  lerr = .true.
    if (ctrl%printl.ge.-1) write (stdout,'(/,a,i5,a,/)') &
         '*** Problem writing residue ', igroup, ' to checkpoint ***'
  end subroutine hdlc_wr_hdlc

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_rd_hdlc
!
! Read all HDLC objects from iunit
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_rd_hdlc (iunit, s_hdlc, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(hdlc_ctrl) :: s_hdlc

! local vars
    integer i, j, igroup
    integer, dimension(1) :: idum
    type(hdlc_obj), pointer :: res

! begin
    lerr = .false.
    igroup = 0

! this should never occur
    if (s_hdlc%ngroups.gt.0) call hdlc_destroy_all (.false., idum)

! read header, return if no HDLC residues
    if (lform) then
       read (iunit,*,err=98) s_hdlc%lhdlc, s_hdlc%internal, s_hdlc%ngroups
    else
       read (iunit,err=98) s_hdlc%lhdlc, s_hdlc%internal, s_hdlc%ngroups
    end if
    if (.not. s_hdlc%lhdlc) return

! loop over all HDLC objects
    do igroup = 1,s_hdlc%ngroups

! allocate new residue and read scalars
       allocate (res)
       if (lform) then
          read (iunit,'(8I8)',err=98) res%name, res%natom, res%start, &
               res%np, res%nconn, res%nbend, res%nrots, res%ncons
          read (iunit,*,err=98) res%err_cnt, res%lgmatok          
       else
          read (iunit,'(8I8)',err=98) res%name, res%natom, res%start, &
               res%np, res%nconn, res%nbend, res%nrots, res%ncons
          read (iunit,*,err=98) res%err_cnt, res%lgmatok
       end if

! allocate contained arrays
       allocate (res%x (res%natom))
       allocate (res%y (res%natom))
       allocate (res%z (res%natom))
       allocate (res%iconn (2,res%nconn))
       allocate (res%ibend (4,res%nbend))
       allocate (res%irots (4,res%nrots))
       allocate (res%icons (6,res%ncons))
       allocate (res%vcons (res%ncons))

! read contained arrays
       if (lform) then
          read (iunit,'(5F16.10)',err=98) (res%x(i), i=1,res%natom)
          read (iunit,'(5F16.10)',err=98) (res%y(i), i=1,res%natom)
          read (iunit,'(5F16.10)',err=98) (res%z(i), i=1,res%natom)
          read (iunit,'(2I8)',err=98) ((res%iconn(i,j), i=1,2), j=1,res%nconn)
          read (iunit,'(4I8)',err=98) ((res%ibend(i,j), i=1,4), j=1,res%nbend)
          read (iunit,'(4I8)',err=98) ((res%irots(i,j), i=1,4), j=1,res%nrots)
          read (iunit,'(6I8)',err=98) ((res%icons(i,j), i=1,6), j=1,res%ncons)
          read (iunit,'(5F16.10)',err=98) (res%vcons(i), i=1,res%ncons)
       else
          read (iunit,err=98) (res%x(i), i=1,res%natom)
          read (iunit,err=98) (res%y(i), i=1,res%natom)
          read (iunit,err=98) (res%z(i), i=1,res%natom)
          read (iunit,err=98) ((res%iconn(i,j), i=1,2), j=1,res%nconn)
          read (iunit,err=98) ((res%ibend(i,j), i=1,4), j=1,res%nbend)
          read (iunit,err=98) ((res%irots(i,j), i=1,4), j=1,res%nrots)
          read (iunit,err=98) ((res%icons(i,j), i=1,6), j=1,res%ncons)
          read (iunit,err=98) (res%vcons(i), i=1,res%ncons)
       end if

! allocate and read UT matrix - see matrixlib
       call hdlc_rd_matrix (iunit, res%ut, lform, lerr)

! linked list stuff
       if (igroup .eq. 1) then
          s_hdlc%first => res
       else 
          s_hdlc%last%next => res
          res%prev => s_hdlc%last
       end if
       s_hdlc%last => res
    end do

! I/O failure jump point - destroy already read objects if necessary
    return
98  lerr = .true.
    s_hdlc%lhdlc = .false.
    if (igroup .gt. 0) call hdlc_destroy_all (.false., idum)
    if (ctrl%printl.ge.-1) write (stdout,'(/,a,i5,a,/)') &
         '*** Problem reading residue ', igroup, ' ***'
  end subroutine hdlc_rd_hdlc

!==============================================================================
! general helper routines follow
!==============================================================================

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_destroy_all
!
! Destroys all HDLC residues and saves the fields residue%err_cnt to the array
! err_cnt if lsave is set
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_destroy_all (lsave, err_cnt)

! args
    logical lsave
    integer, dimension(*) :: err_cnt

! local vars
    integer i
    type(hdlc_obj), pointer :: residue

! begin
    if (ctrl%printl.ge.1) then
       write (stdout, '(a,/)') 'Destroying all HDLC residues'
    end if
    i = 1
    do while (hdlc%ngroups.ne.0)
       residue => hdlc%first
       if (lsave .and. (residue%err_cnt.lt.2000 .or. residue%ncons.gt.0)) then
          err_cnt(i) = residue%err_cnt
          i = i+1
       end if
       call hdlc_destroy (residue, hdlc%first, hdlc%last)
    end do

! end
  end subroutine hdlc_destroy_all

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_linear_checkin
!
! Fill in cartesian components into separate arrays
! On input: mat is a (3,natom) matrix: ((x1,y1,z1),(x2,..),...)
!           or an (3*natom,1) matrix: (x1,y1,z1,x2,...)
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_linear_checkin (mat, x, y, z)

! args
    type(matrix), pointer :: mat
    real(kind=8) x(*), y(*), z(*)

! local vars
    integer i, j, size
    real(kind=8), dimension(:,:), pointer :: xyz_data

! begin
    xyz_data => matrix_data_array (mat)
    if (matrix_dimension(mat,2) .eq. 1) then
       size = matrix_dimension (mat,1)/3
       j = 0
       do i = 1,size
          x(i) = xyz_data(j+1,1)
          y(i) = xyz_data(j+2,1)
          z(i) = xyz_data(j+3,1)
          j = j + 3
       end do
    else
       size = matrix_dimension (mat,2)
       do i = 1,size
          x(i) = xyz_data(1,i)
          y(i) = xyz_data(2,i)
          z(i) = xyz_data(3,i)
       end do
    end if
  end subroutine hdlc_linear_checkin

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_con_checkin
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_con_checkin (imat, iconn)

! args
    type(int_matrix), pointer :: imat
    integer, dimension(2,*) :: iconn

! local vars
    integer i, nconn

! begin
    nconn = int_matrix_dimension (imat,2)
    do i = 1,nconn
       iconn(1,i) = imat%data(1,i)
       iconn(2,i) = imat%data(2,i)
    end do
  end subroutine hdlc_con_checkin

!//////////////////////////////////////////////////////////////////////////////
! subroutine hdlc_report_failure
!//////////////////////////////////////////////////////////////////////////////

  subroutine hdlc_report_failure (residue, failed, location)

! args
    logical failed
    character*(*) location
    type(hdlc_obj) :: residue

! begin
    residue%lgmatok = (.not. failed)
    if (failed) then
       if (ctrl%printl.ge.0) write (stdout,'(/,A,I3,A,A,/)') &
            'G matrix linear dependent in residue ', residue%name, &
            ', routine: ', location
    end if
  end subroutine hdlc_report_failure

!//////////////////////////////////////////////////////////////////////////////
! subroutine vector_initialise
!//////////////////////////////////////////////////////////////////////////////

  subroutine vector_initialise (vector,  number, value)

! args
    integer number
    real(kind=8) value
    real(kind=8), dimension(number) :: vector

! local vars
    integer i

! begin
    do i = 1,number
       vector(i) = value
    end do
  end subroutine vector_initialise

end module hdlclib
