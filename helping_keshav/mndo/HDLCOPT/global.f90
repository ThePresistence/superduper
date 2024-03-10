module global
  implicit none

!------------------------------------------------------------------------------
! Optimisation control data structure
! - reading: hdlc_get_params of the interface library
! - checkpointing: hdlc_wr_ctrl / hdlc_rd_ctrl
!
! mstart   startup procedure:
!          0 -> start up from scratch, no checkpoint used at all
!          1 -> start up using coordinates from checkpoint file
!          2 -> start up using coordinates, resn and ctrl from checkpoint file
!          3 -> continue using the complete optimiser status from checkpoint
!          4 -> as mstart=3, use checkpoint for funct. / gradient code (mdump)
! nbig     number of degrees of freedom of the total system (incl. constraints)
! nvar     number of degrees of freedom of the reaction core
! nrem     number of remembered LBFGS steps
! nrestart periodic restart: number of LBFGS and PRFO steps (0: no restart)
! attypes  nuclear charge of each atom
! toler    tolerance
! rmin     minimum ratio actual / predicted step
! rmax     maximum ratio actual / predicted step
! dmax     initial step size
! ddmax    maximum step size
! osmin    minimum PRFO step size
! omin     minimum overlap criterion
! mode     eigenmode to be followed by PRFO (0: minimise)
! irecalc  recompute Hessian every irecalc PRFO cycles
!          (-1: never, 0: only at init),
! maxcyc   maximum number of cycles (function / gradient evaluations)
! iupd     Hessian update scheme (1: Powell, 2: BFGS)
! mhess    Hessian calculation (0: finite difference, 1: external, if available
! idump    dump the memory every idump cycles (0: never)
! mdump    mode of dump (0: optimiser only, 1: for funct. / gradient code too)
! iprj     project out R/T modes out if no internal coordinates (1: on)
! ilock    mode lock (0: mode switching allowed)
! precon   preconditioning vector for gradient
! contyp   connection type (0: primitives, 1: total connection)
! ctfirst  contyp of the first residue (0: prim, 1: total, 2: cartesian)
! nincon   number of extra connections read in
! incon    extra connections
! mconstr  mode of constraints (0: orthogonalisation)
! nconstr  number of constraints
! iconstr  constraints as specified in constraint.f90, format a
! vconstr  target value of constraints (only if (mconstr.gt.0))
! cfact    scaling factor for cartesians
!          0  -> 1/M as usual
!          <0 -> cfact/M
!          >0 -> cfact
! printl   print level
! nrad     number of covalent radii provided by the user
! irad     nuclear charge of the atoms with user covalent radius
! vrad     user covalent radii
! iccode   calling code (see search)
! icconv   convergence criterion
!          0: step and gradient in the coordinate system(s) of the optimisation
!          1: gradient only in Cartesian coordinates
!------------------------------------------------------------------------------

  type opt_ctrl
     integer mstart, nbig, nvar, nrem
     integer nrestart
     integer, pointer, dimension(:) :: attypes
     real (kind=8) :: toler, rmin, rmax, dmax, ddmax, osmin, omin
     integer mode, irecalc, maxcyc, iupd, mhess, idump, mdump, iprj, ilock
     real (kind=8), pointer, dimension(:) :: precon
     integer nincon, mconstr, nconstr, contyp, ctfirst
     integer, pointer, dimension(:,:) :: incon
     integer, pointer, dimension(:,:) :: iconstr
     real (kind=8), pointer, dimension(:) :: vconstr
     real (kind=8) :: cfact
     integer printl, nrad
     integer, pointer, dimension(:) :: irad
     real(kind=8), pointer, dimension(:) :: vrad
     integer iccode, icconv
  end type opt_ctrl

  type(opt_ctrl) :: ctrl

!------------------------------------------------------------------------------
! Optimiser status 'registers'
! - checkpointing: hdlc_wr_stat / hdlc_rd_stat
!
! ncyc       number of energy / gradient evaluations so far
! nstep      number of tried P-RFO steps
! nsprfo     number of performed P-RFO steps
! nstried    number of tried L-BFGS steps
! nslbfg     number of performed L-BFGS steps
! nshess     number of steps for the Hessian calculation
! ihess      number of Hessian updates
! jump       controls the flow of the optimiser's operations (see flwmsg below)
! lgeook     .true. if the G matrix of all residues is not singular
! lval_hess  .true. if hess(:,:) contains a valid Hessian (or none is required)
! lcart_hess .true. if hess(:,:) contains a Cartesian Hessian
!------------------------------------------------------------------------------

  type opt_stat
     integer ncyc
     integer nstep
     integer nsprfo
     integer nstried
     integer nslbfg
     integer nshess
     integer ihess
     integer jump
     logical lgeook
     logical lval_hess
     logical lcart_hess
  end type opt_stat

  type(opt_stat) :: stat

!------------------------------------------------------------------------------
! Optimiser flow control messages according to stat%jump
!------------------------------------------------------------------------------

  character*44, dimension(-1:7) :: flwmsg

  data flwmsg(-1) / 'restart the optimisation' /
  data flwmsg(0) / 'start the optimisation' /
  data flwmsg(1) / 'start the optimisation, initial point tested' /
  data flwmsg(2) / 'form LBFGS step' /
  data flwmsg(3) / 'perform LBFGS step' /
  data flwmsg(4) / 'test new LBFGS point' /
  data flwmsg(5) / 'form PRFO step' /
  data flwmsg(6) / 'test PRFO step' /
  data flwmsg(7) / 'compute initial hessian' /

!------------------------------------------------------------------------------
! Other global variables
!
! stdout unit number of standard output
! stderr unit number of standard error
! chkpnt unit number of checkpoint file
! chknam name of checkpoint file
! lform  .true. if formatted I/O of checkpoint file
!------------------------------------------------------------------------------

  logical lform
  integer stdout, stderr, chkpnt
  character*12 chknam
  integer, dimension(:), allocatable :: dummyint
  real(kind=8), dimension(:), allocatable :: dummyarg
  data lform / .true. /
  data stdout, stderr, chkpnt / 6, 0, 41 /
  data chknam / 'hdlcopt.chk' /

  interface dummyarg_checkin
     module procedure dummyrel_checkin, dummyint_checkin
  end interface

contains

!------------------------------------------------------------------------------
! dummyarg_checkin / dummyarg_clear
! dummyarg_alloc / dummyarg_checkout
!
! Work around a limitation of F90:
! Assumed-shape arrays may not be contiguous. Passing them to a routine
! with an array as formal argument, they must be copied.
! To avoid unnecessary copies, elements of assumed-shape array may not be
! passed as arrays.
!------------------------------------------------------------------------------

  subroutine dummyrel_checkin (arr, element, number)
    integer, intent(in) :: element, number
    real(kind=8), dimension(*), intent(in) :: arr
    integer i, ia
    if (allocated(dummyarg)) then
       call hdlc_errflag ('dummyarg_checkin: already allocated', 'abort')
    else
       allocate (dummyarg(number))
    end if
    ia = element
    do i = 1,number
       dummyarg(i) = arr(ia)
       ia = ia + 1
    end do
  end subroutine dummyrel_checkin

  subroutine dummyint_checkin (arr, element, number)
    integer, intent(in) :: element, number
    integer, dimension(*), intent(in) :: arr
    integer i, ia
    if (allocated(dummyint)) then
       call hdlc_errflag ('dummyarg_checkin: already allocated', 'abort')
    else
       allocate (dummyint(number))
    end if
    ia = element
    do i = 1,number
       dummyint(i) = arr(ia)
       ia = ia + 1
    end do
  end subroutine dummyint_checkin

  subroutine dummyarg_clear
    if (allocated(dummyarg)) then
       deallocate (dummyarg)
    else if (allocated(dummyint)) then
       deallocate (dummyint)
    else
       call hdlc_errflag ('dummyarg_clear: not allocated', 'abort')
    end if
  end subroutine dummyarg_clear

  subroutine dummyarg_alloc (number)
    integer number
    if (allocated(dummyarg)) then
       call hdlc_errflag ('dummyarg_alloc: already allocated', 'abort')
    else
       allocate (dummyarg(number))
    end if
  end subroutine dummyarg_alloc

  subroutine dummyarg_checkout (arr, element, number)
    integer element, number
    real(kind=8), dimension(*), intent(out) :: arr
    integer i, ia
    if (.not.allocated(dummyarg)) then
       call hdlc_errflag ('dummyarg_checkout: not allocated', 'abort')
    end if
    ia = element
    do i = 1,number
       arr(ia) = dummyarg(i)
       ia = ia + 1
    end do
    deallocate (dummyarg)
  end subroutine dummyarg_checkout

!------------------------------------------------------------------------------
! subroutine hdlc_wr_ctrl
!
! Dump the optimisation control data structure s_ctrl to iunit
!------------------------------------------------------------------------------

  subroutine hdlc_wr_ctrl (iunit, s_ctrl, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    type(opt_ctrl) :: s_ctrl

! local vars
    integer k, m

! begin
    lerr = .false.
    if (lform) then
       write (iunit,'(5I8)',err=98) s_ctrl%mstart, s_ctrl%nbig, &
            s_ctrl%nvar, s_ctrl%nrem, s_ctrl%nrestart
       write (iunit,'(40I2)',err=98) &
            (s_ctrl%attypes(m), m=1,(s_ctrl%nbig/3))
       write (iunit,'(7F11.6)',err=98) s_ctrl%toler, s_ctrl%rmin, &
            s_ctrl%rmax, s_ctrl%dmax, s_ctrl%ddmax, s_ctrl%osmin, s_ctrl%omin
       write (iunit,'(10I4)',err=98) s_ctrl%mode, s_ctrl%irecalc, &
            s_ctrl%maxcyc, s_ctrl%iupd, s_ctrl%mhess, s_ctrl%idump, &
            s_ctrl%mdump, s_ctrl%iprj, s_ctrl%ilock
       write (iunit,'(8F10.5)',err=98) &
            (s_ctrl%precon(m), m=1,s_ctrl%nbig)
       write (iunit,'(5I8)',err=98) s_ctrl%nincon, s_ctrl%mconstr, &
            s_ctrl%nconstr, s_ctrl%contyp, s_ctrl%ctfirst
       write (iunit,'(10I8)',err=98) &
            ((s_ctrl%incon(k,m), k=1,2), m=1,s_ctrl%nincon)
       write (iunit,'(10I8)',err=98) &
            ((s_ctrl%iconstr(k,m), k=1,4), m=1,s_ctrl%nconstr)
       write (iunit,'(8F10.5)',err=98) &
            (s_ctrl%vconstr(m), m=1,s_ctrl%nconstr)
       write (iunit,'(F10.5)',err=98) s_ctrl%cfact
       write (iunit,'(I4)',err=98) s_ctrl%printl
       write (iunit,'(I8)',err=98) s_ctrl%nrad
       write (iunit,'(8I8)',err=98) (s_ctrl%irad(m), m=1,s_ctrl%nrad)
       write (iunit,'(8F10.5)',err=98) (s_ctrl%vrad(m), m=1,s_ctrl%nrad)
       write (iunit,'(I2)',err=98) s_ctrl%iccode
       write (iunit,'(I2)',err=98) s_ctrl%icconv
    else
       write (iunit,err=98) s_ctrl%mstart, s_ctrl%nbig, &
            s_ctrl%nvar, s_ctrl%nrem, s_ctrl%nrestart
       write (iunit,err=98) &
            (s_ctrl%attypes(m), m=1,(s_ctrl%nbig/3))
       write (iunit,err=98) s_ctrl%toler, s_ctrl%rmin, &
            s_ctrl%rmax, s_ctrl%dmax, s_ctrl%ddmax, s_ctrl%osmin, s_ctrl%omin
       write (iunit,err=98) s_ctrl%mode, s_ctrl%irecalc, &
            s_ctrl%maxcyc, s_ctrl%iupd, s_ctrl%mhess, s_ctrl%idump, &
            s_ctrl%mdump, s_ctrl%iprj, s_ctrl%ilock
       write (iunit,err=98) &
            (s_ctrl%precon(m), m=1,s_ctrl%nbig)
       write (iunit,err=98) s_ctrl%nincon, s_ctrl%mconstr, &
            s_ctrl%nconstr, s_ctrl%contyp, s_ctrl%ctfirst
       write (iunit,err=98) &
            ((s_ctrl%incon(k,m), k=1,2), m=1,s_ctrl%nincon)
       write (iunit,err=98) &
            ((s_ctrl%iconstr(k,m), k=1,4), m=1,s_ctrl%nconstr)
       write (iunit,err=98) &
            (s_ctrl%vconstr(m), m=1,s_ctrl%nconstr)
       write (iunit,err=98) s_ctrl%cfact
       write (iunit,err=98) s_ctrl%printl
       write (iunit,err=98) s_ctrl%nrad
       write (iunit,err=98) (s_ctrl%irad(m), m=1,s_ctrl%nrad)
       write (iunit,err=98) (s_ctrl%vrad(m), m=1,s_ctrl%nrad)
       write (iunit,err=98) s_ctrl%iccode
       write (iunit,err=98) s_ctrl%icconv
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_wr_ctrl

!------------------------------------------------------------------------------
! subroutine hdlc_rd_ctrl
!
! Read the optimisation control data structure s_ctrl from iunit and allocate
! the pointer arrays
!------------------------------------------------------------------------------

  subroutine hdlc_rd_ctrl (iunit, s_ctrl, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(opt_ctrl) :: s_ctrl

! local vars
    integer k, m

! begin
    lerr = .false.
    if (lform) then
       read (iunit,'(5I8)',err=98) s_ctrl%mstart, s_ctrl%nbig, &
            s_ctrl%nvar, s_ctrl%nrem, s_ctrl%nrestart
       if (associated (s_ctrl%attypes)) deallocate (s_ctrl%attypes)
       allocate (s_ctrl%attypes (s_ctrl%nbig/3))
       read (iunit,'(40I2)',err=98) &
            (s_ctrl%attypes(m), m=1,(s_ctrl%nbig/3))
       read (iunit,'(7F11.6)',err=98) s_ctrl%toler, s_ctrl%rmin, &
            s_ctrl%rmax, s_ctrl%dmax, s_ctrl%ddmax, s_ctrl%osmin, s_ctrl%omin
       read (iunit,'(10I4)',err=98) s_ctrl%mode, s_ctrl%irecalc, &
            s_ctrl%maxcyc, s_ctrl%iupd, s_ctrl%mhess, s_ctrl%idump, &
            s_ctrl%mdump, s_ctrl%iprj, s_ctrl%ilock
       if (associated (s_ctrl%precon)) deallocate (s_ctrl%precon)
       allocate (s_ctrl%precon (s_ctrl%nbig))
       read (iunit,'(8F10.5)',err=98) &
            (s_ctrl%precon(m), m=1,s_ctrl%nbig)
       read (iunit,'(5I8)',err=98) s_ctrl%nincon, s_ctrl%mconstr, &
            s_ctrl%nconstr, s_ctrl%contyp, s_ctrl%ctfirst
       if (associated (s_ctrl%incon)) deallocate (s_ctrl%incon)
       allocate (s_ctrl%incon (2,s_ctrl%nincon))
       read (iunit,'(10I8)',err=98) &
            ((s_ctrl%incon(k,m), k=1,2), m=1,s_ctrl%nincon)
       if (associated (s_ctrl%iconstr)) deallocate (s_ctrl%iconstr)
       allocate (s_ctrl%iconstr (4,s_ctrl%nconstr))
       read (iunit,'(10I8)',err=98) &
            ((s_ctrl%iconstr(k,m), k=1,4), m=1,s_ctrl%nconstr)
       if (associated (s_ctrl%vconstr)) deallocate (s_ctrl%vconstr)
       allocate (s_ctrl%vconstr (s_ctrl%nconstr))
       read (iunit,'(8F10.5)',err=98) &
            (s_ctrl%vconstr(m), m=1,s_ctrl%nconstr)
       read (iunit,'(F10.5)',err=98) s_ctrl%cfact
       read (iunit,'(I4)',err=98) s_ctrl%printl
       read (iunit,'(I8)',err=98) s_ctrl%nrad
       read (iunit,'(8I8)',err=98) (s_ctrl%irad(m), m=1,s_ctrl%nrad)
       read (iunit,'(8F10.5)',err=98) (s_ctrl%vrad(m), m=1,s_ctrl%nrad)
       read (iunit,'(I2)',err=98) s_ctrl%iccode
       read (iunit,'(I2)',err=98) s_ctrl%icconv
    else
       read (iunit,err=98) s_ctrl%mstart, s_ctrl%nbig, &
            s_ctrl%nvar, s_ctrl%nrem, s_ctrl%nrestart
       if (associated (s_ctrl%attypes)) deallocate (s_ctrl%attypes)
       allocate (s_ctrl%attypes (s_ctrl%nbig/3))
       read (iunit,err=98) &
            (s_ctrl%attypes(m), m=1,(s_ctrl%nbig/3))
       read (iunit,err=98) s_ctrl%toler, s_ctrl%rmin, &
            s_ctrl%rmax, s_ctrl%dmax, s_ctrl%ddmax, s_ctrl%osmin, s_ctrl%omin
       read (iunit,err=98) s_ctrl%mode, s_ctrl%irecalc, &
            s_ctrl%maxcyc, s_ctrl%iupd, s_ctrl%mhess, s_ctrl%idump, &
            s_ctrl%mdump, s_ctrl%iprj, s_ctrl%ilock
       if (associated (s_ctrl%precon)) deallocate (s_ctrl%precon)
       allocate (s_ctrl%precon (s_ctrl%nbig))
       read (iunit,err=98) &
            (s_ctrl%precon(m), m=1,s_ctrl%nbig)
       read (iunit,err=98) s_ctrl%nincon, s_ctrl%mconstr, &
            s_ctrl%nconstr, s_ctrl%contyp, s_ctrl%ctfirst
       if (associated (s_ctrl%incon)) deallocate (s_ctrl%incon)
       allocate (s_ctrl%incon (2,s_ctrl%nincon))
       read (iunit,err=98) &
            ((s_ctrl%incon(k,m), k=1,2), m=1,s_ctrl%nincon)
       if (associated (s_ctrl%iconstr)) deallocate (s_ctrl%iconstr)
       allocate (s_ctrl%iconstr (4,s_ctrl%nconstr))
       read (iunit,err=98) &
            ((s_ctrl%iconstr(k,m), k=1,4), m=1,s_ctrl%nconstr)
       if (associated (s_ctrl%vconstr)) deallocate (s_ctrl%vconstr)
       allocate (s_ctrl%vconstr (s_ctrl%nconstr))
       read (iunit,err=98) &
            (s_ctrl%vconstr(m), m=1,s_ctrl%nconstr)
       read (iunit,err=98) s_ctrl%cfact
       read (iunit,err=98) s_ctrl%printl
       read (iunit,err=98) s_ctrl%nrad
       read (iunit,err=98) (s_ctrl%irad(m), m=1,s_ctrl%nrad)
       read (iunit,err=98) (s_ctrl%vrad(m), m=1,s_ctrl%nrad)
       read (iunit,err=98) s_ctrl%iccode
       read (iunit,err=98) s_ctrl%icconv
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_rd_ctrl

!------------------------------------------------------------------------------
! subroutine hdlc_wr_stat
!
! Dump the optimiser status register structure s_stat to iunit
!------------------------------------------------------------------------------

  subroutine hdlc_wr_stat (iunit, s_stat, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    type(opt_stat) :: s_stat

! begin
    lerr = .false.
    if (lform) then
       write (iunit,*,err=98) s_stat
    else
       write (iunit,err=98) s_stat
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_wr_stat

!------------------------------------------------------------------------------
! subroutine hdlc_rd_stat
!
! Read the optimiser status register structure s_stat from iunit
!------------------------------------------------------------------------------

  subroutine hdlc_rd_stat (iunit, s_stat, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(opt_stat) :: s_stat

! local vars
    integer m

! begin
    lerr = .false.
    if (lform) then
       read (iunit,*,err=98) s_stat
    else
       read (iunit,err=98) s_stat
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_rd_stat

!------------------------------------------------------------------------------
! subroutine flushout
!
! Flush the I/O buffers of standard output - this is just a dummy procedure
!------------------------------------------------------------------------------

  subroutine flushout ()
  end subroutine flushout

!------------------------------------------------------------------------------
! high level error handler - see search.f90
!------------------------------------------------------------------------------

  subroutine hdlc_errflag (message, action)
    character*(*) message, action
    if (action .eq. 'warn') then
       write (stderr,'(A,A)') 'Warning: ', message
    else if (action .eq. 'stop') then
       write (stderr,'(A,A)') 'Stopped: ', message
    else
       write (stderr,'(A,A)') 'Aborting: ', message
    end if
    call search_stat

    write (stderr,*) 'ERROR CODE', ctrl%iccode

    call hdlc_error (ctrl%iccode, message, len(message), &
         action, len(action))

    if (action .ne. 'warn') then
        stop 999
    endif

  end subroutine hdlc_errflag

!------------------------------------------------------------------------------
! subroutine search_stat: write some statistics
!------------------------------------------------------------------------------

  subroutine search_stat
    if (ctrl%printl.ge.0) then
       write (stdout,'(/,A)') 'Summary of the optimisation'
       write (stdout,'(3X,A,I8)') 'Function / gradient evaluations:  ', &
            stat%ncyc
       write (stdout,'(3X,A,I8)') 'Number of tried P-RFO steps:      ', &
            max(0,stat%nstep-1)
       write (stdout,'(3X,A,I8)') 'Number of performed P-RFO steps:  ', &
            stat%nsprfo
       write (stdout,'(3X,A,I8)') 'Number of tried L-BFGS steps:     ', &
            stat%nstried
       write (stdout,'(3X,A,I8)') 'Number of performed L-BFGS steps: ', &
            stat%nslbfg
       write (stdout,'(3X,A,I8,/)') 'Number of Hessian updates:        ', &
            stat%ihess
    end if
  end subroutine search_stat

end module global
