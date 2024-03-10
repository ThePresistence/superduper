!------------------------------------------------------------------------------
! Interface between MNDO and HDLCopt. Although in Fortran90, no interface
! descriptions are used nor are they passed to any side.
! HDLOPT passes control to hdlcopt () and stores all relevant pointers to
! mndogl.
! The relevant common blocks from MNDO are accessible to all subroutines in
! this file and are never directly used by HDLCopt.
! These common blocks from MNDO are used in this file:
!   ATOMC, ATOMS, CGRAD, CONSTF, CYCLES, DFP, ERG, INOPT2, INOPT4, NBFILE,
!   PARM1, PARM3.
!
! These subroutines are provided to MNDO:
!   HDLOPT
!
! These subroutines are provided to HDLCopt (see search.f90):
!   hdlc_get_params
!   hdlc_get_resn
!   hdlc_get_coords
!   hdlc_get_hess
!   hdlc_update
!   hdlc_put_coords
!   hdlc_wr_fandg / hdlc_rd_fandg
!   hdlc_error
!------------------------------------------------------------------------------
! Notation for some common variables:
! - numat   : total number of atoms
! - natopt  : number of optimized atoms (numat-natopt atoms are frozen)
! - nbig    : number of optimized variables (3*natopt)
! - nconstr : number of constrained variables
! - nvar    : number of optimized variables in reaction core
!
!------------------------------------------------------------------------------
! Specification of the interface MNDO <-> HDLCopt
!
! - Defaults: Cartesians, no residues, no constraints, no frozen atoms,
!   tolerances as in Gaussian.
! - No matter whether Z-matrix or Cartesians are used as input in MNDO,
!   only Cartesians are passed between MNDO and HDLCopt.
! - Constraints and extra parameters can be specified in three ways:
!   * The data are read from a dedicated Fortran unit. No limits apply.
!     The input is keyword-driven, additional data are read in free format.
!     MNDO does not 'see' the specifications.
!     Keywords recognised by the routine reading the unit 42 (rdhdlc):
!     . NFCART=<n>: to specify frozen atoms;
!       followed by line(s) with corresponding sequence numbers
!     . NFCOMP=<n>: to specify frozen individual Cartesian coordinates;
!       followed by line(s) with corresponding variable numbers
!     . NFBOND=<n>: to specify n bond length constraints;
!       followed by n lines with relevant sequence numbers and target values
!     . NFANGL=<n>: to specify n bond angle constraints;
!       followed by n lines with relevant sequence numbers and target values
!     . NFDIHE=<n>: to specify n dihedral angle constraints;
!       followed by n lines with relevant sequence numbers and target values
!     . NMRESI=<n>: to specify one residue with n atoms;
!       followed by line(s) with corresponding sequence numbers
!     . NRCORE=<n>: to specify the n-th residue as reaction core;
!     . NDFCOR=<i>: to specify the number of degrees of freedom in the core;
!       nvar=nbig forces a TS search,
!       nbig>nvar>0 forces a microiterative TS search,
!       default value assigned according to NRCORE.
!     . XYZFAC=<value>: to specify the Cartesian scale factor for HDLC
!     . MDCRAD=<i>: to specify covalent radii used to determine connectivity;
!       0: all atoms are treated as H, 1: radii from HDLC code are used;
!       2: covalent radii are defined via input.
!     . NCRAD=<n>: to specify covalent radii for n elements via input;
!       followed by n lines with relevant nuclear charges and radii
!     . NREMST=<n>: to specify the number of remembered L-BFGS steps
!     . NRSHDL=<n>: to force HDLC restarts every n cycles
!   * The data are read from an additional MNDO input block.
!     The conventions are the same as in the preceding case.
!   * Constraints are derived from the input Z-matrix or the input Cartesians.
!     . Restriction: no dummy atoms are allowed.
!     . Restriction: no symmetry is imposed.
! - The atoms are reordered using mndogl%map:
!   i_hdlc = mndogl%map (i_mndo)
! - The residue memberships are stored in mndogl%resn:
!   mndogl%resn (i_hdlc) specifies the residue number of atom i_hdlc.
!   mndogl%resn (i_hdlc) = 0 if atom i_hdlc does not belong to any residue.
!   mndogl%resn (i_hdlc) = 0 implies that Cartesians are used for i_hdlc.
!   mndogl%resn (i_hdlc) =-1 if atom i_hdlc is frozen.
! - The residue membership can be specified in two ways (see above):
!   * by input from a dedicated Fortran unit,
!   * by input from a new MNDO input block.
! - Checkpointing:
!   * HDLCopt reads/writes its own checkpoint to/from an own Fortran unit
!     according to MIDDLE.
!   * HDLCopt does not use DFPSAV.
!   * HDLCopt calls DENSAV according to MIDDLE.
! - Time limit:
!   * hdlc_get_coords() estimates the time for the next cycle and aborts
!     if required according to LIMIT.
! - HDLOPT is called from GEOOPT.
! - hdlc_get_coords() calls SCF directly.
! - hdlcopt() carries along pointers for SCFCAL.
!
!------------------------------------------------------------------------------
! inter-subroutine communication
! type mndoglobal        !                                                init
!    sequence
!    logical lguess      ! initialization flag, true only on first call   true
!    logical lgtres      ! true if residues are defined                  false
!    logical lgtcns      ! true if constraints are defined               false
!    logical lgtfrz      ! true if frozen atoms are defined              false
!    integer icall       ! control flag for CALL SCF                        31
!    integer hdlcsp      ! file number for HDLC input, nb5 or nb42          42
!    integer nrem        ! number of remembered L-BFGS steps                 0
!    integer matt        ! option for covalent radii                         0
!    integer nvar        ! number of degrees of freedom in reaction core     0
!    integer nrestart    ! number of cycles between HDLC restarts            0
!    real :: cfac        ! weighting factor for Cartesians in HDLC         0.0
!    real :: tsum        ! cpu time                                        0.0
!    integer, pointer :: resn (:)    ! residue number, -1 if frozen          0 
!    integer, pointer :: map  (:)    ! i-hdlc = mndogl%map(i-mndo)           0
!    integer, pointer :: irad (:)    ! covalent radii, atomic numbers    undef
!    integer, pointer :: icons (4,:) ! constraints, labels               undef
!    real, pointer :: vcons (:)      ! constraints, target values        undef
!    real, pointer :: vrad (:)       ! covalent radii, values            undef
! end type mndoglobal
! type(mndoglobal) :: mndogl
! common /mndohdlc/ mndogl
!------------------------------------------------------------------------------
!
! -> Files used by HDLOPT:
!    --------------------
! Name   Type   Passing method    Description
! ----   ----   --------------    -----------
! hdlcsp int    in mndolib.f90    Fortran unit number to read the constraints
!                                 and residues from: 42
! chkpnt int    in global.f90     Fortran unit number for restart file: 41
!
! -> Existing MNDO variables read / modified by HDLOPT (from stdin etc.)
!    ------------------------------------------------------------------
! Name   Type   Passing method  HDLOPT Description
! ----   ----   --------------  ------------------
! JOP    int    common /INOPT2/ in     optimisation mode: IN2(3)
! IEF    int    common /INOPT2/ in     choice of optimiser: IN2(8)
! NATOMS int    common /PARM1 / ---    Number of atoms read in input
! NAT    int    common /ATOMS / in     Atomic numbers (no dummies)
! NUMAT  int    common /ATOMS / in     Number of atoms (excluding dummies)
! COORD  real*8 common /ATOMC / intern internally used for GMETRY / SCF
! COORD  real*8 common /ATOMC / in/out Cartesian coordinates used by HDLCopt
! CG     real*8 common /CGRAD / intern Cartesian gradient used by HDLCopt
! IGEOM  int    common /INOPT2/ in     type of coordinates (out: cart): IN2(4)
! LOC    int    common /PARM3 / in     mapping active <-> all variables (w/dummies)
! 
! -> Values of existing and new MNDO parameters for HDLOPT
!    -----------------------------------------------------
! Name   Value  IN2(xx) Description
! ----   -----  ------- -----------
! JOP    *      3       Optimisation mode (as in MNDO)
! IGEOM  0,1    4       Type of input coordinates
! IEF    <0     8       Use HDLCopt
! MAXEND *      32      Maximum number of energy / gradient evals (maxcyc)
! MIDDLE -1     37      No restart, no dump (mdump = 0, idump = 0, mstart = 0)
! MIDDLE 0      37      No restart, but dump (mdump = 0, idump = 1, mstart = 0)
! MIDDLE 1      37      Restart and dump (mdump = 0, idump = 1, mstart = 3)
! IPRINT *      38      Print level (decremented by one for HDLCopt)
! IPREC  *      43      Divisor for tolerances (0.000450)
! ICONV  <>-3   44      Default Gaussian convergence criteria
! ICONV  -3     44      Cartesian RMS/max gradient convergence criteria
! IHDLC1 0,1    51      Coordinates used by HDLC optimizer
! IHDLC2 0-4    52      Connection type of 1st residue
! IHDLC3 0,1,2  53      Mode of special HDLC input
! KTRIAL *      67      KTRIAL=0 reset to KTRIAL=11 upon restart (MIDDLE=1)
! MODE   *      80      Mode to be followed
! IRECLC 0      81      Core Hessian is computed once (JOP=1,4) or never (else)
! IRECLC >0     81      Interval for the calculation of the Hessian
! IUPD   *      82      Mode of Hessian update (0: none, 1: Powell, 2: BFGS)
! NRST   999    191     Take the HDLCopt defaults for irecalc (core Hessian)
! NRST   <>999  191     Regular recalculation of the Hessian (irecalc=nrst)
!
! -> Real parameters and default translation
! Name   MNDOopt HDLCopt XN4(xx) Description
! ----   ------- ------- ------- -----------
! DMAX   0.2      0.2    1       Initial trust radius
! DDMIN  0.001    0.0002 2       Minimum trust radius
! DDMAX  0.5/0.3  1.0    3       Maximum trust radius
! RMIN   0.0      0.0    4       Minimum ratio predicted/actual energy change
! RMAX   4.0     10.0    5       Maximum ratio predicted/actual energy change
! OMIN   0.8      0.6    6       Minimum overlap between the eigenmodes
!
! -> Potential conflicts
!    -------------------
! - The Fortran unit of the HDLCopt checkpoint file is defined in global.f90
! - Upon restart (MIDDLE=1), MIDDLE is set to 0
! - Upon restart (MIDDLE=1), KTRIAL=0 has previously been reset to KTRIAL=11;
!   this has now been commented out since subroutine INPUT sets KTRIAL=-1 for
!   MIDDLE=1, the previous density matrix is read during restart via DENSAV.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Pass control to hdlcopt
!------------------------------------------------------------------------------

subroutine HDLOPT (array, lm5, icall, scfcal)
  use limit, only: LM1
  implicit none

! args
  integer, target :: lm5
  integer icall
  real (kind=8), dimension(lm5), target :: array
  external scfcal

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! MNDO common blocks
  integer NAT, NUMAT, IN2, NBF, NB6
  common /ATOMS / NUMAT,NAT(LM1,3)
  common /INOPT2/ IN2(300)
  common /NBFILE/ NBF(20)

! local vars
  logical linit
  integer ihdlc3, iprint, iret, ktrial, middle, natopt, nconstr, nrad
  real (kind=8) tx1, tx2

! init local data
1 continue
  ihdlc3 = IN2(53)
  if (ihdlc3.eq.1) then
     mndogl%hdlcsp = NBF(5)
  else
     mndogl%hdlcsp = 42
  endif
  middle = IN2(37)
  iprint = IN2(38)
! ktrial = IN2(67)
  NB6 = NBF(6)

! begin
  mndogl%lguess = .true.
  linit = .true.
  if (iprint.ge.0) then
     call CPUSEC (tx1)
  endif

! restart from checkpoint
  if (middle .eq. 1) then
     call densav (0,array,lm5,1)
  end if

! prepare args
  nconstr = 0; nrad = 0
  natopt = NUMAT
  call rdhdlc (mndogl%hdlcsp, nconstr, NUMAT, natopt, nrad)

! pass control to HDLCopt
  call hdlcopt (2, NUMAT, natopt, 0, nconstr, nrad, iret, linit, NB6, &
       scfcal, array, lm5)

! convert return codes from hdlcopt to return codes for MNDO (see search.f90)
  select case (iret)
  case (1)
     icall = 10
  case (-1, -3)
     icall = -3
  case (-4)
     icall = -1
  case (-5)
     icall = -8
  case (-9)
     if (in2(51).eq.0) then
        write (NB6,'(A)') ' Cannot optimize in HDLC coordinates '
        write (NB6,'(A)') ' Switch to Cartesian coordinates '
        in2(51) = 1
        go to 1
     else
        write (NB6,'(A)') ' Failure to optimize in Cartesian coordinates '
        stop 'HDLOPT'
     endif
  case default
     icall = 0
  end select
  if (icall .eq. -1) return

! print cpu times
  if (iprint.ge.0) then
     call CPUSEC (tx2)
     write (NB6,'(/)')
     write (NB6,'(A)') ' Cpu time for HDLC optimization '
     write (NB6,'(A,F16.5,A)') ' Total:', tx2-tx1, ' seconds '
     write (NB6,'(A,F16.5,A)') ' SCF  :', mndogl%tsum, ' seconds '
     write (NB6,'(A,F16.5,A)') ' HDLC :', tx2-tx1-mndogl%tsum, ' seconds '
  endif 

! clean up, due to a bug in SGI F90, these arrays are always allocated
  deallocate (mndogl%icons); deallocate (mndogl%vcons)
  deallocate (mndogl%irad); deallocate (mndogl%vrad)
  deallocate (mndogl%resn); deallocate (mndogl%map)

end subroutine HDLOPT

!------------------------------------------------------------------------------
! Read optimisation control parameters and residue membership information
! from MNDO
!
! Arguments:
! - declaration: fields definition of opt_ctrl
! - exception:   pointers to arrays are replaced by the preallocated arrays
! - linit:       set mstart=4 and exit immediately if .false.
!------------------------------------------------------------------------------

subroutine hdlc_get_params (iccode, linit, mstart, nbig, nvar, nrem, &
       nrestart, attypes, &
       toler, rmin, rmax, dmax, ddmax, ddmin, omin, &
       mode, irecalc, maxcyc, iupd, mhess, &
       idump, mdump, iprj, ilock, precon, &
       nincon, mconstr, nconstr, contyp, ctfirst, &
       incon, iconstr, vconstr, cfact, printl, nrad, irad, vrad, icconv)
  use limit, only: LM1
  implicit none

! args
  logical linit
  integer iccode, mstart, nbig, nvar, nrem
  integer nrestart
  integer, dimension(nbig/3) :: attypes
  real (kind=8) toler, rmin, rmax, dmax, ddmax, ddmin, omin
  integer mode, irecalc, maxcyc, iupd, mhess, idump, mdump, iprj, ilock
  real (kind=8), dimension(nbig) :: precon
  integer nincon, mconstr, nconstr, contyp, ctfirst
  integer, dimension(2,nincon) :: incon
  integer, dimension(4,nconstr) :: iconstr
  real (kind=8), dimension(nconstr) :: vconstr
  real (kind=8) cfact
  integer printl, nrad, icconv
  integer, dimension(nrad) :: irad
  real (kind=8), dimension(nrad) :: vrad

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! MNDO common blocks
  integer IN2, NUMAT, NAT, NFIRST, NLAST, NBF, NB6
  real (kind=8) XN4
  common /ATOMS / NUMAT,NAT(LM1),NFIRST(LM1),NLAST(LM1)
  common /INOPT2/ IN2(300)
  common /INOPT4/ XN4(50)
  common /NBFILE/ NBF(20)

! local vars
  integer i, iconv, ihdlc2, iprec, iprint, ireclc, j, jop, ktrial, &
       maxend, middle, mmode, natom, nrst

! begin, evaluate MNDO integer variables
  jop = IN2(3)
  maxend = IN2(32)
  middle = IN2(37)
  iprint = IN2(38)
  iprec = IN2(43)
  iconv = IN2(44)
  ihdlc2 = IN2(52)
! ktrial = IN2(67)
  mmode = IN2(80)
  ireclc = IN2(81)
  iupd = IN2(82)
  nrst = IN2(191)
  NB6 = NBF(6)

! evaluate MNDO real variables and translate defaults
  toler = 0.000450_8 / real (iprec, 8)
  dmax  = XN4(1)  ! default 0.02_8   ! 0.2_8 for ief=1
  ddmin = XN4(2)  ! default 0.0002_8 ! 0.001 for ief=1
  ddmax = XN4(3)  ! default 1.0_8    ! 0.3_8 for ief=1
  rmin  = XN4(4)  ! default 0.0_8    ! same  for ief=1
  rmax  = XN4(5)  ! default 10.0_8   ! 4.0_8 for ief=1
  omin  = XN4(6)  ! default 0.6_8    ! 0.8_8 for ief=1

! set defaults and parameters
  if (jop.eq.1 .or. jop.eq.4) then
     if (ihdlc2.eq.3) then
        if(numat.gt.2) then
           nvar = nbig-6-nconstr
        else
           nvar = nbig-5-nconstr
        end if
     else
        nvar = nbig-nconstr
     end if
  else
     nvar = 0
  end if
  if (mndogl%nvar.ne.0) then
     nvar = mndogl%nvar
     if (iprint.ge.2) write (NB6,'(1X,A,I5)') &
          ' Forcing TS search by setting nvar=', nvar
  end if
  natom = nbig/3
  if (mndogl%nrem.eq.0) then
     if (nbig .le. 40) then
        nrem = nbig
     else
        nrem = 40
     end if
  else
     nrem = mndogl%nrem
  end if
  do i = 1,natom
     if (mndogl%matt.eq.0) then
        attypes(i) = 1
     else
        attypes(mndogl%map(i)) = NAT(i)
     end if
  end do
  nrestart = mndogl%nrestart
  mode = mmode
  printl = iprint-1
  irecalc = ireclc
  if (nrst .ne. 999) irecalc = nrst
  maxcyc = maxend ! was 1000000
  if (printl.ge.-1) then
     if (mode .eq. 0) then
        if (iupd .ne. 2) &
             write (NB6,'(A,I1,A,I3)') ' Warning: iupd=', iupd, ', mode=', mode
     else
        if (iupd .ne. 1) &
             write (NB6,'(A,I1,A,I3)') ' Warning: iupd=', iupd, ', mode=', mode
     end if
  end if
  mhess = 0
  if (middle .lt. 0) then
     idump = 0
  else
     idump = 1
  end if
  mdump = 0
  if (nvar .eq. nbig) then
     iprj = 0 ! set to one if R/T fit is done
  else
     iprj = 0
  end if
  do i = 1,nbig
     precon(i) = 1.0_8
  end do
  contyp = 0
  ctfirst = ihdlc2 ! 0
  cfact = mndogl%cfac ! 0.0_8
  mconstr = 0

! constraints
  do i = 1,nconstr
     do j = 1,4
        iconstr(j,i) = mndogl%icons(j,i)
     end do
     vconstr(i) = mndogl%vcons(i)
  end do

! all the values above can be superseded on restart from checkpoint
  if (middle .eq. 1) then
     mstart = 3
     middle = 0
!    if (ktrial .eq. 0) ktrial = 11
  end if
  if (.not. linit) mstart = 4

! atomic radii
  if (nrad.gt.0) then
     if (printl.ge.-1 .and. mndogl%matt.eq.0) write (NB6,'(A,I1,A,I2)') &
          ' Warning: nrad=', nrad, ' but matt=', mndogl%matt
     do i = 1,nrad
        irad(i) = mndogl%irad(i)
        vrad(i) = mndogl%vrad(i)
        if (printl.ge.1) write (NB6,'(A,I2,A,F10.4)') &
             ' Covalent radius of atom ', irad(i), ' is ', vrad(i)
     end do
  end if

! convergence criteria
  if (iconv.eq.-3 .or. iconv.eq.4) then
     icconv = 1
  else
     icconv = 0
  end if

! pass possibly changed parameters to MNDO
  IN2(37) = middle
! IN2(67) = ktrial

end subroutine hdlc_get_params

!------------------------------------------------------------------------------
! Get residue memberships
!------------------------------------------------------------------------------

subroutine hdlc_get_resn (iccode, resn, natom)
  implicit none

! args
  integer iccode, natom
  integer, dimension(natom) :: resn

! MNDO common blocks
  integer IN2, NBF, NB6
  common /INOPT2/ IN2(300)
  common /NBFILE/ NBF(20)

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! local vars
  integer i, ihdlc1

  NB6 = NBF(6)

! begin
  ihdlc1 = IN2(51)
  if (mndogl%lgtres) then
     do i = 1,natom
        resn(i) = 0
        if (mndogl%resn(i).ne.-1) resn(i) = mndogl%resn(i)
     end do
     if (ihdlc1.eq.1) then
        write (NB6,'(A,I2)') ' Warning: residues specified but ihdlc1=', ihdlc1
     end if
  else
     if (ihdlc1.eq.0) then
        do i = 1,natom
           resn(i) = 1
        end do
     else
        do i = 1,natom
           resn(i) = 0
        end do
     end if
  end if

end subroutine hdlc_get_resn

!------------------------------------------------------------------------------
! Get coordinates, energy and gradient from MNDO
!------------------------------------------------------------------------------

subroutine hdlc_get_coords (iccode, coords, funct, grad, nbig, ierr, &
     SCFCAL, ARRAY, LM5, mjump)
  use limit, only: LM1, LMV, LM1M
  implicit none

! args
  integer iccode, nbig, ierr, LM5, mjump
  real (kind=8) funct
  real (kind=8), dimension(nbig) :: coords, grad
  real (kind=8), dimension(LM5) :: ARRAY
  external SCFCAL

! local params (CODATA 1986 recommended) ! used in MNDO
  real (kind=8) kcalm, angst
! parameter (kcalm = 1.59360127E-3_8)    ! kcalm = 1.0 / (EV*EVCAL)
! parameter (angst = 1.889725949E0_8)    ! angst = 1.0 / A0

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! MNDO common blocks
  integer ICYC, IN2, N, NAT, NCOUNT, NUMAT, NBF, NB6
  real (kind=8) :: CG, CNORM, COORD, ENERGY, GG, GNORM, X
  real (kind=8) :: A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
  common /ATOMC / COORD(LMV)
  common /ATOMS / NUMAT,NAT(LM1,3)
  common /CGRAD / CG(3*(LM1+LM1M))
  common /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
  common /CYCLES/ ICYC,NCOUNT
  common /DFP   / X(LMV),N
  common /ERG   / ENERGY,GG(LMV),GNORM,CNORM
  common /INOPT2/ IN2(300)
  common /NBFILE/ NBF(20)

! local vars
  integer i, iprint, i_hdlc, i_mndo, j, middle, natopt
  real (kind=8) kcalma, tx1, tx2 
  character*7  currmsg
  character*7, dimension(-1:7) :: mflwmsg
  data mflwmsg(-1) / 'restart' /
  data mflwmsg(0)  / 'start  ' /
  data mflwmsg(1)  / 'start  ' /
  data mflwmsg(2)  / 'LBFGS  ' /
  data mflwmsg(3)  / 'LBFGS  ' /
  data mflwmsg(4)  / 'LBFGS  ' /
  data mflwmsg(5)  / 'PRFO   ' /
  data mflwmsg(6)  / 'PRFO   ' /
  data mflwmsg(7)  / 'Hessian' /

  NB6 = NBF(6)

! begin
  middle = IN2(37)
  iprint = IN2(38)
  angst  = 1.0_8 / A0
  kcalm  = 1.0_8 / (EV*EVCAL)
  kcalma = kcalm / angst

! use the correct coordinates and separate the frozen atoms
  if (mndogl%lguess) then
     N = 3*NUMAT
     i = 0
     natopt = nbig/3
     do i_mndo = 1,NUMAT
        i_hdlc = mndogl%map(i_mndo)
        if (i_hdlc.le.natopt) then
           i_hdlc = 3*(i_hdlc-1)
           do j = 1,3
              coords(i_hdlc+j) = COORD(i+j) * angst
           end do
        end if
        i = i+3
     end do
     mndogl%lguess = .false.
     mndogl%tsum = 0.0_8
     NCOUNT = 0
  endif

! compute energy and gradient 
! mndogl%icall = 31 will always yield the Cartesian gradient
  call CPUSEC (tx1)
  mndogl%icall = 31
  call SCF (ARRAY,LM5,mndogl%icall,SCFCAL)
  NCOUNT = NCOUNT+1
  call CPUSEC (tx2)
  mndogl%tsum = mndogl%tsum + tx2 - tx1
  if (iprint.ge.-1) then
     if (NCOUNT.EQ.1 .or. iprint.gt.0) write (NB6,'(1X)')
     currmsg = ' '
     if (mjump.ge.-1 .and. mjump.le.7) currmsg = mflwmsg(mjump)
!    write (NB6,'( 1X,A,I6,A,F12.3,A,F12.5,A,F15.5,A,I2)') &
     write (NB6,'( 1X,A,I6,A,F12.3,A,F12.5,A,F15.5,3X,A)') &
     'SCF CALL', NCOUNT, ' , TOTAL TIME = ', tx2, &
     ' , HEAT =', ENERGY, ' , CNORM =', CNORM, currmsg
!    ' , HEAT =', ENERGY, ' , CNORM =', CNORM, ' , TASK =', mjump
  endif

! save current density matrix
  if (middle.ge.0) then
     call densav (0,ARRAY,LM5,0)
  endif

! check for error and copy the gradient    
  if (mndogl%icall.eq.-1) then
     ierr = -1
  else
     ierr = 0
     i = 0
     natopt = nbig/3
     do i_mndo = 1,NUMAT
        i_hdlc = mndogl%map(i_mndo)
        if (i_hdlc.le.natopt) then
           i_hdlc = 3*(i_hdlc-1)
           do j = 1,3
!             coords(i_hdlc+j) = angst * X(i+j)
              grad(i_hdlc+j) = kcalma * CG(i+j)
           end do
        end if
        i = i+3
     end do
     funct = kcalm * ENERGY
  end if

end subroutine hdlc_get_coords

!------------------------------------------------------------------------------
! Get analytical Hessian for the reaction core: a dummy routine
!------------------------------------------------------------------------------

subroutine hdlc_get_hess (iccode, hess, nvar, lvalid, lcart)
  implicit none
  logical lvalid, lcart
  integer iccode, nvar
  real (kind=8), dimension(nvar,nvar) :: hess
end subroutine hdlc_get_hess

!------------------------------------------------------------------------------
! Update the pair list: not used for MNDO
!------------------------------------------------------------------------------

subroutine hdlc_update (iccode)
  implicit none
  integer iccode
end subroutine hdlc_update

!------------------------------------------------------------------------------
! Pass new coordinates after each step to MNDO common blocks
! Convert from atomic units to Angstrom
!------------------------------------------------------------------------------

subroutine hdlc_put_coords (iccode, coords, nbig)
  use limit, only: LM1, LMV
  implicit none

! args
  integer iccode, nbig
  real (kind=8), dimension(nbig) :: coords

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! MNDO common blocks
  integer NAT, NUMAT
  real (kind=8) :: COORD
  real (kind=8) :: A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
  common /ATOMC / COORD(LMV)
  common /ATOMS / NUMAT,NAT(LM1,3)
  common /CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP

! local vars
  integer i, i_hdlc, i_mndo, j, natopt

! begin
  i = 0
  natopt = nbig/3
  do i_mndo = 1,NUMAT
     i_hdlc = mndogl%map(i_mndo)
     if (i_hdlc.le.natopt) then
        i_hdlc = 3*(i_hdlc-1)
        do j = 1,3
           COORD(i+j) = A0 * coords(i_hdlc+j)
        end do
     end if
     i = i+3
  end do

end subroutine hdlc_put_coords

!------------------------------------------------------------------------------
! No need to tell MNDO to dump its memory to a checkpoint file
!------------------------------------------------------------------------------

subroutine hdlc_wr_fandg (iccode, chkpnt, lform, lerr)
  implicit none

! args
  logical lerr, lform
  integer chkpnt, iccode

! begin, evaluate MNDO int variables

end subroutine hdlc_wr_fandg

!------------------------------------------------------------------------------
! No need to tell MNDO to restore its memory from a checkpoint file
!------------------------------------------------------------------------------

subroutine hdlc_rd_fandg (iccode, chkpnt, lform, lerr)
  implicit none

! args
  logical lerr, lform
  integer chkpnt, iccode

! begin

end subroutine hdlc_rd_fandg

!------------------------------------------------------------------------------
! Handle HDLCopt exceptions smoothly: 'warn', 'stop', 'abort'
!------------------------------------------------------------------------------

subroutine hdlc_error (iccode, message, lmessage, action, laction)
  implicit none

! args
  character*(*) message, action
  integer iccode, lmessage, laction

! begin
  if (action .eq. 'stop' .or. action .eq. 'abort') then
     stop
  end if

end subroutine hdlc_error

!------------------------------------------------------------------------------
! Helper routine for HDLOPT:
!
! Read the HDLC specifications file into common blk /mndohdlc/
!
! Generate the mapping between atoms in MNDO and atoms in HDLCopt according to
! the residue memberships and the list of atoms to be frozen, satisfying the
! specifications of HDLCopt
!------------------------------------------------------------------------------

subroutine rdhdlc (hdlcsp, nconstr, numat, natopt, nrad)
  use limit, only: LM1, LMV
  implicit none

! args
  integer hdlcsp, nconstr, numat, natopt, nrad

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! MNDO common blocks
  integer IN2, LABELS, LOC, NA, NATOMS, NB, NC, NVAR3, NBF, NB5, NB6
  real (kind=8) GEO
  common /INOPT2/ IN2(300)
  common /NBFILE/ NBF(20)
  common /PARM1 / GEO(3,LM1),NC(LM1),NB(LM1),NA(LM1),LABELS(LM1),NATOMS
  common /PARM3 / LOC(LMV),NVAR3

! MNDO externals
  real (kind=8) READA
  external READA

! local vars
  character*20 str
  logical lfrz
  integer i, i_hdlc, i_mndo, ibottom, ielem, igeom, ihdlc1, ihdlc2, ihdlc3
  integer iprint, ires, iret, j, jop, jprint, k, ksym, l 
  integer nelem, ngeom2, nrcore, nres, nvcore
  integer, dimension(:), allocatable :: imembr, resn_tmp, map_tmp
  integer, dimension(:), pointer :: irad_tmp
  real (kind=8) v
  real (kind=8), dimension(:), pointer :: vrad_tmp

  NB5 = NBF(5)
  NB6 = NBF(6)

! begin, init
  jop = IN2(3)
  igeom  = IN2(4)
  iprint = IN2(38)
  jprint = IN2(42)
  ihdlc1 = IN2(51)
  ihdlc2 = IN2(52)
  ihdlc3 = IN2(53)
  ksym = IN2(75)
  ires = 0
  nres = 0
  nconstr = 0
  ibottom = 1
  nrcore = 0
  nvcore = 0
  natopt = numat
  mndogl%lgtres = .false.; mndogl%lgtcns = .false.; mndogl%lgtfrz = .false.
  mndogl%matt = 0; mndogl%nrem = 0; mndogl%nvar = 0; mndogl%nrestart = 0
  mndogl%cfac = 0.0_8

! one entry is always allocated due to a bug in SGI f90
  allocate (mndogl%icons(4,1)); allocate (mndogl%vcons(1))
  allocate (mndogl%irad(1)); allocate (mndogl%vrad(1))

! init the default membership list and the map
  allocate (mndogl%resn(numat))
  allocate (mndogl%map(numat))
  do i = 1,numat
     mndogl%resn(i) = 0
     mndogl%map(i) = 0
  end do

! default options, constraints derived from MNDO input
  if (ihdlc3.eq.2) goto 2

! loop over entries in the file (not lines)
  if (hdlcsp.ne.NB5) then
     open (unit=hdlcsp, status='OLD', err=2)
  endif
  do while (.true.)
     read (unit=hdlcsp,fmt=*,err=8,end=1) str
     select case (str(1:6))
     case ('NFCART')
        nelem = int (READA (str, 1, 1))
        allocate (imembr(nelem))
        read (unit=hdlcsp,fmt=*,err=8,end=8) (imembr(i), i=1,nelem)
        mndogl%lgtfrz = .true.
        do ielem = 1,nelem
           i = imembr(ielem)
           mndogl%map(i) = natopt
           mndogl%resn(natopt) = -1
           natopt = natopt - 1
        end do
        deallocate (imembr)
     case ('NFCOMP')
        nelem = int (READA (str, 1, 1))
        do ielem = 1,nelem
           call get_cns_vals (hdlcsp, i, j, k, l, v, 2, iret)
           if (iret.eq.-1) goto 8
           k = 0; l = 0
           i = 3*(i-1) + j 
           j = 0
           call cns_ci (nconstr, i, j, k, l, v)
        end do
     case ('NFBOND')
        nelem = int (READA (str, 1, 1))
        do ielem = 1,nelem
           call get_cns_vals (hdlcsp, i, j, k, l, v, 2, iret)
           if (iret.eq.-1) goto 8
           k = 0; l = 0
           call cns_ci (nconstr, i, j, k, l, v)
        end do
     case ('NFANGL')
        nelem = int (READA (str, 1, 1))
        do ielem = 1,nelem
           call get_cns_vals (hdlcsp, i, j, k, l, v, 3, iret)
           if (iret.eq.-1) goto 8
           l = 0
           call cns_ci (nconstr, i, j, k, l, v)
        end do
     case ('NFDIHE')
        nelem = int (READA (str, 1, 1))
        do ielem = 1,nelem
           call get_cns_vals (hdlcsp, i, j, k, l, v, 4, iret)
           if (iret.eq.-1) goto 8
           call cns_ci (nconstr, i, j, k, l, v)
        end do
     case ('NMRESI')
        mndogl%lgtres = .true.
        ires = ires + 1
        nelem = int (READA (str, 1, 1))
        allocate (imembr(nelem))
        read (unit=hdlcsp,fmt=*,err=8,end=8) (imembr(i), i=1,nelem)
        do i = 1,nelem
           mndogl%map(imembr(i)) = ibottom
           mndogl%resn(ibottom) = ires
           ibottom = ibottom + 1
        end do
        deallocate (imembr)
     case ('NRCORE')
        nrcore = int (READA (str, 1, 1))
     case ('NDFCOR')
        mndogl%nvar = int (READA (str, 1, 1))
     case ('XYZFAC')
        mndogl%cfac = READA (str, 1, 1)
     case ('NREMST')
        mndogl%nrem = int (READA (str, 1, 1))
     case ('MDCRAD')
        mndogl%matt = int (READA (str, 1, 1))
     case ('NCRAD')
        nelem = int (READA (str, 1, 1))
        do ielem = 1,nelem
           nrad = nrad + 1
           if (nrad .ne. 1) then
              irad_tmp => mndogl%irad; vrad_tmp => mndogl%vrad
              allocate (mndogl%irad(nrad)); allocate (mndogl%vrad(nrad))
              do l = 1,nrad-1
                 mndogl%irad(l) = irad_tmp(l); mndogl%vrad(l) = vrad_tmp(l)
              end do
              deallocate (irad_tmp); deallocate (vrad_tmp)
           end if
           read (unit=hdlcsp,fmt=*,err=8,end=8) &
                mndogl%irad(nrad), mndogl%vrad(nrad)
        end do
     case ('NRSHDL')
        mndogl%nrestart = int (READA (str, 1, 1))
     case ('      ')
        goto 1
     case default
        goto 8
     end select
  end do

! jump point on end of successful reading
1 continue
  if (hdlcsp.ne.NB5) then
     close (hdlcsp,err=8)
  endif
  ngeom2 = 1
  nres = ires
  mndogl%lgtcns = (nconstr .gt. 0)

! generate the mapping for the non-frozen atoms in Cartesians
  nelem = ibottom - 1
  do i = 1,numat
     if (mndogl%map(i).eq.0) then
        mndogl%map(i) = ibottom
        ibottom = ibottom + 1
     end if
  end do
  ibottom = ibottom - 1

! re-number atom and residue maps if core is not the first residue
  if (nrcore.gt.1) then
     allocate (map_tmp(numat))
     allocate (resn_tmp(numat))

! re-number: check in core first
     ibottom = 1
     do i_mndo = 1,numat ! loop over MNDO atoms
        i_hdlc = mndogl%map(i_mndo)
        ires = mndogl%resn(i_hdlc)
        if (ires.eq.nrcore) then
           map_tmp(i_mndo) = ibottom
           resn_tmp(ibottom) = 1
           ibottom = ibottom + 1
           nvcore = nvcore+3
        end if
     end do

! re-number: check in all other atoms and renumber residue if applicable
     do i_hdlc = 1,numat ! loop over MNDO atoms
        i_mndo = 1
        do while (mndogl%map(i_mndo).ne.i_hdlc .and. i_mndo.le.numat)
           i_mndo = i_mndo + 1  ! back transform of map
        end do
        ires = mndogl%resn(i_hdlc)
        if (ires.ne.nrcore) then
           map_tmp(i_mndo) = ibottom
           if (ires.gt.0) then
              resn_tmp(ibottom) = ires + 1
           else
              resn_tmp(ibottom) = ires
           end if
           ibottom = ibottom + 1
        end if
     end do

! re-number: check in and free temporary space
     ibottom = ibottom - 1
     do i = 1,numat
        mndogl%resn(i) = resn_tmp(i)
        mndogl%map(i) = map_tmp(i)
     end do
     deallocate (resn_tmp)
     deallocate (map_tmp)
  end if

! set default number of variables for reaction core
  if (mndogl%nvar.eq.0 .and. nrcore.gt.0) then
     mndogl%nvar = nvcore
  else
     nvcore = -1
  end if

! stop on duplicate entries
  if (ibottom.ne.natopt) then
     write (NB6,'(A,I5)') ' Number of atoms in total:         ', numat
     write (NB6,'(A,I5)') ' Number of atoms to be optimised:  ', natopt
     write (NB6,'(A,I5)') ' Number of atoms in any residue:   ', nelem
     write (NB6,'(A,I5)') ' Number of duplicate entries:      ', ibottom - natopt
     call hdlc_error (2, ' Overlapping residues', 21, 'abort', 5)
  end if

! jump to printing section
  go to 3

! jump point for default options or for error on opening input file
2 continue
  ngeom2 = 2
! set the default map
  do i = 1,numat
     mndogl%map(i) = i
  end do

! impose constraints from input geometry (k: all, l: active)
! conditions: no symmetry, no dummy atoms, same coords for input and opt
! otherwise: unconstrained optimization
  if (ksym.eq.0 .and. numat.eq.natoms) then
  k = 1; l = 1
! impose constraints from input geometry in internal coords
     if (ihdlc1.eq.0 .and. igeom.eq.0) then
        do i = 1,numat
           do j = 1,3
              if (LOC(l).eq.k) then
                 l = l + 1
              else
                 if (i.gt.1 .and. j.eq.1) then ! bond length
                    call cns_ci (nconstr, i, na(i), 0, 0, 0.0_8)
                 else if (i.gt.2 .and. j.eq.2) then ! bond angle
                    call cns_ci (nconstr, i, na(i), nb(i), 0, 0.0_8)
                 else if (i.gt.3 .and. j.eq.3) then ! dihedral angle
                    call cns_ci (nconstr, i, na(i), nb(i), nc(i), 0.0_8)
                 end if
              end if
              k = k + 1
           end do
        end do
! impose constraints from input geometry in Cartesian coords
     else if (ihdlc1.eq.1 .and. igeom.eq.1) then
        do i = 1,numat
           do j = 1,3
              if (LOC(l).eq.k) then
                 l = l + 1
              else !
                 call cns_ci (nconstr, 3*(i-1)+j, 0, 0, 0, 0.0_8)
              end if
              k = k + 1
           end do
        end do     
     end if
  end if

! printing section - general
3 continue
  if (iprint.ge.-5 .or. jprint.ge.-5) then
     write (NB6,'(//,A,/)') ' Input options for HDLC optimizer.'
     if (igeom.eq.0) then
        write (NB6,'(A)') ' Original MNDO input in internal coordinates.'
     else if (igeom.eq.1) then
        write (NB6,'(A)') ' Original MNDO input in Cartesian coordinates.'
     end if
     if (nres.gt.1 .or. (jop.ne.1 .and. jop.ne.4 .and. ihdlc2.lt.3)) then
        if (ihdlc1.eq.0) then
           write (NB6,'(A)') ' Energy minimization using HDLC coordinates.'
        else if (ihdlc1.eq.1) then
           write (NB6,'(A)') ' Energy minimization using Cartesian coordinates.'
        end if
     end if
     if (jop.eq.1 .or. jop.eq.4 .or. mndogl%nvar.gt.0) then
        if (ihdlc1.eq.1 .or. ihdlc2.eq.2) then
           write (NB6,'(A)') ' Transition state search using Cartesian coordinates.'
        else if (ihdlc2.eq.0) then
           write (NB6,'(A)') ' Transition state search using HDLC primitive coordinates.'
        else if (ihdlc2.eq.1) then
           write (NB6,'(A)') ' Transition state search using HDLC total connection scheme.'
        end if
     end if
     if (nres.le.1 .and. ihdlc2.eq.3 .and. ihdlc1.eq.0) then
        write (NB6,'(A)') ' One fragment: Optimization using DLC primitive coordinates.'
     end if
     if (nres.le.1 .and. ihdlc2.eq.4 .and. ihdlc1.eq.0) then
        write (NB6,'(A)') ' One fragment: Optimization using DLC connection scheme.'
     end if
     if (ksym.ne.0) then
        write (NB6,'(A)') ' Symmetry constraints are not imposed.'
        write (NB6,'(A)') ' Please use other optimizer if symmetry is essential.'
     end if
     if (numat.ne.natoms) then
        write (NB6,'(A)') ' Dummy atoms are removed.'
     end if
     write (NB6,'(1X)')
  end if

! printing section - HDLC data read from file
  if ((iprint.ge.-5 .or. jprint.ge.-5) .and. ngeom2.eq.1) then
     write (NB6,'(A,I5)') ' Specific HDLC data read from file:', hdlcsp
     write (NB6,'(A,I5)') ' Number of atoms in total:         ', numat
     write (NB6,'(A,I5)') ' Number of atoms to be optimised:  ', natopt
     if (numat.ne.natopt) then
        write (NB6,'(A,I5)') ' Number of atoms to be frozen:     ', numat-natopt
     end if
     if (nres.gt.0) then
        write (NB6,'(A,I5)') ' Number of residues defined:       ', nres
        write (NB6,'(A,I5)') ' Number of atoms in residues:      ', nelem
     end if
     if (nrcore.gt.1) then
        write (NB6,'(A,I5)') ' Reaction core defined as residue: ', nrcore
     end if 
     if (mndogl%nvar.gt.0) then
        write (NB6,'(A,I5)') ' Number of variables in core:      ', mndogl%nvar 
        if (nvcore.gt.0 .and. nconstr.gt.0) then
           write (NB6,'(A)') ' Preceding value assigned by default.'
           write (NB6,'(A)') ' Please check for constraints in core.'
        end if
     end if
     if (nconstr.gt.0) then
        write (NB6,'(A,I5)') ' Number of constraints imposed:    ', nconstr
     end if
     if (mndogl%matt.gt.0) then
        write (NB6,'(A,I5)') ' Connectivity option for radii:    ', mndogl%matt 
     end if
     if (nrad.gt.0) then
        write (NB6,'(A,I5)') ' Number of input covalent radii:   ', nrad
     end if
     if (mndogl%nrem.gt.0) then
        write (NB6,'(A,I5)') ' Number of remembered L-BFGS steps:', mndogl%nrem 
     else
        i = min(40,3*natopt)
        write (NB6,'(A,I5)') ' Number of remembered L-BFGS steps:', i
     end if
     if (mndogl%nrestart.gt.0) then
        write (NB6,'(A,I5)') ' Number of cycles between restarts:', mndogl%nrestart
     end if
     if (ihdlc1.eq.0) then
        if (mndogl%cfac.gt.0.0_8) then
           v = mndogl%cfac
        else if (mndogl%cfac.gt.0.0_8) then
           v = - mndogl%cfac / (3*natopt)
        else
           v =  1.0_8 / (3*natopt)
        end if
     write (NB6,'(A,F11.5)') ' Weighting factor for Cartesians:  ', v
     end if
  end if
! printing section - mappings, residues, constraints, radii
  if ((iprint.ge.2 .or. jprint.ge.2) .and. ngeom2.eq.1) then
     write (NB6,'(/,A)') ' Mapping between the atoms in MNDO and HDLCopt:'
     do i = 1,numat
        write (NB6,'(I5,A,I5)') i, ' - ', mndogl%map(i)
     end do
     write (NB6,'(/,A)') ' Residue membership of the HDLCopt atoms:'
     do i = 1,numat
        write (NB6,'(I5,A,I5)') i, ' - ', mndogl%resn(i)
     end do
     if (nconstr.gt.0) then
        write (NB6,'(/,A)') ' Constraints (Angstrom or degree):'
        do i = 1,nconstr
           write (NB6,'(A,I3,A,4I5,A,F10.4)') ' Constraint ', i, ': ', &
                (mndogl%icons(j,i), j=1,4), ', v: ', mndogl%vcons(i)
        end do
     end if
     if (nrad.gt.0) then
        write (NB6,'(/,A)') ' Input covalent radii (Angstrom):'
        do i = 1,nrad
           write (NB6,'(A,I3,A,F10.4)') ' Atomic number', &
                 mndogl%irad(i), ', radius: ', mndogl%vrad(i)
        end do
     end if
  end if

! printing section - default options, no explicit input
  if ((iprint.ge.0 .or. jprint.ge.0) .and. ngeom2.eq.2) then
     write (NB6,'(A)') ' Using default options for HDLC optimizer.'
     if (ihdlc3.ne.2) then
        write (NB6,'(A,I2,A)') ' No HDLC input found on file ', hdlcsp,'.'
     else
        write (NB6,'(A)') ' No explicit HDLC input requested.'
     end if
     write (NB6,'(A)') ' No residues are defined.'
     write (NB6,'(A)') ' No reaction core is specified.'
     write (NB6,'(A)') ' No atoms are frozen.'
     if (nconstr.eq.0) then
        write (NB6,'(A)') ' No constraints are defined.'
     end if
     if (jop.ne.1 .and. jop.ne.4) then
        write (NB6,'(A,I5)') &
             ' Number of remembered L-BFGS steps:', min(40,3*numat)
     end if
     if (ihdlc1.eq.0 .and. ihdlc2.lt.3) then
        v = 1.0_8 / (3*natopt)
        write (NB6,'(A,F11.5)') ' Weighting factor for Cartesians:  ', v
     end if
     write (NB6,'(1X)')
! printing section - constraints
     if (nconstr.gt.0) then
      write (NB6,'(A)') ' Constraints derived from initial input geometry.'
      if (ihdlc1.eq.0 .and. igeom.eq.0) then
        write (NB6,'(A)') ' Fixed internal coordinates found.'
        write (NB6,'(A,I4)') ' Number of such implicit constraints: ', nconstr
        write (NB6,'(//,A)') ' List of constraints (labels of atoms involved):'
        do i = 1,nconstr
           write (NB6,'(A,I3,A,4I5,A)') ' Constraint ', i, ': ', &
                (mndogl%icons(j,i), j=1,4), ', initial input value used'
        end do
      else if (ihdlc1.eq.1 .and. igeom.eq.1) then
        write (NB6,'(A)') ' Fixed Cartesian coordinates found.'
        write (NB6,'(A,I4)') ' Number of such implicit constraints: ', nconstr
        write (NB6,'(//,A)') ' List of constraints:'
        do i = 1,nconstr
           k = (mndogl%icons(1,i)-1)/3 + 1
           l = mndogl%icons(1,i) - 3*(k-1)
           write (NB6,'(A,I3,A,I5,A,I2,A)') ' Constraint ', i, ':  atom', k, &
                ', coordinate', l,', initial input value used'
        end do
      end if
     end if
  end if

! return on successful reading and/or initialization
  return

! jump point on error
8 continue
  write (NB6,'(A,/,A)') ' Error when reading HDLC input, current string: ', str
  call hdlc_error (2, ' Invalid input', 14, 'abort', 5)

end subroutine rdhdlc

!------------------------------------------------------------------------------
! Helper routine for HDLOPT / rdhdlc:
! check in a constraint into /mndohdlc/ and do dynamic memory stuff
!------------------------------------------------------------------------------

subroutine cns_ci (nconstr, i, j, k, l, v)
  implicit none

! args
  integer i, j, k, l, nconstr
  integer, dimension(:,:), pointer :: icons
  real (kind=8) v
  real (kind=8), dimension(:), pointer :: vcons

! inter-subroutine communication
  type mndoglobal
     sequence
     logical lguess, lgtres, lgtcns, lgtfrz
     integer icall, hdlcsp
     integer nrem, matt, nvar, nrestart
     real (kind=8) :: cfac, tsum
     integer, dimension(:), pointer :: resn, map, irad
     integer, dimension(:,:), pointer :: icons
     real (kind=8), dimension(:), pointer :: vcons, vrad
  end type mndoglobal
  type(mndoglobal) :: mndogl
  common /mndohdlc/ mndogl

! local vars
  integer ii, jj
  integer, dimension(:,:), pointer :: icons_tmp
  real (kind=8), dimension(:), pointer :: vcons_tmp

! begin, grow space and work around a problem of SGI F90
  nconstr = nconstr + 1
  if (nconstr .ne. 1) then
     icons_tmp => mndogl%icons
     vcons_tmp => mndogl%vcons
     allocate(mndogl%icons(4,nconstr))
     allocate(mndogl%vcons(nconstr))
     do ii = 1,nconstr-1
        do jj = 1,4
           mndogl%icons(jj,ii) = icons_tmp(jj,ii)
        end do
        mndogl%vcons(ii) = vcons_tmp(ii)
     end do
     deallocate (icons_tmp)
     deallocate (vcons_tmp)
  end if

! check in new entry
  mndogl%icons(1,nconstr) = i; mndogl%icons(2,nconstr) = j
  mndogl%icons(3,nconstr) = k; mndogl%icons(4,nconstr) = l; 
  mndogl%vcons(nconstr) = v
end subroutine cns_ci

!------------------------------------------------------------------------------
! Helper routine for HDLOPT / rdhdlc:
! read a line with variable number of fields
!------------------------------------------------------------------------------

subroutine get_cns_vals (iunit, i, j, k, l, v, num, iret)
  implicit none

! args
  integer iunit, i, j, k, l, num, iret
  real (kind=8) v

! local params
  integer length
  parameter (length = 60)

! local vars
  character*(length) str, substr
  integer pos

! begin, read a line from input
  read (unit=iunit,fmt='(A60)',err=9,end=9) str
  pos = 1

! cycle over mandatory fields
  if (num.ge.1) then
     call get_next_val (str, substr, iret, pos, length)
     if (iret.lt.0) goto 9
     read (substr,err=9,fmt='(I5)') i
  end if
  if (num.ge.2) then
     call get_next_val (str, substr, iret, pos, length)
     if (iret.lt.0) goto 9
     read (substr,err=9,fmt='(I5)') j
  end if
  if (num.ge.3) then
     call get_next_val (str, substr, iret, pos, length)
     if (iret.lt.0) goto 9
     read (substr,err=9,fmt='(I5)') k
  end if
  if (num.ge.4) then
     call get_next_val (str, substr, iret, pos, length)
     if (iret.lt.0) goto 9
     read (substr,err=9,fmt='(I5)') l
  end if

! evaluate optional field and return
  call get_next_val (str, substr, iret, pos, length)
  if (iret.lt.0) then
     v = 0.0_8
     iret = 0
  else
     read (substr,err=9,fmt='(F10.6)') v
  end if
  return

! error condition
9 continue
  iret = -1
end subroutine get_cns_vals

!------------------------------------------------------------------------------
! Helper routine for HDLOPT / rdhdlc / get_cns_vals:
! return a field delimited by blanks or tabs
!------------------------------------------------------------------------------

subroutine get_next_val (str, substr, iret, pos, length)
  implicit none

! args
  character*(*) str, substr
  integer iret, pos, length

! local params: tab character
  integer itab, npos
  parameter (itab = 9)

! begin, advance until next entry
  do while (pos.le.length .and. &
            (str(pos:pos).eq.' '.or.str(pos:pos).eq.char(itab)))
     pos = pos + 1
  end do
  if (pos.gt.length) then
     iret = -1
     return
  else
     iret = 0
  end if

! advance until end of entry
  npos = pos
  do while (npos.le.length .and. &
            (str(npos:npos).ne.' '.and.str(npos:npos).ne.char(itab)))
     npos = npos + 1
  end do

! return the entry and new position
  substr = str(pos:npos-1)
  pos = npos
end subroutine get_next_val
