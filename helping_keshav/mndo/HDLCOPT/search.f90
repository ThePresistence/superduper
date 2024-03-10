!------------------------------------------------------------------------------
!
! THE HDLC OPTIMISER
! ==================
! A geometry optimiser for linear scaling minimum and transition state search
! in hybrid delocalised internal coordindates or Cartesian coordinates
!
! Contact:
! ========
! Salomon Billeter: salomon@chem.psu.edu
! Paul Sherwood:    P.Sherwood@dl.ac.uk
! Walter Thiel:     thiel@mpi-muelheim.mpg.de
!
! Documentation:
! ==============
! - the QUASI project:    http://www.dl.ac.uk/TCSC/QuantumChem/quasi_www/manual/hdlcopt.html
! - printed:              salomon@chem.psu.edu
!
! History:
! ========
! Initial version: 1997, 1998 A.J. Turner, Zurich
! Completed / rewritten: 1998, 1999 S.R. Billeter, Zurich
!
!------------------------------------------------------------------------------
!
! hdlcopt
! =======
!
! Main routine of HDLCopt, not contained in a module in order to be able to be
! called either from Fortran77, Fortran90 or C.
!
! Arguments:
! ==========
! hdlcopt (iccode, nattot, natoms, nincon, nconstr, nrad, iret, linit, nb6, &
!         mndoscf, mndohdc, mndoarr, mndolm)
!
! iccode  (in)  calling code (1 -> CHEMSH, 2 -> MNDO, 3 -> CHARMM)
! nattot  (in)  total number of atoms (including frozen), only for reporting
! natoms  (in)  number of atoms to be optimised
! nincon  (in)  number of user connections (bonds between atoms)
! nconstr (in)  number of constraints (bonds, angles, dihedrals)
! nrad    (in)  number of user covalent radii
! iret    (out) return code (see below)
! linit   (in)  .false. if restart from checkpoint is forced - sets mstart=4
! nb6     (in)  Fortran I/O unit of standard output (usually 6)
! mndoscf (in)  external subroutine used for energy evaluation (MNDO only)
! mndohdc (in)  external subroutine used for gradient evaluation (MNDO only)
! mndoarr (in)  scratch array (MNDO only)
! mndolm5 (in)  limit of scratch array (MNDO only)
!
! Return codes of hdlcopt:
! =======================
!  1: converged and optimised
!  0: not (yet) converged, obsolete
! -1: error: general optimisation or HDLC problem
! -2: converged but not yet optimised (microiterative optimisation), obsolete
! -3: maximum number of cycles reached
! -4: SCF error (from hdlc_get_coords)
! -5: time limit reached (detected by hdlc_get_coords if possible)
!
! Interface:
! ==========
! The following routines must be provided externally (see code of search.f90):
!
! hdlc_get_params 
!   - purpose: pass all optimisation control parameters to hdlcopt, set the
!     defaults
!   - syntax:
!     hdlc_get_params (iccode, linit, <all fields of ctrl>)
!     iccode:  (in)  calling code (see above)
!     linit:   (in)  parameters are to be read only if linit is .true.
!     ctrl:    (out) optimisation control parameters (see global.f90)
!   - note on atom sequence: the atom sequence must match the atom sequence of
!     the Cartesian coordinates, the Cartesian gradient, the Cartesian Hessian
!     and the residue numbers
!   - note on limits: the following limits are rather passed to than from
!     hdlc_get_params: nbig, nconns, nconstr
!
! hdlc_get_resn
!   - purpose: pass the residue membership of the atoms to hdlcopt
!   - syntax:
!     hdlc_get_resn (iccode, resn, natom)
!     iccode:  (in)  calling code (see above)
!     resn:    (out) integer array(natom) of the residue numbers of all atoms
!     natom:   (in)  number of atoms
!   - rules:
!     * atoms not being part of a residue (to be represented in Cartesians)
!       have the residue number 0
!     * the residues are ordered ascending
!     * the residues must only span subsequent atoms (exception: residue 0)
!   - example: 0, 0, 1, 1, 1, 0, 0, 2, 2, 2, 2, 3, 3, 3, 0, 0
!
! hdlc_get_coords
!   - purpose: calculate energy and gradient - pass coordinates, energy and
!     gradient to hdlcopt
!   - syntax:
!     hdlc_get_coords (iccode, coords, funct, grad, nbig, ierr, &
!                      mndoscf, mndoarr, mndolm5)
!     iccode:  (in)  calling code (see above)
!     coords:  (out) array(nbig) of the all Cartesian coordinates
!     funct:   (out) energy function value
!     grad:    (out) array(nbig) of the Cartesian gradient
!     nbig:    (in)  number of degrees of freedom
!     ierr:    (out) return code (see below)
!     mndoscf: (in)  external SCF subroutine - MNDO builtin version only
!     mndoarr: (in)  scratch array - MNDO builtin version only
!     mndolm5: (in)  limit of scratch array - MNDO builtin version only
!     jump:    (in)  flag to specify current task - for printing only
!   - order of the components: x(atom1), y(atom1), z(atom1), x(atom2), ...
!   - return codes ierr from hdlc_get_coords (<0: MNDO style):
!        1 CHEMSH: generic error
!        0 no error
!       -1 no SCF convergence
!       -8 time limit reached
!
! hdlc_get_hess
!   - purpose: pass the analytic Cartesian Hessian of the first nvar degrees of
!     freedom to hdlcopt
!   - syntax:
!     hdlc_get_hess (iccode, hess, nvar, lvalid, lcart)
!     iccode:  (in)  calling code (see above)
!     hess:    (out) array(nvar,nvar) of the Cartesian Hessian
!     nvar:    (in)  number of degrees of freedom
!     lvalid:  (out) .true. if a valid Hessian is returned in hess(:,:)
!     lcart:   (out) .true. in any case
!   - note: this routine may be a dummy routine - it is only called if nvar>0
!     and mhess=1
!
! hdlc_update
!   - purpose: if energy / gradient comes from a force field, the MM code may
!     update the pair list on hdlc_update
!   - syntax:
!     hdlc_update (iccode)
!     iccode:  (in)  calling code (see above)
!
! hdlc_put_coords
!   - purpose: pass the new geometry to the energy / gradient code
!   - syntax:
!     hdlc_put_coords (iccode, coords, nbig)
!     iccode:  (in)  calling code (see above)
!     coords:  (in)  array(nbig) of the Cartesian coordinates
!     nbig:    (in)  number of degrees of freedom
!   - note: the configuration is expected to be passed back to hdlcopt by the
!     subsequent call of hdlc_get_coords
!
! hdlc_wr_fandg / hdlc_rd_fandg
!   - purpose: tell the function / gradient code to dump / undump its memory
!   - syntax:
!     hdlc_wr_fandg (iccode, chkpnt, lform, lerr)
!     hdlc_rd_fandg (iccode, chkpnt, lform, lerr)
!     iccode:  (in)  calling code (see above)
!     chkpnt:  (in)  Fortran unit of the already opened checkpoint file
!     lform:   (in)  .true. if formatted I/O
!     lerr:    (out) .true. if I/O failed
!   - notes:
!     * chkpnt and lform do not need to be observed, lerr is mandatory
!     * the energy / gradient code is expected to handle I/O error on its own
!       for its internal purpose
!
! hdlc_error
!   - purpose: exit procedure on failure
!   - syntax:
!     hdlc_error (iccode, message, action)
!     iccode:  (in)  calling code (see above)
!     message: (in)  message, can be ignored
!     action:  (in)  action: 'warn', 'stop' or 'abort'
!------------------------------------------------------------------------------

subroutine hdlcopt (iccode, nattot, natoms, nincon, nconstr, nrad, &
     iret, linit, nb6, &
     mndoscf, mndoarr, mndolm5)

! includes
  use global
  use hdlclib
  use matrixlib
  use searchlib
  implicit none

! args
  logical linit
  integer iccode, iret, mndolm5, natoms, nattot, nb6, nconstr, nincon, nrad
  real(kind=8), dimension(mndolm5) :: mndoarr
  external mndoscf

! local vars
  logical lstart
  logical interror
  integer group, i, iatom, idum, ierr, ifin, ipcart, iphdlc, istart, j, &
       k, length, m, ncart, ndfhdlc, nconn, ndfcons, ndfopt, &
       ngroups, ngroupsdrop, nmin
  integer, pointer, dimension(:,:) :: iconn
  integer, allocatable, dimension(:) :: resn, err_cnt
  real (kind=8) :: funct
  real (kind=8), dimension(:), pointer :: coords, grad, coords_tmp, grad_tmp, &
       cgrad, prim_tmp
  real (kind=8), dimension(:,:), pointer :: hess
  type(int_matrix), pointer :: con
  type(matrix), pointer :: cns, xyz, cxyz, chdlc, gxyz, ghdlc
  type(hdlc_obj), pointer :: residue

! init the control registers - see also below hdlc_get_params
  if (ctrl%printl.ge.1) then
     write (stdout,'(/,a)') '*---------------------------------------*'
     write (stdout,'(a)')   '| Entering the HDLC optimisation driver |'
     write (stdout,'(a)')   '|          Ver. 1.0 - AJT / SB          |'
     write (stdout,'(a)')   '|          Dec 1997 - Oct 1999          |'
     write (stdout,'(a,/)') '*---------------------------------------*'
  end if
  stat%lgeook = .true.
  stat%ncyc = 0; stat%nstep = 0; stat%nsprfo = 0; stat%nstried = 0
  stat%nslbfg = 0; stat%nshess = 0
  stat%jump = 0
  hdlc%internal = .false.
  hdlc%ngroups  = 0
  nullify (hdlc%first); nullify (hdlc%last)
  nullify (con); nullify (cns)
  ost%lconv = .false.
  ost%lend = .false.

!//////////////////////////////////////////////////////////////////////////////
! begin, initialise the optimisation
!//////////////////////////////////////////////////////////////////////////////

! get optimisation parameters
  ctrl%nbig = 3 * natoms
  ctrl%nconstr = nconstr
  ctrl%nincon = nincon
  ctrl%nrad = nrad
  ctrl%iccode = iccode
  allocate (ctrl%attypes(natoms))
  allocate (ctrl%precon(ctrl%nbig))
  allocate (ctrl%incon(2,nincon))
  allocate (ctrl%iconstr(4,nconstr))
  allocate (ctrl%vconstr(nconstr))
  allocate (ctrl%irad(nrad))
  allocate (ctrl%vrad(nrad))
  call hdlc_get_params (ctrl%iccode, linit, &
       ctrl%mstart, ctrl%nbig, ctrl%nvar, ctrl%nrem, &
       ctrl%nrestart, ctrl%attypes, &
       ctrl%toler, ctrl%rmin, ctrl%rmax, ctrl%dmax, ctrl%ddmax, ctrl%osmin, &
       ctrl%omin, &
       ctrl%mode, ctrl%irecalc, ctrl%maxcyc, ctrl%iupd, ctrl%mhess, &
       ctrl%idump, ctrl%mdump, ctrl%iprj, ctrl%ilock, &
       ctrl%precon, &
       ctrl%nincon, ctrl%mconstr, ctrl%nconstr, ctrl%contyp, ctrl%ctfirst, &
       ctrl%incon, ctrl%iconstr, ctrl%vconstr, ctrl%cfact, ctrl%printl, &
       ctrl%nrad, ctrl%irad, ctrl%vrad, ctrl%icconv)
  stat%lcart_hess = .false.
  stat%lval_hess = (ctrl%nvar .eq. 0)
  interror = .false.

! check for some inconsistency in the input
  if (ctrl%nvar .gt. ctrl%nbig) then
     call hdlc_errflag ('nvar may not be larger than nbig', 'stop')
     goto 9999
  end if

! allocate memory
  allocate (coords(ctrl%nbig)); allocate (grad(ctrl%nbig))
  allocate (coords_tmp(ctrl%nbig)); allocate (grad_tmp(ctrl%nbig))
  allocate (cgrad(ctrl%nbig))
  allocate (hess(ctrl%nvar,ctrl%nvar))
  allocate (resn(natoms))
  call allocate_work (owr, ctrl%nbig, ctrl%nvar, ctrl%nrem)

! read in optimiser state from checkpoint file according to ctrl%mstart
  call hdlc_rd_chk (ctrl%mstart, ctrl%nbig, ctrl%nvar, ctrl%nrem, &
       ctrl, stat, hdlc, ost, owr, coords, hess, resn)
  lstart = (ctrl%mstart .le. 2)
  if (.not. linit) stat%ncyc = 0

! read in residue membership information
  if (ctrl%mstart .le. 1) then
     call hdlc_get_resn (ctrl%iccode, resn, natoms)
     if (ctrl%printl .ge. 2) then
        write (stdout,'(a)') 'Residue memberships:'
        write (stdout,'(3x,15i5)') (resn(i), i=1,natoms)
     end if

! remove the first residue if ctfirst.eq.2
     j = 0
     ngroups = 0
     if (ctrl%ctfirst .eq. 2) then
        if (ctrl%printl.ge.1) then
           write (stdout,'(a,/)') 'Deleting first residue (ctfirst=2)'
        end if
        do i = 1,natoms
           if (resn(i).ne.0) then
              if (resn(i).ne.j) then
                 ngroups = ngroups + 1
                 j = resn(i)
              endif
              if (ngroups.eq.1) then
                 resn(i) = 0
              end if
           end if
        end do
     end if ! (ctrl%ctfirst .eq. 2)
  end if ! (ctrl%mstart .le. 1)

! initialise the Hessian: unit or external (if not yet there and PRFO)
  if (lstart .and. .not.stat%lval_hess) then
     if (ctrl%mhess .eq. 1) then
        call hdlc_get_hess (ctrl%iccode, hess, ctrl%nvar, stat%lval_hess, &
             stat%lcart_hess)
     end if
  end if
  if (.not. stat%lval_hess) then
     do i = 1,ctrl%nvar
        do j = 1,ctrl%nvar
           hess(i,j) = 0.0_8
        end do
        hess(i,i) = 1.0_8
     end do
  end if

! set the initial value of the state machine of the optimiser
  if (lstart) then
     stat%jump = 0
  end if

!//////////////////////////////////////////////////////////////////////////////
! main loop over the iterations
!//////////////////////////////////////////////////////////////////////////////

  iret = 0
  do while (stat%ncyc.lt.ctrl%maxcyc .and. iret.eq.0 .and. (.not.ost%lend))
     stat%ncyc = stat%ncyc + 1

! avoid L-BFGS if PRFO only
     if (ctrl%nbig .eq. ctrl%nvar .and. stat%jump .eq. 2) stat%jump = 5

! print what the code is going to do this cycle (see module global)
     if (ctrl%printl.ge.0) then
        write (stdout,'(/,a,i6,a,i1,a,a)') 'Cycle ', stat%ncyc, ', jump=', &
             stat%jump, ': ', flwmsg(stat%jump)
     end if

! ngroups: number of HDLC residues in resn (hdlc%ngroups: number of objects)
     j = 0
     ncart = 0
     ngroups = 0
     do i = 1,natoms
        if (resn(i).le.0) then
           ncart = ncart + 1
        elseif (resn(i).ne.j) then
           ngroups = ngroups + 1
           if (resn(i).lt.j) then
              call hdlc_errflag ('Residue names must be in ascending order', &
                   'stop')
              goto 9999
           end if
           j = resn(i)
        end if
     end do

! switch hdlc%internal is set if internals only (DLC)
    hdlc%internal = (ngroups.eq.1 .and. ncart.eq.0 .and. ctrl%ctfirst.ge.3 &
         .and. natoms.eq.nattot)
    if (hdlc%internal) then
       nmin = 6
    else
       nmin = 0
    end if
     hdlc%lhdlc = (ngroups .ge. 1)

! allocate space for the HDLC failure counters
     if (stat%ncyc.eq.1) allocate (err_cnt(ngroups))
     do group = 1,ngroups
        err_cnt(group) = 0
     end do
     if (ctrl%printl.ge.2 .or. (ctrl%printl.ge.0 .and. stat%ncyc.eq.1)) then
        if (ngroups.eq.0) then
           write (stdout,'(a,/)') &
                'The entire system is represented in Cartesian coordinates'
        elseif (hdlc%internal) then
           write (stdout,'(a,/)') &
                'The system is represented in pure internal coordinates (DLC)'
        else
           write (stdout,'(a,i4,a,/)') &
                'The system has ', ngroups, ' HDLC residues'
        end if
     end if

! restart the optimiser and update the pair list if required
     if (ctrl%nrestart .ne. 0) then
        if (mod(stat%ncyc,ctrl%nrestart) .eq. 0) then
           if (ctrl%printl.ge.0) write (stdout,'(a,/)') &
                'Restarting the optimiser by user request'
           call hdlc_update (ctrl%iccode)
           stat%jump = -1
           call hdlc_destroy_all (.true., err_cnt)
           stat%lgeook = .false.
        end if
     end if

! get the coordinates and calculate energy / gradient
     call hdlc_get_coords (ctrl%iccode, coords, funct, grad, ctrl%nbig, &
          ierr, mndoscf, mndoarr, mndolm5, stat%jump)
     if (ierr.ne.0) goto 9999

! precondition the gradient
     do i = 1,ctrl%nbig
        grad(i) = grad(i)*ctrl%precon(i)
        cgrad(i) = grad(i)
     end do

!//////////////////////////////////////////////////////////////////////////////
! generate HDLC
!//////////////////////////////////////////////////////////////////////////////

! jump point in case of damaged HDLC
     ndfopt = ctrl%nbig
10   if (ngroups.gt.0) then
        if ((.not.stat%lgeook) .or. lstart) then
           stat%lgeook = .true.
           if (ctrl%printl.ge.1) write (stdout,'(/,a,i4,a)') &
                'Generating new HDLC for ', ngroups, ' residues'

! loop over all atoms to check residue memberships
           iatom = 0
           group = 0
50         do while (iatom.lt.natoms)
              iatom = iatom + 1
              if (resn(iatom).eq.0) then
                 goto 50
              else if (resn(iatom).ne.group) then
                 group = resn(iatom)
                 istart = iatom
              else
                 goto 50
              end if
              ifin = istart
              do while (resn(ifin).eq.group .and. ifin.lt.natoms)
                 ifin = ifin + 1
              end do
              if (ifin.ne.natoms) ifin = ifin - 1
              iatom = ifin
              if (group.eq.-1) goto 50
              length = ifin - istart + 1

! now we have group, istart, ifin
              if (ctrl%printl.ge.1) then
                 write (stdout,'(/,a)') 'Located a new residue for HDLC'
                 write (stdout,'(3x,a,i4)') 'Residue number is ', group
                 write (stdout,'(3x,a,i5)') 'It starts at atom ', istart
                 write (stdout,'(3x,a,i5,/)') 'It ends at atom ', ifin
              end if
              if (length.lt.2) then
                 call hdlc_errflag &
                      ('No residue can contain less than two atoms', 'stop')
                 goto 9999
              end if
              xyz => matrix_create (3*length, 1, 'XYZ')
              call dummyarg_checkin (coords, 3*(istart-1)+1, 3*length)
              idum = matrix_set (xyz, dummyarg)
              call dummyarg_clear
              
! warn if the residue crosses boundary between core and environment
              if (ctrl%nvar.ne.0 .and. &
                   istart.le.ctrl%nvar/3 .and. ifin.gt.ctrl%nvar/3) then
                 if (ctrl%printl.ge.-1) then
                    write (stdout,'(/,A,I4,A)') 'Warning: residue ', group, &
                         ' crosses boundary between core and environment!'
                    write (stdout,'(A,I4,A,I6,A,I6)') 'nvar: ', ctrl%nvar, &
                         ', start: ', istart, ', end: ', ifin
                 end if
              end if

! get connections for primitive internals
              if (ctrl%contyp.eq.0 .and. &
                   ((group.eq.1.and.(ctrl%ctfirst.eq.0.or.ctrl%ctfirst.eq.3)) &
                   .or. group.ne.1)) then
                 call dummyarg_checkin (ctrl%attypes, istart, length)
                 call connect_prim (length, dummyint, nconn, iconn, xyz)
                 call dummyarg_clear

! check in user connections and create connections matrix
                 call ci_conn (con, nconn, iconn, ctrl%nincon, ctrl%incon, &
                      istart, ifin)
                 deallocate (iconn)

! check in constraints
                 call assign_cons (cns, ctrl%iconstr, ctrl%nconstr, &
                      ctrl%vconstr, ctrl%mconstr, istart, ifin, hdlc%internal)

! create the HDLC - primitives
                 call hdlc_create (residue, xyz, con, cns, group, istart)
                 idum = int_matrix_destroy (con)
                 idum = matrix_destroy (cns)

! create the HDLC - total connection scheme - no check if only stretch constr.
              else ! (ctrl%contyp.eq.0 ...)
                 call assign_cons (cns, ctrl%iconstr, ctrl%nconstr, &
                      ctrl%vconstr, ctrl%mconstr, istart, ifin, hdlc%internal)
                 call hdlc_create (residue, xyz, con, cns, group, istart)
                 idum = matrix_destroy (cns)
              end if ! (ctrl%contyp.eq.0 ...)

! restore the error counter
              residue%err_cnt = err_cnt(group)
              idum = matrix_destroy (xyz)
           end do ! while (iatom.lt.natoms)
        end if ! ((.not.stat%lgeook) .or. lstart)

!//////////////////////////////////////////////////////////////////////////////
! convert Cartesian -> HDLC
!//////////////////////////////////////////////////////////////////////////////

        if (ctrl%printl.ge.2) write (stdout, '(/,A)') &
             'Converting Cartesians to HDLC'
        residue => hdlc%first
        ndfcons = 0
        do group = 1,ngroups
           istart = residue%start
           iphdlc = 3*(istart-1) + 1 - ndfcons
           ndfcons = ndfcons + residue%ncons
           length = residue%natom
           ifin = istart + length - 1
           ndfhdlc = 3*length - residue%ncons

! now we have group, istart, ifin, iphdlc, ndfhdlc; convert coords to HDLC
           cxyz => matrix_create (3*length, 1, 'CXYZ')
           chdlc => matrix_create (3*length-nmin, 1, 'CHDLC')
           call dummyarg_checkin (coords, 3*(istart-1)+1, 3*length)
           idum = matrix_set (cxyz, dummyarg)
           call dummyarg_clear
           call coord_cart_to_hdlc (residue, cxyz, chdlc, prim_tmp, .false.)

! convert gradient to HDLC
           gxyz => matrix_create (3*length, 1, 'GXYZ')
           ghdlc => matrix_create (3*length-nmin, 1, 'GHDLC')
           call dummyarg_checkin (grad, 3*(istart-1)+1, 3*length)
           idum = matrix_set (gxyz, dummyarg)
           call dummyarg_clear
           call grad_cart_to_hdlc (residue, cxyz, gxyz, ghdlc)

! separate between active space and constraints - resize CHDLC and GHDLC
           if (residue%ncons .ne. 0) then
              call hdlc_split_cons (residue, chdlc, (ctrl%mconstr.eq.0))
              call hdlc_split_cons (residue, ghdlc, .false.)
           end if

! check in coordinates and gradient, size of chdlc/ghdlc now: ndfhdlc
           call dummyarg_alloc (ndfhdlc)
           idum = matrix_get (chdlc, dummyarg)
           call dummyarg_checkout (coords_tmp, iphdlc, ndfhdlc-nmin)
           call dummyarg_alloc (ndfhdlc)
           idum = matrix_get (ghdlc, dummyarg)
           call dummyarg_checkout (grad_tmp, iphdlc, ndfhdlc-nmin)

! prepare for the next group
           residue => residue%next
           idum = matrix_destroy (cxyz); idum = matrix_destroy (gxyz)
           idum = matrix_destroy (chdlc); idum = matrix_destroy (ghdlc)
        end do ! group = 1,ngroups

! check in Cartesians
        ipcart = 1
        iphdlc = 1
        j = 0
        residue => hdlc%first
        do iatom = 1,natoms
           if (resn(iatom).le.0) then
              do m = 0,2
                 coords_tmp(iphdlc+m) = coords(ipcart+m)
                 grad_tmp(iphdlc+m) = grad(ipcart+m)
              end do
           else if (resn(iatom).ne.j) then
              j = resn(iatom)
              iphdlc = iphdlc - residue%ncons
              if (associated(residue%next)) residue => residue%next
           end if ! (resn(iatom).le.0) ... else ...
           iphdlc = iphdlc + 3
           ipcart = ipcart + 3
        end do ! iatom = 1,natoms

! store the HDLC in the work arrays
        ndfopt = ctrl%nbig - ctrl%nconstr - nmin
        iphdlc = iphdlc - 1 - nmin
        if (ndfopt .ne. iphdlc) then
           write (stdout,'(A,I5,A,I5)') 'NDFOPT: ', ndfopt, ', IPHDLC: ', &
                iphdlc
           call hdlc_errflag &
                ('Wrong number of degrees of freedom checked in', 'stop')
        end if
        do m = 1,ndfopt
           coords(m) = coords_tmp(m)
           grad(m) = grad_tmp(m)
        end do

! Cartesian only: constrain Cartesian components if requested
     else if (ctrl%nconstr.gt.0) then ! (ngroups.gt.0)
        j = 1
        if (ctrl%iconstr(2,j).ne.0) then
           write (stdout,'(A,I5,A,4I5)') 'Constraint ', j, ': ', &
                (ctrl%iconstr(m,j), m=1,4)
           write (stdout,'(A,A)') 'No internal constraints allowed ', &
                'if optimisation is in Cartesians'
           call hdlc_errflag &
                ('Constraints input error', 'stop')
        end if
        ipcart = 1
        iphdlc = 1
        do i = 1,natoms
           do m = 0,2
              if (ctrl%iconstr(1,j).eq.ipcart) then ! constrained
                 coords_tmp(ctrl%nbig-j+1) = coords(ipcart)
                 grad_tmp(ctrl%nbig-j+1) = grad(ipcart)
                 if (j.lt.ctrl%nconstr) then
                    j = j + 1 ! next constraint
                 else
                    j = 1     ! all constraints done
                 end if
                 if (ctrl%iconstr(2,j).ne.0) then
                    write (stdout,'(A,I5,A,4I5)') 'Constraint ', j, ': ', &
                         (ctrl%iconstr(k,j), k=1,4)
                    write (stdout,'(A,A)') 'No internal constraints allowed ', &
                         'if optimisation is in Cartesians'
                    call hdlc_errflag &
                         ('Constraints input error', 'stop')
                 end if
              else ! not constrained
                 coords_tmp(iphdlc) = coords(ipcart)
                 grad_tmp(iphdlc) = grad(ipcart)
                 iphdlc = iphdlc + 1
              end if
              ipcart = ipcart + 1
           end do
        end do
        ndfopt = ctrl%nbig - ctrl%nconstr
        do m = 1,ctrl%nbig
           coords(m) = coords_tmp(m)
           grad(m) = grad_tmp(m)
        end do
     end if ! (ngroups.gt.0) ... else if (ctrl%nconstr.gt.0) ...

!//////////////////////////////////////////////////////////////////////////////
! Hessian if required
!//////////////////////////////////////////////////////////////////////////////

! calculate the Hessian if a restart happens
     if (stat%jump.eq.-1 .and. ctrl%nvar.gt.0) then
        if (ctrl%printl.ge.1) then
           write (stdout,'(A,/)') 'Deleting Hessian due to restart'
        end if
        stat%nshess = 0
        stat%lval_hess = .false.
     end if

! calculate the Hessian if requested
     if (ctrl%irecalc.gt.0 .and. (.not. lstart)) then
        if (mod(stat%ncyc,ctrl%irecalc) .eq. 1) then
           if (ctrl%printl.ge.1) then
              write (stdout,'(A,/)') 'Deleting Hessian by user request'
           end if
           stat%nshess = 0
           stat%lval_hess = .false.
        end if
     end if

! the Hessian is required but not ready so far
     if (stat%jump .eq. 5 .and. (.not. stat%lval_hess)) then
        if (ctrl%printl.ge.1) then
           write (stdout,'(/,a,i1,a,a)') &
                'Hessian not ready for P-RFO, jump=', &
                stat%jump, ': ', flwmsg(stat%jump)
        end if
        stat%nshess = 0
        stat%jump = 7
     end if

! Hessian: finite difference or externally calculated
     if (stat%jump .eq. 7) then
        if (ctrl%mhess .eq. 1) then
           call hdlc_get_hess (ctrl%iccode, hess, ctrl%nvar, stat%lval_hess, &
                stat%lcart_hess)
           stat%jump = 5
        end if
        if (.not. stat%lval_hess) then
           call hess_computer (coords, funct, grad, hess, ctrl%nvar, &
                stat%jump, stat%lval_hess, stat%nshess)
        end if
        if (stat%lval_hess) then              
           stat%nshess = 0
           stat%jump = 5 ! was 5 before Hessian calculation
        end if
     end if

! convert Cartesian Hessian to HDLC: not yet implemented
     if (stat%jump.ne.7 .and. stat%lval_hess .and. stat%lcart_hess) then
        call hdlc_errflag &
             ('Conversion of Cartesian Hessian to HDLC not yet implemented', &
             'stop')
        goto 9999
     end if

!//////////////////////////////////////////////////////////////////////////////
! do the optimisation step
!//////////////////////////////////////////////////////////////////////////////

     if (stat%jump .ne. 7) then
        call prfoss (coords, funct, grad, hess, cgrad, ndfopt, ctrl%nvar, &
             ctrl%nrem, stat%jump, ost%lconv, ost%lend, stat%lval_hess)
     end if

!//////////////////////////////////////////////////////////////////////////////
! convert HDLC -> Cartesian
!//////////////////////////////////////////////////////////////////////////////

     if (ngroups.gt.0) then
        if (ctrl%printl.ge.2) write (stdout, '(/,A)') &
             'Converting HDLC to Cartesians'

! store the HDLC in the work arrays
        do m = 1,ndfopt
           coords_tmp(m) = coords(m)
           grad_tmp(m) = grad(m)
        end do

! loop over all residues
        ngroupsdrop = 0
        residue => hdlc%first
        ndfcons = 0
        do group = 1,ngroups
           istart = residue%start
           iphdlc = 3*(istart-1) + 1 - ndfcons
           ndfcons = ndfcons + residue%ncons
           length = residue%natom
           ifin = istart + length - 1
           ndfhdlc = 3*length - residue%ncons

! now we have group, istart, ifin, iphdlc, ndfhdlc; check out HDLC coordinates
           chdlc => matrix_create (ndfhdlc-nmin, 1, 'CHDLC')
           cxyz => matrix_create (3*length, 1, 'CXYZ')
           call dummyarg_checkin (coords_tmp, iphdlc, ndfhdlc-nmin)
           idum = matrix_set (chdlc, dummyarg)
           call dummyarg_clear

! restore values of constrained variables if required
           if (residue%ncons .ne. 0) then
              call hdlc_rest_cons (residue, chdlc)
           end if

! fit Cartesian coordinates to HDLC coordinates; size of chdlc now: 3*length
           call coord_hdlc_to_cart (residue, cxyz, chdlc)

! conversion HDLC -> Cartesian failed due to singular G matrix
           if (.not. residue%lgmatok) then
              stat%lgeook = .false.
              stat%jump = -1
              residue%err_cnt = residue%err_cnt + 1000
              if (ctrl%printl.ge.0) write (stdout,'(3X,A,I4,A,I4,/)') &
                   'Conversion of residue ', residue%name, &
                   ' failed , HDLC failure gauge: ', residue%err_cnt

! persistent HDLC failure - remove residue if there are no constraints
              if (residue%err_cnt .ge. 2000) then
                 interror = .true.
                 if (residue%ncons .gt. 0) then
                    if (ctrl%printl.ge.-1) write (stdout,'(A,I4,A,I4)') &
                         'Warning: could not remove residue', residue%name, &
                         ', number of constraints: ', residue%ncons
                 else
                    ngroupsdrop = ngroupsdrop + 1
                    if (ctrl%printl.ge.0) write (stdout,'(A,I4,A)') &
                         'Cyclic failure - removing residue ', residue%name, &
                         ' from list'
                    do i = istart,ifin
                       resn(i) = -2
                    end do
                 end if
              end if

! conversion HDLC -> Cartesian was successful
           else
              if (stat%jump .ne. 1) then
                 residue%err_cnt = residue%err_cnt/2
                 if (ctrl%printl.ge.3) write (stdout,'(5X,A,I5,A,I3,/)') &
                      'Residue ', residue%name, ', HDLC failure gauge: ', &
                      residue%err_cnt
              end if
           end if

! check in coordinates
           call dummyarg_alloc (3*length)
           idum = matrix_get (cxyz, dummyarg)
           call dummyarg_checkout (coords, 3*(istart-1)+1, 3*length)

! prepare for the next group
           residue => residue%next
           iphdlc = iphdlc + ndfhdlc
           idum = matrix_destroy (cxyz)
           idum = matrix_destroy (chdlc)
        end do ! (group = 1,ngroups)

! check in Cartesians
        ipcart = 1
        iphdlc = 1
        j = 0
        residue => hdlc%first
        do iatom = 1,natoms
           if (resn(iatom).le.0) then
              do m = 0,2
                 coords(ipcart+m) = coords_tmp(iphdlc+m)
              end do
           else if (resn(iatom).ne.j) then
              j = resn(iatom)
              iphdlc = iphdlc - residue%ncons
              if (associated(residue%next)) residue => residue%next
           end if ! (resn(iatom).le.0) ... else ...
           iphdlc = iphdlc + 3
           ipcart = ipcart + 3
        end do ! iatom = 1,natoms

! Cartesian only: move back constrained components if requested
     else if (ctrl%nconstr.gt.0) then ! (ngroups.gt.0)
        j = 1
        ipcart = 1
        iphdlc = 1
        do i = 1,natoms
           do m = 0,2
              if (ctrl%iconstr(1,j).eq.ipcart) then ! constrained
                 coords_tmp(ipcart) = coords(ctrl%nbig-j+1)
                 if (j.lt.ctrl%nconstr) then
                    j = j + 1 ! next constraint
                 else
                    j = 1     ! all constraints done
                 end if
              else ! not constrained
                 coords_tmp(ipcart) = coords(iphdlc)
                 iphdlc = iphdlc + 1
              end if
              ipcart = ipcart + 1
           end do
        end do
        ndfopt = ctrl%nbig - ctrl%nconstr
        do m = 1,ctrl%nbig
           coords(m) = coords_tmp(m)
           grad(m) = grad_tmp(m)
        end do
     end if ! (ngroups.gt.0) ... else if (ctrl%nconstr.gt.0) ...

! if geometry died, save err_cnt and destroy all HDLC residues and restart
     if (.not. stat%lgeook) then
        call hdlc_destroy_all (.true., err_cnt)
        ngroups = ngroups - ngroupsdrop

! if a residue has been removed, its number is set to -2 to prevent checkin
        do iatom = 1,natoms
           if (resn(iatom).eq.-2) resn(iatom) = -1
        end do
        if (interror) then
           ierr = -9
           go to 9999
        endif
        goto 10
     end if

! remove rotation / translation
     if (stat%jump .ne. 7) then
        if (nattot.eq.natoms) then
!          if (ctrl%printl.ge.2) &
!               write (stdout,'(/,A,/)') 'Removing R/T motions' ! no R/T fit
        else
           if (ctrl%iprj.eq.1 .and. .not.hdlc%lhdlc) then
              call hdlc_errflag &
                   ('Frozen atoms cannot be used with Hessian projection', &
                   'abort')
           else
!             if (ctrl%printl.ge.2) &
!                  write (stdout,'(/,A,/)') 'Skipping R/T motions removal' ! do.
           end if
        end if
     end if

! send new geometry to the calling program
     call hdlc_put_coords (ctrl%iccode, coords, ctrl%nbig)

! dump memory if requested
     if (ctrl%idump .gt. 0) then
        if (mod (stat%ncyc, ctrl%idump) .eq. 0) &
             call hdlc_wr_chk (ctrl%mdump, ctrl%nbig, ctrl%nvar, ctrl%nrem, &
             ctrl, stat, hdlc, ost, owr, coords, hess, resn)
     end if

! prepare for the next cycle
     linit = .false. ! for calling routine
     lstart = .false.
  end do ! while (stat%ncyc.lt.ctrl%maxcyc.and.iret.eq.0.and.(.not.ost%lconv))

!//////////////////////////////////////////////////////////////////////////////
! finalise the optimisation
!//////////////////////////////////////////////////////////////////////////////

! interpret the return code - jump point for errors
9999 continue
  select case (ierr)
  case (0)
     if (ost%lconv) then
        if (ost%lend) then
           iret = 1
        else
           iret = -2
        end if
     else if (stat%ncyc .ge. ctrl%maxcyc) then
        iret = -3
        if (ctrl%printl.ge.-1) write (stdout,'(/,A,I7,/)') &
             'Maximum number of cycles reached: ', ctrl%maxcyc
        call search_stat
     else
        iret = -1
     end if
  case (-1)
     iret = -4
     if (ctrl%printl.ge.-1) write (stdout,'(/,A,/)') &
          'SCF convergence error in MNDO'
     call search_stat
  case (-8)
     iret = -5
     if (ctrl%printl.ge.-1) write (stdout,'(/,A,/)') 'Time limit reached'
     call search_stat
  case (-9)
     iret = -9
  case default 
     if (ctrl%printl.ge.-1) write (stdout,'(/,A,I3,/)') 'Error code', ierr
     iret = -1
  end select

! dump memory if sensible
  if (ctrl%idump.gt.0 .and. iret.ne.-1 .and. iret.ne.-4) then
     call hdlc_wr_chk (ctrl%mdump, ctrl%nbig, ctrl%nvar, ctrl%nrem, ctrl, &
          stat, hdlc, ost, owr, coords, hess, resn)
  end if

! free memory, local and control structures
  call hdlc_destroy_all (.false., err_cnt)
  deallocate (coords); deallocate (grad); deallocate (hess)
  deallocate (coords_tmp); deallocate (grad_tmp)
  deallocate (cgrad)
  deallocate (resn)
  deallocate (err_cnt)
  deallocate (ctrl%attypes); deallocate (ctrl%precon)
  deallocate (ctrl%incon); deallocate (ctrl%iconstr)
  deallocate (ctrl%vconstr)
  deallocate (ctrl%irad); deallocate (ctrl%vrad)
  call deallocate_work (owr)

end subroutine hdlcopt
