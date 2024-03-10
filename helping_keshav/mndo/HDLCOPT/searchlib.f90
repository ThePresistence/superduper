module searchlib
  use global
  use hdlclib
  use matrixlib
  implicit none

!==============================================================================
! Important common variables of search:
! =====================================
! ost:           The 'stack', a collection of scalars important to be saved
! owr:           The work arrays for L-BFGS and P-RFO
!
! Top subroutines of search:
! ==========================
! hess_computer: Calculate a finite difference Hessian
! prfoss:        Implements L-BFGS, P-RFO and microiterative optimisation
!   hdlc_updhes: Update the Hessian
!   hdlc_prjhes: Project rot/trans degress of freedom of Cartesian Hessian
!
! I / O subroutines of search:
! ============================
! hdlc_rd_chk:   Reads all control data and the status from a checkpoint file
! hdlc_wr_chk:   Writes all control data and the status from a checkpoint file
! The routines to dump / undump the 'stack' and the work arrays are
!    hdlc_rd_stack, hdlc_wr_stack, hdlc_rd_work, hdlc_wr_work
! The routines to allocate / deallocate the work arrays are:
!    allocate_work / deallocate_work
!==============================================================================

!------------------------------------------------------------------------------
! Optimiser internal 'registers': saved variables which can be dumped to disk,
! formerly residing on a 'stack' (could be moved to searchlib)
! - checkpointing: hdlc_wr_stack / hdlc_rd_stack
!
! These fields are obsolete and removed:
! ======================================
! toler:   ctrl%toler
! funct:   (argument only)
! rmin:    ctrl%rmin
! rmax:    ctrl%rmax
! dmax:    ctrl%dmax - although it can be changed by the optimiser
! ddmax:   ctrl%ddmax
! osmin:   ctrl%osmin
! omin:    ctrl%omin
! tmnine:  parameter cuteig
!
! nbig:    ctrl%nbig
! nvar:    ctrl%nvar
! nrem:    ctrl%nrem
! mode:    ctrl%mode
! ihess:   stat%ihess
! iprint:  ctrl%printl
! irecalc: ctrl%irecalc
! maxstep: ctrl%maxcyc
! iupd:    ctrl%iupd
! idump:   ctrl%idump
! iprj:    ctrl%iprj
! nstep:   stat%nstep, stat%nsprfo
! nfg:     stat%ncyc
! i, j, k, l: (local)
! dlc:     hdlc%lhdlc=.true. <=> dlc=1
! jump:    ctrl%jump passed as argument
! ilock:   ctrl%ilock
! dlcNam:  hdlc%first -> residue%name -> ... (not used in prfoss())
! lconv:   ost%lconv passed as argument
! lend:    ost%lend passed as argument
!
! These fields have been added:
! =============================
! xlamd:   was there but not saved
! xlamd0:  do.
! skal:    do.
! lcreset: .true. between reset of L-BFGS curvature and first successful step
! bfgsbmb: absolutely unacceptably small L-BFGS step
! bfgscut: cutoff for step etc. considered to be zero
! dmax:    initial value of ctrl%dmax (used by L-BFGS, but changed by P-RFO)
!------------------------------------------------------------------------------

  type opt_stack
     real(kind=8) :: tole1, tole2, tole3, tole4, lbtole1, lbtole2, lbtole3, &
          lbtole4, scale, func, demin, bfgsmin, step, ss, dd, olde, &
          ys, yy, diag, sq, yr, beta, dmag, ddmag, odmax, odd, oolde, depre, &
          xtmp, frodo, grad1, func1, deact, ratio, grad2, tol2, &
          formdstep, oldener, xlamd, xlamd0, skal, bfgsbmb, bfgscut, dmax
     integer negreq, mxstep, ir, iw, ireclc, icalcn, gh, &
          ispt, iypt, ipoint, ibound, iter, ij, neg, npt, icp, iycn, inmc, &
          iscn, ipt, nopt, imode, overlpit
     logical lupd, lts, lrjk, lorjk, rrscal, donr, gnmin, ts, lrdmx, lprj, &
          lcore, lconv, lend, lcreset
  end type opt_stack

  type(opt_stack) :: ost

!------------------------------------------------------------------------------
! Optimiser work arrays (could be moved to searchlib)
! - memory management: allocate_work / deallocate_work
! - checkpointing: hdlc_wr_work / hdlc_rd_work
!
! for L-BFGS:
! ===========
! d     (nbig)
! wlbfg (nbig*(2*nrem+1)+2*nrem) - organisation: see prfoss() in searchlib
!
! for P-RFO:
! =========
! fx     (nvar)
! oldfx  (nvar)
! oldf   (nvar*nvar)
! ooldf  (nvar)
! eigval (nvar)
! oldeig (nvar)
! oldhss (nvar,nvar)
! u      (nvar,nvar)
! oldu   (nvar,nvar)
! svec   (nvar)
! tvec   (nvar)
! vmode  (nvar)
! hessc  (nvar*nvar)
! uc     (nvar,nvar) - points to u - is interpreted as (nvar*nvar)
!------------------------------------------------------------------------------

  type opt_work
     real(kind=8), dimension(:), pointer :: d, wlbfg, fx, oldfx, oldf, ooldf, &
          eigval, oldeig, svec, tvec, vmode, hessc
     real(kind=8), dimension(:,:), pointer :: oldhss, u, uc, oldu
  end type opt_work

  type(opt_work) :: owr

contains

!------------------------------------------------------------------------------
! subroutine hess_computer
!
! Calculate a finite difference Hessian
!
! Input:
! ======
! nvar               Number of degrees of freedom of the reaction core
! nstep              Number of f&g evaluations for Hessian so far (nshess)
! jump               'State machine', see global.f90
! coords (nvar)      The first nvar components of the coordinates
! funct              Energy
! grad   (nvar)      The first nvar components of the gradient
!
! Output:
! =======
! hess   (nvar,nvar) Finite difference Hessian of the core
! lval_hess          .true. if the Hessian is valid
! nstep              equal stat%nshess: (incremented)
! jump               'State machine', see global.f90
!------------------------------------------------------------------------------

  subroutine hess_computer (coords, funct, grad, hess, nvar, jump, lval_hess, &
       nstep)

! args
    logical lval_hess
    integer nvar, nstep, jump
    real(kind=8) funct
    real(kind=8), dimension(nvar) :: coords, grad
    real(kind=8), dimension(nvar,nvar) :: hess

! local params
    real(kind=8) epsilon
    parameter (epsilon = 0.0005_8)

! local vars
    integer i, j, place

! begin, init call
    if (nstep .eq. 0) then
       call hess_wr_ini (nvar)
       nstep = nstep + 1
       coords(1) = coords(1) - epsilon
       return

! final call - get differences
    else
       if (nstep .gt. nvar*2) then
          do i = 1,nvar
             do j = 1,i
                hess(i,j) = (hess(i,j)+hess(j,i)) * 0.5_8
                hess(j,i) = hess(i,j)
             end do
          end do
          nstep = 0
          jump = 0
          lval_hess = .true.
          return
       end if

! neither init nor final call
       place = nstep/2 + mod(nstep,2)
       if (ctrl%printl.ge.2) then
          write (stdout,'(3x,a,i5)') 'Computing step number ', nstep
          write (stdout,'(3x,a,i5)') 'Computing row  number ', place
       end if
    end if

! do finite difference step
    if (mod(nstep,2) .eq. 1) then
       do i = 1,nvar
          hess(place,i) = grad(i)
       end do
       coords(place) = coords(place) + 2*epsilon
    else
       do i=1,nvar
          hess(place,i) = (grad(i)-hess(place,i)) / (2*epsilon)
       enddo
       coords(place) = coords(place) - epsilon
       place = place+1
       if (place .le. nvar) then
          coords(place) = coords(place) - epsilon
       end if
    end if
    nstep = nstep+1

! contained routine: print out initialisation info
  contains
    subroutine hess_wr_ini (nvar)
      integer nvar
      if (ctrl%printl.ge.1) then
         write(stdout,'(/,a)') '##############################################'
         write(stdout,'(a)') '#    ENTERING "SEARCH" HESSIAN COMPUTER      #'
         write(stdout,'(a)') '#                                            #'
         write(stdout,'(a)') '#              Version 1.0                   #'
         write(stdout,'(a)') '#                                            #'
         write(stdout,'(a)') '#     AJT / SB        Dec 1997 / Aug 1999    #'
         write(stdout,'(a)') '##############################################'
         write(stdout,'(/,a,i5,/)') 'Total number of steps ', nvar*2
      end if
    end subroutine hess_wr_ini
  end subroutine hess_computer

!------------------------------------------------------------------------------
! subroutine prfoss
!
! Performs one L-BFGS or P-RFO step.
!
! P-RFO part: it is a quasi Newton-Raphson optimisation routine based on Jack
! Simons P-RFO algorithm as implemented by Jon Baker (J. Comp. Chem. 7, 385).
! The algorithm for step scaling to keep the length within trust radius is
! taken from Culot et al. (Theo. Chim. Acta 82, 189).
! The trust radius can be updated dynamically according to Fletcher safeguards
! on valid step for TS searches based on actual/predicted function change. This
! and change in TS mode are own modifications of AJT.
!
! This version treats the first nvar variables in the normal way and uses a low
! memory BFGS minimisation step on the rest of the variables (microiterative
! optimisation).
! The L-BFGS part minimises the bulk of the system to zero gradient for each
! step of the core part, where zero gradient is defined as 1/3 of the
! respective target values for the core of the RMS and maximum component of
! the gradient.
!
! The L-BFGS uses a trust radius approach, not a line search - to improve
! stability of the algorithm on non-exact surfaces.
!
! Input:
! ======
! nbig   number of degrees of freedom of the entire system (ctrl%nbig-nconstr)
! nvar   number of parameters to be optimised using P-RFO (ctrl%nvar)
! nrem   number of remembered L-BFGS steps (ctrl%nrem)
! coords values of parameters to be optimised
!        dimension: (nbig)
! funct  energy of the system
! grad   gradient of funct w/r to coords
!        dimension: (nbig)
! hess   Hessian of funct w/r to the first nvar coordinates of coords
!        dimension: (nvar,nvar)
! cgrad  Cartesian gradient
!        dimension: (nbig)
! jump   (stat%jump) state machine of the entire optimiser (see global.f90)
!
! Output:
! ======
! coords optimised parameters
! jump   (stat%jump) state machine of the entire optimiser (see global.f90)
! lconv  (stat%lconv)
! lend   (stat%lend)
!
! History:
! ========
! Implementation combined NR, P-RFO and QA algorithm together with trust radius
!   update and step rejection was made October 1992 by F. Jensen, Odense, DK
! Microiterative approach (P-RFO / L-BFGS): 1997, A.J. Turner, Zurich
! Rewritten: 1999, S.R. Billeter, Zurich
!------------------------------------------------------------------------------

  subroutine prfoss (coords, funct, grad, hess, cgrad, nbig, nvar, &
       nrem, jump, lconv, lend, lval_hess)

! args
    logical lconv, lend, lval_hess
    integer nbig, nvar, nrem, jump
    real(kind=8) funct
    real(kind=8), dimension(nbig) :: coords, grad, cgrad
    real(kind=8), dimension(nvar,nvar) :: hess

! externals
    real(kind=8) ddot
    external ddot

! local params
    real(kind=8) cuteig
    parameter (cuteig = 1.0E-9_8)

! local vars
    integer i, j

! begin
    lend = .false.

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump < 0: this section of code is performed on startup +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (jump.le.0) then

! get all initialisation data
       stat%ihess = 0
       stat%nstep = 1
       ost%formdstep = 5.0e-02_8
       ost%bfgsmin = 1.0e-7_8
       ost%bfgsbmb = 1.0e-10_8
       ost%lcreset = .false.
       ost%bfgscut = 1.0e-11_8
       ost%demin = 2.0e-5_8
       ost%icalcn = 0
       ost%donr = .false.

! pointers (new)
       ost%npt = 0

! set up exit tolerances (G94 standard)
       ost%tole3 = ctrl%toler
       ost%tole4 = ctrl%toler * 4.0_8/6.0_8
       ost%tole1 = ctrl%toler * 4.0_8
       ost%tole2 = ost%tole1 * 2.0_8/3.0_8

! set lbfgs to one third of the values above unless nvar=0
      ost%lbtole1 = ost%tole1 
      ost%lbtole2 = ost%tole2 
      if (nvar .ne. 0) then
         ost%scale = 1.0_8/3.0_8
         ost%lbtole1 = 1.0_8
         ost%lbtole2 = 1.0_8
      else
         ost%scale = 1.0_8
      end if
      ost%lbtole3 = ost%tole3 * ost%scale 
      ost%lbtole4 = ost%tole4 * ost%scale 

! init only - not restart, subcondition of: (jump.le.0)
      if (jump.eq.0) then
         ost%lts = .false.
         ost%gnmin = .false.

! PRFO used for TS search
         if (ctrl%mode .ne. 0) then
            ost%lts = .true.
            ost%negreq = 1

! PRFO used for minimisation
         else
            ost%lts = .false.
            ost%negreq = 0
         end if

!//////////////////////////////////////////////////////////////////////////////
! print out parameters
!//////////////////////////////////////////////////////////////////////////////
         
         call prfoss_wr_ini (nbig)

! flag a restart - subcondition of: (jump.le.0)
      else
         if (ctrl%printl.ge.1) then
            write (stdout,'(a)') '############################################'
            write (stdout,'(a)') '#     RESTARTING "SEARCH" OPTIMISER        #'
            write (stdout,'(a)') '############################################'
         else if (ctrl%printl.ge.-1) then
            write (stdout,'(a,/)') 'Restarting the optimisation'
         end if
      end if
      if (ctrl%printl.ge.1) write (stdout,'(/)')

      ost%lorjk = .false.
      ost%lrdmx = .false.
      ost%lprj = .false.
      if ((ctrl%iprj .eq. 1) .and. (.not. hdlc%lhdlc)) ost%lprj = .true.

! for now, allow dynamic trust radius all the time
      ost%lupd = .true.

!//////////////////////////////////////////////////////////////////////////////
! Set up L-BFGS stuff
!
! The work vector owr%wlbfg is organised as follows:
! ------------------------------------------
!
! nbig elements:      (1) ... (nbig)
!                   - store the gradient or some temporary information
! nrem elements:      (nbig+1) ... (nbig+nrem)
!                   - store the scalars rho
! nrem elements:      (nbig+nrem+1) ... (nbig+2*nrem) 
!                   - store the numbers alpha used for the computation of H*g
! nbig*nrem elements: (nbig+2*nrem+1) ... (nbig+2*nrem+nbig*nrem) 
!                   - store the last nrem search steps
! nbig*nrem elements: (nbig+2*nrem+nbig*nrem+1) ... (nbig+2*nrem+2*nbig*nrem) 
!                   - store the last nrem gradient differences
!
! The search steps and gradient differences are stored in a circular order
! controlled by the parameter ipoint
!//////////////////////////////////////////////////////////////////////////////

! set initial L-BFGS trust radius on start up (not on HDLC failure restart)
      if (jump .eq. 0) then
         ost%step = ctrl%dmax
         ost%dmax = ctrl%dmax
      end if

! next action: get initial gradient and do initial step
      jump = 1
   end if ! no need to return here

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 1: start optimisation, initial point tested +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (jump .eq. 1) then
      if (ctrl%irecalc .eq. 0) ctrl%irecalc = -1

! set oldener to greater than energy to avoid error condition (exit)
      ost%oldener = funct + 1.0_8
      do i = 1,nbig
         owr%d(i) = 0.0_8
      end do

! check if start geometry is already optimised
      if (ctrl%printl.ge.0) write (stdout,'(a)') 'Initial configuration'
      if (nvar .lt. nbig) then
         call progress (stat%nstep, stat%ihess, stat%ncyc, funct, owr%d, &
              grad(nvar+1), cgrad(nvar+1), nbig-nvar, nbig-nvar, &
              ost%lbtole1, ost%lbtole2, ost%lbtole3, ost%lbtole4, &
              .false., lconv, ost%oldener, hdlc%lhdlc, ctrl%icconv)
      else
         call progress (stat%nstep, stat%ihess, stat%ncyc, funct, owr%d, &
              grad, cgrad, nbig, nbig, &
              ost%tole1, ost%tole2, ost%tole3, ost%tole4, &
              .true., lconv, ost%oldener, hdlc%lhdlc, ctrl%icconv)
         if (lconv) goto 280
      end if
      ost%oldener = funct

!//////////////////////////////////////////////////////////////////////////////
! initialise L-BFGS with unit diagonal Hessian
!
! - set first record of step in lbfgs to -grad
! - set up two pointers into the array owr%wlbfg
!//////////////////////////////////////////////////////////////////////////////

      do i = 1, nbig*(2*nrem+1)+2*nrem
         owr%wlbfg(i) = 0.0_8
      end do

! set pointer to and reset step area
      ost%ispt = nbig + 2*nrem
      do i = nvar+1, nbig
         owr%wlbfg(ost%ispt+i) = -grad(i)
      end do

! pointer to gradient area
      ost%iypt = ost%ispt + nbig*nrem

! counter for circlar storage of area
      ost%ipoint = 0

! number of stored steps
      ost%ibound= 0

! iteration counter for lbfgs stuff
      ost%iter=0

!//////////////////////////////////////////////////////////////////////////////
! end jump = 1 section, determine next action
!//////////////////////////////////////////////////////////////////////////////

      if (nvar.eq.nbig) then
         jump = 5
      else if (nvar.gt.0 .and. lconv) then
         if (ctrl%printl.ge.1) then
            write (stdout,'(A,A)') 'Start geometry of environment relaxed, ', &
                 'starting with core'
         endif
         jump = 5
      else
         jump = 2
      end if
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 2: form a geometry step using L-BFGS +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

10 if (jump.eq.2) then
      ost%iter = ost%iter + 1
      if (ost%iter .ne. 1) then
         if (ost%iter .gt. nrem) then
            ost%ibound = nrem
         else
            ost%ibound = ost%iter - 1
         end if
         ost%ys = ddot (nbig, owr%wlbfg(ost%iypt+ost%npt+1), 1, &
              owr%wlbfg(ost%ispt+ost%npt+1), 1)
         ost%yy = ddot (nbig, owr%wlbfg(ost%iypt+ost%npt+1), 1, &
              owr%wlbfg(ost%iypt+ost%npt+1), 1)

! diag is a single number here because a unit Hessian is used
         ost%diag = ost%ys / ost%yy
         if (ost%ipoint .eq. 0) then
            owr%wlbfg(nbig+nrem) = 1.0_8 / ost%ys
         else
            owr%wlbfg(nbig+ost%ipoint) = 1.0_8 / ost%ys
         end if
         do i = 1,nvar
            owr%wlbfg(i) = 0.0_8
         end do
         do i = nvar+1,nbig
            owr%wlbfg(i) = -grad(i)
         end do
         ost%icp = ost%ipoint
         do i = 1,ost%ibound
            ost%icp = ost%icp - 1
            if (ost%icp .eq. -1) ost%icp = nrem - 1
            ost%sq = ddot (nbig, owr%wlbfg(ost%ispt+ost%icp*nbig+1), 1, &
                 owr%wlbfg, 1)
            ost%inmc = nbig + nrem + ost%icp + 1
            ost%iycn = ost%iypt + ost%icp*nbig
            owr%wlbfg(ost%inmc) = owr%wlbfg(nbig+ost%icp+1)*ost%sq
            call daxpy (nbig, -owr%wlbfg(ost%inmc), owr%wlbfg(ost%iycn+1), &
                 1, owr%wlbfg, 1)
         end do
         do i = 1,nbig
            owr%wlbfg(i) = owr%wlbfg(i) * ost%diag
         end do
         do i = 1,ost%ibound
            ost%yr = ddot (nbig, owr%wlbfg(ost%iypt+ost%icp*nbig+1), 1, &
                 owr%wlbfg, 1)
            ost%beta = owr%wlbfg(nbig+ost%icp+1) * ost%yr
            ost%inmc = nbig + nrem + ost%icp + 1
            ost%beta = owr%wlbfg(ost%inmc) - ost%beta
            ost%iscn = ost%ispt + ost%icp*nbig
            call daxpy (nbig, ost%beta, owr%wlbfg(ost%iscn+1), 1, owr%wlbfg, 1)
            ost%icp = ost%icp + 1
            if (ost%icp .eq. nrem) ost%icp = 0
         end do
         ost%ipt = ost%ispt + ost%ipoint*nbig
         do i = 1,nbig
            owr%wlbfg(ost%ipt+i) = owr%wlbfg(i)
         end do

! end if (ost%iter .ne. 1)
      end if

! form an L-BFGS step only in non-core coordinates
      do i = nvar+1,nbig
         owr%wlbfg(i) = grad(i)
      end do
      do i = 1,nvar
         owr%wlbfg(i) = 0.0_8
      end do

! update pointers
      ost%npt = ost%ipoint * nbig
      ost%ipoint = ost%ipoint + 1
      if (ost%ipoint.ge.nrem) ost%ipoint = 0

!//////////////////////////////////////////////////////////////////////////////
! end jump = 2 section (L-BFGS step formation), determine next action
!//////////////////////////////////////////////////////////////////////////////

      jump = 3
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 3: perform L-BFGS step +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

15 if (jump.eq.3) then 

! store old energy and one dimensional gradient and scale step to trust radius
      ost%ipt = ost%ispt + ost%npt
      ost%dmag = sqrt (ddot (nbig-nvar,owr%wlbfg(ost%ipt+1+nvar),1, &
           owr%wlbfg(ost%ipt+1+nvar),1))
      if (ost%dmag .gt. ost%step) then
         ost%ddmag = ost%step / ost%dmag
         do i = nvar+1,nbig
            owr%wlbfg(ost%ipt+i) = owr%wlbfg(ost%ipt+i) * ost%ddmag
         end do
      end if
      do i = nvar+1,nbig
         coords(i) = coords(i) + owr%wlbfg(i+ost%ipt)
      end do
      ost%ipt = ost%ispt + ost%npt
      ost%grad1 = ddot (nbig-nvar,grad(nvar+1),1,owr%wlbfg(ost%ipt+1+nvar),1)
      ost%func1 = funct
      if (ost%grad1 .gt. 0.0_8) then
         do i = nvar+1,nbig
            owr%wlbfg(ost%ipt+i) = -owr%wlbfg(ost%ipt+i)
            coords(i) = coords(i) + 2.0_8*owr%wlbfg(ost%ipt+i)
         end do
         if (ctrl%printl.ge.1) write (stdout,'(a,/)') &
              'Inverting L-BFGS step direction'
      end if

!//////////////////////////////////////////////////////////////////////////////
! end jump = 3 section, determine action after energy and gradient evaluation
!//////////////////////////////////////////////////////////////////////////////

      stat%nstried = stat%nstried + 1
      jump = 4
      return
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 4: test if L-BFGS step OK +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (jump.eq.4) then
      if (ctrl%printl.ge.2) &
           write (stdout,'(a,/)') 'Testing quality of L-BFGS step'

! if energy increases - reduce trust radius
      if (funct .gt. ost%func1) then
         if (ctrl%printl.ge.1) &
              write (stdout,'(5x,a)') 'Step failed to reduce function'
         ost%step = min (ost%dmag, ost%step)
         ost%step = ost%step * 0.5_8
         if (ctrl%printl.ge.1) write (stdout,'(5x,a,g13.5)') &
              'Current L-BFGS trust radius ', ost%step
         ost%ipt = ost%ispt + ost%npt
         do i = nvar+1,nbig
            coords(i) = coords(i) - owr%wlbfg(i+ost%ipt)
         end do

! restart L-BFGS if step is too small
         if (ost%step.lt.ost%bfgsmin .and. .not.ost%lcreset) then
            ost%lcreset = .true.
            ost%ipoint = 0
            ost%ibound = 0
            ost%iter = 0
            ost%step = ost%dmax
            if (ctrl%printl.ge.0) write (stdout,'(5x,a,i5,/)') &
                 'Resetting L-BFGS curvature information at step number ', &
                 stat%ncyc
            do i = nvar+1,nbig
               owr%wlbfg(ost%ispt+i) = -grad(i)
            end do
            do i=1,nvar
               owr%wlbfg(ost%ispt+i) = 0.0_8
            end do

!//////////////////////////////////////////////////////////////////////////////
! jump to reform L-BFGS step if step failed and was too small (jump = 2)
!//////////////////////////////////////////////////////////////////////////////

            jump = 2
            goto 10

! step failed but is not too small
         else
            if (ost%step.lt.ost%bfgsbmb) then
               write (stdout,'(A,E10.4)') &
                    'Minimum acceptable L-BFGS trust radius is ', ost%bfgsbmb
               call hdlc_errflag ('L-BFGS trust radius too small', 'stop')
            end if
            funct = ost%func1
            if (ctrl%printl.ge.0) write (stdout,'(5x,a,i5,/)') &
                 'Rejecting L-BFGS step at cycle number ', stat%ncyc

!//////////////////////////////////////////////////////////////////////////////
! get new energy and gradient after smaller step (jump = 3)
!//////////////////////////////////////////////////////////////////////////////

            jump = 3
            goto 15   
         end if

! step is OK, store it for L-BFGS
      else ! if (funct .gt. ost%func1) then ...
         ost%lcreset = .false.
         stat%nslbfg = stat%nslbfg + 1
         do i = nvar+1,nbig
            owr%wlbfg(ost%iypt+ost%npt+i) = grad(i) - owr%wlbfg(i)
         end do
         do i=1,nvar
            owr%wlbfg(ost%iypt+ost%npt+i) = 0.0_8
         end do

! if Wolfe conditions satisfied - increase trust radius
         ost%grad2 = ddot (nbig-nvar, grad(nvar+1), 1, &
              owr%wlbfg(ost%ipt+1+nvar), 1)
         if (funct .le. ost%func1+ost%grad1*1.0e-4_8 .and. &
              abs(ost%grad2) .le. abs(ost%grad1)*0.9_8) then
            ost%step = ost%step * 2.0_8
            if (ctrl%printl.ge.1) write (stdout,'(5x,a)') &
                 'Step has passed the Wolfe conditions'
         else
            if (ctrl%printl.ge.1) write (stdout,'(5x,a)') &
                 'Step has failed the Wolfe conditions'
            ost%step = min (ost%step, ost%dmag)
         end if

! increase trust radius if small but function decreasing
         if (funct .le. ost%func1+ost%grad1*1.0e-4_8 .and. &
              ost%dmag .gt. ost%step) then
            ost%step = ost%step * 1.25_8
         end if
         ost%step = min (ctrl%ddmax,ost%step)
         if (ctrl%printl.ge.1) then
            write (stdout,'(5x,a,g13.5)') &
                 'L-BFGS trust radius is   ', ost%step
            write (stdout,'(5x,a,g13.5,/)') &
                 'L-BFGS predicted step is ', ost%dmag
         end if

         call dummyarg_checkin (owr%wlbfg, ost%ipt+nvar+1, nbig-nvar)
         call progress (stat%nstep, stat%ihess, stat%ncyc, funct, &
              dummyarg, grad(nvar+1), cgrad(nvar+1), nbig-nvar, nbig-nvar, &
              ost%lbtole1, ost%lbtole2, ost%lbtole3, ost%lbtole4, .false., &
              lconv, ost%oldener, hdlc%lhdlc, ctrl%icconv)
         call dummyarg_clear

!//////////////////////////////////////////////////////////////////////////////
! form the next L-BFGS step if not converged yet (jump = 2)
!//////////////////////////////////////////////////////////////////////////////

         if (.not. lconv) then
            ost%oldener = funct
            jump = 2
            goto 10
         end if

!//////////////////////////////////////////////////////////////////////////////
! end optimisation if converged and L-BFGS only (no P-RFO part)
!//////////////////////////////////////////////////////////////////////////////

         if (nvar.eq.0) goto 280

!//////////////////////////////////////////////////////////////////////////////
! perform P-RFO step if converged and nvar>0 (jump = 5)
!//////////////////////////////////////////////////////////////////////////////

         lconv = .false.
         jump = 5
         if (ctrl%printl.ge.0) &
              write (stdout,'(/,a,/)') 'L-BFGS converged, going for P-RFO'
         goto 20
      end if

! end of jump = 4 section
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 5: form P-RFO step +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

20 if (jump.eq.5) then

!//////////////////////////////////////////////////////////////////////////////
! Start of P-RFO loop
!
! We now have gradients and a Hessian. If this is the first time through,
! don't update the Hessian. Hessian recalculation is controlled in search
! using ctrl%irecalc and is coupled to a complete restart of the optimiser
!//////////////////////////////////////////////////////////////////////////////

! calculate the Hessian first if not yet done
      if (.not. lval_hess) return

! store various things for possible omin rejection
      do i = 1,nvar
         owr%oldfx(i) = owr%fx(i)
         owr%ooldf(i) = owr%oldf(i)
         owr%oldeig(i) = owr%eigval(i)
         do j = 1,nvar
            owr%oldhss(i,j) = hess(i,j)
            owr%oldu(i,j) = owr%u(i,j)
         end do
      end do

! update Hessian if required
      if (stat%ihess.gt.0) then
         call hdlc_updhes (owr%svec, owr%tvec, nvar, ctrl%iupd, hess, grad, &
              owr%oldf, owr%d, owr%vmode, owr%u, ost%dd, ctrl%rmin, &
              ctrl%rmax, ctrl%omin, ost%xlamd, ost%xlamd0, ost%skal, &
              ctrl%mode, stat%nstep, ost%negreq, ctrl%printl, ctrl%ilock)
      end if

! project out the rot/trans deg. of freedom from the Hessian if cartesian only
      if (ost%lprj) call hdlc_prjhes (hess, coords, nvar, owr%u)

! test for convergence (remeber to include L-BFGS step)
      do i = nvar+1,nbig
         owr%d(i) = owr%wlbfg(ost%ipt+i)
      end do
      call progress (stat%nstep, stat%ihess, stat%ncyc, funct, owr%d, grad, &
           cgrad, nbig, nbig, ost%tole1, ost%tole2, ost%tole3, ost%tole4, &
           .true., lconv, ost%oldener, hdlc%lhdlc, ctrl%icconv)
      if(lconv) goto 280

      do i = nvar+1,nbig
         owr%d(i) = 0.0_8
      end do
      ost%oldener = funct

! back up stuff before diagonalisation
      ost%olde = funct
      do i = 1,nvar
         owr%oldf(i) = grad(i)
      end do

! pack hessian into lower triangle for compatibility of checkpoint files
      ost%ij = 0
      do i = 1,nvar
         do j = 1,i
            ost%ij = ost%ij+1
            owr%hessc(ost%ij) = hess(j,i)
         end do
      end do

! diagonalise Hessian and get off eigenvectors - use the non-object entry
      i = array_diagonalise (hess, owr%uc, owr%eigval, &
           nvar, nvar, nvar, .true.)
      if (i .ne. 0) then
         call hdlc_errflag ('Hessian could not be diagonalised', 'stop')
      end if

! set small eigenvalues to zero
      do i = 1,nvar
         if (abs(owr%eigval(i)) .lt. cuteig) owr%eigval(i) = 0.0_8
      end do

! this piece of code is obsolete since uc -> u; uncomment whenever changed
!     ost%ij = nvar**2
!     do i = nvar,1,-1
!        do j = nvar,1,-1
!           owr%u(j,i) = owr%uc(ost%ij)
!           ost%ij = ost%ij-1
!        end do
!     end do

! print out number of negative eigenvalues
      ost%neg = 0
      do i = 1,nvar
         if (owr%eigval(i) .lt. 0.0_8) ost%neg = ost%neg + 1
      end do
      if (ctrl%printl .ge. 1) write (stdout,'(3X,A,I3,A)') 'Hessian has ', &
           ost%neg, ' negative eigenvalue(s)'
      if (ctrl%printl .ge. 3) then
         do i = 1,nvar,6
            if (nvar-i .gt. 6) then
               write (stdout,'(5x,6(g13.5,1X))') (owr%eigval(j), j=i,i+5)
            else
               write (stdout,'(5x,6(g13.5,1X))') (owr%eigval(j), j=i,nvar)
            end if
         end do
      end if
      call flushout

! if projection of r/t, remove contribs from components with near zero eigenv.
      do i = 1,nvar
         owr%fx(i) = ddot (nvar, owr%u(1,i), 1, grad, 1)
         if (ost%lprj) then
            if (abs(owr%eigval(i)) .lt. cuteig) owr%fx(i) = 0.0_8
         end if
      end do

   end if ! (jump.eq.5)

!//////////////////////////////////////////////////////////////////////////////
! still jump = 5: in case of step rejection go here
!//////////////////////////////////////////////////////////////////////////////
 
130 if (jump.eq.5) then
      call hdlc_formd (owr%eigval, owr%fx, nvar, ctrl%dmax, ctrl%osmin, &
           ost%lts, ost%lrjk, ost%lorjk, ost%rrscal, ost%donr, &
           owr%oldf, owr%d, owr%vmode, owr%u, ost%dd, ctrl%rmin, ctrl%rmax, &
           ctrl%omin, ost%xlamd, ost%xlamd0, ost%skal, &
           ctrl%mode, stat%nstep, ost%negreq, ctrl%printl, ctrl%ilock, &
           ost%formdstep, ost%overlpit)
      call flushout

! reject previous step if TS mode overlap is less than omin (see hdlc_formd)
      if (ost%lorjk) then
         if (ctrl%printl.ge.1) write (stdout,'(3X,A)') &
              'Now undoing previous P-RFO step'
         ctrl%dmax = ost%odmax
         ost%dd = ost%odd
         ost%olde = ost%oolde

! Hessian eigenvalues / eigenvectors etc.
         do i = 1,nvar
            owr%fx(i) = owr%oldfx(i)
            owr%oldf(i) = owr%ooldf(i)
            owr%eigval(i) = owr%oldeig(i)
            do j = 1,nvar
               hess(i,j) = owr%oldhss(i,j)
               owr%u(i,j) = owr%oldu(i,j)
            end do
         end do

! position
         if (stat%nstep.ne.1) then
            do i = 1,nvar
               coords(i) = coords(i) - owr%d(i)
            end do
         end if

! reduce trust radius
         ctrl%dmax = min(ctrl%dmax,ost%dd)/2.0_8
         ost%odmax = ctrl%dmax
         ost%odd = ost%dd
!        stat%nstep = stat%nstep - 1
         if (ctrl%printl.ge.1) write (stdout,'(3x,a)') &
              'Finished undoing, now trying new step'

! new mode is now set to old mode, so overlap is 1 and a new step is formed
         goto 130
      end if ! (lorjk)

! store new trial geometry from owr%d for convergence test
      do i = nvar+1,nbig
         owr%d(i) = 0.0_8
      end do

! perform the P-RFO step to compare predicted energy change with actual
      do i = 1,nvar
         coords(i) = coords(i) + owr%d(i)
      end do
      ost%depre = 0.0_8
      ost%imode = 1
      if (ctrl%mode .ne. 0) ost%imode = ctrl%mode
      do i = 1,nvar
         ost%xtmp = ost%xlamd
         if (ost%lts .and. i.eq.ost%imode) ost%xtmp = ost%xlamd0
         if (abs(ost%xtmp-owr%eigval(i)) .lt. 1.0e-6_8) then
            ost%ss = 0.0_8
         else
            ost%ss = ost%skal*owr%fx(i)/(ost%xtmp-owr%eigval(i))
         end if
         ost%frodo = ost%ss*owr%fx(i) + 0.5_8*ost%ss*ost%ss*owr%eigval(i)
         ost%depre = ost%depre + ost%frodo
      end do

! counters are incremented here
      stat%ihess = stat%ihess + 1
      stat%nstep = stat%nstep + 1
      stat%nsprfo = stat%nsprfo +1

!//////////////////////////////////////////////////////////////////////////////
! end jump = 5 section, determine action after energy and gradient evaluation
!//////////////////////////////////////////////////////////////////////////////

      jump = 6
      return
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! + jump = 6: continue with core optimisation +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (jump.eq.6) then
      if (ctrl%printl.ge.1) write (stdout,'(a,/)') 'P-RFO trust radius report:'
      ost%deact = funct - ost%olde
      ost%ratio = ost%deact/ost%depre
      if (ctrl%printl.ge.1) &
           write (stdout,'(3X,A,/,5X,A,F14.7,/,5X,A,F14.7,/,5X,A,F14.7)') &
           'Energy change:', 'actual:    ', ost%deact, 'predicted: ', &
           ost%depre, 'ratio:     ', ost%ratio

      ost%lrjk = .false.

! default rmin=0 does not allow the energy to raise; the user can set rmin<>0
      if (.not.ost%lts .and. ost%ratio.lt.ctrl%rmin) then
         if (ctrl%dmax.ge.ctrl%osmin) then
            if (ctrl%printl.ge.1) write (stdout,'(3X,A,F7.4)') &
                 'Ratio below rmin, rejecting step, reducing dmax to ', &
                 min(ctrl%dmax,ost%dd)/2.0_8
            ost%lrjk = .true.

! allow min step for relaxed environment stuff
         else
            if (ctrl%printl.ge.1) write (stdout,'(3X,A,F7.4)') &
                 'Ratio below rmin, but dmax < osmin ', ctrl%osmin
         end if
      end if

! reject step on energy ratio grounds
      if (ost%lts .and. &
           (ost%ratio.lt.ctrl%rmin .or. ost%ratio.gt.ctrl%rmax) .and. &
           (abs(ost%depre).gt.ost%demin .or. abs(ost%deact).gt.ost%demin)) then
         if (ctrl%dmax .ge. ctrl%osmin) then
            if (ctrl%printl.ge.1) write (stdout,'(3X,A,F7.4)') &
                 'Unacceptable ratio, rejecting step, reducing dmax to ', &
                 min(ctrl%dmax,ost%dd)/2.0_8
            ost%lrjk = .true.

! allow min step for relaxed environment stuff
         else
            if (ctrl%printl.ge.1) write (stdout,'(3X,A,F7.4)') &
                 'Unacceptable ratio, but dmax < osmin ', ctrl%osmin
         end if
      end if

! step rejection - with update
      if (ost%lrjk) then
         stat%nsprfo = stat%nsprfo - 1
         ost%ipt = ost%ispt + ost%npt
         do i = nvar+1,nbig
            coords(i) = coords(i) - owr%wlbfg(i+ost%ipt)
         end do
         do i = 1,nvar
            coords(i) = coords(i) - owr%d(i)
         end do
         ctrl%dmax = min(ctrl%dmax,ost%dd)/2.0_8

!//////////////////////////////////////////////////////////////////////////////
! P-RFO step rejected - jump to reform P-RFO step (jump = 5)
!//////////////////////////////////////////////////////////////////////////////

         jump = 5
         goto 130
      end if

!//////////////////////////////////////////////////////////////////////////////
! step successful
!//////////////////////////////////////////////////////////////////////////////

      ost%odmax = ctrl%dmax
      ost%odd = ost%dd
      ost%oolde = ost%olde

!//////////////////////////////////////////////////////////////////////////////
! Dynamical trust radius of core based on energy prediction
! =========================================================
!
! Fletcher recommends dmax=dmax/4 and dmax=dmax*2
! the factors below are a little more conservative since Hessian is being
! updated don't reduce trust radius due to ratio for min searches
!
! Special case (lrdmax=.true. and dmax>osmin):
! scale dmax down to avoid rejected steps
! trying to improve the integration of the overlap/steplength relationship 
!//////////////////////////////////////////////////////////////////////////////

! scale down
      if (ost%lrdmx) then
         if (ctrl%dmax.ge.ctrl%osmin) then
            ctrl%dmax = ctrl%dmax/2.0_8
         else
            if (ctrl%printl.ge.1) write (stdout,'(3X,A)') &
                 'Not reducing dmax: dmax < minstep'
         end if

! dynamic trust radius
      else
         if (ost%lts) then
            if (ost%ratio.ge.0.75_8 .and. ost%ratio.le.(4.0_8/3.0_8) .and. &
                 ost%dd.gt.(ctrl%dmax-1.0E-6_8)) then
               ctrl%dmax = ctrl%dmax * sqrt(2.0_8)
            end if

! allow wider limits for increasing trust radius for minimum searches
         else
            if (ost%ratio.ge.0.5_8 .and. ost%dd.gt.(ctrl%dmax-1.0E-6_8)) then
               ctrl%dmax = ctrl%dmax * sqrt(2.0_8)
            end if
         end if

! be brave if 0.90 < ratio < 1.10 ...
         if (abs(ost%ratio-1.0_8).lt.0.1_8) ctrl%dmax = ctrl%dmax * sqrt(2.0_8)
      end if

! step minimum step length - precaution, this should be ok already
      ctrl%dmax = max (ctrl%dmax, ctrl%osmin)
      ctrl%dmax = min (ctrl%dmax, ctrl%ddmax)
      if (ctrl%printl.ge.1) then
         write (stdout,'(3X,A,F7.5)') 'Current trust radius = ', ctrl%dmax
         write (stdout,'(3X,A,F9.5)') 'Step size used is ', ost%dd
      end if

!//////////////////////////////////////////////////////////////////////////////
! end jump = 6 section, determine the next action
!//////////////////////////////////////////////////////////////////////////////

      if (nvar.eq.nbig) then
         jump = 5
      else
         jump = 2
      end if
      goto 10
   end if

!//////////////////////////////////////////////////////////////////////////////
! optimisation termination
!//////////////////////////////////////////////////////////////////////////////

280 continue
   if (ctrl%printl.ge.0) &
        write (stdout,'(/,A,/)') ' **************** OPTIMISED ****************'
   call search_stat
   lconv = .true.
   lend = .true.
   return

!//////////////////////////////////////////////////////////////////////////////
! end of code of subroutine prfoss
!//////////////////////////////////////////////////////////////////////////////

  contains

!//////////////////////////////////////////////////////////////////////////////
! contained in prfoss: subroutine prfoss_wr_ini
!
! flag parameters on start up
!//////////////////////////////////////////////////////////////////////////////

    subroutine prfoss_wr_ini (nbig)

! args
      integer nbig

! write a banner
      if (ctrl%printl.ge.1) then
         write (stdout,'(/,a)')'##############################################'
         write (stdout,'(a)') '#       ENTERING "SEARCH" OPTIMISER          #'
         write (stdout,'(a)') '#                                            #'
         write (stdout,'(a)') '#              Version 1.0                   #'
         write (stdout,'(a)') '#                                            #'
         write (stdout,'(a)') '#    AJT / SB      -   Dec 1997 / Aug 1999   #'
         write (stdout,'(a)') '##############################################'
      end if
      if (ctrl%printl.ge.0) then
         write (stdout,'(/,a,/)') 'Optimisation setup report:'

! generic
         write (stdout,'(3X,A,I5,A)') &
              'First ', ctrl%nvar, ' variables by (P)RFO '
         write (stdout,'(3X,A,I5,A)') &
              'Rest  ', nbig-ctrl%nvar, ' variables by L-BFGS'

! no core
         if (ctrl%nvar .eq. 0) then
            write (stdout,'(3x,a)') 'Minimising using L-BFGS only'

! there is a core - P-RFO params
         else
            if (ctrl%mode .ne. 0) then
               write (stdout,'(3X,A,I4)') 'Searching mode ', ctrl%mode
            else
               write (stdout,'(3X,A)') 'Minimising'
            end if
            if (ctrl%irecalc .ne. -1) then
               write (stdout,'(3X,A,I5,A)') 'Recompute Hessian every ', &
                    ctrl%irecalc, ' cycles'
            else
               write (stdout,'(3x,a)')    'Use Hessian update only' 
            end if
            if (ctrl%iupd .eq. 2) then
               write (stdout,'(3x,a)') &
                    'Hessian update scheme              BFGS'
            else
               write (stdout,'(3x,a)') &
                    'Hessian update scheme              Powell'
            end if
            if (ctrl%ilock .eq. 0) then
               write (stdout,'(3x,a)') 'Mode switching allowed'
            else
               write (stdout,'(3x,a)') 'Mode locked'
            end if
            if ((ctrl%iprj .eq. 0) .and. (.not. hdlc%lhdlc)) then
               write (stdout,'(3x,a)') 'R/T modes projected out of Hessian'
            else
               write (stdout,'(3x,a)') 'R/T modes left in Hessian'
            end if
         end if

! params
         write (stdout,'(3x,a,g13.5)') &
              'Minimum ratio of act and pred DE ', ctrl%rmin
         write (stdout,'(3x,a,g13.5)') &
              'Maximum ratio of act and pred DE ', ctrl%rmax
         write (stdout,'(3x,a,g13.5)') &
              'Minimum overlap criterion        ', ctrl%omin
         write (stdout,'(3x,a,g13.5)') &
              'Reduce T/R on overlap <          ', sqrt (ctrl%omin)
         write (stdout,'(3x,a,g13.5)') &
              'Minimum DE for dynamic t/r       ', ost%demin
         write (stdout,'(3x,a,g13.5)') &
              'Maximum step size                ', ctrl%ddmax
         write (stdout,'(3x,a,g13.5)') &
              'Initial step size                ', ctrl%dmax
         write (stdout,'(3x,a,g13.5)') &
              'Minimum LBFGS step size          ', ost%bfgsmin
         write (stdout,'(3x,a,g13.5)') &
              'Minimum P-RFO step size          ', ctrl%osmin

! P-RFO targets
         if (ctrl%nvar .gt. 0) then
            write (stdout,'(/,a)') &
                 'P-RFO targets'
            write (stdout,'(3x,a,g13.5)') &
                 'Target maximum step component    ', ost%tole1
            write (stdout,'(3x,a,g13.5)') &
                 'Target maximum step R.M.S.       ', ost%tole2
            write (stdout,'(3x,a,g13.5)') &
                 'Target maximum gradient component', ost%tole3
            write (stdout,'(3x,a,g13.5)') &
                 'Target maximum gradient R.M.S.   ', ost%tole4
         end if

! L-BFGS targets
         write (stdout,'(/,a)') &
              'L-BFGS targets'
         write (stdout,'(3x,a,g13.5)') &
              'Target maximum step component    ', ost%lbtole1
         write (stdout,'(3x,a,g13.5)') &
              'Target maximum step R.M.S.       ', ost%lbtole2
         write (stdout,'(3x,a,g13.5)') &
              'Target maximum gradient component', ost%lbtole3
         write (stdout,'(3x,a,g13.5)') &
              'Target maximum gradient R.M.S.   ', ost%lbtole4
         write (stdout,'(3x,a,i5)') &
              'Number remebered steps in L-BFGS ', ctrl%nrem
      end if
    end subroutine prfoss_wr_ini

!//////////////////////////////////////////////////////////////////////////////
! contained in prfoss: subroutine progress
!
! report optimisation progress to stdout
!
! Parameters:
! ===========
! cutstep:   if step is less than cutstep, don't apply step criterion (accept)
! fgradonly: if gradient is less than tolerance*fgradonly, accept step anyway
!
! Arguments:
! ==========
! grad = gradient
!    d = step
! cgrd = Cartesian gradient
!
! size = size to get RMS
! nbig = actual size
!
! tole1 = tolerance of max step
! tole2 = tolerance of rms step
! tole3 = tolerance of max gradient
! tole4 = tolerance of rms gradient
!
! lcore = .true. for (P)RFO step, .false. for L-BFGS
! lconv = .true. if converged
! lrcrt = .true. if Cartesian gradient is printed out too
!
! iccnv = convergence criterion (0: Cart/HDLC step/grad, 1: Cart grad)
!//////////////////////////////////////////////////////////////////////////////

    subroutine progress (nstep, ihess, nfg, funct, d, grad, cgrd, size, &
         nbig, tole1, tole2, tole3, tole4, lcore, lconv, oldener, &
         lrcrt, iccnv)

! args
      logical lcore, lconv, lrcrt
      integer nstep, ihess, nfg, nbig, size, iccnv
      real(kind=8) tole1, tole2, tole3, tole4, funct, oldener
      real(kind=8), dimension(nbig) :: d, grad, cgrd

! externals
      real(kind=8) ddot
      external ddot

! local params
      real(kind=8) cutstep, fgradonly
      parameter (cutstep = 1.0E-6_8)
      parameter (fgradonly = 0.01_8)

! local vars
      logical lrmsna, lmaxna
      integer i, imaxgr, imaxd
      real(kind=8) rmxgrad, rmxsgrad, rmxd, rmsd, diff, rmsgrad

! begin, get step in energy terms
      diff = funct - oldener

! get RMS gradient for whole system
      rmsgrad = sqrt (ddot(nbig,grad,1,grad,1) / size)

! get maximum component of gradient
      rmxgrad = 0.0_8
      imaxgr = 0
      do i = 1,nbig
         if (abs(grad(i)) .gt. rmxgrad) then
            rmxgrad = abs(grad(i))
            imaxgr = i
         end if
      end do

! get rms step size
      rmsd = sqrt (ddot(nbig,d,1,d,1) / size)

! get maximum step size
      rmxd = 0.0_8
      imaxd = 0
      do i = 1,nbig
         if (abs(d(i)) .gt. rmxd) then
            rmxd = abs(d(i))
            imaxd = i
         end if
      end do

! test for applicability of RMS step criterion
      lrmsna = .false.
      if (abs(diff).le.cutstep .and. diff.le.0.0_8) lrmsna=.true.
      if (rmsgrad .lt. tole2*fgradonly)             lrmsna=.true.

! test for applicability of max step criterion
      lmaxna=.false.
      if (abs(diff).le.cutstep .and. diff.le.0.0_8) lmaxna=.true.
      if (rmxgrad .lt. tole1*fgradonly)             lmaxna=.true.

! report and test
      if (ctrl%printl.ge.1) &
           write (stdout,'(a)') 'Optimisation progress report:'
      if (ctrl%printl.ge.0) &
           write (stdout,500) nstep, ihess, nfg, funct, diff
      if (ctrl%printl.ge.1 .and. lrcrt) write (stdout,1600)

      lconv=.true.

! detect = step 0 condition and flag as unconverged
! this is caused by an HDLC restart
      if ( (abs(rmsd) .lt. ost%bfgscut .and. abs(rmxd) .lt. ost%bfgscut) .and. &
           (abs(rmxgrad).gt.ost%bfgscut .or. abs(rmsgrad).gt.ost%bfgscut) ) then
         if( (lmaxna .and. lrmsna) .or. iccnv.eq.1) then
            if (ctrl%printl.ge.0) write (stdout,1001) 
         else
            if (ctrl%printl.ge.0) write (stdout,1002) 
            lconv = .false.
         endif

! regular case, subsequent steps
      else

         if (rmxd.le.tole1 .and. iccnv.eq.0) then
            if (ctrl%printl.ge.0) write (stdout,1000) rmxd, tole1, 'yes', imaxd
         else if (lmaxna .or. iccnv.eq.1) then
            if (ctrl%printl.ge.0) write (stdout,1000) rmxd, tole1, 'n/a', imaxd
         else
            if (ctrl%printl.ge.0) write (stdout,1000) rmxd, tole1, 'no ', imaxd
            lconv = .false.
         endif

         if (rmsd.le.tole2 .and. iccnv.eq.0) then
            if (ctrl%printl.ge.0) write (stdout,1100) rmsd, tole2, 'yes'
         else if (lrmsna .or. iccnv.eq.1) then
            if (ctrl%printl.ge.0) write (stdout,1100) rmsd, tole2, 'n/a'
         else
            if (ctrl%printl.ge.0) write (stdout,1100) rmsd, tole2, 'no '
            lconv = .false.
         endif

      endif

      if (rmxgrad.le.tole3 .and. iccnv.eq.0) then
         if (ctrl%printl.ge.0) &
              write (stdout,1300) rmxgrad, tole3, 'yes', imaxgr
      else if (iccnv.eq.1) then
         if (ctrl%printl.ge.0) &
              write (stdout,1300) rmxgrad, tole3, 'n/a', imaxgr
      else
         if (ctrl%printl.ge.0) &
              write (stdout,1300) rmxgrad, tole3, 'no ', imaxgr
         lconv = .false.
      endif

      if (rmsgrad.le.tole4 .and. iccnv.eq.0) then
         if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'yes'
      else if (iccnv.eq.1) then
         if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'n/a'
      else
         if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'no '
         lconv = .false.
      end if

! get RMS Cartesian gradient for whole system
      if (lrcrt .or. iccnv.eq.1) then
         rmsgrad = sqrt (ddot(nbig,cgrd,1,cgrd,1) / size)

! get RMS maximum component of gradient
         rmxgrad = 0.0_8
         imaxgr = 0
         do i = 1,nbig
            if (abs(cgrd(i)) .gt. rmxgrad) then
               rmxgrad = abs(cgrd(i))
               imaxgr = i
            end if
         end do
         if (ctrl%printl.ge.0) write (stdout,1500)

         if (rmxgrad.le.tole3 .and. iccnv.eq.1) then
            if (ctrl%printl.ge.0) write (stdout,1300) rmxgrad, tole3, 'yes', &
                 imaxgr
         else if (iccnv.eq.0) then
            if (ctrl%printl.ge.0) write (stdout,1300) rmxgrad, tole3, 'n/a', &
                 imaxgr
         else
            if (ctrl%printl.ge.0) write (stdout,1300) rmxgrad, tole3, 'no ', &
                 imaxgr
            lconv = .false.
         end if

         if (rmsgrad.le.tole4 .and. iccnv.eq.1) then
            if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'yes'
         else if (iccnv.eq.0) then
            if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'n/a'
         else
            if (ctrl%printl.ge.0) write (stdout,1400) rmsgrad, tole4, 'no '
            lconv = .false.
         end if
      end if ! (lcrcrt)

! return
      call flushout
      return

! format statements
 500  format (/,5x,'Steps: ',i5,' Hessian updates: ',i5, &
           ' Function/gradient evals: ',i5, &
           /,5x,'Function: ',g20.12,' Difference: ',g13.5)
1000  format(5x,'MaxStep: ',g12.5,' Target: ',g12.5,' converged? ',a3, &
           ' Component: ',i5)
1001  format(5x,'Step tests not available',19x,'converged? n/a ')
1002  format(5x,'Step tests not available',19x,'converged? no ')
1100  format(5x,'RMSStep: ',g12.5,' Target: ',g12.5,' converged? ',a3)
1300  format(5x,'MaxGrad: ',g12.5,' Target: ',g12.5,' converged? ',a3, &
           ' Component: ',i5)
1400  format(5x,'RMSGrad: ',g12.5,' Target: ',g12.5,' converged? ',a3)
1500  format(/,3x,'Cartesian gradient:')
1600  format(/,3x,'HDLC step and gradient:')
    end subroutine progress

!//////////////////////////////////////////////////////////////////////////////
! contained in prfoss: subroutine hdlc_formd
!
! Form a core optimisation step.
! This version forms geometry step by either pure NR, P-RFO or QA algorithm,
! under the condition that the steplength is less than dm.
! Although working, pure NR is currently disabled.
!
! Note about adaptation (SB): did only remove part of goto / label pairs and
! inserted very few comments so far
!//////////////////////////////////////////////////////////////////////////////

    subroutine hdlc_formd (eigval, fx, nvar, dmax, osmin, &
         lts, lrjk, lorjk, rrscal, donr, &
         oldf, d, vmode, u, dd, rmin, rmax, omin, xlamd, xlamd0, skal, &
         mode, nstep, negreq, printl, ilock, step, overlpit)

! args
      logical lts, lrjk, lorjk, rrscal, donr
      integer nvar, mode, nstep, negreq, printl, ilock, overlpit
      real(kind=8) dmax, osmin, dd, rmin, rmax, omin, xlamd, xlamd0, skal, step
      real(kind=8), dimension(nvar) :: eigval, fx, d, vmode
      real(kind=8), dimension(nvar*nvar) :: oldf
      real(kind=8), dimension(nvar,nvar) :: u

! externals
      real(kind=8) ddot
      external ddot

! local params
      real(kind=8) big, eps, eps2, eps6, sfix, toll, zero
      parameter (big = 1000.0_8, &
           eps = 1.0E-12_8, &
           eps2 = 1.0E-2_8, &
           eps6 = 1.0E-6_8, &
           sfix = 10.0_8, &
           toll = 1.0E-8_8, &
           zero = 0.0_8)

! local vars
      logical rscal, frodo1, frodo2, lrdmx
      integer it, jt, newmod, i, j, k, ncnt
      real(kind=8) lambda, lambda0, eigit, eone, ssmin, ssmax, sstoll,&
           bl, bu, fl, fu, xlambda, fm, temp, d2max,sstep

! begin
      skal = 1.0_8
      rscal = rrscal
      it = 0
      jt = 1

! calculate overlap between Hessian eigenmodes if TS search
      if (lts) then
         if (mode.ne.0) then
            call overlp (dmax, osmin, newmod, nvar, lorjk, lrdmx, &
                oldf, d, vmode, u, dd, rmin, rmax, omin, xlamd, xlamd0, skal, &
                mode, nstep, negreq, printl, ilock, overlpit)
            if (lorjk) return

! on return from overlp, newmod is the TS mode
            if (ilock.eq.0) then
               if (newmod.ne.mode .and. printl.ge.1) then
                  write (stdout,'(A,I3,A,I3)') &
                       'Warning! Switching mode to be followed from ', &
                       mode, ' to ', newmod
               end if
               mode = newmod
            end if
            it = mode
         else
            it = 1
         end if
         eigit = eigval(it)
         if (printl.ge.1) then
            write (stdout,'(/,3X,A,I3,A,g12.5)') 'TS mode is number ', it, &
                 ' with eigenvalue ', eigit
         end if
      end if

! we now have mode it and its eigenvalue eigit if TS search is requested
      if (it.eq.1) jt = 2
      eone = eigval(jt)
      ssmin = max (abs(eone)*eps, (10.0_8*eps))
      ssmax = max (big, abs(eone))
      ssmax = ssmax*big
      sstoll = toll
      d2max = dmax*dmax

! start by bracketing root, then hunt it down with brute force bisection
      frodo1 = .false.
      frodo2 = .false.
      lambda = zero
      lambda0 = zero

! test if pure NR step
      if (lts .and. eigit.lt.zero .and. eone.ge.zero .and. donr) then
         if (printl.ge.1) then
            write (stdout,'(a)') &
                 'TS search, correct Hessian, trying pure NR step'
         end if
         goto 100
      end if
      if (.not.lts .and. eone.ge.zero .and. donr) then
         if (printl.ge.1) then
            write (stdout,'(a)') &
                 'Minimum search, correct Hessian, trying pure NR step'
         end if
         goto 100
      end if

! find lambda / lambda0 if not pure NR step
40    continue
      if (lts) then
         lambda0 = eigval(it) + sqrt(eigval(it)**2 + 4.0_8*fx(it)**2)
         lambda0 = 0.5_8*lambda0
         if (printl.ge.1) write (stdout,'(3X,A,F15.5)') &
              'Lambda that maximises along TS modes  = ', lambda0
      end if
      sstep = step
      if (eone.le.zero) lambda = eone - sstep
      if (eone.gt.zero) sstep = eone
      bl = lambda - sstep
      bu = lambda + sstep*0.5_8
50    fl = zero
      fu = zero
      do i = 1,nvar
         if (i.ne.it) then
            fl = fl + (fx(i)*fx(i))/(bl-eigval(i))
            fu = fu + (fx(i)*fx(i))/(bu-eigval(i))
         end if
      end do
      fl = fl - bl
      fu = fu - bu
      if (fl*fu .lt. zero) goto 70
      bl = bl - (eone-bl)
      bu = bu + 0.5_8*(eone-bu)
      if (bl.le.-ssmax) then
         bl = -ssmax
         frodo1 = .true.
      end if
      if (abs(eone-bu).le.ssmin) then
         bu = eone-ssmin
         frodo2=.true.
      end if
      if (frodo1 .and. frodo2) then
         if (printl.ge.0) then
            write (stdout,'(a)') 'Numerical problems bracketing lambda'
            write (stdout,'(a)')' Going for fixed step size....'
         end if
         goto 180
      end if
      goto 50
70    continue
      ncnt = 0
      xlambda = zero
80    continue
      fl = zero
      fu = zero
      fm = zero
      lambda = 0.5_8*(bl+bu)
      do i = 1,nvar
         if (i.ne.it) then
            fl = fl + (fx(i)*fx(i))/(bl-eigval(i))
            fu = fu + (fx(i)*fx(i))/(bu-eigval(i))
            fm = fm + (fx(i)*fx(i))/(lambda-eigval(i))
         end if
      end do
      fl = fl - bl
      fu = fu - bu
      fm = fm - lambda
      if (abs(xlambda-lambda) .lt. sstoll) goto 100
      ncnt = ncnt + 1
      if (ncnt.gt.1000) then
         call hdlc_errflag &
              ('Too many iterations in lambda bisect', 'stop')
      end if
      xlambda = lambda
      if (fm*fu.lt.zero) bl = lambda
      if (fm*fl.lt.zero) bu = lambda
      goto 80

! jumped here if pure NR or xlambda near lambda or problems in lambda bisect
100   continue
      if (printl.ge.1) write (stdout,'(3X,A,F15.5)') &
           'Lambda that minimises along all modes = ', lambda

! calculate the step
      do i = 1,nvar
         d(i) = zero
      end do
      do i = 1,nvar
         if (lambda.eq.zero .and. abs(eigval(i)).lt.eps2) then
            temp = zero
         else
            temp = fx(i)/(lambda-eigval(i))
         end if
         if (i.eq.it) then
            temp = fx(it)/(lambda0-eigval(it))
         end if
         if (printl.ge.5) write (stdout,'(5x,a)') 'Formd, delta step ', i ,temp
         do j = 1,nvar
            d(j) = d(j) + temp*u(j,i)
         end do
      end do
      dd = sqrt (ddot (nvar, d, 1, d, 1))
      if (printl.ge.1) then
         if (lambda.eq.zero .and. lambda0.eq.zero) &
              write (stdout,'(3X,A,F10.5)') 'Pure NR step has length ', dd
         if (lambda.ne.zero .and. lambda0.ne.-lambda) &
              write (stdout,'(3X,A,F10.5)') 'P-RFO step has length ', dd
      end if
      if (dd .lt. (dmax+eps6)) then
         xlamd = lambda
         xlamd0 = lambda0
         return
      end if
      if (lambda.eq.zero .and. lambda0.eq.zero) goto 40

! scale the calculated step
      if (rscal) then
         skal = dmax/dd
         do i = 1,nvar
            d(i) = d(i)*skal
         end do
         dd = sqrt (ddot (nvar, d, 1, d, 1))
         if (printl.ge.1) write (stdout,'(3X,A,F9.5)') &
              'Calculated step size too large, scaled with ', skal
         xlamd = lambda
         xlamd0 = lambda0
         return
      end if

! fixed step size if bracketing failed
  180 lambda = zero
      frodo1 = .false.
      frodo2 = .false.
      sstep = step
      if (eone.le.zero) lambda = eone - sstep
      if (lts .and. -eigit.lt.eone) lambda = -eigit - sstep
      if (eone.gt.zero) sstep = eone
      bl = lambda - sstep
      bu = lambda + sstep*0.5_8
  190 fl = zero
      fu = zero
      do i = 1,nvar
         if (i.ne.it) then
            fl = fl + (fx(i)/(bl-eigval(i)))**2
            fu = fu + (fx(i)/(bu-eigval(i)))**2
         end if
      end do
      if (lts) then
         fl = fl + (fx(it)/(bl+eigval(it)))**2
         fu = fu + (fx(it)/(bu+eigval(it)))**2
      end if
      fl = fl - d2max
      fu = fu - d2max
      if (fl*fu .ge. zero) then
         bl = bl - (eone-bl)
         bu = bu + 0.5_8*(eone-bu)
         if (bl.le.-ssmax) then
            bl = -ssmax
            frodo1 = .true.
         end if
         if (abs(eone-bu) .le. ssmin) then
            bu = eone - ssmin
            frodo2 = .true.
         end if
         if (frodo1.and.frodo2) then
            if (printl.ge.0) then
               write (stdout,'(a)') 'Numerical problems bracketing lambda'
               write (stdout,'(a)') 'Going for fixed level shifted NR step...'
            end if

! both lambda searches failed, go for fixed level shifted NR
! this is unlikely to produce anything useful, but maybe
            lambda = eone - sfix
            lambda0 = eigit + sfix
            rscal = .true.
            goto 100
         end if
         goto 190
      end if

! go on...
      ncnt = 0
      xlambda = zero
  220 continue
      fl = zero
      fu = zero
      fm = zero
      lambda = 0.5_8*(bl+bu)
      do i = 1,nvar
         if (i.ne.it) then
            fl = fl + (fx(i)/(bl-eigval(i)))**2
            fu = fu + (fx(i)/(bu-eigval(i)))**2
            fm = fm + (fx(i)/(lambda-eigval(i)))**2
         end if
      end do
      if (lts) then
         fl = fl + (fx(it)/(bl+eigval(it)))**2
         fu = fu + (fx(it)/(bu+eigval(it)))**2
         fm = fm + (fx(it)/(lambda+eigval(it)))**2
      end if
      fl = fl - d2max
      fu = fu - d2max
      fm = fm - d2max
      if (abs(xlambda-lambda) .ge. sstoll) then
         ncnt = ncnt + 1
         if (ncnt.gt.1000) then
            call hdlc_errflag ('Too many iterations in lambda bisect', 'stop')
         end if
         xlambda = lambda
         if (fm*fu.lt.zero) bl = lambda
         if (fm*fl.lt.zero) bu = lambda
         goto 220
      end if

! try scaling
      lambda0 = -lambda
      rscal = .true.
      goto 100
    end subroutine hdlc_formd

!//////////////////////////////////////////////////////////////////////////////
! contained in prfoss: subroutine overlp
!
! Calculates the overlap between two Hessian modes
!//////////////////////////////////////////////////////////////////////////////

    subroutine overlp (dmax, osmin, newmod, nvar, lorjk, lrdmx, &
         oldf, d, vmode, u, dd, rmin, rmax, omin, xlamd, xlamd0, skal, &
         mode, nstep, negreq, printl, ilock, it)

! args
      logical lorjk, lrdmx
      integer newmod, nvar, mode, nstep, negreq, printl, ilock, it
      real(kind=8) dmax, osmin, dd, rmin, rmax, omin, xlamd, xlamd0, skal
      real(kind=8), dimension(nvar) :: d, vmode
      real(kind=8), dimension(nvar*nvar) :: oldf
      real(kind=8), dimension(nvar,nvar) :: u

! externals
      real(kind=8) ddot
      external ddot

! local vars
      integer i, j, k, l, icalcn
      real(kind=8) tovlp, ovlp

! init local vars
      save icalcn
      data icalcn/0/

! on the first step simply determine which mode to follow
      if (nstep.eq.1) then
         if (printl.ge.1) write (stdout,'(/,3X,A,I3,/)')  &
              'Hessian mode following switched on, following mode ', mode
         if (mode.gt.nvar) then
            write (stdout,'(A,I3,A,I3)') 'Mode to be followed: ', mode, &
                 ' , number of degrees of freedom (nvar): ', nvar
            call hdlc_errflag ('Error: mode is larger than nvar','stop')
         end if
         it = mode

! on subsequent steps determine which Hessian eigenvector has the greatest
! overlap with the mode we are following
      else
         it = mode
         lorjk = .false.
         tovlp = ddot (nvar, u(1,1), 1, vmode, 1)
         tovlp = abs (tovlp)
         i = 1
         if (printl.ge.4) write (stdout,'(A,I3,A,G12.5)') &
              'Overlap of mode ', i, ' = ', tovlp
         if (ilock .ne. 1) then
            do i = 2,nvar
               ovlp = ddot (nvar, u(1,i), 1, vmode, 1)
               ovlp = abs (ovlp)
               if (printl.ge.4) write (stdout,'(A,I3,A,G12.5)') &
                    'Overlap of mode ', i, ' = ', tovlp
               if (ovlp.gt.tovlp) then
                  tovlp = ovlp
                  it = i
               end if
            end do
         end if

! lrdmx controls if dmax is fixed (dmax is not increased) due to poor overlap
! (overlap < sqrt(omin)), even though the step is not rejected 
         lrdmx = .false.
         if (printl.ge.1) write (stdout,'(3X,A,I3,A,F6.3)') &
              'Overlap of current mode ', it, ' with previous mode is ', tovlp
         if (tovlp .lt. omin) then
            if (dmax .gt. osmin) then
               lorjk=.true.
               if (printl.ge.1) write (stdout,'(3X,A,F6.3,A)') &
                    'Overlap is less than omin = ', omin, &
                    ', rejecting previous step'
               return
            else
               if (printl.ge.1) write (stdout,'(3X,A,F6.3,A,F6.3,A,F6.3,A)') &
                    'Overlap is less than omin = ', omin, &
                    ', but trust radius ', dmax, ' is less than ', osmin, &
                    ', accepting step'
            end if
         else if (tovlp .lt. sqrt(omin)) then
            lrdmx = .true.
            if (printl.ge.1) write (stdout,'(3X,A,A)') &
                 'Overlap is less than sqrt(omin), accepting step, ', &
                 'reducing trust radius'
         end if
      end if

! save the eigenvector in vmode
      do i = 1,nvar
         vmode(i) = u(i,it)
      end do
      newmod = it
    end subroutine overlp

  end subroutine prfoss

!------------------------------------------------------------------------------
! subroutine hdlc_updhes
!
! Update the Hessian of the first nvar degrees of freedom using the BFGS or
! using the Powell formula.
! The new Hessian depends on the old Hessian, the current gradient, the old
! gradient and the correction vector used on the last cycle.
! svec & tvec are used as temporary storage
!
! (i)   The Powell update
!       The Powell formula preserves the symmetric character of the Hessian
!       while allowing its eigenvalue structure to change.
!       It is the default update for a transition state search
! (ii)  The BFGS update
!       The BFGS formula has the important characteristic of retaining
!       positive definiteness (note: this is not rigorously guaranteed, but 
!       can be checked for by the program).
!       It is the default update for a minimum search
!
! switch : iupd
!       iupd = 0  :  skip update
!       iupd = 1  :  Powell
!       iupd = 2  :  BFGS
!
! Note on positive definiteness: see code commented out below
!------------------------------------------------------------------------------

  subroutine hdlc_updhes (svec, tvec, nvar, iupd, hess, grad, &
       oldf, d, vmode, u, dd, rmin, rmax, omin, xlamd, xlamd0, skal, &
       mode, nstep, negreq, printl, ilock)

! args
    integer mode, nstep, negreq, printl, ilock, iupd, nvar
    real(kind=8) dd, rmin, rmax, omin, xlamd, xlamd0, skal
    real(kind=8), dimension(nvar) :: d, vmode, svec, tvec, grad
    real(kind=8), dimension(nvar*nvar) :: oldf
    real(kind=8), dimension(nvar,nvar) :: u, hess

! externals
    real(kind=8) ddot
    external ddot

! local vars
    integer i, icalcn, j
    real(kind=8) zero, dds, ddtd, temp

! init local vars
    save zero, icalcn 
    data icalcn /0/
    data zero /0.0_8/

! begin
    if (iupd.eq.0 .and. printl.ge.1) write (stdout,'(5X,A)') &
         'Hessian is not being updated'
    if (iupd.eq.1 .and. printl.ge.1) write (stdout,'(5X,A,/)') &
         'Hessian is being updated using the Powell update'
    if (iupd.eq.2 .and. printl.ge.1) write (stdout,'(5X,A,/)') &
         'Hessian is being updated using the BFGS update'
    if (iupd.gt.2 .or. iupd.lt.0) then
       write (stdout,'(A,I1)') 'Error ==> iupd out of bounds: ', iupd
       call flushout
       stop
    end if
    call flushout
    if (iupd.eq.0) return

! init temporary storage
    do i = 1,nvar
       tvec(i) = zero
    end do
    do j = 1,nvar
       do i = 1,nvar
          tvec(i) = tvec(i) + hess(i,j)*d(j)
       end do
    end do

! (i) Powell update
    if (iupd.eq.1) then
       do i = 1,nvar
          tvec(i) = grad(i)-oldf(i)-tvec(i)
          svec(i) = grad(i)-oldf(i)
       end do
       dds = dd*dd
       ddtd = ddot (nvar, tvec, 1, d, 1)
       ddtd = ddtd/dds
       do i = 2,nvar
          do j = 1,i-1
             temp = tvec(i)*d(j) + d(i)*tvec(j) - d(i)*ddtd*d(j)
             hess(i,j) = hess(i,j) + temp/dds
             hess(j,i) = hess(i,j)
          end do
       end do
       do i = 1,nvar
          temp = d(i) * (2.0_8*tvec(i)-d(i)*ddtd)
          hess(i,i) = hess(i,i) + temp/dds
       end do
    end if

! (ii) BFGS update
    if (iupd.eq.2) then
       do i = 1,nvar
          svec(i) = grad(i)-oldf(i)
       end do
       dds = ddot (nvar, svec, 1, d, 1)

! If dds is negative, retention of positive definiteness is not
! guaranteed. print a warning and skip update this cycle.
!
!frj With the current level shift technique I think the Hessian should
!frj be allowed to acquire negative eigenvalues. Without updating the
!frj optimization has the potential of stalling
!frj   if (dds.lt.zero) then
!frj      write(iw,100)
!frj      write(iw,110)
!frj      return
!frj   endif
       ddtd = ddot (nvar, d, 1, tvec, 1)
       do i = 2,nvar
          do j = 1,i-1
             temp = (svec(i)*svec(j))/dds - (tvec(i)*tvec(j)) / ddtd
             hess(i,j) = hess(i,j) + temp
             hess(j,i) = hess(i,j)
          end do
       end do
       do i = 1,nvar
          temp = (svec(i)*svec(i))/dds - (tvec(i)*tvec(i)) / ddtd
          hess(i,i) = hess(i,i) + temp
       end do
    end if

  end subroutine hdlc_updhes

!------------------------------------------------------------------------------
! subroutine hdlc_prjhes
!
! Project out roto/translational components from a cartesian Hessian.
!
! Most of this routine is the work of others.
! Overhauled ajt 1996 / sb 1999
!
! This was originally to work with mass weighted coordinates - but that is now
! changed and everything is done non-weighted.
! The projection of the gradient has also been removed.
!
! Arguments:
!   x      : cartesian coordinates
!   f      : cartesian force constants matrix, these are in
!             nvar(x,y,z) format, not in
!             nvar(x),nvar(y),nvar(z) format
!   dx     : this would be the normalised gradient vector but it is set to zero
!   p, cof : buffers
!------------------------------------------------------------------------------

  subroutine hdlc_prjhes (f, x, nvar, cof)

! args
    integer nvar
    real(kind=8), dimension(nvar) :: x
    real(kind=8), dimension(nvar,nvar) :: f, cof

! local params
    real(kind=8) zero, one, eps, cut5, cut8
    parameter (zero = 0.0_8)
    parameter (one = 1.0_8)
    parameter (eps=1.0e-14_8)
    parameter (cut5=1.0e-5_8)
    parameter (cut8=1.0e-8_8)

! local vars
    integer natm, nc1, i, j, k, l, m, n, info, ip, indx, jndx, jend, &
         ia, ib, ic, jp, ja, jb, jc, ii, jj
    real(kind=8) totm, chk, det, trp, sum
    real(kind=8) iscr(6), rot(3,3), scr(3,3), cmass(3)
    real(kind=8), dimension(:), allocatable :: dx
    real(kind=8), dimension(:,:), allocatable :: p

! init local vars
    real(kind=8), dimension(3,3,3) :: tens
    save tens

! totally asymmetric Cartesian tensor
    data tens / &
         0.0_8,  0.0_8,  0.0_8, &
         0.0_8,  0.0_8, -1.0_8, &
         0.0_8,  1.0_8,  0.0_8, &
         0.0_8,  0.0_8,  1.0_8, &
         0.0_8,  0.0_8,  0.0_8, &
        -1.0_8,  0.0_8,  0.0_8, &
         0.0_8, -1.0_8,  0.0_8, &
         1.0_8,  0.0_8,  0.0_8, &
         0.0_8,  0.0_8,  0.0_8  /

! begin
    allocate (dx(nvar))
    allocate (p(nvar,nvar))
    natm = nvar/3
    nc1 = nvar

! set the gradient to zero; a more sophisticated implementation would allow a
! gradient - but for optimiser work we want the gradient in there
! for stationary point g=0 anyway
    do i = 1,nc1
       dx(i) = zero
    end do

! cmass = centre of geometry rather than centre of mass
    cmass(1) = zero
    cmass(2) = zero
    cmass(3) = zero
    l = 0
    do i = 1,natm
       do j = 1,3
          l = l+1
          cmass(j) = cmass(j) + x(l)
       end do
    end do

! total 'mass' equals number of atoms
    totm = natm
    if (ctrl%printl.ge.1) write (stdout,'(a)') &
         'Projecting out rot/trans components from Hessian'

! compute inertia tensor
    do i = 1,3
       do j = 1,3
          rot(i,j) = zero
       end do
    end do
    do i = 1,natm
       l = 3*(i-1) + 1
       rot(1,1) = rot(1,1) + x(l+1)**2 + x(l+2)**2
       rot(1,2) = rot(1,2) - x(l)*x(l+1)
       rot(1,3) = rot(1,3) - x(l)*x(l+2)
       rot(2,2) = rot(2,2) + x(l)**2 + x(l+2)**2
       rot(2,3) = rot(2,3) - x(l+1)*x(l+2)
       rot(3,3) = rot(3,3) + x(l)**2 + x(l+1)**2
    end do
    rot(2,1) = rot(1,2)
    rot(3,1) = rot(1,3)
    rot(3,2) = rot(2,3)

! check the inertia tensor
    chk = rot(1,1)*rot(2,2)*rot(3,3)

! catch special cases of inertia tensor
    if (abs(chk) .le. cut8) then

! x.ne.0
       if (abs(rot(1,1)) .gt. cut8) then

! x,y.ne.0 but z.eq.0
          if (abs(rot(2,2)) .gt. cut8) then
             det = rot(1,1)*rot(2,2) - rot(1,2)*rot(2,1)
             trp = rot(1,1)
             rot(1,1) = rot(2,2)/det
             rot(2,2) = trp/det
             rot(1,2) = -rot(1,2)/det
             rot(2,1) = -rot(2,1)/det
             goto 100

! x,z.ne.0 but y.eq.0
          else if (abs(rot(3,3)) .gt. cut8) then
             det = rot(1,1)*rot(3,3) - rot(1,3)*rot(3,1)
             trp = rot(1,1)
             rot(1,1) = rot(3,3)/det
             rot(3,3) = trp/det
             rot(1,3) = -rot(1,3)/det
             rot(3,1) = -rot(3,1)/det
             goto 100

! x.ne.0 but y,z.eq.0
          else
             rot(1,1) = one/rot(1,1)
             goto 100
          end if

! y.ne.0 but x.eq.0
       else if (abs(rot(2,2)) .gt. cut8) then

! y,z.ne.0 but x.eq.0
          if (abs(rot(3,3)) .gt. cut8) then
             det = rot(3,3)*rot(2,2) - rot(3,2)*rot(2,3)
             trp = rot(3,3)
             rot(3,3) = rot(2,2)/det
             rot(2,2) = trp/det
             rot(3,2) = -rot(3,2)/det
             rot(2,3) = -rot(2,3)/det
             goto 100

! y.ne.0 but x,z.eq.0
          else
             rot(2,2) = one/rot(2,2)
             goto 100
          end if

! z.ne.0 but x,y.eq.0
       else if (abs(rot(3,3)) .gt. cut8) then
          rot(3,3) = one/rot(3,3)
          goto 100

! x,y,z near zero: cannot invert
       else
          if (ctrl%printl.ge.0) then
             write (stdout,'(A)') 'Warnig from hdlc_prjhes'
             write (stdout,'(3x,3g8.3)') &
                  'Every diagonal element of Hessian near zero: ', &
                  rot(1,1), rot(2,2), rot(3,3)
             call flushout
          end if
          return
       end if

! generic case of the inertia tensor
    else

! compute inverse matrix of rot
       info = array_invert (rot, det, .false., 3)
       if (ctrl%printl.ge.0 .and. info.ne.0) then
          write (stdout,'(/,a)') &
               'Warning: inversion of rotation matrix in hdlc_prjhes failed'
          write(stdout,'(a,/)') &
               '          --- no projection performed ----'
          call flushout
          return
       end if
    end if

! the inverse of the inertia tensor is available at this point
100 continue

! compute p matrix
    do ip = 1,natm
       indx = 3*(ip-1)
       do jp = 1,ip
          jndx = 3*(jp-1)
          do ic = 1,3
             jend = 3
             if (jp.eq.ip) jend = ic
             do jc = 1,jend
                sum = zero

! accumulate terms for elements of p 
               do ia = 1,3
                   do ib = 1,3
                      if (tens(ia,ib,ic) .ne. zero) then
                         do ja = 1,3
                            do jb = 1,3
                               if (tens(ja,jb,jc) .ne. zero) then
                                  sum = sum + tens(ia,ib,ic)*tens(ja,jb,jc)* &
                                       rot(ia,ja)*x(indx+ib)*x(jndx+jb)
                               end if ! (tens(ja,jb,jc) .ne. zero)
                            end do ! (ja = 1,3)
                         end do ! (jb = 1,3)
                      end if ! (tens(ia,ib,ic) .ne. zero)
                   end do ! (ib = 1,3)
                end do ! (ia = 1,3)

! end accumulate
                ii = indx + ic
                jj = jndx + jc
                p(ii,jj) = sum + dx(ii)*dx(jj)
                if (ic.eq.jc) p(ii,jj) = p(ii,jj) + one/totm
             end do ! (jc = 1,jend)
          end do ! (ic = 1,3)
       end do ! (jp = 1,ip)
    end do ! (ip = 1,natm)

! compute delta(i,j) - p(i,j)
    do i = 1,nc1
       do j = 1,i
          p(i,j) = -p(i,j)
          if (i.eq.j) p(i,j) = one + p(i,j)
       end do
    end do

! neglect smaller values than 10**-8
    do i = 1,nc1
       do j = 1,i
          if (abs(p(i,j)).lt.cut8) p(i,j) = zero
          p(j,i) = p(i,j)
       end do
    end do

! post and premultiply f by p, use cof for scratch
    do i = 1,nc1
       do j = 1,nc1
          sum = zero
          do k = 1,nc1
             sum = sum + f(i,k)*p(k,j)
          end do
          cof(i,j)=sum
       end do
    end do

! compute p * f * p
    do i = 1,nc1
       do j = 1,nc1
          sum = zero
          do k = 1,nc1
             sum = sum + p(i,k)*cof(k,j)
          end do
          f(i,j)=sum
       end do
    end do

! clean up
    deallocate (p)
    deallocate (dx)

  end subroutine hdlc_prjhes

!==============================================================================
! end of optimisation subroutines - R / T removal routines follow
!==============================================================================

!------------------------------------------------------------------------------
! subroutine rtfit: remove overall rotational / translational motions
!
! L-BFGS minimization with FD gradients to a least squares fit of the result of
! a step to the initial coordinates via rotation and translation
!------------------------------------------------------------------------------

  subroutine rtfit (xinit, xpresent, nvar)

! args
    integer nvar
    real(kind=8), dimension(nvar) :: xinit, xpresent

! local vars

! begin

  end subroutine rtfit

!==============================================================================
! end of R / T removal subroutines - I / O routines follow
!==============================================================================

!------------------------------------------------------------------------------
! subroutine hdlc_rd_chk: dump everything to a checkpoint file
!------------------------------------------------------------------------------

  subroutine hdlc_rd_chk (mstart, nbig, nvar, nrem, s_ctrl, s_stat, s_hdlc, &
       s_ost, s_owr, s_coords, s_hess, s_resn)

! args
    integer mstart, nbig, nvar, nrem
    integer, dimension(nbig/3) :: s_resn
    real(kind=8), dimension(nbig) :: s_coords
    real(kind=8), dimension(nvar,nvar) :: s_hess
    type(opt_ctrl) :: s_ctrl
    type(hdlc_ctrl) :: s_hdlc
    type(opt_stat) :: s_stat
    type(opt_stack) :: s_ost
    type(opt_work) :: s_owr

! local vars
    logical lerr, lval_hess
    character*25 secthead
    integer, save :: w_nbig, w_nvar, w_nrem
    integer l, m, merr, mbak, natom

! begin
    if (mstart .le. 0) return
    if (ctrl%printl.ge.0) write (stdout,'(/,A,/)') 'Reading checkpoint file'

! open checkpoint file and read data for consistency check
    merr = 0
    if (lform) then
       open (unit=chkpnt, file=chknam, form='FORMATTED', status='OLD', err=99)
       read (chkpnt,'(A)',err=98) secthead
       read (chkpnt,'(3I8)',err=98) w_nbig, w_nvar, w_nrem
    else
       open (unit=chkpnt, file=chknam, form='UNFORMATTED', status='OLD', &
            err=99)
       write (chkpnt,err=98) secthead
       write (chkpnt,err=98) w_nbig, w_nvar, w_nrem
    end if

! consistency check
    if (nbig .ne. w_nbig) then
       if (ctrl%printl.ge.-1) then
          write (stdout,'(/,A,I1,A)') &
              '*** Error reading the checkpoint file, mstart: ', mstart, ' ***'
          write (stdout,'(A34,I6)') 'Degrees of freedom in input: ', nbig
          write (stdout,'(A34,I6)') 'Degrees of freedom in checkpoint: ', &
               w_nbig
          write (stdout,'(A,/)') 'Trying the default startup procedure'
       end if
       mstart = 0
       return
    end if
    if (nvar .ne. w_nvar) then
       if (ctrl%printl.ge.-1) then
          write (stdout,'(/,A,I1,A)') &
            '*** Warning reading the checkpoint file, mstart: ', mstart, ' ***'
          write (stdout,'(A46,I3)') &
               'Degrees of freedom of the core in input', nvar
          write (stdout,'(A46,I3)') &
               'Degrees of freedom of the core in checkpoint: ', w_nvar
          write (stdout,'(A,/)') 'Ignoring the control data in checkpoint'
       end if
       if (mstart .ge. 2) then
          mstart = 1
          if (ctrl%printl.ge.0) write (stdout,'(A)') 'Trying mstart = 1'
       end if
    end if
    if (nrem .ne. w_nrem) then
       if (ctrl%printl.ge.-1) then
          write (stdout,'(/,A,I1,A)') &
            '*** Warning reading the checkpoint file, mstart: ', mstart, ' ***'
          write (stdout,'(A41,I5)') 'Remembered steps (memory) in input: ', &
               nrem
          write (stdout,'(A41,I5)') &
               'Remembered steps (memory) in checkpoint: ', w_nrem
          write (stdout,'(A,/)') 'Trying the default startup procedure'
          write (stdout,'(A,/)') 'Ignoring the control data in checkpoint'
       end if
       if (mstart .ge. 2) then
          mstart = 1
          if (ctrl%printl.ge.0) write (stdout,'(A)') 'Trying mstart = 1'
       end if
    end if
    if (nbig.eq.w_nbig .and. nvar.eq.w_nvar .and. nrem.eq.w_nrem) then
       if (ctrl%printl.ge.1) &
            write (stdout,'(A)') 'Consistency check successful'
    end if

! read coordinates section (mstart .ge. 1)
    merr = 0
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
       read (chkpnt,'(3F15.10)',err=98) (s_coords(m), m=1,nbig)
    else
       read (chkpnt,err=98) secthead
       read (chkpnt,err=98) (s_coords(m), m=1,nbig)
    end if
    if (ctrl%printl.ge.1) write (stdout,100) secthead
    call hdlc_put_coords (ctrl%iccode, s_coords, nbig) ! ensure consistency
    if (mstart.le.1) goto 96

! read residue names section (mstart .ge. 2)
    merr = 1
    natom = nbig/3
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
       read (chkpnt,'(16I5)',err=98) (s_resn(m), m=1,natom)
    else
       read (chkpnt,err=98) secthead
       read (chkpnt,err=98) (s_resn(m), m=1,natom)
    end if
    if (ctrl%printl.ge.1) write (stdout,100) secthead

! read control data section (mstart .ge. 2) - see global
    mbak = mstart
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_ctrl (chkpnt, s_ctrl, lform, lerr)
    mstart = mbak
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead
    if (mstart.le.2) return    

! read HDLC section (mstart .ge. 3) - see hdlclib
    merr = 2
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_hdlc (chkpnt, s_hdlc, lform, lerr)
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead

! read optimiser status register section (mstart .ge. 3) - see global
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_stat (chkpnt, s_stat, lform, lerr)
    stat%lval_hess = (s_ctrl%nvar .eq. 0) ! true only if successfully read
    lval_hess = stat%lval_hess
    stat%ncyc = 0
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead

! read optimiser 'stack' section (mstart .ge. 3) - see global
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_stack (chkpnt, s_ost, lform, lerr)
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead

! read Hessian for P-RFO (mstart .ge. 3) - only if a Hessian has been stored
    merr = 2
    if (lval_hess .and. nvar.gt.0) then
       if (lform) then
          read (chkpnt,'(A)',err=98) secthead
          read (chkpnt,'(4E20.12)',err=98) ((s_hess(m,l), m=1,nvar), l=1,nvar)
       else
          read (chkpnt,err=98) secthead
          read (chkpnt,err=98) ((s_hess(m,l), m=1,nvar), l=1,nvar)
       end if
       if (ctrl%printl.ge.1) write (stdout,100) secthead
    end if
    stat%lval_hess = .true.

! read and allocate work arrays and L-BFGS history (mstart .ge. 3) - see global
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_work (chkpnt, s_owr, lform, lerr)
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead
    if (mstart .le. 3) return

! undump function / gradient code memory (mstart .ge. 4) - see interface lib.
    merr = 3
    if (lform) then
       read (chkpnt,'(A)',err=98) secthead
    else
       read (chkpnt,err=98) secthead
    end if
    call hdlc_rd_fandg (ctrl%iccode, chkpnt, lform, lerr)
    if (lerr) goto 98
    if (ctrl%printl.ge.1) write (stdout,100) secthead

! close checkpoint file
96  close (chkpnt,err=97)

! I/O failure jump points
97  return
98  continue
    mstart = merr
    if (ctrl%printl.ge.-1) then
       write (stdout,'(/,A)') 'Warning: problem reading from checkpoint file'
       write (stdout,'(A)') secthead
       write (stdout,'(A,I1,/)') 'Trying fallback procedure mstart = ', mstart
    end if
    close (chkpnt,err=97)
    return
99  continue
    if (ctrl%printl.ge.-1) then
       write (stdout,'(/,A)') &
            '*** Error: could not open checkpoint file for reading ***'
       write (stdout,'(A,/)') 'Trying the default startup procedure mstart = 0'
    end if
    mstart = 0
    return

! format statements
100 format ('Successfully read ',A)

  end subroutine hdlc_rd_chk

!------------------------------------------------------------------------------
! subroutine hdlc_wr_chk: dump everything to a checkpoint file
!------------------------------------------------------------------------------

  subroutine hdlc_wr_chk (mdump, nbig, nvar, nrem, s_ctrl, s_stat, s_hdlc, &
       s_ost, s_owr, s_coords, s_hess, s_resn)

! args
    integer mdump, nbig, nvar, nrem
    integer, dimension(nbig/3) :: s_resn
    real(kind=8), dimension(nbig) :: s_coords
    real(kind=8), dimension(nvar,nvar) :: s_hess
    type(opt_ctrl) :: s_ctrl
    type(hdlc_ctrl) :: s_hdlc
    type(opt_stat) :: s_stat
    type(opt_stack) :: s_ost
    type(opt_work) :: s_owr

! local vars
    logical lerr
    character*25 secthead
    integer i, l, m, natom

! begin, open and write data for consistency check
    natom = nbig/3
    secthead = 'section: nbig, nvar, nrem'
    if (lform) then
       open (unit=chkpnt, file=chknam, form='FORMATTED', status='REPLACE', &
            err=99)
       write (chkpnt,'(A)',err=98) secthead
       write (chkpnt,'(3I8)',err=98) nbig, nvar, nrem
    else
       open (unit=chkpnt, file=chknam, form='UNFORMATTED', status='REPLACE', &
            err=99)
       write (chkpnt,err=98) secthead
       write (chkpnt,err=98) nbig, nvar, nrem
    end if

! write coordinates
    secthead = 'section: coordinates'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
       write (chkpnt,'(3F15.10)',err=98) (s_coords(m), m=1,nbig)
    else
       write (chkpnt,err=98) secthead
       write (chkpnt,err=98) (s_coords(m), m=1,nbig)
    end if
       
! write residue names section
    secthead = 'section: residue numbers'
    natom = nbig/3
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
       write (chkpnt,'(16I5)',err=98) (s_resn(m), m=1,natom)
    else
       write (chkpnt,err=98) secthead
       write (chkpnt,err=98) (s_resn(m), m=1,natom)
    end if

! write control data section - see global
    secthead = 'section: control data'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_ctrl (chkpnt, s_ctrl, lform, lerr)
    if (lerr) goto 98

! write HDLC section - see hdlclib
    secthead = 'section: hdlc'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_hdlc (chkpnt, s_hdlc, lform, lerr)
    if (lerr) goto 98

! write optimiser status registers section - see global
    secthead = 'section: status registers'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_stat (chkpnt, s_stat, lform, lerr)
    if (lerr) goto 98

! write optimiser 'stack' section - see global
    secthead = 'section: optimiser stack'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_stack (chkpnt, s_ost, lform, lerr)
    if (lerr) goto 98

! write Hessian for P-RFO - only if a valid Hessian is there
    if (stat%lval_hess .and. nvar.gt.0) then
       secthead = 'section: Hessian'
       if (lform) then
          write (chkpnt,'(A)',err=98) secthead
          write (chkpnt,'(4E20.12)',err=98) ((s_hess(m,l), m=1,nvar), l=1,nvar)
       else
          write (chkpnt,err=98) secthead
          write (chkpnt,err=98) ((s_hess(m,l), m=1,nvar), l=1,nvar)
       end if
    end if

! write optimiser work arrays and L-BFGS history - see global
    secthead = 'section: work arrays'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_work (chkpnt, s_owr, lform, lerr)
    if (lerr) goto 98

! dump function / gradient code memory (mdump .eq. 1) - see interface library
    if (mdump .ne. 1) goto 96
    secthead = 'section: function / grad'
    if (lform) then
       write (chkpnt,'(A)',err=98) secthead
    else
       write (chkpnt,err=98) secthead
    end if
    call hdlc_wr_fandg (ctrl%iccode, chkpnt, lform, lerr)
    if (lerr) goto 98

! close checkpoint file
96  close (chkpnt,err=97)

! I/O failure jump points
97  return
98  continue
    if (ctrl%printl.ge.-1) then
       write (stdout,'(/,A)') 'Warning: problem writing to checkpoint file'
       write (stdout,'(A,/)') secthead
    end if
    close (chkpnt,err=97)
    return
99  continue
    if (ctrl%printl.ge.-1) then
       write (stdout,'(/,A,/)') &
            'Warning: could not open checkpoint file for writing'
    end if
    return
  end subroutine hdlc_wr_chk

!------------------------------------------------------------------------------
! subroutine hdlc_wr_stack
!
! Dump the optimisers 'stack' structure s_ost to iunit
!------------------------------------------------------------------------------

  subroutine hdlc_wr_stack (iunit, s_ost, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    type(opt_stack) :: s_ost

! begin
    lerr = .false.
    if (lform) then
       write (iunit,*,err=98) s_ost
    else
       write (iunit,err=98) s_ost
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_wr_stack

!------------------------------------------------------------------------------
! subroutine hdlc_rd_stack
!
! Read the optimisers 'stack' structure s_ost from iunit
!------------------------------------------------------------------------------

  subroutine hdlc_rd_stack (iunit, s_ost, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(opt_stack) :: s_ost

! local vars
    integer m

! begin
    lerr = .false.
    if (lform) then
       read (iunit,*,err=98) s_ost
    else
       read (iunit,err=98) s_ost
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_rd_stack

!------------------------------------------------------------------------------
! subroutine hdlc_wr_work
!
! Write the work arrays and L-BFGS history s_owr to iunit
!------------------------------------------------------------------------------

  subroutine hdlc_wr_work (iunit, s_owr, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    character*25 secthead
    type(opt_work) :: s_owr

! begin
    lerr = .false.
    if (lform) then
       write (iunit,*,err=98) s_owr%d
       write (iunit,*,err=98) s_owr%wlbfg
       write (iunit,*,err=98) s_owr%fx
       write (iunit,*,err=98) s_owr%oldfx
       write (iunit,*,err=98) s_owr%oldf
       write (iunit,*,err=98) s_owr%ooldf
       write (iunit,*,err=98) s_owr%eigval
       write (iunit,*,err=98) s_owr%oldeig
       write (iunit,*,err=98) s_owr%oldhss
       write (iunit,*,err=98) s_owr%u
       write (iunit,*,err=98) s_owr%oldu
       write (iunit,*,err=98) s_owr%svec
       write (iunit,*,err=98) s_owr%tvec
       write (iunit,*,err=98) s_owr%vmode
       write (iunit,*,err=98) s_owr%hessc
    else
       write (iunit,err=98) s_owr%d
       write (iunit,err=98) s_owr%wlbfg
       write (iunit,err=98) s_owr%fx
       write (iunit,err=98) s_owr%oldfx
       write (iunit,err=98) s_owr%oldf
       write (iunit,err=98) s_owr%ooldf
       write (iunit,err=98) s_owr%eigval
       write (iunit,err=98) s_owr%oldeig
       write (iunit,err=98) s_owr%oldhss
       write (iunit,err=98) s_owr%u
       write (iunit,err=98) s_owr%oldu
       write (iunit,err=98) s_owr%svec
       write (iunit,err=98) s_owr%tvec
       write (iunit,err=98) s_owr%vmode
       write (iunit,err=98) s_owr%hessc
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_wr_work

!------------------------------------------------------------------------------
! subroutine hdlc_rd_work
!
! Write the work arrays and L-BFGS history s_owr to iunit
!------------------------------------------------------------------------------

  subroutine hdlc_rd_work (iunit, s_owr, lform, lerr)

! args
    logical lform, lerr
    integer iunit
    type(opt_work) :: s_owr

! begin
! begin
    lerr = .false.
    if (lform) then
       read (iunit,*,err=98) s_owr%d
       read (iunit,*,err=98) s_owr%wlbfg
       read (iunit,*,err=98) s_owr%fx
       read (iunit,*,err=98) s_owr%oldfx
       read (iunit,*,err=98) s_owr%oldf
       read (iunit,*,err=98) s_owr%ooldf
       read (iunit,*,err=98) s_owr%eigval
       read (iunit,*,err=98) s_owr%oldeig
       read (iunit,*,err=98) s_owr%oldhss
       read (iunit,*,err=98) s_owr%u
       read (iunit,*,err=98) s_owr%oldu
       read (iunit,*,err=98) s_owr%svec
       read (iunit,*,err=98) s_owr%tvec
       read (iunit,*,err=98) s_owr%vmode
       read (iunit,*,err=98) s_owr%hessc
    else
       read (iunit,err=98) s_owr%d
       read (iunit,err=98) s_owr%wlbfg
       read (iunit,err=98) s_owr%fx
       read (iunit,err=98) s_owr%oldfx
       read (iunit,err=98) s_owr%oldf
       read (iunit,err=98) s_owr%ooldf
       read (iunit,err=98) s_owr%eigval
       read (iunit,err=98) s_owr%oldeig
       read (iunit,err=98) s_owr%oldhss
       read (iunit,err=98) s_owr%u
       read (iunit,err=98) s_owr%oldu
       read (iunit,err=98) s_owr%svec
       read (iunit,err=98) s_owr%tvec
       read (iunit,err=98) s_owr%vmode
       read (iunit,err=98) s_owr%hessc
    end if

! I/O failure jump point
    return
98  lerr = .true.
  end subroutine hdlc_rd_work

!------------------------------------------------------------------------------
! subroutine allocate_work
!
! Allocate space for the work arrays and the L-BFGS history s_owr
!------------------------------------------------------------------------------

  subroutine allocate_work (s_owr, nbig, nvar, nrem)

! args
    integer nbig, nvar, nrem
    type(opt_work) :: s_owr


    integer ierr

! begin
    if (ctrl%printl .ge. 2) write (stdout,'(A,I6,A,I4,A,I4,A)') &
         'Allocating memory, nbig = ', nbig, ', nvar = ', nvar, &
         ', nrem = ', nrem, ' ...'
    allocate (s_owr%d(nbig))
    allocate (s_owr%wlbfg(nbig*(2*nrem+1)+2*nrem), stat=ierr)
    if(ierr.ne.0)stop 'allocate in searchlib'

    allocate (s_owr%fx(nvar))
    allocate (s_owr%oldfx(nvar))
    allocate (s_owr%oldf(nvar*nvar))
    allocate (s_owr%ooldf(nvar))
    allocate (s_owr%eigval(nvar))
    allocate (s_owr%oldeig(nvar))
    allocate (s_owr%oldhss(nvar,nvar))
    allocate (s_owr%u(nvar,nvar))
    allocate (s_owr%oldu(nvar,nvar))
    allocate (s_owr%svec(nvar))
    allocate (s_owr%tvec(nvar))
    allocate (s_owr%vmode(nvar))
    allocate (s_owr%hessc(nvar*nvar))
    s_owr%uc => s_owr%u
    if (ctrl%printl .ge. 2) write (stdout,'(A)') 'Allocation successful'
  end subroutine allocate_work

!------------------------------------------------------------------------------
! subroutine deallocate_work
!
! Deallocate space of the work arrays and the L-BFGS history s_owr
!------------------------------------------------------------------------------

  subroutine deallocate_work (s_owr)

! args
    type(opt_work) :: s_owr

! begin
    nullify (s_owr%uc)
    deallocate (s_owr%d)
    deallocate (s_owr%wlbfg)
    deallocate (s_owr%fx)
    deallocate (s_owr%oldfx)
    deallocate (s_owr%oldf)
    deallocate (s_owr%ooldf)
    deallocate (s_owr%eigval)
    deallocate (s_owr%oldeig)
    deallocate (s_owr%oldhss)
    deallocate (s_owr%u)
    deallocate (s_owr%oldu)
    deallocate (s_owr%svec)
    deallocate (s_owr%tvec)
    deallocate (s_owr%vmode)
    deallocate (s_owr%hessc)
  end subroutine deallocate_work

end module searchlib
