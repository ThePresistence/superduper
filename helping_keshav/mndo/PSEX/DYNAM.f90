      SUBROUTINE DYNAM(ARRAY,LM5,ICALL,SCFCAL)

!==========================================================================C
! This subroutine implements the velocity verlet algorithm.                C
! It is divided into two main parts:                                       C
!  (1) Initialization                                                      C
!       all variables are set up, output files are initialized and         C
!       starting velocities are assigned                                   C
!  (2) MD loop                                                             C
!       the velocity verlet time propagation is performed                  C
!                                                                          C
! Units:                                                                   C
!   energy        kcal/mol                                                 C
!   time          ps                                                       C
!   distance      A                                                        C
!   velocity      A/ps                                                     C
!   temperature   K                                                        C
!   lin. moment.  ???????????                                              C
!   ang. moment.  ???????????                                              C
!--------------------------------------------------------------------------C
!  VARIABLES  (in alphabetical order)                                      C
!                                                                          C
! AMS           Atomic masses                                              C
! an_CC         Flag to activate analytical non-adiabatic coupling         C
!                computation                                               C
! ARRAY         Scratch array (for use in SCFCAL)                          C
! avstat        Flag to turn on statistics averaging                       C
!                If fstat = n > 1 the average value over n steps will be   C
!                saved                                                     C
! C             Electronic population (used in Tully's hopping)            C
! CICOMP        CI components                                              C
! CGX           Cartesian gradients of computed states                     C
! CNORMX        Cartesian norms of computed states                         C
! COORD         Coordinates                                                C
! de            Lexical indexes of CI  components                          C
! delta         Scratch variable for velocity update                       C
! dt            Time step [in picoseconds]                                 C
! E_av          Average total energy                                       C
! Einteg        Integration algorithm to use for electronic equation in    C
!                Tully's hopping                                           C
! E_ini         Total energy from first MD step                            C
! E_old         Total energy from previous MD step                         C
! E_tot         Total energy                                               C
! ELMNT         Chemical symbols for elements                              C
! ENERGX        Energies of computed states                                C
! adapt_ene     Threshold for relative variation of total energy         C
! ehrenfest     Flag to enable Ehrenfest dynamics procedure                C
! fhop          frequency for saving hopping data                          C
!                 = 1  every step                                          C
!                 = n  every n steps                                       C
! fsav          frequency for saving in trajectory and velocity file       C
!                 = 1  every step                                          C
!                 = n  every n steps                                       C
! fstat         frequency for saving in statistics file                    C
!                 = 1  every step                                          C
!                 = n  every n steps                                       C
! fvs           frequency for performing velocity scaling                  C
!                 = 1  every step                                          C
!                 = n  every n steps                                       C
! Hold          Effective hamiltonian for tully's hopping at previous      C
!                MD step                                                   C
! hopfile       Hopping data file name                                     C
! hopunit       Fortran unit for hopping data file                         C
! i             Counter                                                    C
! ic2           scratch variable for initialization of nran2 and nran3     C
! ICALL         control variable for SCF calculations and error flag       C
! ICICAL        Control flag for SCFMUL                                    C
! IFMAP         Flag for mapping failure (when IMOMAP=2)                   C
! ihop          Control flag for hopping. If .TRUE. hopping happend        C
! IGRST         States for which the gradient is requested                 C
! init_stat     Initial state                                              C
! io_disk       Error flag for I/O                                         C
! iout          Level of output (=0 none; =1 minimal; =2 big; =3 debug)    C
! istat         Current state to propagate                                 C
! IN2           MNDO options                                               C
! K             Kinetic energy                                             C
! K_av          Average kinetic energy                                     C
! L2_av         Average square angular momentum                            C
! Lx_av,Ly_av,Lz_av  Average cartesian components of angular momentum      C
! LM1           Maximum number of atoms (see LIMIT.f90)                    C
! LM1M          Maximum number of external points (see LIMIT.f90)          C
! LM5           Buffer lenght (for use in SCFCAL)                          C
! LMGRD         Maximum number of GUGA-CI gradient vectors (see LIMIT.f90) C
! LMV           Muximum number of geometrical variables (see LIMIT.f90)    C
! mdstep        current MD step                                            C
! naco          number of active orbitals in CI                            C
! NAT           Atomic numbers                                             C
! NCICONF       Number of CI configurations                                C
! NE            Number of steps for electronic integration (in Tully's hop)C
! nf            Number of degrees of freedom                               C
! norot         Flag to select initial velocities such as angular momentum C
!                is zero                                                   C
! notra         Flag to select initial velocities such as linear momentum  C
!                is zero                                                   C
! nstep         Max number of MD steps                                     C
! NUMAT         Number of atoms                                            C
! num_CC        Flag to activate numerical non-adiabatic coupling          C
!                 computation                                              C
! old_CICOMP    CI components in the previous MD step                      C
! old_de        Lexical indexes in the previous MD step                    C
! old_ENERGX    Excited states energies in the previous MD step            C
! P2_av         Average square linear momentum                             C
! Px_av,Py_av,Pz_av  Average cartesian componets of linear momentum        C
! rdold         Dot product of velocity and non-adiabatric vectors         C
!                at previosu MD step                                       C 
! rdummy        scratch variable for initialization of nran2 and nran3     C
! restart       Flag to allow restarting form file                         C
! rfile         Restart file name                                          C
! rnd_gen       Algorithm for the generation of random numbers in tully's  C
!               hopping:                                                   C
!                 standard  -> standard FORTRAN algorithm [RANDOM_NUMBER]  C
!                 PM_BD     -> Park and Miller algorithm with Bays-Durham  C
!                                shuffle                                   C
!                 knuth     -> Knuth subtractive algorithm                 C
! rnd_seed                                                                 C
! sca           conversion factor from kcal/(mol*A*AMU) to A/ps^2          C
! SCFCAL        External routine for SCF calculations                      C
! sfile         Statistics file name                                       C
! simhop        Variable to control simple hopping                         C
!                  = 0  No simple hopping                                  C
!                  = 1,2,3,4,5,6   Simple hopping                          C
! sunit         Fortran unit for statistics file                           C
! T_av          Average temperature                                        C
! Tfix          Temperature to be reached with velocity scaling            C
! T_ist         Instantaneous temperature                                  C
! temp0         Initial temperature                                        C
! tully_hop     Flag to enable Tully's hopping procedure                   C
! tully_rest    Flag to activate restart of tully hopping                  C
! v_av          Average potential energy                                   C
! vfile         Velocity file name                                         C
! vunit         Fortran unit for velocity file                             C
! vs            Flag to enable velocity scaling                            C
! vs_Emax       Maximum kinetic energy correction allowed in velocity      C
!                scaling procedure. Negative number means unlimited        C
! vx,vy,vz      Cartesian components of velocity                           C
! xfile         Trajectory file name                                       C
! xunit         Fortran unit for trajectory file                           C
! write_hop     Flag to allow writing hopping data to file                 C
! write rest    Flag to allow writing the restart file                     C
! write_stats   Flag to allow writing of statistics                        C
! write_traj    Flag to allow writing trajectory file                      C
! write_vel     Flag to allow writing of velocities                        C
! dec_cor       Decoherence correction by Persico & Granucci               C
!               (0.0=off, standard value=0.1)                              C
! nac_ph        threshold for nonadiabatic coupling phase following        C
!               algorithm (0.0=off, standard=0.1)                          C
! cuthop        cut-off Energy difference to avoid unphysical hops         C
!               >0 : simple cut at cuthop value                            C
!               <0 : cut by gaussian scale of probability cenetered at     C
!               cuthop value                                               C
!--------------------------------------------------------------------------C
!  CALLINGS                                                                C
!                                                                          C
!  nmlst                     Manadgment of name list variables             C
!  read_rest                 Read restart file                             C
!  v_rnd                     create random velocities                      C
!  v_remove_tra              Removes global translation from velocity set  C
!  v_remove_rot              Removes global rotations fom velocity set     C
!  term_scal                 perform termal scaling                        C
!  SCFMUL                    calculates energy and gradients               C
!  Kene                      computes kinetic energy                       C
!  Tist                      computes intantaneous temperature             C
!  outpr                     prints initial statistics                     C
!  RANDOM_SEED               initialize FORTRAN random number generator    C
!  system_clock              reads system clock (used for initialization   C
!                             of nran2 and nran3)                          C
!  nran2                     Random number generator                       C
!  nran3                     Random number generator                       C
!  SURHOP                    Routine for Tully's hopping                   C
!  ehrendyn                  Routine for Ehrenfest dynamics                C
!  estat_fol                 Mapping of MD active state                    C
!  vel_scal                  Performs velocity scaling                     C
!  outpr2                    prepares and prints statistics                C
!  write_res                 writes restart file                           C
!                                                                          C
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       C
! MAIL: efabiano@mpi-muelheim.mpg.de                                       C
! DATE: March 2007                                                         C
!==========================================================================C

      USE LIMIT, ONLY: LM1, LMGRD, LM1M, LMV, LMCONF

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: LM5
      REAL*8, DIMENSION(LM5) :: ARRAY
      INTEGER :: ICALL

!  COMMON
      REAL*8 :: AMS
      REAL*8 :: COORD
      INTEGER :: NUMAT
      INTEGER :: IN2
      INTEGER :: IGRST
      INTEGER :: NAT
      CHARACTER(2) :: ELEMNT
      INTEGER :: NCICONF
      REAL*8 :: CICOMP
      INTEGER :: de
      INTEGER :: naco
      INTEGER :: IFMAP
      LOGICAL :: NEWMAP
      INTEGER :: IMOMAP
      REAL*8 :: ctime_v
      INTEGER :: wtime_v
      REAL*8 :: ctime_g
      INTEGER :: wtime_g
      REAL*8 :: ctime_sf
      INTEGER :: wtime_sf
      REAL*8 :: ctime_acc
      INTEGER :: wtime_acc
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON /AMASS/ AMS(LM1)
      COMMON /ATOMC / COORD(3,LM1)
      COMMON /ATOMS / NUMAT,NAT(LMV)
      COMMON /INOPT2/ IN2(300)
      COMMON /GRDORG/ IGRST(LMGRD)
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /DYNVA/ CICOMP(LMCONF,6),NCICONF
      COMMON /DYNVA2/ de(40,LMCONF),naco
      COMMON /CIMAP / IFMAP,NEWMAP
      COMMON /ctimemd/ ctime_v,ctime_g,ctime_sf,ctime_acc
      COMMON /wtimemd/ wtime_v,wtime_g,wtime_sf,wtime_acc
      COMMON /ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
      COMMON /wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

      INTEGER :: wait_steps
!  LOCAL VARIABLES
      INTEGER :: i,j
      INTEGER :: io_disk
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz
      REAL*8, DIMENSION(NUMAT) :: vx_temp,vy_temp,vz_temp
      REAL*8, DIMENSION(NUMAT) :: hop_vx,hop_vy,hop_vz
      REAL*8, DIMENSION(3,LM1) :: COORD_temp
      INTEGER :: nf,adapt_try,mapthr_temp,step_tot,step_corr,ICALL_temp
      INTEGER :: IDIIS,KTRIAL,NSTART,IMOMAP_old
      REAL*8 :: K
      REAL*8 :: T_ist,t_tot,t_goal,dt_temp
      REAL*8 :: t_tot_temp,t_goal_temp
      REAL*8, DIMENSION(LMGRD) :: ENERGX
      REAL*8, DIMENSION(3,LM1+LM1M,LMGRD) :: CGX
      REAL*8, DIMENSION(3,LM1+LM1M,LMGRD) :: CGX_temp
      REAL*8, DIMENSION(LMGRD) :: CNORMX
      INTEGER :: ICICAL
      REAL*8 :: E_ini, E_tot, E_therm
      REAL*8 :: E_old, E_K_old, E_P_old, E_therm_old
      REAL*8 :: E_old_temp, E_K_old_temp, E_P_old_temp, E_therm_old_temp
      REAL*8 :: delta
      REAL*8 :: E_av
      REAL*8 :: K_av
      REAL*8 :: V_av
      REAL*8 :: T_av
      REAL*8 :: P2_av
      REAL*8 :: Px_av
      REAL*8 :: Py_av
      REAL*8 :: Pz_av
      REAL*8 :: L2_av
      REAL*8 :: Lx_av
      REAL*8 :: Ly_av
      REAL*8 :: Lz_av
      LOGICAL :: ihop
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: Hold
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: Hold_temp
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rdold
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rdold_temp
      COMPLEX*16, DIMENSION(IN2(159)) :: C
      COMPLEX*16, DIMENSION(IN2(159)) :: C_temp
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP_temp
      REAL*8, DIMENSION(LMGRD,3) :: old_ENERGX
      REAL*8, DIMENSION(LMGRD,3) :: old_ENERGX_temp
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      INTEGER, DIMENSION(40,LMCONF) :: old_de_temp
      LOGICAL :: tully_rest, step_fail, adaptive, THERE
      REAL*8, DIMENSION(3,LM1+LM1M) :: ehrgr
      REAL*8, DIMENSION(3,LM1+LM1M) :: temp_grad
      INTEGER, DIMENSION(IN2(159)) :: sw_CC
      INTEGER, DIMENSION(IN2(159)) :: inde
      INTEGER, DIMENSION(IN2(159)) :: old_inde
      INTEGER, DIMENSION(IN2(159)) :: old_inde_temp
      INTEGER :: istat, istatt, mdstep, startStep
      INTEGER :: istat_temp,istatt_temp
      INTEGER, ALLOCATABLE :: rnd_seed_array(:)
      REAL*8, DIMENSION(:), ALLOCATABLE :: thermo_eta, thermo_p
      REAL*8, DIMENSION(:), ALLOCATABLE :: thermo_eta_temp, thermo_p_temp
      INTEGER :: rsa_size
      INTEGER :: ic2
      REAL*8 :: rdummy
      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir
      REAL*8 :: ctime3, ctime4
      INTEGER :: wtime3, wtime4
      CHARACTER(150) :: str

!  NAMELIST VARIABLES
      INTEGER :: iout
      INTEGER :: nstep
      INTEGER :: jst
      REAL*8 :: dt
      REAL*8 :: temp0
      LOGICAL :: norot
      LOGICAL :: notra
      INTEGER :: init_stat
      LOGICAL :: write_traj
      INTEGER :: xunit
      CHARACTER(20) :: xfile
      LOGICAL :: write_vel
      INTEGER :: vunit
      CHARACTER(20) :: vfile
      LOGICAL :: ene_cons
      REAL*8 :: ene_scal
      REAL*8 :: adapt_ene
      INTEGER :: adapt_ene_mode
      INTEGER :: adapt_map
      INTEGER :: adapt_tries
      LOGICAL :: adapt_proceed
      INTEGER :: fsav
      LOGICAL :: write_stats
      INTEGER :: sunit
      CHARACTER(20) :: sfile
      INTEGER :: fstat
      LOGICAL :: restart
      LOGICAL :: write_rest
      CHARACTER(20) :: rfile
      LOGICAL :: avstat
      LOGICAL :: vs
      REAL*8 :: Tfix
      INTEGER :: fvs
      REAL*8 :: vs_Emax
      LOGICAL :: ehrenfest
      LOGICAL :: write_hop
      INTEGER :: hopunit
      CHARACTER(20) :: hopfile
      INTEGER :: fhop
      LOGICAL :: tully_hop
      INTEGER :: NE
      LOGICAL :: num_CC
      LOGICAL :: an_CC
      CHARACTER(5) :: Einteg
      CHARACTER(8) :: rnd_gen
      INTEGER :: rnd_seed
      INTEGER :: rnd_count
      INTEGER :: simhop
      LOGICAL :: fol_stat
      REAL*8  :: dec_cor
      REAL*8  :: nac_ph
      REAL*8  :: cuthop
      REAL*8 :: conv
      REAL*8 :: g0
      INTEGER :: thermostat
      INTEGER :: thermo_len
      INTEGER :: thermo_equi
      REAL*8 :: thermo_tau
      REAL*8 :: thermo_T

      CHARACTER(20) :: NAMEVAR

!  EXTERNAL ROUTINeS
      EXTERNAL SCFCAL

!  DATA
      REAL*8, PARAMETER :: sca = 418.4 ! from kcal/(mol*A*AMU) to A/ps^2

      DATA ICICAL/0/
      SAVE ICICAL

!============================================================================C

      IF(IN2(205).GT.0) THEN
         WRITE(*,'(/1X,"ATOMIC MASSES IN USE FOR MD:")')
         WRITE(*,'(/5X,"I    Z    MASS(AMU)")')
         DO I=1, NUMAT
            WRITE(*,'(1X,2I5,F12.5)') I, NAT(I), AMS(I)
         END DO
      END IF

! TIMING - overall time, start
      CALL CPU_TIME(ctime1)
      CALL SYSTEM_CLOCK(wtime1,check_ir)

!  Initialization
      IMOMAP = IN2(156)

      ihop = .FALSE.
      t_tot = 0.D0
      step_tot = 0
      step_corr = 1
      E_tot = 0.D0
      E_ini = 0.D0
      E_old = 0.D0
      E_K_old = 0.D0
      E_P_old = 0.D0
      E_therm_old = 0.D0
      E_therm = 0.D0
      T_ist = 0.D0
      vx = 0.d0
      vy = 0.d0
      vz = 0.d0
      hop_vx = 0.d0
      hop_vy = 0.d0
      hop_vz = 0.d0
      E_av = 0.D0
      K_av = 0.D0
      V_av = 0.D0
      T_av = 0.D0
      P2_av = 0.D0
      Px_av = 0.D0
      Py_av = 0.D0
      Pz_av = 0.D0
      L2_av = 0.D0
      Lx_av = 0.D0
      Ly_av = 0.D0
      Lz_av = 0.D0
      tully_rest = .FALSE.
      rnd_count = 0
      IFMAP = 0
      NEWMAP = .FALSE.
      ctime_v = 0.D0
      wtime_v = 0
      ctime_g = 0.D0
      wtime_g = 0
      ctime_sf = 0.D0
      wtime_sf = 0
      ctime_acc = 0.D0
      wtime_acc = 0
      ctime_eh = 0.D0
      wtime_eh = 0
      ctime_eint = 0.D0
      wtime_eint = 0
      ctime_ncc = 0.D0
      wtime_ncc = 0

!  Namelist initialization
      CALL nmlst(iout,nstep,dt,temp0,norot,notra,init_stat,write_stats,sunit,sfile,write_traj,xunit,xfile,write_vel,vunit,vfile,adapt_ene,fstat,avstat,fsav,restart,write_rest,rfile,vs,Tfix,fvs,vs_Emax,ehrenfest,tully_hop,Einteg,NE,num_CC,an_CC,rnd_gen,rnd_seed,write_hop,hopunit,hopfile,fhop,simhop,fol_stat,dec_cor,nac_ph,cuthop,ene_cons,ene_scal,adapt_map,adapt_tries,adapt_proceed,adapt_ene_mode,thermostat,thermo_len,thermo_tau,thermo_T,thermo_equi,g0)

      Allocate(thermo_eta(thermo_len), thermo_p(thermo_len))
      Allocate(thermo_eta_temp(thermo_len), thermo_p_temp(thermo_len))
      thermo_eta = 0.D0
      thermo_p = 0.D0
      thermo_eta_temp = 0.D0
      thermo_p_temp = 0.D0

!  Read restart file if requested
      IF (restart) THEN
         CALL read_rest(iout,rfile,NUMAT,COORD(1,:),COORD(2,:),COORD(3,:),vx,vy,vz,tully_hop,IN2(159),init_stat,C,Hold,rdold,old_ENERGX(:,3),old_CIcomp,old_de,tully_rest,rnd_count)

         IF (NUMAT .EQ. 1) THEN
            nf = 3
         ELSE
            IF (norot .AND. notra) THEN
               IF (NUMAT .EQ. 2) THEN
                  nf = 3*NUMAT -5
               ELSE
                  nf = 3*NUMAT - 6
               END IF
            ELSE IF (notra .AND. .NOT.norot) THEN
               nf = 3*NUMAT - 3
            ELSE IF (norot .AND. .NOT. notra) THEN
               IF (NUMAT .EQ. 2) THEN
                  nf = 3*NUMAT - 2
               ELSE
                  nf = 3*NUMAT - 3
               END IF
            ELSE IF (.NOT.norot .AND. .NOT.notra) THEN
               nf = 3*NUMAT
            END IF
         END IF
      END IF

!......................................................................C

! check adaptive input
      adaptive = (adapt_map.GT.0 .OR. adapt_ene.GT.0.D0)
      IF (adaptive) THEN
         CALL SYSTEM("rm -f fort.11*")
         CALL SYSTEM("rm -f hop.tmp*")
         CALL SYSTEM("rm -f imomap.tmp*")
         CALL SYSTEM("touch hop.tmp")
         IMOMAP_old = IMOMAP
         if (IMOMAP .NE. 3) then
            CALL SYSTEM("rm -f imomap.dat")
         end if
         IF (adapt_map.GT.0 .or. IMOMAP.GT.0) THEN
            IF (IN2(77).LE.5) THEN
               IN2(156)=3
               IMOMAP=3
               WRITE(*,*) 
               WRITE(*,'(A)') "Set imomap=3 for adaptive procedure."
               WRITE(*,*) 
            END IF
         END IF
      END IF

! check decoherence input
      IF((dec_cor.gt.0).and.((Einteg.ne."UP1").and.(Einteg.ne."UP2").and.(Einteg.ne."UP3"))) THEN
         write(*,*) 'Decoherence correction is currently onlypossible with the unitary propagators (UP1-3)'
         STOP
      END IF

!     initialize random number generator
!TWK  for both Tully hopping and velocity assignment for all methods
!TWK  Initialise with the user-specified seed if available, otherwise
!TWK  use a random seed as before.
!TWK  (This section moved before velocity assignment so we can use the same 
!TWK  seeding throughout to ensure reproducibility for testing purposes)

      IF (rnd_gen .EQ. "standard") THEN
         CALL RANDOM_SEED(SIZE=rsa_size)
         ALLOCATE(rnd_seed_array(rsa_size))
         IF (rnd_seed .EQ. -1) THEN
            CALL RANDOM_SEED
            CALL RANDOM_SEED(GET=rnd_seed_array(1:rsa_size))
         ELSE
            rnd_seed_array = rnd_seed
            CALL RANDOM_SEED(PUT=rnd_seed_array(1:rsa_size))
         ENDIF
         IF (iout .GE. 0) THEN
            WRITE(*,*)
            WRITE(*,*) 'Random number seeds: ', rnd_seed_array
            WRITE(*,*)
         ENDIF
         DEALLOCATE(rnd_seed_array)
      ELSE IF (rnd_gen .EQ. "PM_BD" .OR. rnd_gen .EQ. "knuth") THEN
         IF (rnd_seed .EQ. -1) THEN
            CALL system_clock(count=ic2)
         ELSE
            ic2 = rnd_seed
         ENDIF
         IF (iout .GE. 0) THEN
            WRITE(*,*)
            WRITE(*,*) 'Random number seed: ', ic2
            WRITE(*,*)
         ENDIF
         IF (rnd_gen .EQ. "PM_BD") THEN
            CALL nran2(rdummy,ic2,.TRUE.)
         ELSE
            CALL nran3(rdummy,ic2,.TRUE.)
         ENDIF
      ELSE
         WRITE(*,*)
         WRITE(*,'("Unknown random numbers generator ",A8)')rnd_gen
         WRITE(*,*)
         STOP 'Abnormal termination in sub. DYNAM'
      END IF

!TWK  If reproducibility required (i.e. rnd_seed != -1), 
!TWK  we must cycle the random number generator until it reaches
!TWK  the same state it was in when the restart file was written.
!TWK  NB: Only necessary for Tully hopping, as the only other random
!TWK  number used is in the velocity assignment, which is not called
!TWK  after a restart
      IF (restart .and. tully_rest .and. rnd_seed .ne. -1) THEN
         DO i = 1, rnd_count
            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(rdummy)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(rdummy,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(rdummy,1,.FALSE.)
            END IF
         END DO
      END IF

!......................................................................C

! TIMING - assignment of velocity, start
      IF (iout .GE. 2 ) THEN
         CALL CPU_TIME(ctime3)
         CALL SYSTEM_CLOCK(wtime3,check_ir)
      END IF

! Assignment of starting velocity
!  If temp0 < 0 stop
      IF (.NOT. restart) THEN
         IF (temp0 .LT. 0.D0) THEN
            WRITE(*,*)
            WRITE(*,*)"ERROR!!"
            WRITE(*,*)"Initial temperature set to a negative value"
            WRITE(*,*)
            STOP 'Abnormal termination in sub. DYNAM'
         END IF
      END IF

!  if we start from restart file we don't need to do anything here
!TWK removed temp0=0.d0 test so nf is set correctly and prevent
!TWK floating point strangeness in the test
      IF (.NOT. RESTART) THEN

         IF (NUMAT .EQ. 1) THEN
            WRITE(*,*)
            WRITE(*,*)"Single atom velocity assignment"
            WRITE(*,*)
            CALL v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)
            nf = 3
         ELSE
            IF (norot .AND. notra) THEN
               CALL v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)
               CALL v_remove_tra(iout,NUMAT,AMS,vx,vy,vz)
               CALL v_remove_rot(iout,NUMAT,AMS,COORD(1,:),COORD(2,:),COORD(3,:),vx,vy,vz)
               IF (NUMAT .EQ. 2) THEN
                  nf = 3*NUMAT -5
               ELSE
                  nf = 3*NUMAT - 6
               END IF
            ELSE IF (notra .AND. .NOT.norot) THEN
               CALL v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)
               CALL v_remove_tra(iout,NUMAT,AMS,vx,vy,vz)
               nf = 3*NUMAT - 3
            ELSE IF (norot .AND. .NOT. notra) THEN
               CALL v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)
               CALL v_remove_rot(iout,NUMAT,AMS,COORD(1,:),COORD(2,:),COORD(3,:),vx,vy,vz)
               IF (NUMAT .EQ. 2) THEN
                  nf = 3*NUMAT - 2
               ELSE
                  nf = 3*NUMAT - 3
               END IF
            ELSE IF (.NOT.norot .AND. .NOT.notra) THEN
               CALL v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)
               nf = 3*NUMAT
            END IF
         END IF
      
!     Perform initial termal scaling

         CALL term_scal(iout,NUMAT,AMS,vx,vy,vz,temp0,nf)

      END IF

!   Output
      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,'("+-----------------------------------------------+")')
         WRITE(*,'("|",8X,"Initial velocities assignement",9X,"|")')
         WRITE(*,'("+-----------------------------------------------+")')
         IF (.NOT. restart) THEN
            WRITE(*,'("| Initial temperature: ",F10.4," K",13X,"|")')temp0
            WRITE(*,'("| Degrees of fredom: ",I3,24X,"|")')nf
         ELSE
           WRITE(*,'("| Data from restart file ",A)')rfile
         END IF 
         WRITE(*,'("+-----------------------------------------------+")')
         WRITE(*,'("|",15X,"Velocities [A/ps]",15X,"|")')
         WRITE(*,'("| Atom",8X,"vx",12X,"vy",12X,"vz",4X,"|")')
         DO i = 1, NUMAT
            WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,'("+-----------------------------------------------+")')
      END IF

      IF (iout .GE. 2) THEN
! TIMING - velocity assignment, end
         CALL CPU_TIME(ctime4)
         CALL SYSTEM_CLOCK(wtime4,check_ir)
         wtime_v = wtime_v + wtime4-wtime3
            IF (ctime4 .GT. REAL(wtime4-wtime3)/REAL(check_ir)+ctime3)ctime4 = REAL(wtime4-wtime3)/REAL(check_ir)+ctime3
         ctime_v = ctime_v + ctime4-ctime3
      END IF

! TIMING - output
      IF (iout .GE. 3) THEN
            WRITE(*,*)
            WRITE(*,'("-----------------------------------")')
            WRITE(*,*)"*** TIMING ***"
            WRITE(*,*)"Assignment of starting velocities"
            WRITE(*,*)
            CALL out_time(ctime3,ctime4,wtime3,wtime4,check_ir)
            WRITE(*,'("-----------------------------------")')
      END IF

!.....................................................................C 

! TIMING - Gradient computation, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime3)
         CALL SYSTEM_CLOCK(wtime3,check_ir)
      END IF

!  Computation of the gradient for inizialization
      CALL SCFMUL(ARRAY,LM5,ICALL,SCFCAL,ENERGX,CGX,CNORMX,NUMAT)
      ICICAL = ICICAL+1
         
!  Check for convergency problems
!TWK Also check the case of mapping failure when IMOMAP=2
      IF (ICALL .EQ. -1) THEN
         IF (IMOMAP .GE. 2 .AND. IFMAP .EQ. 1) THEN
            WRITE(*,*)
            WRITE(*,*) "ACTIVE ORBITAL MAPPING FAILURE"
         ELSE 
            WRITE(*,*)
            WRITE(*,*)"CONVERGENCY PROBLEMS!!!"
         END IF
         WRITE(*,*)"The MD simulation will end here."
         WRITE(*,*)
         STOP 'Abnormal termination in sub. DYNAM'
      END IF
      
      IF (iout .GE. 2) THEN
! TIMING - Gradient computation, end
         CALL CPU_TIME(ctime4)
         CALL SYSTEM_CLOCK(wtime4,check_ir)
         wtime_g = wtime_g + wtime4-wtime3
            IF (ctime4 .GT. REAL(wtime4-wtime3)/REAL(check_ir)+ctime3)ctime4 = REAL(wtime4-wtime3)/REAL(check_ir)+ctime3
         ctime_g = ctime_g + ctime4-ctime3
      END IF

! TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Gradient computation"
         WRITE(*,*)
         CALL out_time(ctime3,ctime4,wtime3,wtime4,check_ir)
         WRITE(*,'("-----------------------------------")')
      END IF

!  Computation of total energy
      CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
      CALL Tist(iout,K,nf,T_ist)

      E_tot = ENERGX(init_stat) + K

      E_K_old = K
      E_P_old = ENERGX(init_stat)
      E_old = E_tot
      E_ini = E_tot

      IF (iout .GE. 1) CALL outpr(NUMAT,nf,AMS,vx,vy,vz,E_tot,K,ENERGX(init_stat))

!....................................................................C

! Initialization of output files

! Open file for output of statistics if requested
      IF (write_stats) THEN
         OPEN(UNIT=sunit,NAME=sfile,STATUS="REPLACE",ACTION="WRITE",IOSTAT=io_disk)
         IF (io_disk .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,'("File ",A," already exist.")')sfile
            WRITE(*,'("Please change the name of output file")')
            WRITE(*,*)
            STOP 'Abnormal termination in sub. DYNAM'
         END IF
         WRITE(sunit,'("#mdstep",2X,"total energy",2X,"kinetic energy",2X,"potent. energy",2X,"temperature",8X,"p**2",12X,"Px",13X,"Py",13X,"Pz",12X,"L**2",12X,"Lx",13X,"Ly",13X,"Lz",7X,"thermostat energy")')
      END IF

!  Open files for output of trajectory if requested
      IF (write_traj) THEN
         OPEN(UNIT=xunit,NAME=xfile,STATUS="REPLACE",ACTION="WRITE",IOSTAT=io_disk)
         IF (io_disk .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,'("File ",A," already exist.")')xfile
            WRITE(*,'("Please change the name of output file")')
            WRITE(*,*)
            STOP 'Abnormal termination in sub. DYNAM'
         END IF
         WRITE(xunit,'(a)')"[Molden Format]"
         WRITE(xunit,'(a)')"[Atoms] Angs"
         DO i = 1, NUMAT
            WRITE(xunit,'(1X,A2,2X,I3,2X,I3,2X,F10.5,2X,F10.5,2X,F10.5)')ELEMNT(NAT(i)),i,NAT(i),COORD(1,i),COORD(2,i),COORD(3,i)
         END DO
         WRITE(xunit,'(a)')"[GEOMETRIES] XYZ"
      END IF

!  Open file for output of velocity if requested
      IF (write_vel) THEN
         OPEN(UNIT=vunit,NAME=vfile,STATUS="REPLACE",ACTION="WRITE",IOSTAT=io_disk)
         IF (io_disk .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,'("File ",A," already exist.")')vfile
            WRITE(*,'("Please change the name of output file")')
            WRITE(*,*)
            STOP 'Abnormal termination in sub. DYNAM'
         END IF
      END IF

!  Open file for output of hopping data if requested
      IF (write_hop) THEN
         OPEN(UNIT=hopunit,NAME=hopfile,STATUS="REPLACE",ACTION="WRITE",IOSTAT=io_disk)
         IF (io_disk .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,'("File ",A," already exist.")')hopfile
            WRITE(*,'("Please change the name of output file")')
            WRITE(*,*)
            STOP 'Abnormal termination in sub. DYNAM'
         END IF
      END IF
!....................................................................C

!  Inizialization for MD loop
      istat = init_stat

      ehrgr(:,:) = CGX(:,:,istat)

!     Output
      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)"***** STARTING MD loop ****"
         WRITE(*,'("Initial state for dynamics: ",I1)')istat
         WRITE(*,*)
         WRITE(*,*)
      END IF

!....................................................................C
      IF (adaptive) THEN
         IF (adapt_map) THEN
            WRITE(*,*)
            WRITE(*,'(A,I3,A)') "Set mapthr=",adapt_map," for adaptive procedure."
            WRITE(*,*)
            IN2(166)=adapt_map
         ELSE
            IF (IMOMAP_old .EQ. 1) IN2(166)=1
         END IF
         mapthr_temp = IN2(166)
      END IF

      istatt = istat

      adapt_try=0
      wait_steps=0
      IF (adaptive) THEN
!           Save all the stats for future restoration
         vx_temp = vx
         vy_temp = vy
         vz_temp = vz
         COORD_temp = COORD
         E_old_temp = E_tot
         E_K_old_temp = K
         E_P_old_temp = ENERGX(init_stat)
         E_therm_old_temp = 0.0D0
         CGX_temp = CGX
         Hold_temp = Hold
         rdold_temp = rdold
         old_CICOMP_temp = old_CICOMP
         old_ENERGX_temp = old_ENERGX
         old_de_temp = old_de
         old_inde_temp = old_inde
         C_temp = C
         istat_temp = istat
         istatt_temp = istatt
         t_tot_temp = t_tot
         t_goal_temp = t_tot + dt
         thermo_eta_temp = thermo_eta
         thermo_p_temp = thermo_p
!           Copy MO mapping to temporary file
         if (IMOMAP .EQ. 3) then
            CALL SYSTEM("cp -f imomap.dat imomap.tmp1")
         end if
         CALL SYSTEM("cp -f fort.11 fort.11.tmp1")
      END IF

!-------------------------C
!  Now starts the MD loop C
!-------------------------C
      dt_temp = dt
      ICALL_temp = ICALL
      IDIIS = IN2(9)
      KTRIAL = IN2(67)
      NSTART = IN2(78)
      
      IF (thermostat .le. 0) thermo_equi=0
      startStep = -thermo_equi+1

      DO mdstep = startStep, nstep

         t_goal = t_tot + dt

         IF (write_traj .AND. MOD(mdstep,fsav).EQ.0 .AND. mdstep .GE. 1) THEN
            WRITE(xunit,*)NUMAT
            WRITE(xunit,*)
            DO i = 1, NUMAT
               WRITE(xunit,'(A2,3X,F13.8,3X,F13.8,3X,F13.8)')ELEMNT(NAT(i)), COORD(1,i),COORD(2,i), COORD(3,i)
            END DO
         END IF
      
!     Velocity
         IF (write_vel .AND. MOD(mdstep,fsav).EQ.0 .AND. mdstep .GE. 1) THEN 
            DO i = 1, NUMAT
                  WRITE(vunit,'(F13.8,3X,F13.8,3X,F13.8)')vx(i),vy(i),vz(i)
            END DO
            WRITE(vunit,*)
         END IF

         jst = 0
         DO WHILE (t_tot .LT. t_goal)

            jst = jst + 1
            step_tot = step_tot + 1
            ICALL = ICALL_temp
            IN2(9) = IDIIS
            IN2(67) = KTRIAL
            IN2(78) = NSTART

            IFMAP=0
            if (wait_steps .LE. 0 .and. adapt_try.gt.0) THEN
               adapt_try = adapt_try-1
            END IF
            dt_temp = dt/(2.0D0**adapt_try)
            IF (t_tot + dt_temp .GT. t_goal) dt_temp = t_goal - t_tot
 

            IF (iout .GE. 0) THEN
               WRITE(*,*)
               WRITE(*,'(92("*"))')
               WRITE(*,*)
               WRITE(*,'("MD step:            ",I5)')mdstep
            END IF

            IF (iout .GE. 1) THEN
               IF (adaptive) THEN
                  WRITE(*,'("Attempted substeps: ",I5)') step_tot
                  WRITE(*,'(" Accepted substeps: ",I5)') step_corr
               END IF
               WRITE(*,*)
               WRITE(*,'(A)') "Simulation time:"
               WRITE(*,'(A,F15.8,A)') "  Act time:   ",t_tot*1000.D0," fs"
               WRITE(*,'(A,F15.8,A)') "  Next goal:  ",t_goal*1000.D0," fs"
               WRITE(*,'(A,F15.8,A)') "  Timestep:   ",dt_temp*1000.D0," fs"
               IF (adaptive) THEN
                  WRITE(*,'(A,F15.8,A)') "  Mean attempt step: ",(t_tot+dt_temp)/step_tot*1000.D0," fs"
                  WRITE(*,'(A,F15.8,A)') "  Mean correct step: ",(t_tot+dt_temp)/step_corr*1000.D0," fs"
               END IF
            END IF

            IF (adaptive .AND. adapt_proceed) THEN
!           Make sure that MO mapping is accepted at the last try
               IF (adapt_try.EQ.adapt_tries) THEN
                  IN2(166) = 1
               ELSE
                  IN2(166) = mapthr_temp
               END IF
            END IF

            IF (g0.gt.0 .and. mdstep.eq.1) THEN
               IF (ehrenfest) THEN
                  call langevin(AMS, vx, vy, vz, g0, temp0, ehrgr, NUMAT, dt, mdstep, jst)
               ELSE
                  call langevin(AMS, vx, vy, vz, g0, temp0, CGX(1,1,istatt), NUMAT, dt, mdstep, jst)
               ENDIF
            ENDIF

!     Update velocity to half step forward
!     and position to one step forward
            IF (thermostat .EQ. 1) THEN
               call noseHoover_advance(iout,dt_temp,vx,vy,vz,AMS,NUMAT,thermo_eta,thermo_p,thermo_tau,thermo_len,thermo_T,nf,E_therm,.false.)
            END IF
            IF (ehrenfest) THEN
               temp_grad = ehrgr(:,:)
            ELSE
               temp_grad = CGX(:,:,istatt)
            END IF
            DO i = 1, NUMAT
               delta = sca*0.5D0*dt_temp/AMS(i)
               vx(i) = vx(i) - temp_grad(1,i)*delta
               vy(i) = vy(i) - temp_grad(2,i)*delta
               vz(i) = vz(i) - temp_grad(3,i)*delta   
            END DO

            DO i = 1, NUMAT
               COORD(1,i) = COORD(1,i) + vx(i)*dt_temp
               COORD(2,i) = COORD(2,i) + vy(i)*dt_temp
               COORD(3,i) = COORD(3,i) + vz(i)*dt_temp
            END DO

            IF (iout .GE. 3) THEN
               WRITE(*,*)
               WRITE(*,'("+-----------------------------------------------+")')
               WRITE(*,'("|",3X,"MD step: ",I4," - First update of velocity",4X,"|")')mdstep
               WRITE(*,'("|",47X,"|")')
               WRITE(*,'("| Atom",8X,"vx",12X,"vy",12X,"vz",4X,"|")')
               DO i = 1, NUMAT
                  WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')i,vx(i),vy(i),vz(i)
               END DO
               WRITE(*,'("+-----------------------------------------------+")')
               WRITE(*,*)istatt
               WRITE(*,'("+-----------------------------------------------+")')
               WRITE(*,'("|",5X,"MD step: ",I4," - Update of coordinates",5X,"|")')mdstep
               WRITE(*,'("|",47X,"|")')
               WRITE(*,'("| Atom",9X,"x",13X,"y",13X,"z",4X,"|")')
               DO i = 1, NUMAT
                  WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')i,COORD(1,i),COORD(2,i),COORD(3,i)
               END DO
               WRITE(*,'("+-----------------------------------------------+")')
            END IF

!  Calculate energy and gradient for all states
            CALL SCFMUL(ARRAY,LM5,ICALL,SCFCAL,ENERGX,CGX,CNORMX,NUMAT)
            ICICAL = ICICAL+1

            ihop = .false.

            IF (ICALL .NE. -1) THEN
               CALL estat_fol(iout,fol_stat,IN2(159),NCICONF,CICOMP,ENERGX,old_CICOMP,old_ENERGX,mdstep,startStep,naco,de,old_de,sw_CC,tully_rest,inde,old_inde)
               istatt = inde(istat)

!   Update ovelocity to half step forward
               IF (ehrenfest) THEN
                  temp_grad = ehrgr(:,:)
               ELSE
                  temp_grad = CGX(:,:,istatt)
               END IF
               DO i = 1, NUMAT
                  delta = sca*0.5D0*dt_temp/AMS(i)
                  vx(i) = vx(i) - temp_grad(1,i)*delta
                  vy(i) = vy(i) - temp_grad(2,i)*delta
                  vz(i) = vz(i) - temp_grad(3,i)*delta   
               END DO
               IF (thermostat .EQ. 1) THEN
                  call noseHoover_advance(iout,dt_temp,vx,vy,vz,AMS,NUMAT,thermo_eta,thermo_p,thermo_tau,thermo_len,thermo_T,nf,E_therm,.true.)
               END IF

               IF (iout .GE. 3) THEN
                  WRITE(*,*)
                  WRITE(*,'("+-----------------------------------------------+")')
                  WRITE(*,'("|",3X,"MD step: ",I4," - Second update of velocity",3X,"|")')mdstep
                  WRITE(*,'("|",47X,"|")')
                  WRITE(*,'("| Atom",8X,"vx",12X,"vy",12X,"vz",4X,"|")')
                  DO i = 1, NUMAT
                     WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')i,vx(i),vy(i),vz(i)
                  END DO
                  WRITE(*,'("+-----------------------------------------------+")')
               END IF


!     Hopping procedures
               hop_vx = vx
               hop_vy = vy
               hop_vz = vz
               IF (tully_hop) THEN
                  CALL SURHOP(iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,dt_temp,vx,vy,vz,AMS,NUMAT,istat,ENERGX,CGX,ihop,write_hop,hopunit,Hold,rdold,C,NE,old_CICOMP,old_ENERGX,old_de,num_CC,an_CC,Einteg,tully_rest,istatt,rnd_gen,rnd_count,fol_stat,dec_cor,nac_ph,cuthop,inde,sw_CC)
               ELSE IF (ehrenfest) THEN
                  CALL ehrendyn(iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,dt_temp,vx,vy,vz,AMS,NUMAT,istat,ENERGX,CGX,write_hop,hopunit,Hold,rdold,C,NE,old_CICOMP,old_ENERGX(:,3),old_de,num_CC,an_CC,Einteg,tully_rest,ehrgr,fol_stat,istatt,inde,sw_CC)
               ELSE IF (simhop .GT. 0) THEN
                  CALL SIMPHOP(iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,dt_temp,NUMAT,istat,ENERGX,CGX,old_CICOMP,istatt,rnd_gen,AMS,vx,vy,vz,ihop,simhop)
               END IF
               IF (.not.ihop) THEN
                  vx = hop_vx
                  vy = hop_vy
                  vz = hop_vz
               END IF

! store data for next step
               DO i = 1, IN2(159)
                  old_CICOMP(:,i) = CICOMP(:,inde(i))
                  old_ENERGX(i,1) = old_ENERGX(i,2)
                  old_ENERGX(i,2) = old_ENERGX(i,3)
                  old_ENERGX(i,3) = ENERGX(inde(i))
               END DO
               old_de = de

! LANGEVIN
               IF (g0.gt.0) THEN
                  IF (ehrenfest) THEN
                     call langevin(AMS, vx, vy, vz, g0, temp0, ehrgr, NUMAT, dt, mdstep, jst)
                  ELSE
                     call langevin(AMS, vx, vy, vz, g0, temp0, CGX(1,1,istatt), NUMAT, dt, mdstep, jst)
                  ENDIF
               ENDIF
 
                  IF (iout .GE. 3) THEN
!     output gradient
                  IF (ehrenfest) THEN
                     temp_grad = ehrgr(:,:)
                  ELSE
                     temp_grad = CGX(:,:,istatt)
                  END IF
                  WRITE(*,*)
                  WRITE(*,'("+-----------------------------------------------+")')
                  WRITE(*,'("|",11X,"MD step: ",I4," - Gradient",12X,"|")')mdstep
                  WRITE(*,'("|",47X,"|")')
                  WRITE(*,'("| Atom",8X,"Gx",12X,"Gy",12X,"Gz",4X,"|")')
                  DO i = 1, NUMAT
                  WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')i,temp_grad(1,i),temp_grad(2,i),temp_grad(3,i)
                  END DO
                  WRITE(*,'("+-----------------------------------------------+")')

               END IF
            END IF

            IF (adaptive) THEN
!  Check for adaptive procedure
               step_fail = .FALSE.
               IF (ICALL .EQ. -1) THEN
                  IF (IFMAP .EQ. 1) THEN
                     WRITE(*,*) "ACTIVE ORBITAL MAPPING FAILURE"
                     IF (adapt_map .EQ. 0) THEN
                        STOP 'Abnormal termination in sub. DYNAM'
                     END IF
                  ELSE 
                     WRITE(*,*)"CONVERGENCY PROBLEMS!!!"
                  END IF
!       Adaptive failure because of Mapping or convergence
                  step_fail = .TRUE.
               ELSE IF (adapt_ene .GT. 0.D0 .AND. ihop .NE. .TRUE.) THEN
                  CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
                  E_tot = ENERGX(istatt) + K + E_therm
                  IF (iout .GE. 1) THEN
                     WRITE(*,*)
                     WRITE(*,*) "Adaptive energy conservation"
                     IF (adapt_ene_mode .EQ. 0) THEN
                        WRITE(*,'(11X,3(A16,3X))') "now","old","total change %"
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_tot",E_tot,E_old,100.D0*(E_tot-E_old)/E_old
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_kin",K,E_K_old,100.D0*(K-E_K_old)/E_K_old
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_pot",ENERGX(istatt),E_P_old,100.D0*(ENERGX(istatt)-E_P_old)/E_P_old
                        IF (thermostat>0) WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_therm",E_therm,E_therm_old,100.D0*(E_therm-E_therm_old)/E_therm_old
                        WRITE(*,'(3X,A8,F9.4,"%")') "thresh",adapt_ene*100.D0
                        IF (DABS((E_tot-E_old)/E_old) .GT.adapt_ene) THEN
      !       Adaptive failure because of total energy change
                           step_fail = .TRUE.
                        END IF
                     ELSE
                        WRITE(*,'(11X,3(A16,3X))') "now","old","total change"
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_tot",E_tot,E_old,(E_tot-E_old)
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_kin",K,E_K_old,(K-E_K_old)
                        WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_pot",ENERGX(istatt),E_P_old,ENERGX(istatt)-E_P_old
                        IF (thermostat>0) WRITE(*,'(3X,A8,2(F16.6,3X),F16.4)') "E_therm",E_therm,E_therm_old,E_therm-E_therm_old
                        WRITE(*,'(3X,A8,F9.4," kcal/mol")') "thresh",adapt_ene
                        IF (DABS(E_tot-E_old).GT.adapt_ene) THEN
      !       Adaptive failure because of total energy change
                           step_fail = .TRUE.
                        END IF
                     END IF
                     WRITE(*,*)
                  END IF
               END IF

!       If adaptive fails, increase number of tries
               IF (step_fail) THEN
                  adapt_try = adapt_try +1
                  IF (adapt_try.GT.adapt_tries) THEN
                     WRITE(*,'(A,I2,A)')"No correct calculation after ",adapt_tries," tries."
                     IF (ICALL .EQ. -1 .OR. .NOT. adapt_proceed) THEN
                        WRITE(*,*) "The MD simulation will end here."
                        STOP 'Adaptive failure in sub. DYNAM'
!       Break whole simulation if convergence errors after adapt_tries
                     ELSE
                        WRITE(*,*) "Proceed with maybe wrong values."
!       If no convergence errors, proceed
                     END IF
                  ELSE
!       If tries < max tries, retry wih restored properties
                     WRITE(*,'("Retry #",I2," with smaller timestep")')adapt_try

                     vx = vx_temp
                     vy = vy_temp
                     vz = vz_temp
                     COORD = COORD_temp
                     E_old = E_old_temp
                     E_K_old = E_K_old_temp
                     E_P_old = E_P_old_temp
                     E_therm_old = E_therm_old_temp
                     CGX = CGX_temp
                     Hold = Hold_temp
                     rdold = rdold_temp
                     old_CICOMP = old_CICOMP_temp
                     old_ENERGX = old_ENERGX_temp
                     old_de = old_de_temp
                     old_inde = old_inde_temp
                     C = C_temp
                     istat = istat_temp
                     istatt = istatt_temp
                     t_goal = t_goal_temp
                     thermo_eta = thermo_eta_temp
                     thermo_p = thermo_p_temp
                     wait_steps = NINT((t_tot+dt_temp-t_tot_temp)/(dt_temp/2.0D0))
                     t_tot = t_tot_temp
!       Restore the MO mapping
                     if (IMOMAP .EQ. 3) then
                        INQUIRE( FILE='imomap.tmp1', EXIST=THERE ) 
                        CALL SYSTEM("cp -f imomap.tmp1 imomap.dat")
                        IF ( THERE ) NEWMAP = .TRUE.
                     end if
                     IF (ICALL .EQ. -1) THEN
                        INQUIRE( FILE='fort.11', EXIST=THERE ) 
                        IF ( THERE ) CALL SYSTEM("rm fort.11")
                     ELSE
                        INQUIRE( FILE='fort.11.tmp1', EXIST=THERE ) 
                        IF ( THERE ) CALL SYSTEM("cp -f fort.11.tmp1 fort.11")
                     END IF
                     INQUIRE( FILE='hop.tmp1', EXIST=THERE ) 
                     IF ( THERE ) CALL SYSTEM("cp -f hop.tmp1 hop.tmp")
!       Jump to the start of adaptive procedure
                     CYCLE
                  END IF
               END IF
               step_corr = step_corr + 1
               NEWMAP = .FALSE.
!           Save all the stats for future restoration

               IF (wait_steps .GT. 0) wait_steps = wait_steps - 1
               COORD_temp = COORD
               CGX_temp = CGX
               Hold_temp = Hold
               rdold_temp = rdold
               old_CICOMP_temp = old_CICOMP
               old_ENERGX_temp = old_ENERGX
               old_de_temp = old_de
               old_inde_temp = old_inde
               C_temp = C
               istat_temp = istat
               istatt_temp = istatt
               t_goal_temp = t_goal
               thermo_eta_temp = thermo_eta
               thermo_p_temp = thermo_p
!           Copy MO mapping to temporary file
               if (IMOMAP .EQ. 3) then
                  INQUIRE( FILE='imomap.dat', EXIST=THERE ) 
                  IF ( THERE ) CALL SYSTEM("cp -f imomap.dat imomap.tmp1")
               end if
               INQUIRE( FILE='fort.11', EXIST=THERE ) 
               IF ( THERE ) CALL SYSTEM("cp -f fort.11 fort.11.tmp1")
               INQUIRE( FILE='hop.tmp', EXIST=THERE ) 
               IF ( THERE ) CALL SYSTEM("cp -f hop.tmp hop.tmp1")
            ELSE
!       If adaptive is turned off, break simulation in case of errors
               IF (ICALL .EQ. -1) THEN
                  IF (IMOMAP .GE. 2 .AND. IFMAP .EQ. 1) THEN
                     WRITE(*,*)
                     WRITE(*,*) "ACTIVE ORBITAL MAPPING FAILURE"
                  ELSE 
                     WRITE(*,*)
                     WRITE(*,*)"CONVERGENCY PROBLEMS!!!"
                  END IF
                  WRITE(*,*)"The MD simulation will end here."
                  WRITE(*,*)
                  STOP 'Abnormal termination in sub. DYNAM'
               END IF
            END IF

!     Check for energy conservation
            IF (ene_cons) THEN
!       Calculate total energy
               CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
               E_tot = ENERGX(istatt) + K + E_therm
               IF (iout .GE. 1) THEN
                  WRITE(*,*)
                  WRITE(*,'("Velocity scaling thermostat")')
                  WRITE(*,'("   E_ini     ",F16.8)') E_ini
                  WRITE(*,*)
                  WRITE(*,'("   E_tot     ",F16.8)')E_tot
                  WRITE(*,'("   E_pot     ",F16.8)')ENERGX(istatt)
                  WRITE(*,'("   E_kin     ",F16.8)')K
               END IF

!       Calculate conversion factor for slow energy adaption to initial value
               E_tot = E_ini+(E_tot-E_ini)*(ene_scal**(dt_temp/1.D-3))
               conv = DSQRT((E_tot-ENERGX(istatt)-E_therm)/K)
               vx = conv * vx
               vy = conv * vy
               vz = conv * vz
!       Calculate new kinetic and total energy
               CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
               E_tot = ENERGX(istatt) + K + E_therm

               IF (iout .GE. 1) THEN
               WRITE(*,*)
                  WRITE(*,'("   Adapt total energy to initial value")')
                  WRITE(*,'("   Scale velocities with ",F12.8)') conv
                  WRITE(*,*)
                  WRITE(*,'("   E_tot     ",F16.8)') E_tot
                  WRITE(*,'("   E_pot     ",F16.8)') ENERGX(istatt)
                  WRITE(*,'("   E_kin     ",F16.8)') K
               END IF
            END IF

!       Increment the actual time counter
            t_tot = t_tot + dt_temp
!       Double the time step, but no larger than the next goal
            IF (DABS(t_tot - t_goal) .LT. 1.0D-12) THEN
               t_tot = t_goal
            END IF

            IF (adaptive) THEN
               vx_temp = vx
               vy_temp = vy
               vz_temp = vz
               CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
               E_tot = ENERGX(istatt) + K + E_therm
               E_old_temp = E_tot
               E_K_old_temp = K
               E_P_old_temp = ENERGX(istatt)
               E_therm_old_temp = E_therm
               E_old = E_tot
               E_K_old = K
               E_P_old = ENERGX(istatt)
               E_therm_old = E_therm
               IF (DABS(t_tot - t_goal) .LT. 1.0D-12) THEN
                  t_goal_temp = t_tot+dt
               END IF
               t_tot_temp = t_tot
!       If no adaptive errors, proceed with the next step
               WRITE(*,'(A)')"Correct calculation in adaptive scheme."
            END IF
!        Adaptive-do
         END DO


         CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)
         E_tot = ENERGX(istatt) + K + E_therm
    
!     Velocity scaling
         IF (vs) THEN
            IF (mod(mdstep,fvs) .EQ. 0) THEN
               IF (.NOT. write_rest) THEN
                  WRITE(*,*)
                  WRITE(*,*)"WARNING!!"
                  WRITE(*,*)"A velocity scaling run is attempted"
                  WRITE(*,*)"with no final saving for restart"
                  WRITE(*,*)"Velocity scaling is highly recomended"
                  WRITE(*,*)"to be used only for equilibration."
                  WRITE(*,*)"Thus, writing a restart file would be"
                  WRITE(*,*)"a good idea"
                  WRITE(*,*)
               END IF
               CALL vel_scal(iout,NUMAT,AMS,vx,vy,vz,Tfix,vs_Emax,nf)
               IF (notra) THEN
                  CALL v_remove_tra(iout,NUMAT,AMS,vx,vy,vz)
               END IF
               IF (norot) THEN
                  CALL v_remove_rot(iout,NUMAT,AMS,COORD(1,:),COORD(2,:),COORD(3,:),vx,vy,vz)
               END IF
            END IF
         END IF

         IF (write_hop) THEN
            INQUIRE( FILE='hop.tmp', EXIST=THERE ) 
            IF ( THERE ) THEN
               open(77,STATUS='OLD',FILE='hop.tmp',IOSTAT=io_disk)
               if (io_disk.eq.0) THEN
                  do
                     read(77,'(A)',IOSTAT=io_disk) str
                     if (io_disk.ne.0) exit
                     write(hopunit,'(A)') TRIM(str)
                  end do
                  close(77,STATUS='DELETE')
               END IF
            END IF
         END IF

         IF (write_stats) THEN
            IF (avstat .OR. MOD(mdstep,fstat).EQ.0) THEN
               CALL outpr2(iout,sunit,mdstep,NUMAT,nf,AMS,vx,vy,vz,COORD,IN2(159),istatt,ENERGX(istatt),avstat,fstat,E_av,K_av,V_av,T_av,P2_av,Px_av,Py_av,Pz_av,L2_av,Lx_av,Ly_av,Lz_av,C,ehrenfest,E_therm)
            END IF
         END IF

!  Write restart file if requested
         IF (write_rest) THEN
            CALL write_res(iout,rfile,NUMAT,COORD(1,:),COORD(2,:),COORD(3,:),vx,vy,vz,tully_hop,istat,IN2(159),C,Hold,rdold,NCICONF,old_CICOMP,naco,old_de,ENERGX,rnd_count)
         END IF


!  End of MD loop
      END DO
         
!....................................................................C


!  CLOSE OUTPUT FILES IF NECESSARY
      IF (write_stats) THEN
         CLOSE (UNIT=sunit)
      END IF

      IF (write_traj) THEN
         CLOSE(UNIT=xunit)
      END IF

      IF (write_vel) THEN
         CLOSE(UNIT=vunit)
      END IF

      IF (write_hop) THEN
         CLOSE(unit=hopunit)
      END IF

      IF (iout .GE. 2) THEN
! TIMING - overall time, end
         CALL CPU_TIME(ctime2)
         CALL SYSTEM_CLOCK(wtime2,check_ir)

! TIMING - output
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING - Final statistics ***"
         WRITE(*,*)
         WRITE(*,*)"Overall time for MD simulation"
         CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
         WRITE(*,*)
         WRITE(*,*)"Time for velocity and coordinates managment"
         CALL out_time(0.D0,ctime_v,0,wtime_v,check_ir)
         WRITE(*,*)
         WRITE(*,*)"Time for gradient computation"
         CALL out_time(0.D0,ctime_g,0,wtime_g,check_ir)
         WRITE(*,*)
         WRITE(*,*)"Time for state tracking"
         CALL out_time(0.D0,ctime_sf,0,wtime_sf,check_ir)
         IF (an_CC) THEN
         WRITE(*,*)
         WRITE(*,*)"Time for analytical nonadiabatic coupling vectors"
         CALL out_time(0.D0,ctime_acc,0,wtime_acc,check_ir)
         END IF
         IF (num_CC) THEN
         WRITE(*,*)
         WRITE(*,*)"Time for numerical nonadiabatic coupling"
         CALL out_time(0.D0,ctime_ncc,0,wtime_ncc,check_ir)
         END IF
         WRITE(*,*)
         WRITE(*,*)"Time for electronic effective Hamiltonian"
         CALL out_time(0.D0,ctime_eh,0,wtime_eh,check_ir)
         WRITE(*,*)
         WRITE(*,*)"Time for electronic integration"
         CALL out_time(0.D0,ctime_eint,0,wtime_eint,check_ir)
         WRITE(*,'("-----------------------------------")')
      END IF

      END SUBROUTINE



!=============================================================================C

      SUBROUTINE nmlst(iout,nstep,dt,temp0,norot,notra,init_stat,write_stats,sunit,sfile,write_traj,xunit,xfile,write_vel,vunit,vfile,adapt_ene,fstat,avstat,fsav,restart,write_rest,rfile,vs,Tfix,fvs,vs_Emax,ehrenfest,tully_hop,Einteg,NE,num_CC,an_CC,rnd_gen,rnd_seed,write_hop,hopunit,hopfile,fhop,simhop,fol_stat,dec_cor,nac_ph,cuthop,ene_cons,ene_scal,adapt_map,adapt_tries,adapt_proceed,adapt_ene_mode,thermostat,thermo_len,thermo_tau,thermo_T,thermo_equi,g0)

!     This subroutine  manages  name list options

!  VARIABLES  (in namelist order)
!
! nstep         Max number of MD steps 
! dt            Time step [in picoseconds]
! temp0         Initial temperature   
! norot         Flag to select initial velocities so that angular 
!                 momentum is zero
! notra         Flag to select initial velocities so that linear 
!                 momentum  is zero
! init_stat     Initial state 
! write_stats   Flag to allow writing of statistics 
! sunit         Fortran unit for statistics file 
! sfile         Statistics file name 
! write_traj    Flag to allow writing trajectory file
! xunit         Fortran unit for trajectory file 
! xfile         Trajectory file name 
! write_vel     Flag to allow writing of velocities 
! vunit         Fortran unit for velocity file
! vfile         Velocity file name
! adapt_ene       Threshold for relative variation of total energy
! fstat         frequency for saving in statistics file                    
!                 = 1  every step                                          
!                 = n  every n steps
! avstat        Flag to turn on statistics averaging                       
!                If fstat = n > 1 the average value over n steps will be   
!                saved   
! fsav          frequency for saving in trajectory and velocity file       
!                 = 1  every step                                          
!                 = n  every n steps
! restart       Flag to allow restarting form file 
! write rest    Flag to allow writing the restart file   
! rfile         Restart file name
! ehrenfest     Flag to enable Ehrenfest dynamics
! tully_hop     Flag to enable Tully's hopping procedure
! simhop        Flag to enable simple hopping
! Einteg        Integration algorithm to use for electronic equation in
!                Tully's hopping:
!                    euler  -> Euler method
!                    RK2    -> Runge-Kutta second-order algorithm
!                    RK4    -> Runge-Kutta fourth-order algorithm
!                    ABM4   -> Adams-Bashforth-Moulton fourth order        
!                              predictor-corrector with error estimation   
!                    ABM2   -> Adams-Bashforth-Moulton second order        
!                              predictor-corrector                         
!                    ABM5   -> Adams-Bashforth-Moulton fifth order         
!                              predictor-corrector with error estimation   
!                    imp    -> implicit middle-point algorithm             
!                    UP1    -> Unitary propagator evaluated at middle-point
!                              built from sylvester theorem
!                    UP2    -> Unitary propagator evaluated at the         
!                              beginning of time interval
!                    UP3    -> Unitary propagator evaluated at middle-point  
!                              bulit from taylor expansion                 
! NE            Number of steps for electronic integration (in Tully's hop)
! num_CC        Flag to enable calculation of numerical non-adiabatic couplings
! an_CC         Flag to enable calcualtion of analytical non-adiabatic couplings
! rnd_gen       Algorithm for the generation of random numbers in tully's hopping:
!                    standard  -> standard FORTRAN algorithm [RANDOM_NUMBER]
!                    PM_BD    -> Park and Miller algorithm with Bays-Durham shuffle
!                    knuth     -> Knuth subtractive algorithm        
! rnd_seed      Integer to seed random number generator.
!               If -1 (default), a seed is generated randomly, e.g. with the system clock
! fol_stat      Flag to enable tracking of electronic states
! write_hop     Flag to allow writing hopping data to file
! hopunit       Fortran unit for hopping data file
! hopfile       Hopping data file name  
! fhop          frequency for saving hopping data   
!                 = 1  every step              
!                 = n  every n steps
! vs            Flag to enable velocity scaling 
! Tfix          Temperature to be reached with velocity scaling
! fvs           frequency for performing velocity scaling                
!                 = 1  every step                                          
!                 = n  every n steps
! vs_Emax       Maximum kinetic energy correction allowed in velocity
! dec_cor       Decoherence correction by Persico & Granucci (standard value=0.1)
! nac_ph        threshold for nonadiabatic coupling phase following        C
!               algorithm (0,0=off, standard=0.1)                          C
! thermostat    Choice of thermostat
! thermo_len    Length of Nose-Hoover chain
! thermo_equi   Number of equilibration steps
! thermo_tau    Characteristic time tau of Nose-Hoover thermostat
! thermo_T      Goal temperature for thermostat
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                         

      IMPLICIT NONE

!  NAMELIST VARIABLES
      INTEGER :: iout
      INTEGER :: nstep
      REAL*8 :: dt
      REAL*8 :: temp0
      LOGICAL :: norot
      LOGICAL :: notra
      INTEGER :: init_stat
      LOGICAL :: write_stats
      INTEGER :: sunit
      CHARACTER(20) :: sfile
      LOGICAL :: write_traj
      INTEGER :: xunit
      CHARACTER(20) :: xfile
      LOGICAL :: write_vel
      INTEGER :: vunit
      CHARACTER(20) :: vfile
      LOGICAL :: ene_cons
      REAL*8 :: ene_tol
      REAL*8 :: ene_scal
      REAL*8 :: adapt_ene
      INTEGER :: adapt_ene_mode
      INTEGER :: adapt_map
      INTEGER :: adapt_tries
      LOGICAL :: adapt_proceed
      INTEGER :: fstat
      LOGICAL :: avstat
      INTEGER :: fsav
      LOGICAL :: restart
      LOGICAL :: write_rest
      CHARACTER(20) :: rfile
      LOGICAL :: ehrenfest
      LOGICAL :: tully_hop
      CHARACTER(5) :: Einteg
      INTEGER :: NE
      LOGICAL :: num_CC
      LOGICAL :: an_CC
      CHARACTER(8) :: rnd_gen
      INTEGER :: rnd_seed
      LOGICAL :: write_hop
      INTEGER :: hopunit
      CHARACTER(20) :: hopfile
      INTEGER :: fhop
      LOGICAL :: vs
      REAL*8 :: Tfix
      INTEGER :: fvs
      REAL*8 :: vs_Emax
      INTEGER :: simhop
      LOGICAL :: fol_stat
      REAL*8 :: dec_cor
      REAL*8 :: nac_ph 
      REAL*8 :: cuthop     
      INTEGER :: thermostat
      INTEGER :: thermo_len
      INTEGER :: thermo_equi
      REAL*8 :: thermo_tau
      REAL*8 :: thermo_T
      REAL*8 :: g0

      NAMELIST / dynvar / iout, nstep, dt, temp0, norot, notra, init_stat, write_stats, sunit, sfile, write_traj, xunit, xfile, write_vel, vunit, vfile, fstat, avstat, ene_tol, fsav, restart, write_rest, rfile, ehrenfest, tully_hop, simhop, Einteg, NE, num_CC, an_CC, rnd_gen, rnd_seed, fol_stat, write_hop, hopunit, hopfile, fhop, vs, Tfix, fvs, vs_Emax, dec_cor, nac_ph, cuthop, adapt_ene, ene_scal, ene_cons, adapt_map, adapt_tries, adapt_proceed, adapt_ene_mode, thermostat, thermo_len, thermo_tau, thermo_T, thermo_equi, g0

!  LOCAL VARIABLES
      INTEGER :: file_exist
      INTEGER :: io_disk
      INTEGER :: tmp_a,tmp_b,tmp_c,tmp_d

!  COMMON
      INTEGER :: IN2

      COMMON /INOPT2/ IN2(300)

! ---------------- C

!  Initialization of variables
      iout = 0
      ehrenfest = .FALSE.
      tully_hop = .FALSE.
      simhop = 0
      nstep = 10
      dt = 5.d-5
      temp0 = 300.d0
      norot = .TRUE.
      notra = .TRUE.
      init_stat = 1
      write_stats = .TRUE.
      sunit = 30
      sfile = "stat.out"
      write_traj = .TRUE.
      xunit = 31
      xfile = "traj.out"
      write_vel = .FALSE.
      vunit = 32
      vfile = "vel.out"
      ene_cons = .FALSE.
      ene_scal = 1.0D0
      adapt_ene = 0.0D0
      adapt_ene_mode = 0
      adapt_map = 0
      adapt_tries = 7
      adapt_proceed = .FALSE.
      fstat = 1
      avstat = .FALSE.
      fsav = 1
      restart = .FALSE.
      write_rest = .TRUE.
      rfile = "dynam.restart"
      Einteg = "UP3"
      NE = 100
      num_CC = .FALSE.
      an_CC = .TRUE.
      rnd_gen = "standard"
      rnd_seed = -1
      fol_stat = .FALSE.
      write_hop = .FALSE.
      hopunit = 33
      hopfile = "hopping.out"
      fhop = 1
      vs = .FALSE.
      Tfix = temp0
      fvs = 1
      vs_Emax = -1.d0    ! negative values mean no limit
      dec_cor = 0.1d0
      nac_ph = 0.1d0
      cuthop= 0.0d0
      ene_tol= 0.0d0
      thermostat = 0
      thermo_len = 1
      thermo_equi = 0
      thermo_tau = 1.D0
      thermo_T = 300.D0
      g0= 0.0d0

!  Check if a namelist file is already there. If don't write one with
!  default values. If it already exists read it
      INQUIRE(FILE="dynvar.in", EXIST=file_exist)

      IF(.NOT. file_exist) THEN
         OPEN(UNIT=31,NAME="dynvar.in",STATUS="NEW",ACTION="WRITE",IOSTAT=io_disk)
         IF (io_disk .EQ. 0) THEN
            WRITE(31,NML=dynvar)
         ELSE
            WRITE(*,*)
            WRITE(*,*)"Unknown disk error"
            WRITE(*,*)
            STOP 'Abnormal temination in sub. DYNAM'
         END IF
      ELSE
         OPEN(UNIT=31,NAME="dynvar.in",STATUS="UNKNOWN",ACTION="READ",IOSTAT=io_disk)
         IF (io_disk .EQ. 0) THEN
            READ(31,NML=dynvar)
         ELSE
            WRITE(*,*)
            WRITE(*,*)"Unknown disk error"
            WRITE(*,*)
            STOP 'Abnormal temination in sub. DYNAM'
         END IF
      END IF

      CLOSE(UNIT=31)

!  Check for constistency of option
      IF (ene_tol .NE. 0.0D0) THEN
         WRITE(*,*)
         WRITE(*,*)"Warning: Deprecated option ene_tol used. Will be ignored."
         WRITE(*,*)
      END IF

      IF (thermo_len .LT. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"thermo_len must be greater equal than 1"
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF
      IF (thermo_tau .LT. 0.0D0) THEN
         WRITE(*,*)
         WRITE(*,*)"thermo_tau must be greater equal than 0.0"
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF
      IF (thermo_equi .LT. 0) THEN
         WRITE(*,*)
         WRITE(*,*)"thermo_equi must be greater equal than 0"
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF
      IF (thermo_T .LT. 0.0D0) THEN
         WRITE(*,*)
         WRITE(*,*)"thermo_T must be greater equal than 0.0"
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF
      IF (simhop .LT. 0) THEN
         WRITE(*,*)
         WRITE(*,*)"SIMHOP must be greater equal than zero"
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF

      tmp_a = 0
      tmp_b = 0
      tmp_c = 0
      
      IF (tully_hop) tmp_a = 1
      IF (ehrenfest) tmp_b = 1
      IF (simhop .GT. 0) tmp_c = 1

      tmp_d = tmp_a + tmp_b + tmp_c
      
      IF (tmp_d .GT. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"Tully hopping,Ehrenfest and simple hopping"
         WRITE(*,*)"are incompatible. Please set only one to .TRUE."
         WRITE(*,*)
         STOP 'Abnormal termination in sub NMLST'
      END IF

      IF (tmp_d.EQ.1 .AND. .NOT.write_hop) THEN
         WRITE(*,*)
         WRITE(*,*)"Hopping selected!"
         WRITE(*,*)"Hopping file writing enabled"
         WRITE(*,*)
         write_hop = .TRUE.
      END IF
      IF (tmp_d.EQ.0 .AND. write_hop) THEN
         WRITE(*,*)
         WRITE(*,*)"No hopping selected!"
         WRITE(*,*)"Hopping file writing disabled"
         WRITE(*,*)
         write_hop = .FALSE.
      END IF

      IF (tully_hop) THEN
         IF (.NOT. num_CC .AND. .NOT. an_CC) THEN
            WRITE(*,*)
            WRITE(*,*)"No method selected for the calculation of"
            WRITE(*,*)"non-adiabatic coupling vectors"
            WRITE(*,*)"Please make a choice"
            WRITE(*,*)
            STOP 'Abnormal termination in sub DYNAM'
         END IF
      END IF

!  Output
      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"       Namelist variables"
         WRITE(*,*)
         WRITE(*,'("+---------------+---------------------+")')
         WRITE(*,'("| iout ",9X,"|",19X,I1," |")')iout
         WRITE(*,'("| nstep ",8X,"|",15X,I5," |")')nstep
         WRITE(*,'("| dt ",11X,"|",11X,F9.6," |")')dt
         WRITE(*,'("| temp0 ",8X,"|",9X,F11.5," |")')temp0
         WRITE(*,'("| norot ",8X,"|",19X,L1," |")')norot
         WRITE(*,'("| notra ",8X,"|",19X,L1," |")')notra
         WRITE(*,'("| init_stat ",4X,"|",19X,I1," |")')init_stat
         WRITE(*,'("| write_stats ",2X,"|",19X,L1," |")')write_stats
         WRITE(*,'("|   sunit ",6X,"|",17X,I3," |")')sunit
         WRITE(*,'("|   sfile ",6X,"|",A20," |")')sfile
         WRITE(*,'("| write_traj ",3X,"|",19X,L1," |")')write_traj
         WRITE(*,'("|   xunit ",6X,"|",17X,I3," |")')xunit
         WRITE(*,'("|   xfile ",6X,"|",A20," |")')xfile
         WRITE(*,'("| write_vel ",4X,"|",19X,L1," |")')write_vel
         WRITE(*,'("|   vunit ",6X,"|",17X,I3," |")')vunit
         WRITE(*,'("|   vfile ",6X,"|",A20," |")')vfile
         WRITE(*,'("| ene_cons ",5X,"|",19X,L1," |")')ene_cons
         WRITE(*,'("| ene_scal ",5X,"|",14X,F6.4," |")')ene_scal
         WRITE(*,'("| adapt_ene ",4X,"|",14X,F6.4," |")')adapt_ene
         WRITE(*,'("| adapt_ene_mode|",17X,I3," |")')adapt_ene_mode
         WRITE(*,'("| adapt_map ",4X,"|",17X,I3," |")')adapt_map
         WRITE(*,'("| adapt_tries ",2X,"|",17X,I3," |")')adapt_tries
         WRITE(*,'("| adapt_proceed |",19X,L1," |")')adapt_proceed
         WRITE(*,'("| fstat ",8X,"|",17X,I3," |")')fstat
         WRITE(*,'("|   avstat ",5X,"|",19X,L1," |")')avstat
         WRITE(*,'("| fsav ",9X,"|",17X,I3," |")')fsav
         WRITE(*,'("| restart ",6X,"|",19X,L1," |")')restart
         WRITE(*,'("| write_rest ",3X,"|",19X,L1," |")')write_rest
         WRITE(*,'("|   rfile ",6X,"|",A20," |")')rfile
         WRITE(*,'("| ehrenfest ",4X,"|",19X,L1," |")')ehrenfest
         WRITE(*,'("| tully_hop ",4X,"|",19X,L1," |")')tully_hop
         WRITE(*,'("| dec_cor ",6X,"|",14X,F6.4," |")')dec_cor
         WRITE(*,'("| nac_ph ",7X,"|",14X,F6.4," |")')nac_ph
         WRITE(*,'("| cuthop ",7X,"|",14X,F6.4," |")')cuthop
         WRITE(*,'("| simhop ",7X,"|",19X,I1," |")')simhop
         WRITE(*,'("|   Einteg ",5X,"|",15X,A5," |")')Einteg
         WRITE(*,'("|   NE ",9X,"|",15X,I5," |")')NE
         WRITE(*,'("|   num_CC ",5X,"|",19X,L1," |")')num_CC
         WRITE(*,'("|   an_CC ",6X,"|",19X,L1," |")')an_CC
         WRITE(*,'("|   rnd_gen ",4X,"|",12X,A8," |")')rnd_gen
         WRITE(*,'("|   rnd_seed ",3X,"|",9X,I11," |")') rnd_seed
         WRITE(*,'("|   fol_stat ",3X,"|",19X,L1," |")') fol_stat
         WRITE(*,'("| write_hop ",4X,"|",19X,L1," |")')write_hop
         WRITE(*,'("|   hopunit ",4X,"|",17X,I3," |")')hopunit
         WRITE(*,'("|   hopfile ",4X,"|",A20," |")')hopfile
         WRITE(*,'("|   fhop ",7X,"|",16X,I4," |")')fhop
         WRITE(*,'("| vs ",11X,"|",19X,L1," |")')vs
         WRITE(*,'("|   Tfix ",7X,"|",9X,F11.5," |")')Tfix
         WRITE(*,'("|   fvs ",8X,"|",16X,I4," |")')fvs
         WRITE(*,'("|   vs_Emax ",4X,"|",8X,F12.6," |")')vs_Emax
         WRITE(*,'("|   g0      ",4X,"|",8X,F12.6," |")')g0 
         WRITE(*,'("| Thermostat ",3X,"|",16X,I4," |")')thermostat
         WRITE(*,'("|   Thermo_T ",3X,"|",8X,F12.6," |")')thermo_T
         WRITE(*,'("|   Thermo_len ",1X,"|",16X,I4," |")')thermo_len
         WRITE(*,'("|   Thermo_tau ",1X,"|",8X,F12.6," |")')thermo_tau
         WRITE(*,'("|   Thermo_equi",1X,"|",12X,I8," |")')thermo_equi
         WRITE(*,'("+---------------+---------------------+")')
      END IF

      END SUBROUTINE





      SUBROUTINE noseHoover_advance(iout,dt,vx,vy,vz,mass,N,eta,p,tau,M,T,nf,E_therm,output)
! This subroutine advances the Nose-Hoover thermostat variables.
! Formulas are taken from J. Phys. A, 39 (2006) p. 5629

! IN/OUT variables
      INTEGER :: iout
      REAL*8 :: dt
      INTEGER :: N
      REAL*8, DIMENSION(N) :: vx,vy,vz
      REAL*8, DIMENSION(N) :: mass
      INTEGER :: M
      REAL*8, DIMENSION(M) :: eta, p
      REAL*8 :: tau
      REAL*8 :: T
      REAL*8 :: E_therm
      logical :: output

! Local variables
      REAL*8 :: K
      REAL*8, DIMENSION(M) :: Q, G
      INTEGER :: i

      REAL*8 :: k2kcal = 0.0019872041  ! k[kcal/mol*K]
 
! Calculation of thermostat masses
      Q = k2kcal*T*tau**2
      Q(1) = nf*Q(1)

! Calculation of forces acting on the thermostat
      CALL Kene(iout,N,mass,vx,vy,vz,K)
      G(1) = K*2.D0 - nf*k2kcal*T
      DO i = 2, M
         G(i) = p(i-1)**2/Q(i-1) - k2kcal*T
      END DO

! First update of friction values
      p(M) = p(M) + G(M) * dt/4.D0
      DO i = M-1, 1, -1
         p(i) = p(i) * exp(-p(i+1)/Q(i+1) * dt/8.D0)
         p(i) = p(i) + G(i) * dt/4.D0
         p(i) = p(i) * exp(-p(i+1)/Q(i+1) * dt/8.D0)
      END DO

! Update of atomic velocities
      DO i = 1, N
         vx(i) = vx(i) * exp(-p(1)/Q(1) * dt/2.D0)
         vy(i) = vy(i) * exp(-p(1)/Q(1) * dt/2.D0)
         vz(i) = vz(i) * exp(-p(1)/Q(1) * dt/2.D0)
      END DO

! Update of thermostat values (only needed for total energy contributions)
      DO i = 1, M
         eta(i) = eta(i) - p(i)/Q(i) * dt/2.D0
      END DO

! Recalculation of forces acting on the thermostat
      CALL Kene(iout,N,mass,vx,vy,vz,K)
      G(1) = K*2.D0 - nf*k2kcal*T
      DO i = 2, M
         G(i) = p(i-1)**2/Q(i-1) - k2kcal*T
      END DO

! Second update of friction values
      DO i = 1, M-1
         p(i) = p(i) * exp(-p(i+1)/Q(i+1) * dt/8.D0)
         p(i) = p(i) + G(i) * dt/4.D0
         p(i) = p(i) * exp(-p(i+1)/Q(i+1) * dt/8.D0)
      END DO
      p(M) = p(M) + G(M) * dt/4.D0

! Contribution of thermostat to total energy
      E_therm = 0.D0
      DO i = 1, M
         E_therm = E_therm + p(i)**2/(2.D0*Q(i))
      END DO
      E_therm = E_therm - nf*k2kcal*T*eta(1)
      DO i = 2, M
         E_therm = E_therm - k2kcal*T*eta(i)
      END DO

! Printing section
      IF (output) THEN
         IF (iout .GE. 1) THEN
            WRITE(*,*)
            WRITE(*,'(A)') "Nose-Hoover chains thermostat:"
            WRITE(*,'(A,F16.6,A,F16.6)') "goal temperature =",T,", actual temperature = ",2.d0*K/(k2kcal*nf)
            WRITE(*,'(A,F16.6,A,F16.6)') "thermostat contribution to total energy =",E_therm
            WRITE(*,'(A,F16.6)') "velocity scaling factor = ",exp(-p(1)/Q(1) * dt/2.D0)
         END IF
         IF (iout .GE. 3) THEN
            WRITE(*,'(A,F16.6)') "new velocities: "
            DO i = 1, N
               WRITE(*,'(A,I,A,3F16.6)') "v(",i,") = ",vx(i),vy(i),vz(i)
            END DO
         END IF
         IF (iout .GE. 3) THEN
            DO i = 1, M
               WRITE(*,'(A,i2,A,F16.6)') "G",i," = ",G(i)
            END DO
            DO i = 1, M
               WRITE(*,'(A,i2,A,E16.6)') "Q",i," = ",Q(i)
            END DO
            DO i = 1, M
               WRITE(*,'(A,i2,A,F16.6)') "eta",i," = ",eta(i)
            END DO
            DO i = 1, M
               WRITE(*,'(A,i2,A,F16.6)') "p",i," = ",p(i)
            END DO
         END IF
      END IF

      END SUBROUTINE




!=============================================================================C


      SUBROUTINE v_rnd(iout,NUMAT,AMS,rnd_gen,rnd_count,vx,vy,vz)

!  This subroutine creates random velocities
!     Starting velocities are assignet as follows:
!     (1) a random number 0=<q<1 is chosen
!     (2) v = (q/(1-q))/mass
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      CHARACTER(8) :: rnd_gen

!  INPUT/OUTPUT
      INTEGER :: rnd_count

!  OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: vx, vy, vz

!  LOCAL VARIABLES
      INTEGER :: i, j
      REAL*8, DIMENSION(3) :: q

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>v_rnd"
         WRITE(*,*)
         WRITE(*,*)"Random velocities assignement"
         WRITE(*,*)
         WRITE(*,*)"Random numbers:"
      END IF

!TWK   Seeded in DYNAM
!      CALL RANDOM_SEED

      DO i = 1, NUMAT
!TWK     Choose generator consistently  
         IF (rnd_gen .EQ. "standard") THEN
            CALL RANDOM_NUMBER(q)
         ELSE IF (rnd_gen .EQ. "PM_BD") THEN
            DO j = 1, 3
               CALL nran2(q(j),1,.FALSE.)
            END DO
         ELSE IF (rnd_gen .EQ. "knuth") THEN
            DO j = 1, 3
               CALL nran3(q(j),1,.FALSE.)
            END DO
         END IF

         rnd_count = rnd_count + 3

         IF (iout .GE. 3) THEN
            WRITE(*,'(16X,F8.6,2X,F8.6,2X,F8.6)')q(1),q(2),q(3)
         ENDIF

         vx(i) = (q(1)/(1-q(1)))/AMS(i)
         vy(i) = (q(2)/(1-q(2)))/AMS(i)
         vz(i) = (q(3)/(1-q(3)))/AMS(i)
      END DO

!   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)        
         WRITE(*,*)"Velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"v_rnd<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE v_remove_tra(iout,NUMAT,AMS,vx,vy,vz)

!     This routine removes any overall translational motion
!     from the molecular system

!     (1) total momentum p=sum(m_i*v_i) is computed
!     (2) M = sum(m_i) 
!     (3) v = v - (p/M)
!         In this way we have zero total linear momentum
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS

!  INPUT/OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: vx, vy, vz

!  LOCAL VARIABLES
      INTEGER :: i
      REAL*8 :: px,py,pz
      REAL*8 :: M
      REAL*8 :: facx,facy,facz

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>v_remove_tra"
      END IF

! (1)(2)
      px = 0.D0
      py = 0.D0
      pz = 0.D0

      DO i = 1, NUMAT
         px = px + AMS(i)*vx(i)
         py = py + AMS(i)*vy(i)
         pz = pz + AMS(i)*vz(i)
         M = M + AMS(i)
      END DO

      facx = px/M
      facy = py/M
      facz = pz/M

! (3)
      DO i = 1, NUMAT
         vx(i) = vx(i) - facx
         vy(i) = vy(i) - facy
         vz(i) = vz(i) - facz
      END DO

!   output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Linear momentum [AMU*A/ps]"
         WRITE(*,'(5X,"px",12X,"py",12X,"pz")')
         WRITE(*,'(F12.5,2X,F12.5,2X,F12.5)')px,py,pz
         WRITE(*,*)
         WRITE(*,'("Total mass = ",F10.5)')M
         WRITE(*,*)
         WRITE(*,*)"Translation corrections"
         WRITE(*,'("facx = ", F11.5)')facx
         WRITE(*,'("facy = ", F11.5)')facy
         WRITE(*,'("facz = ", F11.5)')facz
         WRITE(*,*)
         WRITE(*,*)"New velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"v_remove_tra<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE v_remove_rot(iout,NUMAT,AMS,x,y,z,vx,vy,vz)

!     This subroutine removes any overall rotation from the 
!     molecular system.

!     (1) Calculate total angular momentum vector L
!     (2) Calculate moment of inertia tensor I
!     (3) Find the angular velocity vector w by
!         by solving the linear equation L=Iw
!     (4) Apply the reverse angular velocity to each atom in turn,
!         cancelling out the overall rotation of the system
!--------------------------------------------------------------------------C
! AUTHOR: T. Keal
! MAIL: keal@mpi-muelheim.mpg.de                                       
! DATE: November 2006  

      IMPLICIT NONE

!     INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: x,y,z

!     INPUT/OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: vx, vy, vz

!     LOCAL VARIABLES
      INTEGER :: i, la_return
      REAL*8, DIMENSION(3) :: L, dum1
      REAL*8, DIMENSION(3,3) :: MoI
      REAL*8, DIMENSION(100) :: dum2
      REAL*8, DIMENSION(NUMAT) :: xx, yy, zz

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>v_remove_rot"
      END IF

! (1) Calculate total angular momentum vector L(1:3)
      CALL ang_mom(iout,NUMAT,AMS,x,y,z,vx,vy,vz,L(1),L(2),L(3))

! (2) Calculate moment of inertia tensor I (MoI)
!     with respect to centre of mass
      CALL pos_wrt_com(iout,NUMAT,AMS,x,y,z,xx,yy,zz)

      MoI = 0.D0

      DO i = 1, NUMAT
         MoI(1,1) = MoI(1,1) + (AMS(i) * (yy(i)*yy(i) + zz(i)*zz(i)))
         MoI(2,2) = MoI(2,2) + (AMS(i) * (xx(i)*xx(i) + zz(i)*zz(i)))
         MoI(3,3) = MoI(3,3) + (AMS(i) * (xx(i)*xx(i) + yy(i)*yy(i)))
         MoI(1,2) = MoI(1,2) - (AMS(i) * xx(i) * yy(i))
         MoI(1,3) = MoI(1,3) - (AMS(i) * xx(i) * zz(i))
         MoI(2,3) = MoI(2,3) - (AMS(i) * yy(i) * zz(i))
      END DO

!     The moment of inertia tensor is symmetrical
      MoI(2,1) = MoI(1,2)
      MoI(3,1) = MoI(1,3)
      MoI(3,2) = MoI(2,3)

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Moment of inertia"
         WRITE(*,'(9X,"1",11X,"2",11X,"3")')
         WRITE(*,'("1",4X,F10.4,2X,F10.4,2X,F10.4)')MoI(1,1),MoI(1,2),MoI(1,3)
         WRITE(*,'("2",16X,F10.4,2X,F10.4)')MoI(2,2),MoI(2,3)
         WRITE(*,'("3",28X,F10.4)')MoI(3,3)
      END IF

! (3) Calculate the angular velocity vector w
!     The LAPACK routine overwrites the input array L with this quantity
!     so from this point on L is actually the angular velocity vector
!     MoI is also overwritten
!     See http://www.netlib.org/lapack/double/dsysv.f for more details

      CALL DSYSV('U', 3, 1, MoI, 3, dum1, L, 3, dum2, 100, la_return)

      IF (la_return .NE. 0) THEN
         STOP 'Angular velocities could not be found in v_remove_rot!'
      ENDIF

! (4) Apply the angular velocities in reverse
!     thus cancelling out the overall angular momentum

!     The velocities are corrected according to
!     v -> v - w x r
!     where r are atomic positions relative to the centre of mass

      DO i = 1, NUMAT
         vx(i) = vx(i) - (L(2)*zz(i) - L(3)*yy(i))
         vy(i) = vy(i) - (L(3)*xx(i) - L(1)*zz(i))
         vz(i) = vz(i) - (L(1)*yy(i) - L(2)*xx(i))
      END DO

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Angular velocity vector"
         WRITE(*,'(4X,"wx",10X,"wy",10X,"wz")')
         WRITE(*,'(F10.4,2X,F10.4,2X,F10.4)')L(1),L(2),L(3)
         WRITE(*,*)
         WRITE(*,*)"New velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"v_remove_rot<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE v_remove_rot



!=============================================================================C


      SUBROUTINE term_scal(iout,NUMAT,AMS,vx,vy,vz,temp0,N)

!  This subroutine performs termal scaling
!     Velocities are scaled as follows:
!     (1) the "instantaneous" velocity T=(1/2)*sum(m_i*v_i^2)/(k*N) is computed
!         k is bolzman konstant and N is the number of degrees of freedom
!     (2) scaling factor a = SQRT(T_0/T) with T_0 desired temperature (temp0)
!     (3) v = a*v; in this way we have T=T_0
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8 :: temp0
      INTEGER :: N

!  INPUT/OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: vx, vy, vz

!  LOCAL VARIABLES
      REAL*8 :: K
      REAL*8 :: T_ist
      REAL*8 :: conv
      INTEGER :: i

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>term_scal"
      END IF

!  (1)
      CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)

      CALL Tist(iout,K,N,T_ist)

!  (2)
      conv = DSQRT(temp0/T_ist)

!  (3)
      vx = conv*vx
      vy = conv*vy
      vz = conv*vz

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Conversion factor = ",F10.4)')conv
         WRITE(*,*)
         WRITE(*,*)"New velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"term_scal<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE lin_mom(iout,NUMAT,AMS,vx,vy,vz,Px,Py,Pz)

!  This subroutine computes total linear momentum
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz

!  OUTPUT VARIABLES
      REAL*8 :: Px,Py,Pz

! LOCAL VARIABLES
      INTEGER :: i

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>lin_mom"
      END IF

      Px = 0.D0
      Py = 0.D0
      Pz = 0.D0

      DO i = 1, NUMAT
         Px = Px + AMS(i)*vx(i)
         Py = Py + AMS(i)*vy(i)
         Pz = Pz + AMS(i)*vz(i)
      END DO

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Linear momentum [AMU*A/ps]"
         WRITE(*,'(5X,"px",12X,"py",12X,"pz")')
         WRITE(*,'(F12.5,2X,F12.5,2X,F12.5)')px,py,pz
         WRITE(*,*)
         WRITE(*,*)"lin_mom<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE ang_mom(iout,NUMAT,AMS,x,y,z,vx,vy,vz,Lx,Ly,Lz)

!  This subroutine computes total angular momentum with respect
!  to the center of mass
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: x,y,z
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz

!  OUTPUT VARIABLES
      REAL*8 :: Lx,Ly,Lz

!  LOCAL VARIABLES
      INTEGER :: i
      REAL*8, DIMENSION(NUMAT) :: xx,yy,zz

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>ang_mom"
      END IF

!     Get atomic positions with respect to centre of mass
      CALL pos_wrt_com(iout,NUMAT,AMS,x,y,z,xx,yy,zz)

!     Compute angular momentum
      Lx = 0.D0
      Ly = 0.D0
      Lz = 0.D0

      DO i = 1, NUMAT
         Lx = Lx + AMS(i)*(yy(i)*vz(i) - zz(i)*vy(i))
         Ly = Ly + AMS(i)*(zz(i)*vx(i) - xx(i)*vz(i))
         Lz = Lz + AMS(i)*(xx(i)*vy(i) - yy(i)*vx(i))
      END DO

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Angolar momentum [AMU*A^2/ps]"
         WRITE(*,'(5X,"Lx",12X,"Ly",12X,"Lz")')
         WRITE(*,'(F12.5,2X,F12.5,2X,F12.5)')Lx,Ly,Lz
         WRITE(*,*)
         WRITE(*,*)"ang_mom<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE pos_wrt_com(iout,NUMAT,AMS,x,y,z,xx,yy,zz)

!     Calculates atom positions with respect to centre of mass
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!     INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: x,y,z

!     OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: xx,yy,zz

!     LOCAL VARIABLES
      INTEGER :: i
      REAL*8 :: xcm,ycm,zcm
      REAL*8 :: M

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>pos_wrt_com"
      END IF

!     Find center of mass
      xcm = 0.D0
      ycm = 0.D0
      zcm = 0.D0
      M = 0.D0

      DO i = 1, NUMAT
         xcm = xcm + x(i)*AMS(i)
         ycm = ycm + y(i)*AMS(i)
         zcm = zcm + z(i)*AMS(i)
         M = M + AMS(i)
      END DO

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M

!     Translate origin of coordinates in the center of mass
      DO i = 1, NUMAT
         xx(i) = x(i) - xcm
         yy(i) = y(i) - ycm
         zz(i) = z(i) - zcm
      END DO

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Coordinates of center of mass [A]"
         WRITE(*,'("X = ",F10.5,3X,"Y = ",F10.5,3X,"Z = ",F10.5)')xcm,ycm,zcm
         WRITE(*,*)
         WRITE(*,*)"Coordinates in center of mass reference system [A]"
         WRITE(*,'("Atom",9X,"x",13X,"y",13X,"z")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,xx(i),yy(i),zz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"pos_wrt_com<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE pos_wrt_com



!=============================================================================C


      SUBROUTINE Kene(iout,NUMAT,AMS,vx,vy,vz,K)
!  This subroutine computes total kinetic energy
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007  
                                                      
      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz

!  OUTPUT VARIABLES
      REAL*8 :: K

!  LOCAL VARIABLES
      REAL*8 :: v2
      INTEGER :: i

!  PARAMETERS
      REAL*8 :: conv = 2.390057D-3  ! AMU*A^2/ps^2 -> kcal/mol

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>Kene"
      END IF

      v2 = 0.D0
      K = 0.D0

! bug correction by spoerkel
      DO i = 1, NUMAT
         v2 = vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
         K = K + AMS(i)*v2
      END DO

      K = conv*0.5D0*K

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("v^2 = ",F10.5," A^2/ps^2")')v2
         WRITE(*,*)
         WRITE(*,'("Kinetic energy = ",F10.5," AMU*A^2/ps^2")')K
         WRITE(*,*)
         WRITE(*,*)"Kene<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE Tist(iout,K,N,T)

!  This subroutine computes instantaneous temperature
!      2*K
!  T = ---
!      B*N
!
!   K = kinetik energy
!   B = boltzmann constant
!   N = number of degrees of freedom
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      REAL*8 :: K
      INTEGER :: N

!  OUTPUT VARIABLES
      REAL*8 :: T

!  PARAMETERS
      REAL*8, PARAMETER :: boltz=1.987168981D-3 ! Boltzmann constant 
                                                ! [kcal/(mol*K)]

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>Tist"
      END IF

      T = 2.d0*K/(boltz*REAL(N))

!  output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Instantaneous temperature = ",F10.5," K")')T
         WRITE(*,*)
         WRITE(*,*)"Tist<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE outpr(NUMAT,nf,AMS,vx,vy,vz,E_tot,K,V)

!  This subroutine prints output initial statistics
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      USE LIMIT, ONLY: LM1

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: NUMAT
      INTEGER :: nf
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz
      REAL*8 :: E_tot
      REAL*8 :: K
      REAL*8 :: V
    
      REAL*8 :: COORD

      COMMON /ATOMC / COORD(3,LM1)

!  LOCAL VARIABLES
      REAL*8 :: T_ist
      REAL*8 :: Px,Py,Pz
      REAL*8 :: Lx,Ly,Lz
      REAL*8, DIMENSION(NUMAT) :: x,y,z
      INTEGER :: i

!-----------------------C

      do i = 1, NUMAT
         x(i) = COORD(1,i)
         y(i) = COORD(2,i)
         z(i) = COORD(3,i)
      end do

      WRITE(*,*)
      WRITE(*,'("-------------------------------------------------------------------")')
      WRITE(*,*)
      WRITE(*,*)"                       INITIAL STATISTICS"

      CALL Tist(0,K,nf,T_ist)

      WRITE(*,*)
      WRITE(*,'("total energy",3X,"kinet energy",3X,"poten energy",3X,"temperature")')
      WRITE(*,'(F12.6,3X,F12.6,3X,F12.6,3X,F12.6)')E_tot,K,V,T_ist
      WRITE(*,'("...........................................................")')

      CALL lin_mom(0,NUMAT,AMS,vx,vy,vz,Px,Py,Pz)

      WRITE(*,'(12X,"Px",11X,"Py",11X,"Pz")')
      WRITE(*,'(8X,F10.5,3X,F10.5,3X,F10.5)')Px,Py,Pz

      CALL ang_mom(0,NUMAT,AMS,x,y,z,vx,vy,vz,Lx,Ly,Lz)
      
      WRITE(*,'(12X,"Lx",11X,"Ly",11X,"Lz")')
      WRITE(*,'(8X,F10.5,3X,F10.5,3X,F10.5)')Lx,Ly,Lz
      WRITE(*,'("-------------------------------------------------------------------")')
     
      END SUBROUTINE



!=============================================================================C


      SUBROUTINE outpr2(iout,sunit,mdstep,NUMAT,nf,AMS,vx,vy,vz,COORD,nExStat,istat,EX,avstat,fstat,E_av,K_av,V_av,T_av,P2_av,Px_av,Py_av,Pz_av,L2_av,Lx_av,Ly_av,Lz_av,C,ehrenfest,E_therm)

!  This subroutine processes averages and prints statistics
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      USE LIMIT, ONLY: LM1

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: sunit
      INTEGER :: mdstep
      INTEGER :: NUMAT
      INTEGER :: nf
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz
      INTEGER :: nExStat
      INTEGER :: istat
      REAL*8 :: EX, E_therm
      LOGICAL :: avstat
      INTEGER :: fstat
      COMPLEX*16, DIMENSION(nExStat) :: C
      LOGICAL :: ehrenfest
 
!  INPUT/OUTPUT VARIABLES
      REAL*8 :: E_av
      REAL*8 :: K_av
      REAL*8 :: V_av
      REAL*8 :: T_av
      REAL*8 :: P2_av
      REAL*8 :: Px_av
      REAL*8 :: Py_av
      REAL*8 :: Pz_av
      REAL*8 :: L2_av
      REAL*8 :: Lx_av
      REAL*8 :: Ly_av
      REAL*8 :: Lz_av
      REAL*8 :: E_therm_av

!  COMMON
      REAL*8, DIMENSION(3,LM1) :: COORD

!  LOCAL VARIABLES
      REAL*8 :: T_ist
      REAL*8 :: Px,Py,Pz
      REAL*8 :: Lx,Ly,Lz
      REAL*8, DIMENSION(NUMAT) :: x,y,z
      REAL*8 :: P2
      REAL*8 :: L2
      INTEGER :: i      
      REAL*8 :: E_tot
      REAL*8 :: K
      REAL*8 :: V

!-----------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>outpr2"
      END IF

      do i = 1, NUMAT
         x(i) = COORD(1,i)
         y(i) = COORD(2,i)
         z(i) = COORD(3,i)
      end do

      CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)

      CALL Tist(iout,K,nf,T_ist)

      E_tot = 0.D0
      V = 0.D0
      IF (ehrenfest) THEN
         DO i = 1, nExStat
            V = V + REAL(C(i)*CONJG(C(i)))*EX
         END DO
      ELSE
         V = EX
      END IF
      E_tot = V + K + E_therm

      CALL lin_mom(iout,NUMAT,AMS,vx,vy,vz,Px,Py,Pz)

      CALL ang_mom(iout,NUMAT,AMS,x,y,z,vx,vy,vz,Lx,Ly,Lz)

      P2 = (Px*Px)+(Py*Py)+(Pz*Pz)
      L2 = (Lx*Lx)+(Ly*Ly)+(Lz*Lz)

      IF (avstat) THEN
         E_av = E_av + E_tot
         K_av = K_av + K
         V_av = V_av + V
         T_av = T_av + T_ist
         P2_av = P2_av + P2
         Px_av = Px_av + Px
         Py_av = Py_av + Py
         Pz_av = Pz_av + Pz
         L2_av = L2_av + L2
         Lx_av = Lx_av + Lx
         Ly_av = Ly_av + Ly
         Lz_av = Lz_av + Lz
         E_therm_av = E_therm_av + E_therm
      END IF

      IF (avstat .AND. MOD(mdstep,fstat).EQ. 0) THEN
         E_av = E_av/REAL(fstat)
         K_av = K_av/REAL(fstat)
         V_av = V_av/REAL(fstat)
         T_av = T_av/REAL(fstat)
         P2_av = P2_av/REAL(fstat)
         Px_av = Px_av/REAL(fstat)
         Py_av = Py_av/REAL(fstat)
         Pz_av = Pz_av/REAL(fstat)
         L2_av = L2_av/REAL(fstat)
         Lx_av = Lx_av/REAL(fstat)
         Ly_av = Ly_av/REAL(fstat)
         Lz_av = Lz_av/REAL(fstat)
         E_therm_av = E_therm_av/REAL(fstat)
      END IF

      IF (MOD(mdstep,fstat) .EQ. 0) THEN

         IF (avstat) THEN 
            WRITE(sunit,'(I5,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7)')mdstep,E_av,K_av,V_av,T_av,P2_av,Px_av,Py_av,Pz_av,L2_av,Lx_av,Ly_av,Lz_av,E_therm_av

            E_av = 0.D0
            K_av = 0.D0
            V_av = 0.D0
            T_av = 0.D0
            P2_av = 0.D0
            Px_av = 0.D0
            Py_av = 0.D0
            Pz_av = 0.D0
            L2_av = 0.D0
            Lx_av = 0.D0
            Ly_av = 0.D0
            Lz_av = 0.D0
         ELSE
            WRITE(sunit,'(I5,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7,1X,F14.7)')mdstep,E_tot,K,V,T_ist,P2,Px,Py,Pz,L2,Lx,Ly,Lz,E_therm
         END IF
         
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"outpr2<<<"
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE read_rest(iout,rfile,N,x,y,z,vx,vy,vz,tully_hop,nExStat,init_stat,C,Hold,rdold,energx,CI_comp,old_de,tully_rest,rnd_count)

!  This subroutine reads the restart file
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      USE LIMIT, ONLY: LMCONF
      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      CHARACTER(20) :: rfile
      INTEGER :: N
      LOGICAL :: tully_hop
      INTEGER :: nExStat
      
!  OUTPUT VARIABLES
      REAL*8, DIMENSION(N) :: vx,vy,vz
      REAL*8, DIMENSION(N) :: x,y,z

      INTEGER :: init_stat
      COMPLEX*16, DIMENSION(nExStat) :: C
      COMPLEX*16, DIMENSION(nExStat,nExStat) :: Hold
      REAL*8, DIMENSION(nExStat,nExStat) :: rdold
      REAL*8, DIMENSION(nExStat) :: energx
      REAL*8, DIMENSION(LMCONF,6) :: CI_comp
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      LOGICAL :: tully_rest
      INTEGER :: rnd_count

!  LOCAL VARIABLES
      INTEGER :: i,j,k
      INTEGER :: disk_error
      LOGICAL :: file_exist
      CHARACTER(13) :: dummyc
      CHARACTER(7) :: restype
      INTEGER :: natom
      INTEGER :: nExStat2
      REAL*8 :: c1,c2
      REAL*8 :: dummyr
      INTEGER :: nciconfig
      INTEGER :: n_a_ci_o

! ---------------- C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>read_rest"
      END IF

!  Check if restart file exists
      INQUIRE(FILE=rfile, EXIST=file_exist)
      
      IF (.NOT. file_exist) THEN
         WRITE(*,*)
         WRITE(*,*)"Restart file not found"
         WRITE(*,*)
         STOP 'Abnormal termination in sub. DYNAM'
      END IF

      OPEN(UNIT=29,NAME=rfile,STATUS="OLD", ACTION="READ",IOSTAT=disk_error)
      IF (disk_error .NE. 0) THEN
         WRITE(*,*)
         WRITE(*,*)"Unkown disk error"
         WRITE(*,*)
         STOP 'Abnormal termination in sub. read_rest'
      END IF

! read header
      READ(29,*)dummyc,restype

      IF (restype .EQ. "Tully") THEN
         IF (.NOT. tully_hop) THEN
            WRITE(*,*)
            WRITE(*,'("Restart file is inconsistent!")')
            WRITE(*,'("Tully keyword found in a non-Tully run!")')
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF
         tully_rest = .TRUE.
      ELSE IF (restype .EQ. "MD_only") THEN
         tully_rest = .FALSE.
      ELSE
         WRITE(*,*)
         WRITE(*,'("Restart file is inconsistent!")')
         WRITE(*,'(A, "is an unknown command")')restype
         WRITE(*,*)
         STOP 'Abnormal termination in read_rest'
      END IF

! read number of atoms
      READ(29,*)dummyc,natom

      IF (dummyc .NE. "#N_atoms") THEN
         WRITE(*,*)
         WRITE(*,*)"Restart file is inconsistent!"
         WRITE(*,'(A, "is wrong")')dummyc
         WRITE(*,*)
         STOP 'Abnormal termination in read_rest'
      END IF

      IF (natom .NE. N) THEN
         WRITE(*,*)
         WRITE(*,*)"number of atoms inconsistent in restart file"
         WRITE(*,*)
         STOP 'Abnormal termination in sub. read_rest'
      END IF

! read coordinates
      READ(29,*)dummyc

      IF (dummyc .NE. "#coordinates") THEN
         WRITE(*,*)
         WRITE(*,*)"Restart file is inconsistent!"
         WRITE(*,'(A, "is wrong")')dummyc
         WRITE(*,*)
         STOP 'Abnormal termination in read_rest'
      END IF

      DO i = 1, N
         READ(29,*)x(i),y(i),z(i)
      END DO

! read velocity
      READ(29,*)dummyc

      IF (dummyc .NE. "#velocity") THEN
         WRITE(*,*)
         WRITE(*,*)"Restart file is inconsistent!"
         WRITE(*,'(A, "is wrong")')dummyc
         WRITE(*,*)
         STOP 'Abnormal termination in read_rest'
      END IF

      DO i = 1, N
         READ(29,*)vx(i),vy(i),vz(i)
      END DO

! read the followings only if tully_rest
      IF (tully_rest) THEN

! read init_stat
         READ(29,*)dummyc,init_stat
         IF (dummyc .NE. "#init_stat") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF
 
! read N_ex_stat
         READ(29,*)dummyc,nExStat2
         IF (dummyc .NE. "#N_ex_stat") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         IF (nExStat2 .NE. nExStat) THEN
            WRITE(*,*)
            WRITE(*,*)"Number of excited states inconsistent"
            WRITE(*,*)
            STOP 'Abnormal termination in sub. read_rest'
         END IF

! read population
         READ(29,*)dummyc
         
         IF (dummyc .NE. "#population") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         DO i = 1, nExStat
            READ(29,*)c1,c2
            C(i) = DCMPLX(c1,c2)
         END DO

! read Hold
         READ(29,*)dummyc
         
         IF (dummyc .NE. "#old_H") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         DO i = 1, nExStat
            DO j = i, nExStat
               READ(29,*)dummyr
               IF (i .EQ. j) Hold(i,i) = DCMPLX(0.d0,dummyr)
               IF (i .NE. j) THEN
                  Hold(i,j) = DCMPLX(dummyr,0.D0)
                  Hold(j,i) = Hold(i,j)
               END IF
            END DO
         END DO

! read rdold
         READ(29,*)dummyc
         
         IF (dummyc .NE. "#old_rd") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         rdold = 0.D0
         DO i = 1, nExStat
            DO j = i, nExStat 
               READ(29,*)rdold(i,j)
               rdold(j,i) = rdold(i,j)
            END DO
         END DO

! read energx
         READ(29,*)dummyc

         IF (dummyc .NE. "#energx") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         DO i = 1, nExStat
            READ(29,*)energx(i)
         END DO

! read NCICONF
         READ(29,*)dummyc,nciconfig
         IF (dummyc .NE. "#N_ci_conf") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

! read old_CICOMP
         READ(29,*)dummyc
         
         IF (dummyc .NE. "#old_CI_comp") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         DO i = 1, nExStat
            DO j = 1, nciconfig
               READ(29,*)CI_comp(j,i)
            END DO
         END DO

! read naco (number of active CI orbitals
         READ(29,*)dummyc,n_a_ci_o
         IF (dummyc .NE. "#N_act_CI_orb") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

! read old_de
         READ(29,*)dummyc
         
         IF (dummyc .NE. "#old_de") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

         DO i = 1, nciconfig
            DO j = 1, n_a_ci_o
               READ(29,*)old_de(j,i)
            END DO
         END DO

!TWK read rnd_count
         READ(29,*)dummyc,rnd_count
         IF (dummyc .NE. "#rnd_count") THEN
            WRITE(*,*)
            WRITE(*,*)"Restart file is inconsistent!"
            WRITE(*,'(A, "is wrong")')dummyc
            WRITE(*,*)
            STOP 'Abnormal termination in read_rest'
         END IF

      END IF

      CLOSE(UNIT=29)

!   output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         IF (tully_rest) THEN
            WRITE(*,*)"Restart type: Tully"
         ELSE
            WRITE(*,*)"Restart type: MD"
         END IF
         WRITE(*,*)
         WRITE(*,'(" Number of restart atoms: ",I5)')natom

         WRITE(*,*)
         WRITE(*,*)"Restart coordinates [A]"
         WRITE(*,'("Atom",9X,"x",13X,"y",13X,"z")')
         DO i = 1, N
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,x(i),y(i),z(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"Restart velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, N
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO

         IF (tully_rest) THEN
            WRITE(*,*)
            WRITE(*,'(" Restart initial state: ",I2)')init_stat
            WRITE(*,*)
            WRITE(*,'(" Restart number of excited states: ",I2)')nExStat
            WRITE(*,*)
            WRITE(*,*)"Restart electronic population"
            DO i = 1, nExStat
               IF (AIMAG(C(i)) .GE. 0.D0) THEN
                  WRITE(*,'("C(",I1") = ",F10.6," + i ",F10.6)')i,REAL(C(i)),AIMAG(C(i))
               ELSE
                  WRITE(*,'("C(",I1") = ",F10.6," - i ",F10.6)')i,REAL(C(i)),ABS(AIMAG(C(i)))
               END IF
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart effective electronic Hamiltonian"
            DO i = 1, nExStat
               WRITE(*,'(I1,5X)',ADVANCE='NO')i
               DO j = 1, i-1
                  WRITE(*,'(16X,3X)',ADVANCE='NO')
               END DO
               WRITE(*,'(F16.6,1X,"i",1X)',ADVANCE='NO')AIMAG(Hold(i,i))
               IF (i .LT. nExStat) THEN
                  DO j = i+1, nExStat-1
                     WRITE(*,'(F16.6,3X)',ADVANCE='NO')REAL(Hold(i,j))
                  END DO
                  WRITE(*,'(F16.6)')REAL(Hold(i,nExStat))
               ELSE
                  WRITE(*,*)
               END IF
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart non-adiabatic coupling terms"
            WRITE(*,'(1X,10I15)')(j,j=2,nExStat)
            DO i = 1, nExStat-1
               WRITE(*,'(I1,4X,10F15.6)')i,(rdold(i,j),j=i+1,nExStat)
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart excitation energies"
            DO i = 1, nExStat
               WRITE(*,'(F12.6)')energx(i)
            END DO
            WRITE(*,*)
            WRITE(*,'(" Restart number of CI components: ", I5)')nciconfig
            WRITE(*,*)
            WRITE(*,'(" Restart number of CI active orbitals: ", I3)')n_a_ci_o 
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)"Restart CI components and lexical index"
            DO i = 1, nExStat
               WRITE(*,'(" -- state ",I2," --")')i
               DO j = 1, nciconfig
                  WRITE(*,'(I5,2X,F12.6,20I3)')j,CI_comp(j,i),(old_de(k,j),k=1,n_a_ci_o)
               END DO
            END DO
            WRITE(*,*)
         END IF

         WRITE(*,*)
         WRITE(*,*)"read_rest<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE write_res(iout,rfile,N,x,y,z,vx,vy,vz,tully_hop,init_stat,nExStat,C,Hold,rdold,nciconfig,CI_comp,n_a_ci_o,old_de,energx,rnd_count)

!  This writes the restart file
!
!  Format is (options in square brakets):
!
!  #restart  [MD_only; Tully]
!  #N_atoms .....
!  #coordinates
!  ..........  .......... .........
!  #velocity
!  ..........  .......... .........
!  #init_stat  ..
!  #N_ex_stat  ..
!  #population
!  ..........  ..........
!  #old_H
!  ..........
!  #old_rd
!  ..........
!  #energx
!  ..........
!  #N_ci_conf  .....
!  #old_CI_comp
!  ..........
!  #N_act_CI_orb
!  #old_de
!  ..........
!  #rnd_count
!
!TWK Increased format to F20.10 to hopefully ensure better reproducibility
!
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        
      USE LIMIT, ONLY: LMCONF

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout
      CHARACTER(20) :: rfile
      INTEGER :: N
      REAL*8, DIMENSION(N) :: vx,vy,vz
      REAL*8, DIMENSION(N) :: x,y,z
      LOGICAL :: tully_hop
      INTEGER :: init_stat
      INTEGER :: nExStat
      COMPLEX*16, DIMENSION(nExStat) :: C
      COMPLEX*16, DIMENSION(nExStat,nExStat) :: Hold
      REAL*8, DIMENSION(nExStat,nExStat) :: rdold
      INTEGER :: nciconfig
      REAL*8, DIMENSION(LMCONF,6) :: CI_comp
      INTEGER :: n_a_ci_o
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      REAL*8, DIMENSION(nExStat) :: energx
      INTEGER :: rnd_count

!  LOCAL VARIABLES
      INTEGER :: i,j,k
      INTEGER :: disk_error

!-------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>write_res"
      END IF

      OPEN(UNIT=29,NAME=rfile,STATUS="REPLACE", ACTION="WRITE",IOSTAT=disk_error)
      IF (disk_error .NE. 0) THEN
         WRITE(*,*)
         WRITE(*,*)"Unkown disk error"
         WRITE(*,*)
         STOP 'Abnormal termination in sub. DYNAM'
      END IF

! write header. MD_only indicates that only coordinates and velocities are 
! written. Tully indicates that all data are witten
      IF (.NOT. tully_hop) THEN
         WRITE(29,'("#restart   MD_only")')
      ELSE
         WRITE(29,'("#restart   Tully")')
      END IF

! write number of atoms
      WRITE(29,'("#N_atoms ",I5)')N

! write coordinates
      WRITE(29,'("#coordinates")')
      DO i = 1, N
         WRITE(29,'(3F20.10)')x(i),y(i),z(i)
      END DO

! write velocities
      WRITE(29,'("#velocity")')
      DO i = 1, N
         WRITE(29,'(3F20.10)')vx(i),vy(i),vz(i)
      END DO

! next informations are written only if Tully hopping
      IF (tully_hop) THEN

! write state
         WRITE(29,'("#init_stat ",I2)')init_stat

! write number of excited states
         WRITE(29,'("#N_ex_stat ",I2)')nExStat

! write population
         WRITE(29,'("#population")')
         DO i = 1, nExStat
            WRITE(29,'(2F20.10)')REAL(C(i)),AIMAG(C(i))
         END DO

! write Hold
         WRITE(29,'("#old_H")')
         DO i = 1, nExStat
            DO j = i, nExStat
               IF (i .EQ. j) WRITE(29,'(F20.10)')AIMAG(Hold(i,j))
               IF (i .NE. j) WRITE(29,'(F20.10)')REAL(Hold(i,j))
            END DO
         END DO

! write rdold
         WRITE(29,'("#old_rd")')
         DO i = 1, nExStat
            DO j = i, nExStat
               WRITE(29,'(F20.10)')rdold(i,j)
            END DO
         END DO

! write energx
         WRITE(29,'("#energx")')
         DO i = 1, nExStat
            WRITE(29,'(F20.10)')ENERGX(i)
         END DO

! write NCICONF
         WRITE(29,'("#N_ci_conf ",I5)')nciconfig

! write old_CICOMP
         WRITE(29,'("#old_CI_comp")')
         DO i = 1, nExStat
            DO j = 1, nciconfig
               WRITE(29,'(F20.10)')CI_comp(j,i)
            END DO
         END DO

! write naco (number of active CI orbitals
         WRITE(29,'("#N_act_CI_orb ",I3)')n_a_ci_o

! write old_de
         WRITE(29,'("#old_de")')
         DO i = 1, nciconfig
            DO j = 1, n_a_ci_o
               WRITE(29,'(I2)')old_de(j,i)
            END DO
         END DO

! write rnd_count
         WRITE(29,'("#rnd_count ",I3)') rnd_count

      END IF

      CLOSE(UNIT=29)

!   output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         IF (tully_hop) THEN
            WRITE(*,*)"Restart type: Tully"
         ELSE
            WRITE(*,*)"Restart type: MD"
         END IF
         WRITE(*,*)
         WRITE(*,'(" Number of restart atoms: ",I5)')N

         WRITE(*,*)
         WRITE(*,*)"Restart coordinates [A]"
         WRITE(*,'("Atom",9X,"x",13X,"y",13X,"z")')
         DO i = 1, N
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,x(i),y(i),z(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"Restart velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, N
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO


         IF (tully_hop) THEN

            WRITE(*,*)
            WRITE(*,'(" Restart initial state: ",I2)')init_stat
            WRITE(*,*)
            WRITE(*,'(" Restart number of excited states: ",I2)')nExStat
            WRITE(*,*)
            WRITE(*,*)"Restart electronic population"
            DO i = 1, nExStat
               IF (AIMAG(C(i)) .GE. 0.D0) THEN
                  WRITE(*,'("C(",I1") = ",F10.6," + i ",F10.6)')i,REAL(C(i)),AIMAG(C(i))
               ELSE
                  WRITE(*,'("C(",I1") = ",F10.6," - i ",F10.6)')i,REAL(C(i)),ABS(AIMAG(C(i)))
               END IF
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart effective electronic Hamiltonian"
            DO i = 1, nExStat
               WRITE(*,'(I1,5X)',ADVANCE='NO')i
               DO j = 1, i-1
                  WRITE(*,'(16X,3X)',ADVANCE='NO')
               END DO
               WRITE(*,'(F16.6,1X,"i",1X)',ADVANCE='NO')AIMAG(Hold(i,i))
               IF (i .LT. nExStat) THEN
                  DO j = i+1, nExStat-1
                     WRITE(*,'(F16.6,3X)',ADVANCE='NO')REAL(Hold(i,j))
                  END DO
                  WRITE(*,'(F16.6)')REAL(Hold(i,nExStat))
               ELSE
                  WRITE(*,*)
               END IF
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart non-adiabatic coupling terms"
            WRITE(*,'(1X,10I15)')(j,j=2,nExStat)
            DO i = 1, nExStat-1
               WRITE(*,'(I1,4X,10F15.6)')i,(rdold(i,j),j=i+1,nExStat)
            END DO
            WRITE(*,*)
            WRITE(*,*)"Restart excitation energies"
            DO i = 1, nExStat
               WRITE(*,'(F12.6)')ENERGX(i)
            END DO
            WRITE(*,*)
            WRITE(*,'(" Restart number of CI components: ", I5)')nciconfig
            WRITE(*,*)
            WRITE(*,'(" Restart number of CI active orbitals: ", I3)')n_a_ci_o 
            WRITE(*,*)
            WRITE(*,*)
            WRITE(*,*)"Restart CI components and lexical index"
            DO i = 1, nExStat
               WRITE(*,'(" -- state ",I2," --")')i
               DO j = 1, nciconfig
                  WRITE(*,'(I5,2X,F12.6,20I3)')j,CI_comp(j,i),(old_de(k,j),k=1,n_a_ci_o)
               END DO
            END DO
            WRITE(*,*)
         END IF

         WRITE(*,*)
         WRITE(*,*)"write_res<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!=============================================================================C


      SUBROUTINE vel_scal(iout,NUMAT,AMS,vx,vy,vz,Tfix,vs_Emax,nf)

!  This subroutine performs velocity scaling
!     Velocities are scaled as follows:
!     (1) The kinetic energy corresponding to desidered temperature is calculated
!     (2) Actual kinetic energy is calculated
!     (3) If requested a energy difference limit:
!         (3.1) Check if required energy correction is largen than energy limit, if yes:
!          (3.1.1) change target temperature so that new kinetic energy will
!                   differ from old one exactly of the maximum ammount possible
!     (4) scaling factor conv = SQRT(Tfix/T) with Tfix desired temperature
!     (5) new velocity calculated v = conv*v
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NUMAT
      REAL*8, DIMENSION(NUMAT) :: AMS
      REAL*8 :: Tfix
      REAL*8 :: vs_Emax
      INTEGER :: nf

!  INPUT/OUTPUT VARIABLES
      REAL*8, DIMENSION(NUMAT) :: vx,vy,vz

!  COMMONS
      REAL*8 :: ctime_v
      INTEGER :: wtime_v
      REAL*8 :: ctime_g
      INTEGER :: wtime_g
      REAL*8 :: ctime_sf
      INTEGER :: wtime_sf
      REAL*8 :: ctime_acc
      INTEGER :: wtime_acc

      COMMON /ctimemd/ ctime_v,ctime_g,ctime_sf,ctime_acc
      COMMON /wtimemd/ wtime_v,wtime_g,wtime_sf,wtime_acc

!  LOCAL VARIABLES
      REAL*8 :: KK
      REAL*8 :: K
      REAL*8 :: T_ist
      REAL*8 :: Tfix_old
      REAL*8 :: conv
      INTEGER :: i
      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir

!  PARAMETERS
      REAL*8, PARAMETER :: boltz=1.987168981D-3 ! Boltzmann constant 
                                                ! [kcal/(mol*K)]
      
!-------------C

! TIMING - velocity scaling time, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime1)
         CALL SYSTEM_CLOCK(wtime1,check_ir)
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>vel_scal"
      END IF
   
! (1)
      KK = 0.5D0*Tfix*REAL(nf)*boltz

! (2)
      CALL Kene(iout,NUMAT,AMS,vx,vy,vz,K)

      CALL Tist(iout,K,nf,T_ist)
     
! (3)
      IF (vs_Emax .GE. 0.D0) THEN

! (3.1)
         IF ((K-KK) .GT. vs_Emax) THEN
!(3.1.1)
            Tfix_old = Tfix
            Tfix = (2.D0*(K-vs_Emax))/(nf*boltz)

            IF (iout .GE. 3) THEN
               WRITE(*,*)
               WRITE(*,*)"Change in energy required by velocity"
               WRITE(*,*)"scaling is larger than maximum allowed"
               WRITE(*,'("Temperature will be scaled to ",F10.5," K")')Tfix
            END IF

            ELSE IF((KK-k) .GT. vs_Emax) THEN

               Tfix_old = Tfix
               Tfix = (2.D0*(K+vs_Emax))/(nf*boltz)

               IF (iout .GE. 3) THEN
                  WRITE(*,*)
                  WRITE(*,*)"Change in energy required by velocity"
                  WRITE(*,*)"scaling is larger than maximum allowed"
                  WRITE(*,'("Temperature will be scaled to ",F10.5," K")') Tfix
               END IF

            END IF

         END IF

! (4)
      conv = DSQRT(Tfix/T_ist)

! (5)
      vx = conv*vx
      vy = conv*vy
      vz = conv*vz

!   output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Fixed temperature = ",F10.5," K")')Tfix
         WRITE(*,'("Actual temperature = ",F10.5," K")')T_ist
         WRITE(*,'("Conversion factor = ",F10.5)')conv
         WRITE(*,*)"Scaled velocities [A/ps]"
         WRITE(*,'("Atom",8X,"vx",12X,"vy",12X,"vz")')
         DO i = 1, NUMAT
            WRITE(*,'(I3,2X,F12.5,2X,F12.5,2X,F12.5)')i,vx(i),vy(i),vz(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"vel_scal<<<"
         WRITE(*,*)
      END IF

!  Restore fixed temperature value
      IF (vs_Emax .GE. 0.D0) THEN
         IF (ABS(K-KK) .GT. vs_Emax) THEN
            Tfix = Tfix_old
         END IF
      END IF

      IF (iout .GE. 2) THEN
! TIMING - velocity scaling, end
         CALL CPU_TIME(ctime2)
         CALL SYSTEM_CLOCK(wtime2,check_ir)
         wtime_v = wtime_v + wtime2-wtime1
            IF (ctime2 .GT. REAL(wtime2-wtime1)/REAL(check_ir)+ctime1)ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
         ctime_v = ctime_v + ctime2-ctime1
      END IF

! TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Velocity scaling"
         WRITE(*,*)
         CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
         WRITE(*,'("-----------------------------------")')
      END IF
         
      END SUBROUTINE



!=============================================================================C


      SUBROUTINE estat_fol(iout,fol_stat,nExStat,NCICONF,CICOMP,ENERGX,OLD_CICOMP,old_ENERGX,mdstep,startStep,naco,de,old_de,sw_CC,tully_rest,inde,old_inde)

! This subroutine maps states along MD steps.
! At first an energy criterium is used. If energy variation
! between old and new MD step if enough lower than energy gap 
! between states the "active" state is mapped onto the state
! with closest energy.
! If energy criterium fails the mapping is performed by comparison
! of CI components. Two comparison are made: "overlap" and term by
! term absolute difference. The "active" state is mapped onto the
! state with the higher "overlap" and the lower difference.
! If the two methods yield different results the program stops.
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      USE LIMIT, ONLY: LMGRD, LMCONF 

      IMPLICIT NONE

!  INPUT VARIABLES
      INTEGER :: iout                       ! flag for output verbosity
      INTEGER :: fol_stat                   ! flag to enable state following 
      INTEGER :: nExStat                    ! Number of excited states
      INTEGER :: NCICONF                    ! Number of CI configurations
      REAL*8, DIMENSION(LMCONF,6) :: CICOMP   ! CI components
      REAL*8, DIMENSION(LMGRD) :: ENERGX    ! Energies       
      REAL*8, DIMENSION(LMGRD) :: ENERGX_EXPECT    ! extrapolated Energies       
      INTEGER :: mdstep                       ! actual MD step
      INTEGER :: startStep                       ! First MD step
      INTEGER :: naco                       ! number of active CI orbitals
      INTEGER, DIMENSION(40,LMCONF) :: de     ! lexical indexes
      LOGICAL :: tully_rest                 ! flag to perform tully restart
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP  ! old CI components
      REAL*8, DIMENSION(LMGRD,3) :: old_ENERGX   ! old energies
      INTEGER, DIMENSION(40,LMCONF) :: old_de    ! old lexical indexes
      INTEGER, DIMENSION(nExStat) :: old_inde  ! index mapping old states
  
!  OUTPUT VARIABLES
      INTEGER, DIMENSION(nExStat) :: sw_CC     ! phase of CI components
      INTEGER, DIMENSION(nExStat) :: inde      ! index mapping old states
                                                   ! to new states
!  COMMONS
      REAL*8 :: ctime_v
      INTEGER :: wtime_v
      REAL*8 :: ctime_g
      INTEGER :: wtime_g
      REAL*8 :: ctime_sf
      INTEGER :: wtime_sf
      REAL*8 :: ctime_acc
      INTEGER :: wtime_acc
      INTEGER :: IN2

      COMMON /ctimemd/ ctime_v,ctime_g,ctime_sf,ctime_acc
      COMMON /wtimemd/ wtime_v,wtime_g,wtime_sf,wtime_acc
      COMMON /INOPT2 / IN2(300)

!  LOCAL VARIABLES
      INTEGER :: i,j,l,k,tmp
      REAL*8 :: a,e
      LOGICAL :: chbelow
      REAL*8 :: tmped,maxFactorTot,factorSum
      REAL*8, DIMENSION(nExStat,nExStat) :: a_sign
      REAL*8, DIMENSION(nExStat,nExStat) :: factor
      INTEGER, DIMENSION(nExStat) :: maxPath
      INTEGER, DIMENSION(NCICONF) :: index
      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir

!  PARAMETERS
      REAL*8, PARAMETER :: the2 = 0.1D0

!----------------------C

! TIMING - State following time, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime1)
         CALL SYSTEM_CLOCK(wtime1,check_ir)
      END IF

      IF (iout .GE. 4) THEN
         WRITE(*,*)
         WRITE(*,*)">>>estat_fol"
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)"Energy of excited states [kcal/mol]"
         DO i = 1, nExStat
            WRITE(*,'(I2,2X,F12.6)')i,ENERGX(i)
         END DO
         WRITE(*,*)
         WRITE(*,'("Number of CI components = ",I5)')NCICONF
         WRITE(*,'(6X,"CI components      Lexical index")')
         DO i = 1, nExStat
            WRITE(*,'(" -- state ",I2," --")')i
            DO j = 1, NCICONF
               WRITE(*,'(I5,2X,F12.6,2x,20I3)')j,CICOMP(j,i),(de(k,j),k=1,naco)
            END DO
         END DO
      END IF
         
! Initialization
      DO i = 1, nExStat
         inde(i) = i
      END DO

!  In the first MD step if not tully_restart exit
      IF (mdstep-startStep .LE. 1 .AND. .NOT. tully_rest) THEN
         RETURN        
      END IF   
      
!  Part rewritten by lasse spoerkel
!  call routine to match lexical indexes
      IF(IN2(77).LT.6) THEN
        CALL mindd(iout,NCICONF,naco,de,old_de,index,CICOMP,old_CICOMP,nExStat)
      ELSE
        DO I=1,NCICONF
          index(I) = I
        ENDDO
      ENDIF

      DO j = 2,1,-1
         IF (old_ENERGX(1,j).EQ.0.0D0) old_ENERGX(:,j)=old_ENERGX(:,j+1)
      END DO
!   Calculate the extrapolated energies
      IF (mdstep-startStep .GT. 2) THEN
         DO j = 1, nExStat
            ENERGX_EXPECT(j) = 3.D0*(old_ENERGX(j,3)-old_ENERGX(j,2))+ old_ENERGX(j,1)
            IF (iout .GE. 2) THEN
               WRITE(*,'("Energy expectation value for state ",I2)') j
               WRITE(*,'("Old Energies: ",3F14.6)') old_ENERGX(j,1),old_ENERGX(j,2),old_ENERGX(j,3)
               WRITE(*,'("Extrapolated: ",F14.6)') ENERGX_EXPECT(j)
               WRITE(*,'("Actual      : ",F14.6)') ENERGX(j)
               WRITE(*,*)
            END IF
         END DO
      ELSE
         ENERGX_EXPECT=old_ENERGX(:,3)
      ENDIF

!  Check if more than one state differences are below threshold
      chbelow=.false.
      DO i = 1, nExStat
         tmp=0
         DO j = 1, nExStat
            IF (DABS(ENERGX_EXPECT(i)-ENERGX(j)) .LT. the2) THEN
               tmp=tmp+1
            END IF
         END DO
         IF (tmp .GT. 1) chbelow=.true.
      END DO

!  Start of state mapping
      DO i = 1, nExStat
         factorSum=0.D0
         DO j = 1, nExStat
!  First factor: Energy difference
            tmped = DABS(ENERGX_EXPECT(i)-ENERGX(j))
!       if more are below threshold, set difference to threshold
            IF (chbelow .AND. tmped .LT. the2) tmped=the2

!  Second factor: Overlap of CI components
            a = 0.D0
            DO l = 1, NCICONF
               a = a + old_CICOMP(index(l),i)*CICOMP(l,j)
            END DO
            
            a_sign(i,j) = 1
            IF (a .LT. 0.D0) a_sign(i,j) = -1

!  Total factor
            IF (tmped .GT. 0.D0) THEN
               factor(i,j)=DABS(a)/(tmped)
            ELSE
               factor(i,j)=1.0D0
            END IF
            factorSum = factorSum + factor(i,j)

            IF (iout .GE. 2) THEN
               WRITE(*,'("E_exp(",I1,")-E(",I1,")       ",F15.6)')i,j,DABS(ENERGX_EXPECT(i)-ENERGX(j))
               WRITE(*,'("<",I1,"(t+dt)|",I1,"(t)>      ",F15.6)')i,j,a
               WRITE(*,'("Weighting Factor    ",F15.8)') factor(i,j)
               WRITE(*,*)
            END IF
         END DO

         DO j = 1, nExStat
            factor(i,j)=factor(i,j)/factorSum
         END DO
      END DO

      IF (iout .GE. 2) THEN
         WRITE(*,'("Normalized factors")')
         DO i = 1, nExStat
            DO j = 1, nExStat
               WRITE(*,'(" ",I2," -> ",I2,F15.8)') i,j,factor(i,j)
            END DO
         END DO
         WRITE(*,*)
      END IF

      maxPath=0
      maxFactorTot=0.0D0

!  call the recursive function that finds the best mapping
      CALL find_path(iout,factor,nExStat,maxPath,1,maxPath,maxFactorTot)
      inde=maxPath

      IF (iout .GE. 2) THEN
         write(*,*) 
         write(*,'(A,$)') "Maximum overlap path: "
         do k = 1, nExStat
            write(*,'(I3,$)') maxPath(k)
         end do
         write(*,'(F15.8)') maxFactorTot
      END IF

!  Output
      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"    Excited states MD mapping"
         WRITE(*,'("Old state --> New state")')
         DO i = 1, nExStat
            WRITE(*,'(4X,I1,5X,"-->"5X,I1)')i,inde(i)
         END DO
         WRITE(*,*)
      END IF

!  Check the phase of the CI components
      sw_CC = 1
      DO i = 1, nExStat
         a = 0.D0
         DO j = 1, NCICONF
            a = a + CICOMP(j,inde(i))*old_CICOMP(index(j),i)
         END DO
         IF (a .LT. 0.D0) THEN
            CICOMP(:,inde(i)) = -CICOMP(:,inde(i))
            sw_CC(inde(i)) = -1
            IF (iout .GE. 2) THEN
               write(*,'(A,I2,A,F8.4)')"Inverting CI factors of state ",inde(i),", DOTPRO =",a
               WRITE(*,*)
            END IF
         END IF
      END DO

      if (.not. fol_stat) then
         DO i = 1, nExStat
            inde(i) = i
         END DO
      endif

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"estat_fol<<<"
         WRITE(*,*)
      END IF

      IF (iout .GE. 3) THEN
! TIMING - state following, end
         CALL CPU_TIME(ctime2)
         CALL SYSTEM_CLOCK(wtime2,check_ir)
         wtime_sf = wtime_sf + wtime2-wtime1
            IF (ctime2 .GT. REAL(wtime2-wtime1)/REAL(check_ir)+ctime1)ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
         ctime_sf = ctime_sf + ctime2-ctime1
      END IF

! TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"State following"
         WRITE(*,*)
         CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
         WRITE(*,'("-----------------------------------")')
      END IF

      END SUBROUTINE



      recursive subroutine find_path(iout,factor,nExStat,path,pos,maxPath,maxFactorTot)

      implicit none
!  Recursive subroutine for finding an optimum state mapping.
!  Creates all possible state mappings without a repitition
!  returns mapping with highest overlap factors

!  INPUT VARIABLES
      INTEGER, INTENT(in) :: nExStat,pos,iout
      REAL*8, DIMENSION(nExStat,nExStat), INTENT(in) :: factor
      INTEGER, DIMENSION(nExStat), INTENT(in) :: path

!  INOUTPUT VARIABLES
      REAL*8, INTENT(inout) :: maxFactorTot
      INTEGER, DIMENSION(nExStat), INTENT(inout) :: maxPath

!  LOCAL VARIABLES
      INTEGER, DIMENSION(nExStat) :: pathTemp
      INTEGER :: i,j,k
      REAL*8 :: factorTot
      LOGICAL :: Check

      pathTemp=path

!     Every position can be every state
      do i = 1, nExStat
         pathTemp(pos)=i

         CHECK=.true.
!     Check if the actual state has been assigned to a different position
         do j = 1, pos-1
            if (pathTemp(j) .EQ. i) then
               CHECK=.false.
            end if
         end do

         if (Check) then
!     If not at the end of position, call new find_path with pos+1
            if (pos .lt. nExStat) then
               call find_path(iout,factor,nExStat,pathTemp,pos+1,maxPath,maxFactorTot)
!     Else calculate total factor of mapping and return, if higher than max
            else
               factorTot=0.D0
               do j = 1, pos
                  factorTot=factorTot+factor(j,pathTemp(j))
               end do
               if (factorTot .gt. maxFactorTot) then
                  maxPath=pathTemp
                  maxFactorTot=factorTot
               end if
               if (iout .GE. 2) then
                  write(*,'(A,$)') "Overlap path: "
                  do k = 1, nExStat
                     write(*,'(I3,$)') pathTemp(k)
                  end do
                  write(*,'(F15.8)') factorTot
               end if
            end if
         end if
      end do
      
      return
      end subroutine



!-----------------------------------------------------------

      SUBROUTINE mindd(iout,NCICONF,naco,de,old_de,index,CICOMP,old_CICOMP,nExStat)

!  This subroutine creates an index to match lexical indexes of two
!  subsequent CI runs (de and old_de).
!  Since sometimes the CI components order does not follow the ordering 
!  of lexical index at the end a check is performed by overlap
!  criterium
!--------------------------------------------------------------------------C
! AUTHOR: E. Fabiano                                                       
! MAIL: efabiano@mpi-muelheim.mpg.de                                       
! DATE: March 2007                                                        

      USE LIMIT, ONLY: LMCONF

      IMPLICIT NONE

! INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: NCICONF
      INTEGER :: naco
      INTEGER, DIMENSION(40,LMCONF) :: de
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      REAL*8, DIMENSION(LMCONF,6) :: CICOMP
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP
      INTEGER :: nExStat

! INPUT/OUTPUT VARIABLES
      INTEGER, DIMENSION(NCICONF) :: index

! LOCAL VARIABLES
      INTEGER :: i,j,k
      LOGICAL, DIMENSION(NCICONF) :: done
      LOGICAL :: cheq
      REAL*8 :: a,b,c,d

!-------------C

      IF (iout .GE. 4) THEN
         WRITE(*,*)
         WRITE(*,*)">>>mindd"
         WRITE(*,*)
         WRITE(*,*)"-----------------------------------"
         WRITE(*,*)"         Lexical indexes"
         DO i = 1, NCICONF
            WRITE(*,'(12X,"old:",20I2)')(old_de(j,i),j=1,naco)
            WRITE(*,'("conf: ",I3)')i
            WRITE(*,'(12X,"new:",20I2)')(de(j,i),j=1,naco)
            WRITE(*,*)
         END DO
         WRITE(*,*)"-----------------------------------"
      END IF

      done = .FALSE.

      DO i = 1, NCICONF
         DO j = 1, NCICONF

            IF (.NOT. done(j)) THEN
               cheq = .TRUE.
               DO k = 1, naco
!                 transitions with 1 or 2 are the same. 
!                 Therefore we chose to call the all 1
                  IF (old_de(k,i) .EQ. 2) old_de(k,i) = 1
                  IF (de(k,j) .EQ. 2) de(k,j) = 1

                  IF (old_de(k,i) .NE. de(k,j)) THEN
                     cheq = .FALSE.
                     EXIT
                  END IF
               END DO
            
               IF (cheq) THEN
                  index(i) = j
                  done(j) = .TRUE.
                  EXIT
               END IF

            END IF

         END DO
      END DO

      IF (iout .GE. 4) THEN
         WRITE(*,*)
         WRITE(*,*)"Final index changes:"
         WRITE(*,'("  i",3X,"index(i)")')
         DO i = 1, NCICONF
            IF (i .NE. index(i)) WRITE(*,'(I3," -> ",I3)')i,index(i)
         END DO
         WRITE(*,*)
         WRITE(*,*)"Checking for different ordering of CI components..."
      END IF
      
      b = 0.D0
      d = 0.D0

      DO i = 1, nExStat
         a = 0.D0
         c = 0.D0
         DO j = 1, NCICONF
            a = a + CICOMP(j,i)*old_CICOMP(index(j),i)
            c = c + CICOMP(j,i)*old_CICOMP(j,i)
         END DO
         b = b + DABS(a)
         d = d + DABS(c)
      END DO

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Overall overlap of indexed CI comp:     ",F8.5)')b
         WRITE(*,'("Overall overlap of non-indexed CI comp: ",F8.5)')d
         WRITE(*,'("Diff: ",E15.5)') d-b
      END IF

      IF (DABS(d-b) .GT. 0.0000001D0) THEN
         IF (iout .GE. 3) THEN
            WRITE(*,*)
            WRITE(*,*)"Test failed!!!"
            WRITE(*,*)"CI components have a different ordering"
            WRITE(*,*)"index is set to natural ordering [index(i)=i]"
         END IF
         DO i = 1, NCICONF
            index(i) = i
         END DO
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"mindd<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



!-----------------------------------------------------------

      SUBROUTINE out_time(ct1,ct2,wt1,wt2,ir)

      IMPLICIT NONE

!  INPUT VARIABLES
      REAL*8 :: ct1
      REAL*8 :: ct2
      INTEGER :: wt1
      INTEGER :: wt2
      INTEGER :: ir

!  LOCAL VARIABLES
      REAL*8 :: ctot
      REAL*8 :: wtot
      INTEGER :: c_ore
      INTEGER :: c_min
      REAL*8 :: c_sec
      INTEGER :: w_ore
      INTEGER :: w_min
      REAL*8 :: w_sec

!---------------------C

      ctot = ct2-ct1
      wtot = REAL(wt2-wt1)/REAL(ir)

      c_ore = INT(ctot/3600.D0)
      c_min = INT((ctot-REAL(c_ore)*3600.D0)/60.D0)
      c_sec = ctot - REAL(c_ore)*3600 - REAL(c_min)*60

      w_ore = INT(wtot/3600.D0)
      w_min = INT((wtot-REAL(w_ore)*3600.D0)/60.D0)
      w_sec = wtot - REAL(w_ore)*3600 - REAL(w_min)*60

      WRITE(*,'("cpu-time:  ",I2," h, ",I2," m, ",F7.4," s")')c_ore,c_min,c_sec
      WRITE(*,'("wall-time: ",I2," h, ",I2," m, ",F7.4," s")')w_ore,w_min,w_sec

      END SUBROUTINE




      INTEGER FUNCTION delete(filename)

      IMPLICIT NONE

      INTEGER :: stat
      CHARACTER(150) :: filename

      open(unit=1234, iostat=stat, file=trim(adjustl(filename)),status='old')
      if (stat == 0) then
         close(1234, iostat=stat, status='delete')
         if (stat == 0) then 
            delete = 0
         else
            delete = 2
         end if
      else
         delete = 1
      end if


      END FUNCTION


!-------------------------------------------------------------------
      subroutine langevin(mass,vx,vy,vz,g0,T,grad,natom,dt,nstep,nsstep)
      use random
      implicit none
      integer natom, i, natom3, nstep, nsstep
      real*8  g0,T,dt,a0
      real*8, dimension(natom) :: mass, vx, vy, vz
      real*8, dimension(3,natom) :: grad
      real*8, dimension(:), allocatable, save :: q
      real*8, dimension(:), allocatable:: g
      real*8, parameter :: boltz=1.987168981D-3 ! Boltzmann constant [kcal/(mol*K)]
      real*8, parameter :: sca = 418.4 ! from kcal/(mol*A*AMU) to A/ps^2
      

      allocate(g(natom))
      if(.not. allocated(q)) allocate(q(3*natom))

      if(dt .gt. 1.0d-10) then
        g(:) = g0
      else
        g(:) = mass(:)/dt
      endif

      do i=1,natom
        a0 = g(i)*mass(i)/sca
        grad(1,i) = grad(1,i) + vx(i)*a0
        grad(2,i) = grad(2,i) + vy(i)*a0
        grad(3,i) = grad(3,i) + vz(i)*a0
      enddo

      natom3 = 3*natom
      if(nsstep.eq.1) then
        call ran_num(q,natom3)
      endif
      do i=1,natom
        a0 = sqrt(2.d0*boltz*T*g(i)*mass(i)/dt/sca)
        grad(1,i) = grad(1,i) - q((i-1)*3+1)*a0 
        grad(2,i) = grad(2,i) - q((i-1)*3+2)*a0
        grad(3,i) = grad(3,i) - q((i-1)*3+3)*a0
      enddo

      deallocate(g)
      return
      end subroutine langevin

