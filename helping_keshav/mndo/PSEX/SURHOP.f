      SUBROUTINE SURHOP (iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,
     1                   dt,vx,vy,vz,mass,N,istat,ENERGX,CGX,ihop,
     2                   write_hop,hopunit,Hold,rdold,C,NE,old_CICOMP,
     3                   old_ENERGX,old_de,num_CC,an_CC,Einteg,
     4                   tully_rest,istatt,rnd_gen,rnd_count,fol_stat,
     5                   dec_cor,nac_ph,cuthop,inde,sw_CC)

C==========================================================================C
C This subroutine implements Tully's surface hopping. For a dettailed      C
C description see Ref. [1,2].                                              C
C                                                                          C
C The subroutine is divided into five parts:                               C
C   1- Multiple CI runs to get needed data [SCFMUL]                        C
C   2- Propagate velocity half time step forward: in this way velocity and C
C        coordinates are computed at the same time step                    C
C   3- Computation of nonadiabatic coupling vectors [SCF]                  C
C   4- Evaluation of states' population (C(t) in Ref. [1,2])               C
C     4.1 - create effective hamiltonian for propagation of population     C
C             (see Eq. (26) in Ref. [2]                                    C
C     4.2 - propagate population (solve Eq. (26) of Ref. [2]) and compute  C
C             transition probabilities (Eq. (29) of Ref. [2])              C
C   5- Fewest switch algorithm                                             C
C     5.1 - If hopping happens adjust velocity to get energy consevation   C
C                                                                          C
C  Note that both the effective hamiltonian used for propagation of        C
C  population and the dot product of velocity times nonadiabatic vectors   C
C  used for computation of transition probabilities, are time dependent.   C
C  Since for solution of Eq. (26) and computation of integral (29) the     C
C  value of these two quantities is required at many times intermediates   C
C  between two successive time points used for MD propagation, their value C
C  is linearly interpolated between two successive MD steps.               C
C     H(t) = H((n-1)dt) + (H(ndt) - H((n-1dt)))*((t-(n-1)dt)/dt)           C
C  The same for r*d                                                        C
C  To this end we need the value of each at the previous time step.        C
C  To do this at the first MD step we just compute nonadiabatic vectors    C
C  and Hamiltonian, store them and exit.                                   C
C  In next steps the old value can be used and at the end actual value     C 
C  becomes old value for next step.                                        C
C                                                                          C
C  It is possible to compute non-adiabatic coupling terms both analyticallyC
C  and numerically. The numerical procedure is based on formula (28) of    C
C  Ref. [2]. The analytical procedure instead computes non-adiabatic       C
C  coupling vectors and then the dot product with velocity is performed.   C
C  Since numerical procedure does not produce non-adiabatic coupling       C
C  vectors and these are needed for velocity adjustment if hopping is      C
C  performed, in this case (quite rare) the analytical procedure is called C
C  to get the vector.                                                      C
C  The two treatment can also be used toghether. In this case only results C
C  from the analytical calculations are used. Numerical results are only   C
C  printed out (if iout .GE. 1) for data analysis and debugging purposes.  C
C                                                                          C
C REFERENCES                                                               C
C [1] J.C. Tully, J Chem Phys 93, 1061 (1990)                              C
C [2] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C                                                                          C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C ARRAY         Scratch array (for use in SCFCAL)                          C
C C             Electronic population                                      C
C CCX           Nonadiabatic couplingg vectors for all pairs of states     C
C CG            Cartesian gradient                                         C
C CGX           Cartesian gradients of computed states                     C
C CNORMX        Cartesian norms of computed states                         C
C csi           Random number for fewest switch algorithm                  C
C DE            Energy difference between states involved in hopping       C
C delta         Working variable for propagation of velocity               C
C dt            Time step [in picoseconds]                                 C
C Einteg        Integration algorithm to use for electronic equation in    C
C                Tully's hopping:                                          C
C                    euler  -> Euler method                                C
C                    RK2    -> Runge-Kutta second-order algorithm          C
C                    RK4    -> Runge-Kutta fourth-order algorithm          C
C                    ABM4   -> Adams-Bashforth-Moulton fourth order        C
C                              predictor-corrector with error estimation   C
C                    ABM2   -> Adams-Bashforth-Moulton second order        C
C                              predictor-corrector                         C
C                    ABM5   -> Adams-Bashforth-Moulton fifth order         C
C                              predictor-corrector with error estimation   C
C                    imp    -> implicit middle-point algorithm             C
C                    UP1    -> Unitary propagator evaluated at middle-pointC
C                              built from sylvester theorem                C
C                    UP2    -> Unitary propagator evaluated at the         C
C                              beginning of time interval                  C
C                    UP3    -> Unitary propagator evaluated at middle-pointC  
C                              bulit from taylor expansion                 C
C ENERGX        Energies of computed states                                C
C F2            Vector along velocity must be adjusted (is equal to the    C
C                 nonadiabatic coupling vector of the two states involved  C
C                 in hopping)                                              C
C g             Transition probabilities                                   C
C GAMMA         Correction factor for velocity                             C
C H             Effective Hamiltonian for time propagation of population   C
C hbar          Plank's constant/2*pi in (kcal/mol)ps                      C
C Hold          Effective hamiltonian for tully's hopping at previous      C
C                MD step                                                   C
C hopunit       Fortran unit for hopping data file                         C
C i             Counter                                                    C
C ICALL         control variable for SCF calculations and error flag       C
C ICALL2        control variable and error flag for SCF calculation of     C
C                 nonadiabatic copuling vectors                            C
C IGRST                                                                    C
C ihop          Control flag for hopping. If .TRUE. hopping happend        C
C IN2           MNDO options                                               C
C istat         Current state to propagate                                 C
C ISTATE        First index for nonadiabatic coupling vector to be         C
C                 computed by SCF                                          C
C j             Counter                                                    C
C JSTATE        Second index for nonadiabatic coupling vector to be        C
C                 computed by SCF                                          C
C l             Counter                                                    C
C LM1           Maximum number of atoms (see LIMIT.f90)                    C
C LM1M          Maximum number of external points (see LIMIT.f90)          C
C LM5           Buffer lenght (for use in SCFCAL)                          C
C LMGRD         Maximum number of GUGA-CI gradient vectors (see LIMIT.f90) C
C LMV           Muximum number of geometrical variables (see LIMIT.f90)    C
C mass          Atomic masses                                              C
C mdstep        current MD step                                            C
C startStep     Fisrt MD step                                              C
C N             Number of atoms                                            C
C NE            Number of steps for electronic integration (in RK2)        C
C vx,vy,vz   Cartesian components of velocity at half time step foorward   C
C prob          Total probability of transition                            C
C rd            Dot products of velocity times nonadiabatic vectors        C
C rdold         Dot product of velocity and non-adiabatric vectors         C
C                at previosu MD step                                       C
C sca           conversion factor from kcal/(mol*A*AMU) to A/ps^2          C
C SCFCAL        External routine for SCF calculations                      C
C tmph          Scratch variable                                           C
C tmpin2_160    Temporary scratch variable                                 C
C tmrd          Scratch variable                                           C
C vx,vy,vz      Cartesian components of velocity                           C
C write_hop     Flag to allow writing hopping data to file                 C
C                                                                          C
C dec_cor       empirical decoherence correction by Persico & Granucci     C
C               (standard value = 0.1)                                     C
C nac_ph        threshold for nonadiabatic coupling phase following        C
C               algorithm (0,0=off, standard=0.1)                          C
C cuthop        cutt-off Energy difference to avoid unphysical hops        C
C               >0 : simple cut at cuthop value                            C
C               <0 : cut by gaussian scale of probability cenetered at     C
C               cuthop value                                               C
C--------------------------------------------------------------------------C
C  CALLINGS                                                                C
C                                                                          C
C  SCFMUL                    calculates energy and gradients               C
C  estat_fol                 Mapping of states                             C
C  m_num_CC                  Computes numerical nonadiabatic couplings     C
C  SCF                       Calculates nonadiabatic coupling vectors      C
C  int_euler                 Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_RK2                   Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_RK4                   Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_ABM4                  Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_ABM2                  Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_ABM5                  Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_imp                   Integration of Eq. (26 and (29) of Ref. [2]   C
C  int_UP                    Integration of Eq. (26 and (29) of Ref. [2]   C
C  pop_norm                  Normalization of quantum amplitudes           C
C  RANDOM_NUMBER             Standard FORTRAN random number generator      C
C  nran2                     PM-BD random number generator                 C
C  nran3                     Knuth random number generator                 C
C  ADJEK2                    verifies if hopping is energetically allowed  C
C                              and calculates correction factor for        C
C                              velocity (GAMMA)                            C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  C
C == REVIEW ==                                                             C
C Added callings to ABM4, ABM2, ABM5, imp, UP1, UP2 and UP3                C
C                                                                          C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C

      USE LIMIT, ONLY: LM1, LMGRD, LM1M, LMV, LMCONF

      IMPLICIT NONE

C  COMMON
      INTEGER :: IN2
      INTEGER :: IGRST
      INTEGER :: ISTATE
      INTEGER :: JSTATE
      REAL*8 :: CG
      INTEGER :: NCICONF
      REAL*8 :: CICOMP
      INTEGER :: dee
      INTEGER :: naco
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

      COMMON
     ./INOPT2/ IN2(300)
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./CGRAD / CG(3,LM1+LM1M)
     ./DYNVA/ CICOMP(LMCONF,6),NCICONF
     ./DYNVA2/ dee(40,LMCONF),naco
     ./ctimemd/ ctime_v,ctime_g,ctime_sf,ctime_acc
     ./wtimemd/ wtime_v,wtime_g,wtime_sf,wtime_acc
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: LM5
      REAL*8, DIMENSION(LM5) :: ARRAY
      INTEGER :: ICALL
      INTEGER :: mdstep, startStep
      REAL*8 :: dt
      INTEGER :: N
      REAL*8, DIMENSION(N) :: mass
      LOGICAL :: write_hop
      INTEGER :: hopunit
      INTEGER :: NE
      LOGICAL :: num_CC
      LOGICAL :: an_CC
      CHARACTER(5) :: Einteg
      LOGICAL :: tully_rest
      CHARACTER(8) :: rnd_gen
      LOGICAL :: fol_stat
      REAL*8  :: dec_cor
      REAL*8  :: nac_ph
      REAL*8  :: cuthop
      INTEGER, DIMENSION(IN2(159)) :: sw_CC
      INTEGER, DIMENSION(IN2(159)) :: inde
      
C  INPUT/OUTPUT VARIABLES
      REAL*8, DIMENSION(N) :: vx,vy,vz
      INTEGER :: istat
      REAL*8, DIMENSION(LMGRD) :: ENERGX
      REAL*8, DIMENSION(3,LM1+LM1M,LMGRD) :: CGX
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: Hold
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rdold
      COMPLEX*16, DIMENSION(IN2(159)) :: C
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP
      REAL*8, DIMENSION(LMGRD,3) :: old_ENERGX
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      INTEGER :: rnd_count

C  OUTPUT VARIABLES
      LOGICAL :: ihop
      INTEGER :: istatt

C  LOCAL VARIABLES 
      INTEGER :: i,j,l,k
      REAL*8, DIMENSION(LMGRD) :: CNORMX
      INTEGER :: tmpin2_160
      REAL*8, DIMENSION(3,N,IN2(159),IN2(159)) :: CCX
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: H
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rd
      REAL*8, DIMENSION(IN2(159)) :: g
      REAL*8 :: csi
      REAL*8 :: prob,g_temp
      REAL*8 :: DE
      REAL*8, DIMENSION(3,N) :: F2
      INTEGER :: ICICAL
      REAL*8 :: delta
      INTEGER :: ICALL2
      REAL*8 :: tmph
      REAL*8 :: GAMMA
      REAL*8 :: tmprd
      REAL*8 :: sws
      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir
      REAL*8 :: factor,scaling,E_kin,E_diff,v1_sq,v2_sq

C  EXTERNAL ROUTINES
      EXTERNAL SCFCAL

C  PARAMETERS
      REAL*8, PARAMETER :: hbar = 1.517873021D-2 ! hbar in (kcal/mol)ps
      REAL*8, PARAMETER :: sca = 418.4 ! from kcal/(mol*A*AMU) to A/ps^2

C--------------------------------------------------------------------------C

      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)">>>SURHOP"
         WRITE(*,*)
      END IF

C  (3)

C  Compute nonadiabatic coupling vector

C      NUMERICAL
      IF (num_CC) THEN
         CALL m_num_CC(iout,IN2(159),NCICONF,CICOMP,OLD_CICOMP,dt,
     &                mdstep,startStep,naco,dee,old_de,rd,tully_rest,
     &                inde)
      END IF

C TIMING - Analytical nonadiabatic coupling vectors, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime1)
         CALL SYSTEM_CLOCK(wtime1,check_ir)
      END IF

C      ANALYTICAL
C  Uses calling to SCF with the following special conventions:
C   1- ICALL = ICALL2 = 3  (used in SCFCAL, local)
C   2- (ICROSS) IN2(160) = 16       (used in PSDST1 and PSOMX, reset)
      IF (an_CC) THEN
         CCX = 0.D0
         ICALL2 = 3

         DO i = 2, IN2(159)
            DO j = 1, i-1
               ISTATE = IGRST(inde(i))
               JSTATE = IGRST(inde(j))
               tmpin2_160 = IN2(160) 
               IN2(160) = 16

               CALL SCF(ARRAY,LM5,ICALL2,SCFCAL)
               
               IN2(160) = tmpin2_160

               IF (ICALL2 .EQ. -1) STOP 'ICALL2=-1 in SURHOP'

               sws = 1.D0
C  Setting phase of CCX vectors 
               IF (mdstep-startStep .GE. 2 .or. tully_rest) THEN
                  sws = REAL(sw_CC(inde(i))*sw_CC(inde(j)))
               END IF

C  Make index ranges explicit to avoid array boundary violations (AK).
               CCX(1:3, 1:N, i, j) =  sws * CG(1:3, 1:N)
               CCX(1:3, 1:N, j, i) = -sws * CG(1:3, 1:N)

            END DO
         END DO

         IF (iout .GE. 2) THEN
C TIMING - Analytical nonadiabatic coupling vectors, end
            CALL CPU_TIME(ctime2)
            CALL SYSTEM_CLOCK(wtime2,check_ir)
            wtime_acc = wtime_acc + wtime2-wtime1
               IF (ctime2 .GT. 
     &               REAL(wtime2-wtime1)/REAL(check_ir)+ctime1) 
     &               ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
            ctime_acc = ctime_acc + ctime2-ctime1
         END IF

C TIMING - output
         IF (iout .GE. 3) THEN
            WRITE(*,*)
            WRITE(*,'("-----------------------------------")')
            WRITE(*,*)"*** TIMING ***"
            WRITE(*,*)"Analytical nonadiabatic coupling vectors"
            WRITE(*,*)
            CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
            WRITE(*,'("-----------------------------------")')
         END IF

         IF (iout .GE. 3) THEN
            WRITE(*,*)
            WRITE(*,*)"     Nonadiabatic coupling vectors"
            DO i = 1, IN2(159)-1
               DO j = i+1, IN2(159)
                  WRITE(*,'(" --- CCX(",I1,",",I1,") ---")')i,j
                  WRITE(*,'("Atom",11X,"CCx",11X,"CCy",11X,"CCz")')
                  DO l = 1, N
                     WRITE(*,'(1X,I3,4X,F12.6,2X,F12.6,2X,F12.6)')
     &                       l,CCX(1,l,i,j),CCX(2,l,i,j),CCX(3,l,i,j)
                  END DO
               END DO
               WRITE(*,*)
            END DO
         END IF

      END IF

C  (4)
C
C  (4.1)

C TIMING - Effective Hamiltonian computation, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime1)
         CALL SYSTEM_CLOCK(wtime1,check_ir)
      END IF

C  Compute Effective Hamiltonian and Rd products
      H = DCMPLX(0.D0,0.D0)

C     Diagonal part
      DO i = 1, IN2(159)
         tmph = (ENERGX(inde(i))-ENERGX(inde(1)))/hbar
         H(i,i) = DCMPLX(0.D0,-tmph)
      END DO

C*********************************************************
C Decoherence correction implemented by spoerkel, 06/2013
C Taken from 
C N. Shenvi, J.E. Subotnik, W. Yang, J. Chem. Phys. 135, 024101 (2011)
      if (dec_cor .eq. -1) then
         IF (iout .GE. 2) THEN
            write(*,*) 
            write(*,*) "Decoherence correction"
         END IF
         factor=2.390057361377D-3 ! A^2/ps^2*AMU -> kcal/mol
         do i = 1, IN2(159)
            !Calculate kinetic energy for all atoms
            v1_sq = 0.0D0
            do j=1,N
               v1_sq = v1_sq+
     &            (vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)) * mass(j)
            end do
            
            E_kin = v1_sq*factor/2.0D0

            E_diff = (ENERGX(inde(i))-ENERGX(inde(istat)))

            scaling = dsqrt( abs(1.0D0 - E_diff/E_kin) )
            tmph = (-2.0D0*E_kin*scaling) / hbar

            IF (iout .GE. 2) THEN
               write(*,'(3(A,F14.6))') "Ekin: ",E_kin," Ediff: ",E_diff,
     &                                  " scal: ",scaling

               write(*,'(A,F16.6,SP,F16.6,A)') "   Hold: ",
     &                                  dreal(H(i,i)),dimag(H(i,i)),"*i"
            END IF

            if (E_kin .gt. E_diff) then
               H(i,i) = dcmplx(0.D0,tmph)
            else
               H(i,i) = dcmplx(tmph,0.D0)
            end if

            IF (iout .GE. 2) THEN
               write(*,'(A,F16.6,SP,F16.6,A)') "   Hnew: ",
     &                                  dreal(H(i,i)),dimag(H(i,i)),"*i"
            END IF

         end do
      end if
C*********************************************************

C     Off-diagonal part
      IF (an_CC) THEN

         rd = 0.D0

         DO i = 2, IN2(159)
            DO j = 1, i-1
               DO l = 1, N
                  rd(i,j) = rd(i,j) + CCX(1,l,i,j)*vx(l) + 
     &                 CCX(2,l,i,j)*vy(l) + CCX(3,l,i,j)*vz(l)
               END DO
               rd(j,i) = -rd(i,j)
               H(i,j) = DCMPLX(-rd(i,j),0.D0)
               H(j,i) = -H(i,j)
            END DO
         END DO

      ELSE
         DO i = 2, IN2(159)
            DO j = 1, i-1
               H(i,j) = DCMPLX(-rd(i,j),0.D0)
               H(j,i) = -H(i,j)
            END DO
         END DO
      END IF

      IF (iout .GE. 2) THEN
C TIMING - Effective Hamiltonian computation, end
         CALL CPU_TIME(ctime2)
         CALL SYSTEM_CLOCK(wtime2,check_ir)
         wtime_eh = wtime_v + wtime2-wtime1
            IF (ctime2 .GT. REAL(wtime2-wtime1)/REAL(check_ir)+ctime1) 
     &          ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
         ctime_eh = ctime_v + ctime2-ctime1
      END IF

C TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Effective Hamiltonian computation"
         WRITE(*,*)
         CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
         WRITE(*,'("-----------------------------------")')
      END IF

C OUTPUT
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"       Non-adiabatic coupling terms"
         WRITE(*,*)"     [dot product of velocity and CCX]"
         WRITE(*,'(1X,10I15)')(j,j=2,IN2(159))
         DO i = 1, IN2(159)-1
            WRITE(*,'(I1,4X,10F15.6)')i,(rd(i,j),j=i+1,IN2(159))
         END DO
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"       Effective electronic Hamiltonian"
         DO i = 1, IN2(159)
            WRITE(*,'(I1,5X)',ADVANCE='NO')i
            DO j = 1, i-1
               WRITE(*,'(16X,3X)',ADVANCE='NO')
            END DO
            WRITE(*,'(F16.6,1X,"i",1X)',ADVANCE='NO')AIMAG(H(i,i))
            IF (i .LT. IN2(159)) THEN
               DO j = i+1, IN2(159)-1
                  WRITE(*,'(F16.6,3X)',ADVANCE='NO')REAL(H(i,j))
               END DO
               WRITE(*,'(F16.6)')REAL(H(i,IN2(159)))
            ELSE
               WRITE(*,*)
            END IF
         END DO
      END IF


C  If mdstep .LE. 1 (2 for numerical CC) then store H in Hold and 
C  rd in rdold and exit
      IF (.NOT. tully_rest) THEN
         IF ((.NOT. an_CC .AND. mdstep-startStep .LE. 2) .OR. 
     &             (an_CC .AND. mdstep-startStep .LE. 1)) THEN
            Hold = H
            rdold = rd
            C = DCMPLX(0.D0,0.D0)
            C(istat) = DCMPLX(1.D0,0.D0)
            RETURN
         END IF
      END IF

C   (4.2)

C  Solve equation for electronic population and integrate for probabilities
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"     Starting population coefficients"
         WRITE(*,'("State",6X,"REAL(C)",4X,"IM(C)")')
         DO i = 1, IN2(159)
            WRITE(*,'(2X,I1,5X,2F10.6)')i,C(i)
         END DO
      END IF

      IF (Einteg .EQ. "euler") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"Euler algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_euler(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,istat,
     &                  C,g)
      ELSE IF (Einteg .EQ. "RK2") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"RK2 algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_RK2(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,istat,
     &                C,g)
      ELSE IF (Einteg .EQ. "RK4") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"RK4 algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_RK4(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,istat,
     &                C,g)  

      ELSE IF (Einteg .EQ. "ABM4") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"ABM4 algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_ABM4(iout,mdstep,startStep,N,IN2(159),H,Hold,rd,rdold,
     &                 dt,NE,istat,C,g)
      ELSE IF (Einteg .EQ. "ABM2") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"ABM2 algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_ABM2(iout,mdstep,startStep,N,IN2(159),H,Hold,rd,rdold,
     &                    dt,NE,istat,C,g)
      ELSE IF (Einteg .EQ. "ABM5") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"ABM5 algorithm will be used for integration"
            WRITE(*,*)"of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_ABM5(iout,mdstep,startStep,N,IN2(159),H,Hold,rd,rdold,
     &                    dt,NE,istat,C,g)
      ELSE IF (Einteg .EQ. "imp") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"Implicit midle-point algorithm will be used" 
            WRITE(*,*)"for integration of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_imp(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,istat,C,g)
      ELSE IF (Einteg .EQ. "UP1") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"Unitary propagator algorithm with"
            WRITE(*,*)"midle-point time appoximation will be used" 
            WRITE(*,*)"for integration of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_UP(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,
     &                    istat,C,g,1,vx,vy,vz,mass,ENERGX,dec_cor)
      ELSE IF (Einteg .EQ. "UP2") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"Unitary propagator algorithm with"
            WRITE(*,*)"initial time appoximation will be used" 
            WRITE(*,*)"for integration of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_UP(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,
     &                    istat,C,g,2,vx,vy,vz,mass,ENERGX,dec_cor)
      ELSE IF (Einteg .EQ. "UP3") THEN
         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"Unitary propagator algorithm with"
            WRITE(*,*)"direct construction of exp  will be used" 
            WRITE(*,*)"for integration of electronic equation"
            WRITE(*,*)
         END IF
         CALL int_UP(iout,N,IN2(159),H,Hold,rd,rdold,dt,NE,
     &                    istat,C,g,3,vx,vy,vz,mass,ENERGX,dec_cor)
      ELSE
         WRITE(*,*)Einteg
         WRITE(*,*)"Bad definition of electronic integrator!"
         WRITE(*,*)
         STOP 'Abnormal termination in SURHOP'
      END IF

      DO i = 1, IN2(159)
         IF (g(i) .GT. 1.D0) THEN
            write(*,*)"warning g=",g(i),i
            g(i) = 1.D0
         END IF
      END DO
      
      CALL pop_norm(iout,IN2(159),C)

      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,*)"    Population of excited states"
         WRITE(*,'("State",2X,"population")')
         DO i = 1, IN2(159)
            WRITE(*,'(2X,I1,5X,F8.5)')i,REAL(C(i)*CONJG(C(i)))
         END DO
         WRITE(*,*)
         WRITE(*,*)"    Transition probabilities"
         WRITE(*,'("State",2X,"probability")')
         DO i = 1, IN2(159)
            WRITE(*,'(2X,I1,5X,F9.6)')i,g(i)
         END DO
      END IF

C  update Hold and rdold
      Hold = H
      rdold = rd

C  (5)
C  Fewest switch criterion
C    SELECT RANDOM NUMBER
      IF (rnd_gen .EQ. "standard") THEN
         CALL RANDOM_NUMBER(csi)
      ELSE IF (rnd_gen .EQ. "PM_BD") THEN
         CALL nran2(csi,1,.FALSE.)
      ELSE IF (rnd_gen .EQ. "knuth") THEN
         CALL nran3(csi,1,.FALSE.)
      END IF

      rnd_count = rnd_count + 1

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Random number for hopping: ",F9.6)')csi
      END IF

C   Write hopping data if requested
      IF (write_hop) THEN
         OPEN(77,STATUS='REPLACE',FILE='hop.tmp')
         WRITE(77,'("mdstep: ",I6,3X,"state: ",I1,3X,"ordered: "
     &                      ,I1)')mdstep,istat,istatt
         WRITE(77,'("energy",3F15.5)') (ENERGX(inde(j)),
     &                                       j=1,IN2(159))
         WRITE(77,'("energy differences")')
         DO i = 2, IN2(159)
            WRITE(77,'(4F15.5)') 
     &              (ENERGX(inde(i))-ENERGX(inde(j)),j=1,i-1)
         END DO
         WRITE(77,'("coupling")')
         DO i = 2, IN2(159)
            WRITE(77,'(4F15.5)') (rd(i,j),j=1,i-1)
         END DO
         WRITE(77,'("population   probability")')
         DO j = 1, IN2(159)
            write(77,'(1X,F8.5,5X,F8.5)')
     &                REAL(C(j)*DCONJG(C(j))),g(j)
         END DO
         WRITE(77,'("rnd number",F8.5)')csi
         WRITE(77,*)"---------------------------------"
         CLOSE(77)
      END IF

C       bug correction by spoerkel
      prob = 0.D0
      write(*,*) 
      write(*,*) "CUTHOP = ",cuthop
      write(*,*) "Random number = ",csi
      write(*,*) 

      DO j = 1, IN2(159)

C       can not jump on same state
         if (j .eq. istatt) cycle

c       only work with temporary probabilities
         g_temp=g(j)

         DE = ENERGX(istatt) - ENERGX(inde(j))

         WRITE(*,'(2(A,I3))') "Test for ",istat,"->",j
         WRITE(*,'(A,F12.6)') "   Absolute EDiff=",abs(DE)

C       gaussian damping only applies to the actual probability
         IF (cuthop.lt.0) THEN
            IF ((abs(DE).ge.abs(cuthop))) THEN
              WRITE(*,'(2A)') "   SPURIOUS HOP, damping by",
     &          " function exp(-0.02*(dE+cuthop)**2)"
              WRITE(*,'(A,E12.4,A)') "     Old probability=",g_temp
              g_temp=g_temp*exp(-0.02*(abs(DE)-abs(cuthop))**2)
              WRITE(*,'(A,E12.4,A)') "     New probability=",g_temp
            END IF
         END IF

C       actual hopping probability is set to 0 for spurious hop
         IF ((cuthop.gt.0).and.(abs(DE) .ge. cuthop)) THEN
            WRITE(*,'(A,F9.6,A)') "   SPURIOUS HOP, set probability "
     &          ,g_temp," to: 0.0"
            g_temp=0.0
         END IF

         prob = prob + g_temp

         WRITE(*,'(2(A,F9.6))')   "   g(j)=",g(j)," g_actual=",g_temp
         WRITE(*,'(2(A,F9.6),A)') "   Is ",csi," < ",prob,"?"

         WRITE(*,*)

C       a hop occurs, if random number points to the actual 
C       probability in the full probability catalog
         IF (csi .LE. prob) THEN

            WRITE(*,*)"DeltaE is:",DE

            ihop =.TRUE.        !JUMP

C.....................
C           If only numerical non-adiabatic-coupling is used
C           CCX are not computed. Anyway they are needed here.
C           Therefore only in this cases when hopping happens
C           we use analytical procedure to get CCX

            IF (.NOT. an_CC) THEN
               
C TIMING - Analytical nonadiabatic coupling vectors, start
               IF (iout .GE. 2) THEN
                  CALL CPU_TIME(ctime1)
                  CALL SYSTEM_CLOCK(wtime1,check_ir)
               END IF

               CCX = 0.D0
            
               ICALL2 = 3
            
               ISTATE = IGRST(istatt)
               JSTATE = IGRST(inde(j))

C  Uses calling to SCF with the following special conventions:
C   1- ICALL = ICALL2 = 3  (used in SCFCAL, local)
C   2- (ICROSS) IN2(160) = 16       (used in PSDST1 and PSOMX, reset)
               tmpin2_160 = IN2(160) 

               IN2(160) = 16

               CALL SCF(ARRAY,LM5,ICALL2,SCFCAL)

               IN2(160) = tmpin2_160

               IF (ICALL2 .EQ. -1) STOP 'ICALL2=-1 in SURHOP'

C              Set phase of CCX vector
               sws = REAL(sw_CC(istatt)*sw_CC(inde(j)))
                  
               CCX(:,:,istat,j) = sws*CG

C              output
               IF (iout .GE. 3) THEN
                  WRITE(*,*)
                  WRITE(*,*)"     Non-adiabatic coupling vector"
                  WRITE(*,'(" --- CCX(",I1,",",I1,") ---")')istat,j
                  WRITE(*,'("Atom",11X,"CCx",11X,"CCy",11X,"CCz")')
                  DO l = 1, N
                     WRITE(*,'(1X,I3,4X,F12.6,2X,F12.6,2X,F12.6)')
     &                       l,CCX(1,l,istat,k),CCX(2,l,istat,j),
     &                       CCX(3,l,istat,k)
                  END DO
                  WRITE(*,*)
               END IF

               IF (iout .GE. 2) THEN
C TIMING - Analytical nonadiabatic coupling vectors, end
                  CALL CPU_TIME(ctime2)
                  CALL SYSTEM_CLOCK(wtime2,check_ir)
                  wtime_acc = wtime_acc + wtime2-wtime1
                     IF (ctime2 .GT. REAL(wtime2-wtime1)/
     &                REAL(check_ir)+ctime1)
     &                ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
                  ctime_acc = ctime_acc + ctime2-ctime1
               END IF

C TIMING - output
               IF (iout .GE. 3) THEN
                  WRITE(*,*)
                  WRITE(*,'("-----------------------------------")')
                  WRITE(*,*)"*** TIMING ***"
                  WRITE(*,*)"Analytical nonadiabatic coupling vectors"
                  WRITE(*,*)
                  CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
                  WRITE(*,'("-----------------------------------")')
               END IF

            END IF

C.....................
            DO l = 1, N
               F2(1,l) = CCX(1,l,istat,j)
               F2(2,l) = CCX(2,l,istat,j)
               F2(3,l) = CCX(3,l,istat,j)
            END DO

C  (5.1)

C TIMING - Velocity adjustment after hopping, start
            IF (iout .GE. 2) THEN
               CALL CPU_TIME(ctime1)
               CALL SYSTEM_CLOCK(wtime1,check_ir)
            END IF

            CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,DE,j,istat,GAMMA,
     &                  mdstep)

            istatt = inde(istat)

            DO i = 1, N
               vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
               vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
               vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
            END DO

            IF (iout .GE. 2) THEN
C TIMING -Velocity adjustment after hopping , end
               CALL CPU_TIME(ctime2)
               CALL SYSTEM_CLOCK(wtime2,check_ir)
               wtime_v = wtime_v + wtime2-wtime1
                  IF (ctime2 .GT. REAL(wtime2-wtime1)/
     &              REAL(check_ir)+ctime1) 
     &              ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
               ctime_v = ctime_v + ctime2-ctime1
            END IF

C TIMING - output
            IF (iout .GE. 3) THEN
               WRITE(*,*)
               WRITE(*,'("-----------------------------------")')
               WRITE(*,*)"*** TIMING ***"
               WRITE(*,*)"Velocity adjustment after hopping"
               WRITE(*,*)
               CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
               WRITE(*,'("-----------------------------------")')
            END IF

C OUTPUT            
            IF (iout .GE. 2) THEN
               WRITE(*,*)
               WRITE(*,'("+------------------------------
     &-----------------+")')
               WRITE(*,'("|",1X,"MD step: ",I4," - update of velocity
     & after hop.",1X,"|")')mdstep
               WRITE(*,'("|",47X,"|")')
               WRITE(*,'("| Atom",8X,"vx",12X,"vy",12X,"vz",4X,"|")')
               DO i = 1, N
                  WRITE(*,'("|",1X,I3,2X,F12.5,2X,F12.5,2X,F12.5," |")')
     &               i,vx(i),vy(i),vz(i)
               END DO
               WRITE(*,'("+------------------------------
     &-----------------+")')
            END IF

            EXIT

         END IF

      END DO

      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"SURHOP<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_RK2(iout,N,Ngrd,H,Hold,rd,rdold,dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements second order Runge-Kutta algorithm for the    C
C solution of Eq. (26) of REf. [1] and at the same time it performs the    C
C integral in Eq (29) of Ref. [1].                                         C
C For Runge-Kutta algorithm see for example Ref. [2]                       C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] "Numerical Recipes in FORTRAN 77: The Art of Scientific Computing",  C
C      by William H. Press, Brian P. Flannery, Saul A. Teukolsky,          C
C      and William T. Vetterling                                           C
C     also: http://www.library.cornell.edu/nr/bookfpdf/f16-1.pdf           C
C                                                                          C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Quantum amplitudes                                         C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C iout          Flag for output verbosity                                  C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k1,k2         Intermediates used in Runge-Kutta method                   C
C l             Counter                                                    C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16, DIMENSION(Ngrd) :: k1,k2
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh
      
      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C--------------------------------------------------------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_RK2>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C   Initializatio of variables

      d = dt/REAL(NE)

      g = 0.D0
      
C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C   Start cycle over subintervals
      DO i = 0, NE-1

         old_C = C

C      Computation of quantum amplitudes
         DO j = 1, Ngrd
            k1(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
               k1(j) = k1(j) + hh*c(l)
            END DO
            k1(j) = k1(j)*d
         END DO


         DO j = 1, Ngrd
            k2(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+0.5D0)/REAL(NE)
               k2(j) = k2(j) + hh*(c(l)+0.5D0*k1(l))
            END DO
            k2(j) = k2(j)*d
         END DO

         DO j = 1, Ngrd
            C(j) = C(j) + k2(j)
         END DO

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF
         
C      Computation of probabilities
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO
         
C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over subintervals    
      END DO

C   TIMING - C integration time, end
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for RK2 integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")')
      END IF

C   Normalization of probabilities 
      rh = C(istat)*DCONJG(C(istat))
      
      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_RK2<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE ADJEK2(iout,N,DIR,mass,vx,vy,vz,dt,DE,j,istat,GAMMA,
     &                  mdstep)

C==========================================================================C
C This subroutine checks if there is enough kinetic energy to compensate   C
C the variation of potential energy due to hopping and in this case        C
C adjustes velocity to assure energy conservation                          C
C                                                                          C
C A dettailed description is found in Ref. [1]                             C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C                                                                          C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C A             See definition in THEORY                                   C
C B             See definition in THEORY                                   C
C conv          conversion factor from  Kcal/mol to  AMU*A**2/ps**2        C
C DE            Energy difference between states involved in hopping       C
C DIR           Vector along which velocity is adjusted                    C
C discr         = (B*B) + 4.D0*(DE*A)  (see THEORY)                        C
C dt            Time step [in picoseconds]         OBSOLETE!!              C
C GAMMA         Correction factor for velocity                             C
C i             Counter                                                    C
C istat         current state                                              C
C j             state to which hopping should occur                        C
C mass          Atomic masses                                              C
C N             Number of atoms                                            C
C OLD_DE        Variable to store original value of DE                     C
C OLD_DIR       Variabole to store original value of DIR                   C
C vx,vy,vz      Cartesian components of velocity                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C==========================================================================C

C      USE LIMIT, ONLY : LM1, LMGRD, LM1M
      
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      REAL*8, DIMENSION(3,N) :: DIR
      REAL*8, DIMENSION(N) :: mass
      REAL*8, DIMENSION(N) :: vx,vy,vz
      REAL*8 :: dt
      REAL*8 :: DE
      INTEGER :: j
      INTEGER :: mdstep

C INPUT/OUTPUT VARIABLES
      REAL*8 :: GAMMA
      INTEGER :: istat

C LOCAL VARIABLES
      REAL*8 :: A
      REAL*8 :: B
      REAL*8 :: discr
      INTEGER :: i
      REAL*8, DIMENSION(3,N) :: OLD_DIR
      REAL*8 :: OLD_DE

C  PARAMETERS
      REAL*8, PARAMETER :: conv = 418.4  ! Kcal/mol -> AMU*A**2/ps**2

C--------------------------------------------------------------------------C

      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,*)">>>ADJEK2"
         WRITE(*,*)
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"Velocity adjustment will be performed along"
         WRITE(*,*)"the following vector:"
         WRITE(*,'("Atom",11X,"CCx",11X,"CCy",11X,"CCz")')
         DO i = 1, N
            WRITE(*,'(1X,I3,4X,F12.6,2X,F12.6,2X,F12.6)')
     &           i,DIR(1,i),DIR(2,i),DIR(3,i)
         END DO
      END IF

C  inizialization
C  
C    Note that DE has originally dimensions of
C    kcal/mol. But we need it to be in AMU*A^2/ps^2
      OLD_DE = DE

      DE = conv*DE 

      GAMMA = 0.D0

C  computation of A, B and discr
      A = 0.D0
      B = 0.D0

      DO i = 1, N
         A = A + (DIR(1,i)**2+DIR(2,i)**2+DIR(3,i)**2)/mass(i)

         B = B + DIR(1,i)*vx(i) + DIR(2,i)*vy(i) + DIR(3,i)*vz(i) 
      END DO

      A  = A*0.5D0

      discr = (B*B) + 4.D0*(DE*A)

      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,'("discr = ",F16.5)')discr
         IF (iout .GE. 3) THEN
            WRITE(*,'("A = ",F12.5)')A
            WRITE(*,'("B = ",F12.5)')B
         END IF
      END IF
  
C  If DISCR < 0  --> no hopping and GAMMA = B/A
      IF (discr .LT. 0.D0) THEN
         GAMMA = B/A

         IF (iout .GE. 2) THEN
            WRITE(*,*)
            WRITE(*,*)"There is not enough kinetic energy to"
            WRITE(*,*)"perform hopping!"
            WRITE(*,*)"Velocity components will be reversed"
         END IF
      END IF

C  If DISCR .GE. 0  --> hopping
C     If B < 0     -->  GAMMA = (B+SQRT(DISCR))/(2.0D0*A)
C     If B .GE. 0  -->  GAMMA = (B-SQRT(DISCR))/(2.0D0*A)
      IF (discr .GE. 0.D0) THEN
           IF (iout .GE. 1) THEN
            WRITE(*,*)
            WRITE(*,'("MD step: ",I5,"  -  Hopping: ",I1," -> ",I1)')
     &               mdstep,istat,j
         END IF
       
         istat = j

         IF (B .LT. 0.D0) THEN
            GAMMA = (B+DSQRT(discr))/(2.D0*A)
         ELSE IF (B .GE. 0.D0) THEN
            GAMMA = (B-DSQRT(discr))/(2.D0*A)
         END IF

      END IF

      DE = OLD_DE

      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,*)"Correction factor for velocity
     & scaling after hopping"
         WRITE(*,'(" GAMMA = ",F12.6)')GAMMA
      END IF

      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,*)"ADJEK2<<<"
         WRITE(*,*)
      END IF

      END


C--------------------------------------------------------------------------------C

      SUBROUTINE pop_norm(iout,N,C)

C This subroutine normalizes to one the electronic population C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       
C MAIL: efabiano@mpi-muelheim.mpg.de                                       
C DATE: March 2007                                                         

      IMPLICIT NONE

      INTEGER :: iout
      INTEGER :: N
      COMPLEX*16, DIMENSION(N) :: C

      INTEGER :: i
      REAL*8 :: K
      REAL*8 :: tmp


      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>pop_norm"
         WRITE(*,*)
      END IF

      tmp = 0.D0
      DO i = 1, N
         tmp = tmp + REAL(C(i)*DCONJG(C(i)))
      END DO

      K = DSQRT(1/tmp)
      
      DO i = 1, N
         C(i) = C(i)*K
C         WRITE(*,*)"pop",i,REAL(C(i)*DCONJG(C(i)))
      END DO

      IF (iout .GE. 3) THEN
         WRITE(*,'("Scaling factor for population
     & normalization: ",F10.6)')K
         WRITE(*,*)
         WRITE(*,*)"pop_norm<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C-----------------------------------------------------------------------------C

      SUBROUTINE m_num_CC(iout,nexstat,NCICONF,CICOMP,OLD_CICOMP,dt,
     &                   mdstep,startStep,naco,de,old_de,nCC,tully_rest,
     &                   inde)


C  This subroutine compute numerically non-adiabatic coupling terms
C  using the formula (28) of Ref. [1]:
C      nCC(i,j) = (1/2dt)[<i(t)|j(t+dt)> - <i(t+dt)|j(t)>]
C
C  To do this it it necessary to map excited states of previous MD 
C  step into new ones. This is done according to an energy difference
C  criterion if possible. Otherwise comparing CI components.
C  The subroutine also gives as output the mapping of actual state used
C  MD propagation.
C
C  Reference:
C  [1] S. Hammes-Schiffer, J.C. Tully: J. Chem. Phys. 101, 4657 (1994).
C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       
C MAIL: efabiano@mpi-muelheim.mpg.de                                       
C DATE: March 2007                                                         

      USE LIMIT, ONLY: LMGRD, LMCONF

      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: IN2
      INTEGER :: iout                       ! flag for output verbosity
      INTEGER :: nexstat                    ! Number of excited states
      INTEGER :: NCICONF                    ! Number of CI configurations
      REAL*8, DIMENSION(LMCONF,6) :: CICOMP   ! CI components   
      REAL*8 :: dt                          ! time step
      INTEGER :: mdstep                       ! actual MD step
      INTEGER :: startStep                       ! first MD step
      INTEGER :: naco                       ! number of active CI orbitals
      INTEGER, DIMENSION(40,LMCONF) :: de     ! lexical indexes
      LOGICAL :: tully_rest                 ! flag to perform tully restart
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP  ! old CI components
      INTEGER, DIMENSION(40,LMCONF) :: old_de    ! old lexical indexes
      INTEGER, DIMENSION(nexstat) :: inde      ! mapping of states

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(nexstat,nexstat) :: nCC

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./INOPT2/ IN2(300)
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      INTEGER, DIMENSION(NCICONF) :: index
      INTEGER :: i,j,l

      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir

C----------------------C

C TIMING - overall numerical CC time, start
      IF (iout .GE. 2) THEN
         CALL CPU_TIME(ctime1)
         CALL SYSTEM_CLOCK(wtime1,check_ir)
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)">>>m_num_CC"
         WRITE(*,*)
      END IF

C  INITIALIZATION
      nCC = 0.D0

C  In the first MD step if not tully_restart exit
      IF (mdstep-startStep .LE. 1) THEN
         IF (.NOT. tully_rest) THEN
            RETURN        
         END IF
      END IF   

C     compute index for tracking of MOs mapping
      IF(IN2(77).LT.6) THEN
        CALL mindd(iout,NCICONF,naco,de,old_de,index,CICOMP,old_CICOMP,
     &           nexstat)
      ELSE
        DO I=1,NCICONF
          index(I) = I
        ENDDO
      ENDIF

C  Compute numerical non-adiabatic coupling coefficients
      DO i = 1, nexstat
         Do j = 1, nexstat
            DO l = 1, NCICONF
               nCC(i,j) = nCC(i,j) + old_CICOMP(index(l),i)*
     &                    CICOMP(l,inde(j)) - old_CICOMP(index(l),j)*
     &                    CICOMP(l,inde(i))
            END DO
         END DO
      END DO

      nCC = nCC/(2.D0*dt)

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("Numerical non-adiabatic coupling elements -
     & mdstep: ",I5)')mdstep
         WRITE(*,'(1X,10I15)')(j,j=2,nexstat)
         DO i = 1, nexstat-1
            WRITE(*,'(I1,4X,10F15.6)')i,(nCC(i,j),j=i+1,nexstat)
         END DO
      END IF

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"m_num_CC<<<"
         WRITE(*,*)
      END IF

      IF (iout .GE. 2) THEN
C TIMING - numerical CC, end
         CALL CPU_TIME(ctime2)
         CALL SYSTEM_CLOCK(wtime2,check_ir)
         wtime_ncc = wtime_ncc + wtime2-wtime1
            IF (ctime2 .GT. REAL(wtime2-wtime1)/REAL(check_ir)+ctime1) 
     &          ctime2 = REAL(wtime2-wtime1)/REAL(check_ir)+ctime1
         ctime_ncc = ctime_ncc + ctime2-ctime1
      END IF

      IF (iout .GE. 3) THEN
C TIMING - output
            WRITE(*,*)
            WRITE(*,'("-----------------------------------")')
            WRITE(*,*)"*** TIMING ***"
            WRITE(*,*)"Numerical nonadiabatic couplings"
            WRITE(*,*)
            CALL out_time(ctime1,ctime2,wtime1,wtime2,check_ir)
            WRITE(*,'("-----------------------------------")')
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_euler(iout,N,Ngrd,H,Hold,rd,rdold,dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements euler algorithm for the                       C
C solution of Eq. (26) of REf. [1] and at the same time it performs the    C
C integral in Eq (29) of Ref. [1].                                         C
C                                                                          C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C                                                                          C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C istat         Current state                                              C
C j             Counter                                                    C
C l             Counter                                                    C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C qq            Scratch variable used in Euler method                      C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 previous time step                                       C
C rh            Scratch variable                                           C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh      
      COMPLEX*16, DIMENSION(Ngrd) :: qq

C--------------------------------------------------------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_euler>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

      d = dt/REAL(NE)

      g = 0.D0

      DO i = 0, NE-1
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i)/REAL(NE)
            cb = C(j)*DCONJG(C(istat))*tmprd
            b = -2.D0*REAL(cb)

            g(j) = g(j) + b*d 
         END DO

         DO j = 1, Ngrd
            qq(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
               qq(j) = qq(j) + hh*C(l)
            END DO
         END DO

         DO j = 1, Ngrd
            C(j) = C(j) + qq(j)*d
C  Change C with tmp_C, comment above line and use the following line 
C  if you want just to print C but not to update the value 
C  (print tmp_C in this case)
C             tmp_C(j) = tmp_C(j) + qq(j)*d
         END DO
      END DO

      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_euler<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------C

      SUBROUTINE int_RK4(iout,N,Ngrd,H,Hold,rd,rdold,dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements fourth-order Runge-Kutta algorithm for the    C
C solution of Eq. (26) of REf. [1] and at the same time it performs the    C
C integral in Eq (29) of Ref. [1].                                         C
C For Runge-Kutta algorithm see for example Ref. [2]                       C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] "Numerical Recipes in FORTRAN 77: The Art of Scientific Computing",  C
C      by William H. Press, Brian P. Flannery, Saul A. Teukolsky,          C
C      and William T. Vetterling                                           C
C     also: http://www.library.cornell.edu/nr/bookfpdf/f16-1.pdf           C
C                                                                          C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C iout          Flag for output verbosity                                  C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k1,k2,k3,k4   Intermediates used in Runge-Kutta method                   C
C l             Counter                                                    C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C==========================================================================C

      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16, DIMENSION(Ngrd) :: k1,k2,k3,k4
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh
      
      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C--------------------------------------------------------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_RK4>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C   Initializatio of variables
      d = dt/REAL(NE)

      g = 0.D0

C   TIMING - C integration time, start   
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C   Start cycle over subintervals
      DO i = 0, NE-1
         old_C = C

C        Computation of quantum amplitudes
         DO j = 1, Ngrd
            k1(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
               k1(j) = k1(j) + hh*c(l)
            END DO
            k1(j) = k1(j)*d
         END DO

         DO j = 1, Ngrd
            k2(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+0.5D0)/REAL(NE)
               k2(j) = k2(j) + hh*(c(l)+0.5D0*k1(l))
            END DO
            k2(j) = k2(j)*d
         END DO

         DO j = 1, Ngrd
            k3(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+0.5D0)/REAL(NE)
               k3(j) = k3(j) + hh*(c(l)+0.5D0*k2(l))
            END DO
            k3(j) = k3(j)*d
         END DO

         DO j = 1, Ngrd
            k4(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+1.D0)/REAL(NE)
               k4(j) = k4(j) + hh*(c(l)+k3(l))
            END DO
            k4(j) = k4(j)*d
         END DO

         DO j = 1, Ngrd
            C(j) = C(j) + (k1(j)/6.D0) + (k2(j)/3.D0) + 
     &                    (k3(j)/3.D0) + (k4(j)/6.D0)
         END DO

C   TIMING - probability integration time, start 
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF    

C      Computation of probabilities         
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO

C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over subintervals  
      END DO

C   TIMING - C integration time, end
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C   TIMING - output
         IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for RK4 integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")')
      END IF
     
C   Normalization of probabilities
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_RK4<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE nran2(z,idum0,start)

C This subroutine is a minimal random number generator of
C Park and Miller plus Bays-Durham shuffle [1].
C First call must be used to initialize idum, for example
C       CALL nran1(z,idum,.TRUE.)
C with idum an integer [suggested use of CALL system_clock(count=idum)]
C Successive call must be made with start = .FALSE. [idum is automatically
C updated and must be not changed between two succesive calls]
C With start .FALSE. the value of idum0 is unimportant.
C The subroutine returns a random number between 0 and 1.
C
C REFERENCE
C [1] Numerical recipes in fortran
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       
C MAIL: efabiano@mpi-muelheim.mpg.de                                       
C DATE: March 2007                                                        

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: start ! flag for initialization of idum
      INTEGER, INTENT(IN) :: idum0 ! seed for random number generation

      REAL*8, INTENT(OUT) :: z  ! random number

      INTEGER, PARAMETER :: a = 48271
      INTEGER, PARAMETER :: m = 2147483647
      INTEGER, PARAMETER :: q = 44488
      INTEGER, PARAMETER :: r = 3399
      INTEGER, PARAMETER :: mask = 123459876  
      INTEGER, PARAMETER :: ntab = 32
      INTEGER, PARAMETER :: ndiv = 1+(m-1)/ntab
      REAL*8, PARAMETER :: invm = 1.D0/m

      INTEGER :: idum
      INTEGER :: i
      INTEGER :: k
      INTEGER, DIMENSION(ntab) :: sx
      INTEGER :: y
      INTEGER :: j

C..............C

      SAVE idum
      SAVE sx
      SAVE y

      IF (start) THEN
         idum  = IEOR(idum0,mask)

         DO i = ntab+10, 1, -1
            k = idum/q
            idum = a*(idum-(k*q))-(r*k)
            IF (idum .LT. 0) idum = idum + m
            
            IF (i .LE. ntab) sx(i) = idum

            IF (i .EQ. ntab+10) idum = IEOR(idum,mask)
         END DO

         y = sx(1)
      END IF

      k = idum/q
      idum = a*(idum-(k*q))-(r*k)
      IF (idum .LT. 0) idum = idum + m

      j = 1 + (y/ndiv)


      y = sx(j)

      sx(j) = idum

      z = y*invm

      END SUBROUTINE nran2



C---------------------------------------------------------------------!

      SUBROUTINE nran3(z,idum0,start)
      
C This subroutine random number generator based on subtractive method
C of Knuth [1,2].
C First call must be used to initialize idum, for example
C       CALL nran3(z,idum,.TRUE.)
C with idum an integer [suggested use of CALL system_clock(count=idum)]
C Successive call must be made with start = .FALSE. [idum is automatically
C updated and must be not changed between two succesive calls]
C With start .FALSE. the value of idum0 is unimportant.
C The subroutine returns a random number between 0 and 1.
C
C REFERENCE
C [1] D.E. Knuth, Seminumerical Algorithms in The Art of Computer Programming,
C         Addison-Wesely, 1981.
C [2] Numerical recipes in fortran
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       
C MAIL: efabiano@mpi-muelheim.mpg.de                                       
C DATE: March 2007                                                        

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idum0
      LOGICAL, INTENT(IN) :: start

      REAL*8, INTENT(OUT) :: z

      INTEGER :: mj
      INTEGER, DIMENSION(55) :: ma
      INTEGER :: mk
      INTEGER :: i
      INTEGER :: ii
      INTEGER :: inext
      INTEGER :: inextp
      INTEGER :: idum

      INTEGER, PARAMETER :: mseed = 161803398
      INTEGER, PARAMETER :: mbig = 1000000000
      REAL*8, PARAMETER :: fac = 1.D0/mbig

      SAVE inext
      SAVE inextp
      SAVE ma
      SAVE idum

      IF (start) THEN
         mj = ABS(mseed-ABS(idum0))
         mj = MOD(mj,mbig)
         ma(55) = mj
         mk = 1

         DO i = 1, 54
            ii = MOD(21*i,55)
            ma(ii) = mk
            mk = mj-mk
            IF (mk .LT. 0) mk = mk+mbig
            mj = ma(ii)
         END DO
         
         DO ii = 1, 4
            DO i = 1, 55
               ma(i) = ma(i) - ma(1+MOD(i+30,55))
               IF (ma(i) .LT. 0) ma(i) = ma(i)+mbig
            END DO
         END DO

         inext = 0
         inextp = 31
         idum = 1
      END IF

      inext = inext + 1
      IF (inext .EQ. 56) inext = 1

      inextp = inextp + 1
      IF (inextp .EQ. 56) inextp = 1
      
      mj = ma(inext) - ma(inextp)
      IF (mj .LT. 0) mj = mj + mbig

      ma(inext) = mj

      z = mj*fac

      END SUBROUTINE nran3



C--------------------------------------------------------------------------------C

      SUBROUTINE int_ABM4(iout,mdstep,startStep,N,Ngrd,H,Hold,rd,rdold,
     &                    dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements fourth order Adams-Bashforth-Moulton          C
C predictor-corrector algorithm for the solution of Eq. (26) of REf. [1]   C
C and at the same time it performs the integral in Eq (29) of Ref. [1].    C
C For dettails on the algorithm see for example Ref. [2]                   C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] ???????                                                              C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C erro          Error in the computation of c in the actual step           C
C f             H*c                                                        C
C f1            H*c at previous time step                                  C
C f2            H*c two time steps before                                  C
C f3            H*c three time steps before                                C
C f_pred        H*c predicted in next step                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C iout          Flag for output verbosity                                  C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k1,k2,k3,k4   Intermediates used in Runge-Kutta method                   C
C l             Counter                                                    C
C max_err       Maximum error in the computation of c                      C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C pred          value of c from predictor formula                          C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmp_c         Scratch variable                                           C
C tmp_p         Scratch variable                                           C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: mdstep,startStep
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16, DIMENSION(Ngrd) :: k1,k2,k3,k4
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh
      
      COMPLEX*16, DIMENSION(4) :: f1,f2,f3 ! CHANGE THE DIMENSION IF NEEDED
      COMPLEX*16, DIMENSION(Ngrd) :: f
      COMPLEX*16, DIMENSION(Ngrd) :: pred
      COMPLEX*16, DIMENSION(Ngrd) :: corr
      COMPLEX*16, DIMENSION(Ngrd) :: f_pred
      COMPLEX*16, DIMENSION(Ngrd) :: tmp_p, tmp_c
      REAL*8 :: erro
      REAL*8 :: max_err
      INTEGER :: ite
      REAL*8 :: max_err2
      COMPLEX*16, DIMENSION(Ngrd) :: old_corr
      REAL*8 :: err_tol

      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C--------------------------------------------------------------------------C

      SAVE f3
      SAVE f2
      SAVE f1

C----------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM4>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C Initialitation of time interval

      d = dt/REAL(NE)

C Initialization of variables
      g = 0.D0
      max_err2 = 0.D0

C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C Starting cycle over integration steps
      DO i = 0, NE-1
         old_C = C

C       Computation of quantum amplitudes
         IF (mdstep-startStep .LT. 6) THEN 
C           For the first steps use RK4 to start
            DO j = 1, Ngrd
               k1(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  k1(j) = k1(j) + hh*c(l)
               END DO
               k1(j) = k1(j)*d
            END DO

            DO j = 1, Ngrd
               k2(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                 (REAL(i)+0.5D0)/REAL(NE)
                  k2(j) = k2(j) + hh*(c(l)+0.5D0*k1(l))
               END DO
               k2(j) = k2(j)*d
            END DO
           
            DO j = 1, Ngrd
               k3(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                 (REAL(i)+0.5D0)/REAL(NE)
                  k3(j) = k3(j) + hh*(c(l)+0.5D0*k2(l))
               END DO
               k3(j) = k3(j)*d
            END DO

            DO j = 1, Ngrd
               k4(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                 (REAL(i)+1.D0)/REAL(NE)
                  k4(j) = k4(j) + hh*(c(l)+k3(l))
               END DO
               k4(j) = k4(j)*d
            END DO
            
            DO j = 1, Ngrd
               C(j) = C(j) + (k1(j)/6.D0) + (k2(j)/3.D0) + 
     &              (k3(j)/3.D0) + (k4(j)/6.D0)
            END DO

C           Update stored values of H*C (only in last steps)
            IF (mdstep-startStep .GE. 2) THEN
               IF (i .GE. (NE-6)) THEN! 4
                  f3 = f2
                  f2 = f1
                  DO j = 1, Ngrd
                     f1(j) = DCMPLX(0.D0,0.D0)
                     DO l = 1, Ngrd
                        hh = Hold(j,l)+(H(j,l)-Hold(j,l))*
     &                                          REAL(i)/REAL(NE)
                        f1(j) = f1(j) + hh*C(l)
                     END DO
                  END DO
               END IF
            END IF

C        In successive steps use ABM
         ELSE
            DO j = 1, Ngrd
               f(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  f(j) = f(j) + hh*C(l)
               END DO
            END DO

            max_err = 0.D0
            DO j = 1, Ngrd
               tmp_p(j) = (55.D0*f(j)) - (59.D0*f1(j)) + (37.D0*f2(j)) 
     &                    - (9.D0*f3(j)) 
               pred(j) = C(j) + (d/24.D0)*tmp_p(j)
C               C(j) = pred(j) ! uncomment this and comment everything below
C                              ! to use only Adams-Bashforth
               f_pred(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i+1)/REAL(NE)
                  f_pred(j) = f_pred(j) + hh*pred(l)
               END DO
               tmp_c(j) = (9.D0*f_pred(j)) + (19.D0*f(j)) - (5.D0*f1(j))
     &                     + f2(j)
               corr(j) = c(j) +  (d/24.D0)*tmp_c(j)
               erro = SQRT((REAL(corr(j))-REAL(pred(j)))**2 + 
     &                  (AIMAG(corr(j))-AIMAG(pred(j)))**2) 
               erro = 0.07037037D0*erro ! =(19/270)*erro
               IF (erro .GE. max_err) max_err = erro
            END DO

C         If needed compute iterative solution to
c         required accuracy (err_tol)
            err_tol = 1.D-10

            IF (max_err .GE. err_tol) THEN
               DO ite = 1, 20
                  max_err = 0.D0
                  DO j = 1, Ngrd
                     old_corr(j) = corr(j)
                     f_pred(j) = DCMPLX(0.D0,0.D0)
                     DO l = 1, Ngrd
                        hh = Hold(j,l) + (H(j,l)-Hold(j,l))
     &                                      *REAL(i+1)/REAL(NE)
                        f_pred(j) = f_pred(j) + hh*corr(l)
                     END DO
                     tmp_c(j) = (9.D0*f_pred(j)) + (19.D0*f(j))
     &                            - (5.D0*f1(j)) + f2(j)
                     corr(j) = c(j) +  (d/24.D0)*tmp_c(j)
                     erro = SQRT((REAL(corr(j))-REAL(old_corr(j)))**2 + 
     &                           (AIMAG(corr(j))-AIMAG(old_corr(j)))**2)
                     IF (erro .GE. max_err) max_err = erro
                  END DO
                  IF (max_err .LE. err_tol) EXIT
               END DO
            END IF

            C = corr

            IF (max_err .GE. max_err2) max_err2 = max_err
            

            f3 = f2
            f2 = f1
            f1 = f
         END IF

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF

C        Computation of probabilities
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO

C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over integration steps   
      END DO

C   TIMING - C integration time, end
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1/REAL(ir1))+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C  Output - Write integration error     
      IF (iout.GE.3 .AND. mdstep-startStep.GE.5) THEN
         WRITE(*,*)
         WRITE(*,'("Maximum estimated error in integration: ",E12.4)')
     &             max_err2
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for ABM4 integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")') 
      END IF
     
C  Normalize probabilities 
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM4<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_ABM2(iout,mdstep,startStep,N,Ngrd,H,Hold,rd,rdold,
     &                    dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements second order Adams-Bashforth-Moulton          C
C predictor-corrector algorithm for the solution of Eq. (26) of REf. [1]   C
C and at the same time it performs the integral in Eq (29) of Ref. [1].    C
C For dettails on the algorithm see for example Ref. [2]                   C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] ???????                                                              C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C f             H*c                                                        C
C f1            H*c at previous time step                                  C
C f_pred        H*c predicted at next time step                            C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C iout          Flag for output verbosity                                  C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k1,k2         Intermediates used in Runge-Kutta method                   C
C l             Counter                                                    C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C pred          value of c from predictor formula                          C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmp_c         Scratch variable                                           C
C tmp_p         Scratch variable                                           C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: mdstep,startStep
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16, DIMENSION(Ngrd) :: k1,k2
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh
      
      COMPLEX*16, DIMENSION(4) :: f1 ! CHANGE THE DIMENSION IF NEEDED
      COMPLEX*16, DIMENSION(Ngrd) :: f
      COMPLEX*16, DIMENSION(Ngrd) :: pred
      COMPLEX*16, DIMENSION(Ngrd) :: f_pred
      COMPLEX*16, DIMENSION(Ngrd) :: tmp_p, tmp_c

      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C--------------------------------------------------------------------------C

      SAVE f1

C----------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM2>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C Initialitation of time interval
      d = dt/REAL(NE)

C Initialization of variables
      g = 0.D0
     
C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C Starting cycle over integration steps
      DO i = 0, NE-1
         old_C = C

C       Computation of quantum amplitudes
         IF (mdstep-startStep .LT. 6) THEN 
C           For the first steps use RK2 to start
            DO j = 1, Ngrd
               k1(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  k1(j) = k1(j) + hh*c(l)
               END DO
               k1(j) = k1(j)*d
            END DO

            DO j = 1, Ngrd
               k2(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                 (REAL(i)+0.5D0)/REAL(NE)
                  k2(j) = k2(j) + hh*(c(l)+0.5D0*k1(l))
               END DO
               k2(j) = k2(j)*d
            END DO
            
            DO j = 1, Ngrd
               C(j) = C(j) + k2(j)
            END DO

C           Update stored values of H*C (only in last steps)
            IF (i.GE.(NE-6) .AND. mdstep-startStep .GE. 3) THEN 
               DO j = 1, Ngrd
                  f1(j) = DCMPLX(0.D0,0.D0)
                  DO l = 1, Ngrd
                     hh = Hold(j,l)+(H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                     f1(j) = f1(j) + hh*c(l)
                  END DO
               END DO
            END IF

C        In successive steps use ABM
         ELSE
            DO j = 1, Ngrd
               f(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  f(j) = f(j) + hh*c(l)
               END DO
            END DO

            DO j = 1, Ngrd
               tmp_p(j) = (3.D0*f(j)) - f1(j)
               pred(j) = c(j) + (d/2.D0)*tmp_p(j)
C               c(j) = pred(j) ! uncomment this and comment everything below
C                              ! to use only Adams-Bashforth
               f_pred(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i+1)/REAL(NE)
                  f_pred(j) = f_pred(j) + hh*pred(l)
               END DO
               tmp_c(j) = f_pred(j) + f(j)
               c(j) = c(j) +  (d/2.D0)*tmp_c(j)
            END DO

            f1 = f
         END IF

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF
         
C        Computation of probabilities
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO

C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over integration steps   
      END DO

C   TIMING - C integration time, end
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for ABM2 integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")')
      END IF

C  Normalize probabilities 
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM2<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_ABM5(iout,mdstep,startStep,N,Ngrd,H,Hold,rd,rdold,
     &                    dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements fifth order Adams-Bashforth-Moulton           C
C predictor-corrector algorithm for the solution of Eq. (26) of REf. [1]   C
C and at the same time it performs the integral in Eq (29) of Ref. [1].    C
C For dettails on the algorithm see for example Ref. [2]                   C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] ???????                                                              C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C erro          Error in the computation of c in the actual step           C
C f             H*c                                                        C
C f1            H*c at previous time step                                  C
C f2            H*c two time steps before                                  C
C f3            H*c three time steps before                                C
C f4            H*c four time steps before                                 C
C f_pred        H*c predicted in next step                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C iout          Flag for output verbosity                                  C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k1,k2,k3,k4   Intermediates used in Runge-Kutta method                   C
C l             Counter                                                    C
C max_err       Maximum error in the computation of c                      C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C pred          value of c from predictor formula                          C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmp_c         Scratch variable                                           C
C tmp_p         Scratch variable                                           C
C tmprd         Scratch variable                                           C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: mdstep,startStep
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      COMPLEX*16 :: hh
      INTEGER :: i,j,l
      COMPLEX*16, DIMENSION(Ngrd) :: k1,k2,k3,k4
      COMPLEX*16 :: cb
      REAL*8 :: b
      REAL*8 :: tmprd
      COMPLEX*16 :: rh
      
      COMPLEX*16, DIMENSION(4) :: f1,f2,f3,f4 ! CHANGE THE DIMENSION IF NEEDED
      COMPLEX*16, DIMENSION(Ngrd) :: f
      COMPLEX*16, DIMENSION(Ngrd) :: pred
      COMPLEX*16, DIMENSION(Ngrd) :: corr
      COMPLEX*16, DIMENSION(Ngrd) :: f_pred
      COMPLEX*16, DIMENSION(Ngrd) :: tmp_p, tmp_c
      REAL*8 :: erro
      REAL*8 :: max_err
      INTEGER :: ite
      REAL*8 :: max_err2
      COMPLEX*16, DIMENSION(Ngrd) :: old_corr
      REAL*8 :: err_tol

      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C--------------------------------------------------------------------------C

      SAVE f4
      SAVE f3
      SAVE f2
      SAVE f1

C----------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM5>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C Initialitation of time interval
      d = dt/REAL(NE)

C Initialization of variables
      g = 0.D0
      max_err2 = 0.D0

C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C Starting cycle over integration steps
      DO i = 0, NE-1
         old_C = C

C       Computation of quantum amplitudes
         IF (mdstep-startStep .LT. 6) THEN
C           For the first steps use RK4 to start
            DO j = 1, Ngrd
               k1(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  k1(j) = k1(j) + hh*c(l)
               END DO
               k1(j) = k1(j)*d
            END DO

            DO j = 1, Ngrd
               k2(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+0.5D0)/REAL(NE)
                  k2(j) = k2(j) + hh*(c(l)+0.5D0*k1(l))
               END DO
               k2(j) = k2(j)*d
            END DO

            DO j = 1, Ngrd
               k3(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+0.5D0)/REAL(NE)
                  k3(j) = k3(j) + hh*(c(l)+0.5D0*k2(l))
               END DO
               k3(j) = k3(j)*d
            END DO

            DO j = 1, Ngrd
               k4(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*
     &                     (REAL(i)+1.D0)/REAL(NE)
                  k4(j) = k4(j) + hh*(c(l)+k3(l))
               END DO
               k4(j) = k4(j)*d
            END DO

            DO j = 1, Ngrd
               C(j) = C(j) + (k1(j)/6.D0) + (k2(j)/3.D0) + 
     &                    (k3(j)/3.D0) + (k4(j)/6.D0)
            END DO

C           Update stored values of H*C (only in last steps)
            IF (mdstep-startStep .GE. 2) THEN
               IF (i .GE. (NE-6)) THEN
                  f4 = f3
                  f3 = f2
                  f2 = f1
                  DO j = 1, Ngrd
                     f1(j) = DCMPLX(0.D0,0.D0)
                     DO l = 1, Ngrd
                        hh = Hold(j,l)+(H(j,l)-Hold(j,l))*
     &                                          REAL(i)/REAL(NE)
                        f1(j) = f1(j) + hh*c(l)
                     END DO
                  END DO
               END IF
            END IF

C        In successive steps use ABM
         ELSE
            DO j = 1, Ngrd
               f(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i)/REAL(NE)
                  f(j) = f(j) + hh*c(l)
               END DO
            END DO

            max_err = 0.D0
            DO j = 1, Ngrd
               tmp_p(j) = (1901.D0*f(j)) - (2774.D0*f1(j)) +
     &                    (2616.D0*f2(j)) - (1274.D0*f3(j)) +
     &                    (251.D0*f4(j))
               pred(j) = c(j) + (d/720.D0)*tmp_p(j)
C               c(j) = pred(j) ! uncomment this and comment everything below
C                              ! to use only Adams-Bashforth
               f_pred(j) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i+1)/REAL(NE)
                  f_pred(j) = f_pred(j) + hh*pred(l)
               END DO
               tmp_c(j) = (251.D0*f_pred(j)) + (646.D0*f(j)) -
     &                    (264.D0*f1(j)) + (106.D0*f2(j)) -
     &                    (19.D0*f3(j))
               corr(j) = c(j) +  (d/720.D0)*tmp_c(j)

               erro = SQRT((REAL(c(j))-REAL(pred(j)))**2 + 
     &                  (AIMAG(c(j))-AIMAG(pred(j)))**2) 
               erro = 0.060267857D0*erro ! =(27/448)*erro
               IF (erro .GE. max_err) max_err = erro
            END DO

C         If needed compute iterative solution to
C         required accuracy (err_tol)
            err_tol = 1.D-10

            IF (max_err .GE. err_tol) THEN
               DO ite = 1, 20
                  max_err = 0.D0
                  DO j = 1, Ngrd
                     old_corr(j) = corr(j)
                     f_pred(j) = DCMPLX(0.D0,0.D0)
                     DO l = 1, Ngrd
                        hh = Hold(j,l) + (H(j,l)-Hold(j,l))
     &                                     *REAL(i+1)/REAL(NE)
                        f_pred(j) = f_pred(j) + hh*corr(l)
                     END DO
                     tmp_c(j) = (251.D0*f_pred(j)) + (646.D0*f(j)) -
     &                    (264.D0*f1(j)) + (106.D0*f2(j)) -
     &                    (19.D0*f3(j))
                     corr(j) = c(j) +  (d/720.D0)*tmp_c(j)
                     erro = SQRT((REAL(corr(j))-REAL(old_corr(j)))**2 + 
     &                    (AIMAG(corr(j))-AIMAG(old_corr(j)))**2) 
                     IF (erro .GE. max_err) max_err = erro
                  END DO
                  IF (max_err .LE. err_tol) EXIT
               END DO
            END IF

            C = corr

            IF (max_err .GE. max_err2) max_err2 = max_err

            f4 = f3
            f3 = f2
            f2 = f1
            f1 = f
         END IF

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF
 
C        Computation of probabilities
C        act_C = old_C + (C-old_C)*REAL(i+0.5)/REAL(NE)
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO
  
C   TIMING - probability integration time, end       
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over integration steps   
      END DO

C   TIMING - C integration time, end 
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C  Output - Write integration error     
      IF (iout.GE.3 .AND. mdstep-startStep.GE.5) THEN
         WRITE(*,*)
         WRITE(*,'("Maximum estimated error in integration: ",E12.4)')
     &             max_err2
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for ABM5 integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")') 
      END IF

C  Normalize probabilities 
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_ABM5<<<",istat
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_imp(iout,N,Ngrd,H,Hold,rd,rdold,dt,NE,istat,C,g)

C==========================================================================C
C This subroutine implements implicit midle-point algorithm for the        C
C solution of Eq. (26) of REf. [1]                                         C
C At the same time it performs the integral in Eq (29) of Ref. [1].        C
C For dettails on the algorithm see for example Ref. [2]                   C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] ???????                                                              C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cc            Scratch variable                                           C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C hh            Scratch variable                                           C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C info          Scratch variable for Lapack control                        C
C iout          Flag for output verbosity                                  C
C ipiv          Scratch array for Lapack                                   C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k             Counter                                                    C
C l             Counter                                                    C
C m1            Matrix (I-(d/2)*H(t_n+1)) and then its inverse             C
C m11           Scratch matrix                                             C
C m2            Matrix (I+(d/2)*H(t_n))                                    C
C m3            Scratch matrix                                             C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Values of rd at the end of the previous call               C
C rh            Scratch variable                                           C
C sc1,2,3,4     Scratch variables for timing                               C
C tmprd         Scratch variable                                           C
C work          Scratch variable for Lapack                                C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      REAL*8 :: d
      INTEGER :: i,j,l,k
      REAL*8 :: tmprd
      COMPLEX*16 :: cb
      REAL*8 :: b
      COMPLEX*16 :: hh
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: m1,m2,m3
      INTEGER, DIMENSION(Ngrd) :: ipiv
      COMPLEX*16, DIMENSION(Ngrd*2) :: work
      INTEGER :: info
      COMPLEX*16 :: rh
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: m11
      COMPLEX*16, DIMENSION(Ngrd) :: cc

      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C------------------------------------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_imp>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
      END IF

C Initialitation of time interval
      d = dt/REAL(NE)

C Initialization of variables
       g = 0.D0

C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C Starting cycle over integration steps
      DO i = 0, NE-1
         old_C = C

C       Computation of quantum amplitudes
         DO j = 1, Ngrd
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i+0.5)/REAL(NE)
               m2(j,l) = hh*d*.5D0
            END DO
            m2(j,j) = m2(j,j) + DCMPLX(1.D0,0.D0)
         END DO


         DO j = 1, Ngrd
            DO l = 1, Ngrd
               hh = Hold(j,l) + (H(j,l)-Hold(j,l))*REAL(i+0.5)/REAL(NE)
               m1(j,l) = -hh*d*.5D0
            END DO
            m1(j,j) = m1(j,j) + DCMPLX(1.D0,0.D0)
         END DO

C        Matrix inversion
         CALL ZGETRF(Ngrd,Ngrd,m1,Ngrd,ipiv,info)
         IF (info .NE. 0) THEN
            WRITE(*,*)"LU factorization failed"
            STOP 'Abnormal termination in int_imp'
         END IF

         CALL ZGETRI(Ngrd, m1, Ngrd, ipiv, work, 2*Ngrd, info)
         IF (info .NE. 0) THEN
            WRITE(*,*)"Inversion failed"
            STOP 'Abnormal termination in int_imp'
         END IF

C       Matrix multiplication (1-dH)^-1*(1+dH)
         DO j = 1, Ngrd
            DO k = 1, Ngrd
               m3(j,k) = DCMPLX(0.D0,0.D0)
               DO l = 1, Ngrd
                  m3(j,k) = m3(j,k) + m1(j,l)*m2(l,k)
               END DO
            END DO
         END DO

C       Matrix-vector multiplication [(1-dH)^-1*(1+dH)]*C
         DO j = 1, Ngrd
            cc(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               cc(j) = cc(j) + m3(j,l)*C(l)
            END DO
         END DO
         C = cc

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF

C        Computation of probabilities
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
               tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
               cb = act_C(j)*DCONJG(act_C(istat))*tmprd
               b = -2.D0*REAL(cb)
               g(j) = g(j) + b*d 
         END DO

C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over integration steps
      END DO

C   TIMING - C integration time, end  
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,*)"Time for imp integration"
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")')   
      END IF

C  Normalize probabilities  
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_imp<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



C--------------------------------------------------------------------------------C

      SUBROUTINE int_UP(iout,N,Ngrd,H,Hold,rd,rdold,dt,NE,istat,C,g,
     &                    type,vx,vy,vz,mass,ENERGX,dec_cor)

C     &                    istat,C,g,type)
C  changed for decoherence correction

C==========================================================================C
C This subroutine implements the unitary propagator algorithm for the      C
C solution of Eq. (26) of REf. [1].                                        C
C Two possible choiches can be made for the time at which the  propagator  C
C is computed: middle-point (type = 1) and initial time (type = 2).        C
C At the same time it performs the integral in Eq (29) of Ref. [1].        C
C For dettails on the algorithm see for example Ref. [2]                   C
C                                                                          C
C REFERENCES                                                               C
C [1] S. Hammes-Schiffer, J.C. Tully, J Chem Phys 101, 4657 (1994)         C
C [2] ???????                                                              C
C--------------------------------------------------------------------------C
C                                                                          C
C  VARIABLES  (in alphabetical order)                                      C
C                                                                          C
C act_C         Scratch variable                                           C
C b             function to be integrated in Eq. (29) of Ref. [1]          C
C C             Population                                                 C
C cc            Scratch array                                              C
C cb            Scratch variable                                           C
C d             Time increment used in integration                         C
C dt            Time step [in picoseconds]                                 C
C g             Transition probabilities                                   C
C H             Effective Hamiltonian at current time step                 C
C HH            Effective Hermitian Hamiltonian at the chosen time         C
C                after diagonalization: eigenvectors of HH                 C
C Hold          Effective Hamiltonian at previous time step                C
C i             Counter                                                    C
C IH            (-i)*H this makes the matrix Hermitian (for diagonalization)
C IHold         (-i)*Hold this makes the matrix Hermitian                  C
C info          Variable to check diagonalizzation exit status             C
C inv_HH        Inverse of the matrix of eigenvectors of HH                C
C                Actually is just the transpose complex conjugate          C
C iout          Output level                                               C
C ir1,2,3,4     Scratch variables for timing                               C
C istat         Current state                                              C
C j             Counter                                                    C
C k             Counter                                                    C
C l             Counter                                                    C
C LWORK         Scratch variable for Lapack                                C
C m1            Diagonal matrix of exponentials of eigenvalues             C
C m2            product of eigenvectors matrix times m1                    C
C m3            Product of m2 times inv_HH                                 C
C N             Number of atoms                                            C
C NE            Number of steps for integration                            C
C Ngrd          Number of states                                           C
C old_C         Quantum amplitudes at the previous step                    C
C rd            Dot product of velocity times nonadiabatic vectors at      C
C                 current time step                                        C
C rdold         Dot product of velocity times nonadiabatic vectors at      C
C                 previous time step                                       C
C rh            Scratch variable                                           C
C RWORK         Scratch array for Lapack                                   C
C sc1,2,3,4     Scratch variables for timing                               C
C tmprd         Scratch variable                                           C
C type          Flag to chose at which time to compute the propagator:     C
C                 = 1   computed at [t_n + t_(n+1)]/2                      C
C                 = 2   computed at t_n                                    C
C w             Eigenvalues of HH                                          C
C WORK          Scratch array for Lapack                                   C
C dec_cor       empirical decoeherence correction by Persico & Granucci    C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: June 2007                                                          C
C==========================================================================C
 
      IMPLICIT NONE

C  INPUT VARIABLES
      INTEGER :: iout
      INTEGER :: N
      INTEGER :: Ngrd
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: H
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: Hold
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rd
      REAL*8, DIMENSION(Ngrd,Ngrd) :: rdold
      REAL*8 :: dt
      INTEGER :: NE
      INTEGER :: istat
      INTEGER :: type
      REAL*8  :: dec_cor

C  INPUT/OUTPUT VARIABLES
      COMPLEX*16, DIMENSION(Ngrd) :: C

C  OUTPUT VARIABLES
      REAL*8, DIMENSION(Ngrd) :: g

C  COMMON
      REAL*8 :: ctime_eh
      INTEGER :: wtime_eh
      REAL*8 :: ctime_eint
      INTEGER :: wtime_eint
      REAL*8 :: ctime_ncc
      INTEGER :: wtime_ncc

      COMMON
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc

C  LOCAL VARIABLES
      INTEGER :: i,j,l,k
      REAL*8 :: d
      REAL*8 :: tmprd
      COMPLEX*16 :: cb
      REAL*8 :: b
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: IH
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: IHold
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: HH
      REAL*8, DIMENSION(Ngrd) :: w
      COMPLEX*16, DIMENSION(2*Ngrd) :: WORK
      INTEGER :: LWORK
      REAL*8, DIMENSION(3*Ngrd-2) :: RWORK
      INTEGER :: info
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: inv_HH
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: m1,m2,m3
      COMPLEX*16, DIMENSION(Ngrd) :: cc
      COMPLEX*16 :: rh
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: old_HH
      REAL*8 :: tmp_diff
      COMPLEX*16, DIMENSION(Ngrd,Ngrd) :: He,He1,He2
      INTEGER :: lop
      LOGICAL :: check
      COMPLEX*16 :: fact

      INTEGER :: sc1,sc2,sc3,sc4
      INTEGER :: ir1,ir2,ir3,ir4
      INTEGER :: tttot
      REAL*8 :: cttot
      REAL*8 :: ctime1, ctime2
      REAL*8 :: ctime3, ctime4

      COMPLEX*16, DIMENSION(Ngrd) :: old_C
      COMPLEX*16, DIMENSION(Ngrd) :: act_C

C Decoherence Part
      REAL*8 :: EKIN,ekinnu,sum_c,TAU
      REAL*8, DIMENSION(N) :: vx,vy,vz
      REAL*8, DIMENSION(Ngrd) :: ENERGX
      REAL*8, DIMENSION(N) :: mass

C----------------------------------------------------------C

      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_UP>>>"
         WRITE(*,*)
         WRITE(*,'("Number of integration steps: ",I6)')NE
         WRITE(*,'("Type: ",I1)')type
         WRITE(*,'("value for decoherence correction (0.0=off): " 
     &   ,F4.1)'),dec_cor
      END IF

C Initialitation of time interval
      d = dt/REAL(NE)

C Initialization of variables
      g = 0.D0

C   TIMING - C integration time, start
      IF (iout .GE. 2) THEN
         tttot = 0
         cttot = 0.D0
         call SYSTEM_CLOCK(sc1,ir1)
         CALL CPU_TIME(ctime1)
      END IF

C Starting cycle over integration steps
      DO i = 0, NE-1
         old_C = C

C       Computation of quantum amplitudes
C        Build Hermitian matrices
         IH = H * DCMPLX(0.D0,-1.D0)
         IHold = Hold * DCMPLX(0.D0,-1.D0)

C        Compute Hermittian effective Hamiltonian at given time
         DO j = 1, Ngrd
            DO l = 1, Ngrd
               IF (type.EQ.1 .OR. type.EQ.3 ) THEN
                  HH(j,l) = IHold(j,l) + (IH(j,l)-IHold(j,l))
     &                      *(REAL(i)+0.5D0)/REAL(NE)
               ELSE IF (type .EQ. 2) THEN
                  HH(j,l) = IHold(j,l) + (IH(j,l)-IHold(j,l))
     &                      *REAL(i)/REAL(NE)
               END IF
            END DO
         END DO

         HH = HH*d
         
         old_HH = HH

C      IF eigenvalue expansion then diagonalize
         IF (type .NE. 3) THEN

C        Diagonalization        
            CALL ZHEEV('V','U',Ngrd,HH,Ngrd,w,WORK,2*Ngrd,RWORK,info)
            IF (info .NE. 0) THEN
               WRITE(*,*)"Diagonalization failed"
               STOP 'Abnormal termination in int_UP'
            END IF

C        Build inverse of the matrix of eigenvectors         
            DO j = 1, Ngrd
               DO l = 1, Ngrd
                  inv_HH(j,l) = CONJG(HH(l,j))
               END DO
            END DO
        
C        Build diagonal matrix with exponentials of eigenvalues
            m1 = DCMPLX(0.D0,0.D0)
            DO j = 1, Ngrd
               m1(j,j) = DCMPLX(DCOS(w(j)),DSIN(w(j)))
            END DO

C        Matrix product HH*m1         
            DO j = 1, Ngrd
               DO k = 1, Ngrd
                  m2(j,k) = DCMPLX(0.D0,0.D0)
                  DO l = 1, Ngrd
                     m2(j,k) = m2(j,k) + HH(j,l)*m1(l,k)
                  END DO
               END DO
            END DO

C        Matrix product m2*inv_HH
            DO j = 1, Ngrd
               DO k = 1, Ngrd
                  m3(j,k) = DCMPLX(0.D0,0.D0)
                  DO l = 1, Ngrd
                     m3(j,k) = m3(j,k) + m2(j,l)*inv_HH(l,k)
                  END DO
               END DO
            END DO

         END IF

C      IF taylor expansion then built the series
         IF (type .EQ. 3) THEN

C build directly exponential         
            He = DCMPLX(0.D0,0.D0)
            DO j = 1, Ngrd
               He(j,j) = DCMPLX(1.D0,0.D0)
            END DO

            He1 = He
            fact = DCMPLX(1.D0,0.D0)
            check = .TRUE.
            DO lop = 1, 100
               DO j = 1, Ngrd
                  DO k = 1, Ngrd
                     He2(j,k) = DCMPLX(0.D0,0.D0)
                     DO l = 1, Ngrd
                        He2(j,k) = He2(j,k) + old_HH(l,k)*He1(j,l)
                     END DO
                  END DO
               END DO
               fact = fact*DCMPLX(0.D0,1.D0/REAL(lop))
               He1 = He2
               He2 = He2*fact

               He = He + He2

               DO j = 1, Ngrd
                  DO k = 1, Ngrd
                     tmp_diff = DSQRT(REAL(He2(j,k))**2 +
     &                       DIMAG(He2(j,k))**2)
                     IF (tmp_diff .GE. 1.D-10) THEN
                        check = .FALSE.
                        EXIT
                     END IF
                  END DO
               END DO

               IF (check) EXIT
               
               check = .TRUE.
            END DO

            IF (lop .GE. 100) THEN
               WRITE(*,*)
               WRITE(*,*)"convergency problems for exponential"
               STOP 'Abnormal termination in int_UP'
            END IF

            m3 = He

         END IF

C        Matrix product m3*C         
         DO j = 1, Ngrd
            cc(j) = DCMPLX(0.D0,0.D0)
            DO l = 1, Ngrd
               cc(j) = cc(j) + m3(j,l)*C(l)
            END DO
         END DO
         C = cc

C  *************************************
C  Decoherence correction  OW 06.01.2010
C  *************************************
         IF (dec_cor.gt.0) THEN
C         write(*,*)"************************************"
C         write(*,*)"Decoherence Correction of Amplitudes"
C         write(*,*)"************************************"
C         write(*,*)" "
            sum_C=0.0
            EKIN=0
            ekinnu=0
            do j=1,N
               ekinnu=vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)
               EKIN=EKIN+mass(j)*ekinnu
            end do
            EKIN=EKIN*2.39057D-3*0.5D0
            EKIN=EKIN/627.51
            DO j = 1, Ngrd
               IF (j.ne.istat) THEN
                  tau=627.51/(abs(ENERGX(j)-ENERGX(istat)))*
     &           (1+(dec_cor/EKIN))
                  C(j)=C(j)*DCMPLX(exp(-d*41341.105/tau))  

                  sum_C=sum_C+abs(C(j))**2
               END IF
            END DO
            C(istat)=C(istat)*sqrt((DCMPLX(1)-DCMPLX(sum_C))
     &   /DCMPLX(abs(C(istat)**2)))
         END IF

C   TIMING - probability integration time, start
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc3,ir3)
            CALL CPU_TIME(ctime3)
         END IF

C        Computation of probabilities
         act_C = old_C + (C-old_C)*0.5D0
         DO j = 1, Ngrd
            tmprd = rdold(j,istat) + (rd(j,istat)-rdold(j,istat))
     &                        *REAL(i+0.5)/REAL(NE)
            cb = act_C(j)*DCONJG(act_C(istat))*tmprd
            b = -2.D0*REAL(cb)
            g(j) = g(j) + b*d 
         END DO

C   TIMING - probability integration time, end
         IF (iout .GE. 2) THEN
            call SYSTEM_CLOCK(sc4,ir4)
            CALL CPU_TIME(ctime4)
            tttot = tttot + (sc4-sc3)
            IF (ctime4 .GT. REAL(sc4-sc3)/REAL(ir4)+ctime3) 
     &          ctime4 = REAL(sc4-sc3)/REAL(ir4)+ctime3
            cttot = cttot + (ctime4-ctime3)
         END IF

C   End cycle over integration steps  
      END DO

C   TIMING - C integration time, end      
      IF (iout .GE. 2) THEN
         call SYSTEM_CLOCK(sc2,ir2)
         CALL CPU_TIME(ctime2)
         wtime_eint = wtime_eint + sc2-sc1
            IF (ctime2 .GT. REAL(sc2-sc1)/REAL(ir1)+ctime1) 
     &          ctime2 = REAL(sc2-sc1)/REAL(ir1)+ctime1
         ctime_eint = ctime_eint + ctime2-ctime1
      END IF

C   TIMING - output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,'("-----------------------------------")')
         WRITE(*,*)"*** TIMING ***"
         WRITE(*,'(" Time for UP",I1," integration")')type
         CALL out_time(ctime1,ctime2-cttot,sc1,sc2,ir1)
         WRITE(*,*)
         WRITE(*,*)"Time for probability integration"         
         CALL out_time(0.D0,cttot,0,tttot,ir1)
         WRITE(*,'("-----------------------------------")') 
      END IF

C  Normalize probabilities  
      rh = C(istat)*DCONJG(C(istat))

      DO j = 1, Ngrd
         IF(g(j) .LE. 0.D0) THEN
            g(j) = 0.D0
         ELSE
            g(j) = g(j)/REAL(rh)
         END IF
      END DO

C   Output
      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"int_UP<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE



