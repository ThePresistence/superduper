      SUBROUTINE ehrendyn(iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,
     1                   dt,vx,vy,vz,mass,N,istat,ENERGX,CGX,write_hop,
     2                   hopunit,Hold,rdold,C,NE,old_CICOMP,old_ENERGX,
     3                   old_de,num_CC,an_CC,Einteg,tully_rest,ehrgr,
     4                   fol_stat,istatt,inde,sw_CC)

C==========================================================================C
C This subroutine computes the effective gradient required for Ehrenfest   C
C molecular dynamics [1].                                                  C
C                                                                          C
C The subroutine is divided into five parts:                               C
C                                                                          C
C   1- Multiple CI runs to get needed data [SCFMUL]                        C
C   2- Propagate velocity half time step forward: in this way velocity and C
C        coordinates are computed at the same time step                    C
C   3- Computation of nonadiabatic coupling vectors [SCF]                  C
C   4- Evaluation of states' population (a in Ref. [1], C here)            C
C      population is propagated according to Eq. (26) of Ref. [2]          C
C      this equation is equivalent to Eq. 64 of Ref. [1] except for a      C
C      phase factor                                                        C
C     4.1 - create effective hamiltonian for propagation of population     C
C             (see Eq. (26) in Ref. [2]                                    C
C     4.2 - propagate population (solve Eq. (26) of Ref. [2])              C
C   5- Build effective gradient (see Eq. (83) of Ref. [1])                 C
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
C [1] N.L. Doltsinis, D. Marx, J.Theor.Comput.Chem. 1, 319 (2002)          C
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
C delta         Working variable for propagation of velocity               C
C dt            Time step [in picoseconds]                                 C
C Einteg        Integration algorithm to use for electronic equation in    C
C                Tully's hopping:                                          C
C                    euler  -> Euler method                                C
C                    RK2    -> Runge-Kutta second-order algorithm          C
C                    RK4    -> Runge-Kutta fourth-order algorithm          C
C ENERGX        Energies of computed states                                C
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
C startStep     fisrt MD step                                              C
C N             Number of atoms                                            C
C NE            Number of steps for electronic integration (in RK2)        C
C rd            Dot products of velocity times nonadiabatic vectors        C
C rdold         Dot product of velocity and non-adiabatric vectors         C
C                at previosu MD step                                       C
C sca           conversion factor from kcal/(mol*A*AMU) to A/ps^2          C
C SCFCAL        External routine for SCF calculations                      C
C tmph          Scratch variable                                           C
C tmpin2_160    Temporary scratch variable                                 C
C vx,vy,vz      Cartesian components of velocity                           C
C write_hop     Flag to allow writing hopping data to file                 C
C                                                                          C
C--------------------------------------------------------------------------C
C AUTHOR: E. Fabiano                                                       C
C MAIL: efabiano@mpi-muelheim.mpg.de                                       C
C DATE: March 2007                                                         C
C--------------------------------------------------------------------------C

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

      COMMON
     ./INOPT2/ IN2(300)
     ./GRDORG/ IGRST(LMGRD),ISTATE,JSTATE
     ./CGRAD / CG(3,LM1+LM1M)
     ./DYNVA/ CICOMP(LMCONF,6),NCICONF
     ./DYNVA2/ dee(40,LMCONF),naco


C  INPUT VARIABLES

      INTEGER :: iout
      INTEGER :: LM5
      REAL*8, DIMENSION(LM5) :: ARRAY
      INTEGER :: ICALL
      INTEGER :: mdstep,startStep
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
      LOGICAL :: fol_stat


C  INPUT/OUTPUT VARIABLES

      REAL*8, DIMENSION(N) :: vx,vy,vz
      INTEGER :: istat
      REAL*8, DIMENSION(LMGRD) :: ENERGX
      REAL*8, DIMENSION(3,LM1+LM1M,LMGRD) :: CGX
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: Hold
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rdold
      COMPLEX*16, DIMENSION(IN2(159)) :: C
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP
      REAL*8, DIMENSION(LMGRD) :: old_ENERGX
      INTEGER, DIMENSION(40,LMCONF) :: old_de
      REAL*8, DIMENSION(3,LM1+LM1M) :: ehrgr
      REAL*8, DIMENSION(N,IN2(159),IN2(159)) :: nacx_old
      REAL*8, DIMENSION(N,IN2(159),IN2(159)) :: nacy_old
      REAL*8, DIMENSION(N,IN2(159),IN2(159)) :: nacz_old


C  LOCAL VARIABLES 

      INTEGER :: i,j,l,k
      REAL*8, DIMENSION(LMGRD) :: CNORMX
      INTEGER :: tmpin2_160
      REAL*8, DIMENSION(3,N,IN2(159),IN2(159)) :: CCX
      COMPLEX*16, DIMENSION(IN2(159),IN2(159)) :: H
      REAL*8, DIMENSION(IN2(159),IN2(159)) :: rd
      INTEGER :: ICICAL
      REAL*8 :: delta
      INTEGER :: ICALL2
      REAL*8 :: tmph
      REAL*8 :: sws
      INTEGER, DIMENSION(IN2(159)) :: sw_CC
      REAL*8, DIMENSION(IN2(159)) :: g
      REAL*8 :: eh1,eh2
      INTEGER, DIMENSION(IN2(159)) :: inde
      INTEGER :: istatt


C  EXTERNAL ROUTINES

      EXTERNAL SCFCAL


C  PARAMETERS

      REAL*8, PARAMETER :: hbar = 1.517873021D-2 ! hbar in (kcal/mol)ps
      REAL*8, PARAMETER :: sca = 418.4 ! from kcal/(mol*A*AMU) to A/ps^2


C--------------------------------------------------------------------------C


      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)">>>EHRENFEST"
         WRITE(*,*)
      END IF

C  (3)

C  Compute nonadiabatic coupling vector

C      NUMERICAL
      IF (num_CC) THEN
         CALL m_num_CC(iout,IN2(159),NCICONF,CICOMP,OLD_CICOMP,
     &                dt,mdstep,startStep,naco,dee,old_de,rd,tully_rest,
     &                inde)
      END IF


C                ANALITICAL
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
C  First method for setting phase of CCX vectors 
C  (only possible with state following)
                  IF (mdstep-startStep .GE. 2 .or. tully_rest) THEN
                     sws = REAL(sw_CC(inde(i))*sw_CC(inde(j)))
                  END IF
                  CCX(:,:,i,j)=sws*CG/(ENERGX(ISTATE)-ENERGX(JSTATE))
                  CCX(:,:,j,i)=sws*CG/(ENERGX(JSTATE)-ENERGX(ISTATE))
                  
               END DO
            END DO

            IF (iout .GE. 3) THEN
               WRITE(*,*)
               WRITE(*,*)"     Non-adiabatic coupling vectors"
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

C  Compute Effective Hamiltonian and Rd products

      H = DCMPLX(0.D0,0.D0)


C     Diagonal part
      
      DO i = 1, IN2(159)
         tmph = (ENERGX(inde(i))-ENERGX(inde(1)))/hbar
         H(i,i) = DCMPLX(0.D0,-tmph)
      END DO


C     Off-diagonal part

      IF (an_CC) THEN

         rd = 0.D0

         DO i = 2, IN2(159)
            DO j = 1, i-1
               DO l = 1, N
                  rd(i,j) = rd(i,j) + CCX(1,l,i,j)*vx(l) + 
     &                 CCX(2,l,i,j)*vy(l) + CCX(3,l,i,j)*vz(l)
               END DO
               rd(j,i) = rd(i,j)
               H(i,j) = DCMPLX(-rd(i,j),0.D0)
               H(j,i) = H(i,j)
            END DO
         END DO

      ELSE
         DO i = 2, IN2(159)
            DO j = 1, i-1
               H(i,j) = DCMPLX(-rd(i,j),0.D0)
               H(j,i) = H(i,j)
            END DO
         END DO
      END IF



      IF (iout .GE. 3) THEN
         WRITE(*,*)
         WRITE(*,*)"       Non-adiabatic coupling terms"
         WRITE(*,*)"     [dot product of velocity and CCX]"
         WRITE(*,'(1X,10I15)')(j,j=2,IN2(159))
         DO i = 1, IN2(159)-1
            WRITE(*,'(I1,4X,10F15.6)')i,(rd(i,j),j=i+1,IN2(159))
         END DO
      END IF


      IF (iout .GE. 2) THEN
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


            


C store data for next step
      
      DO i = 1, IN2(159)
         old_CICOMP(:,i) = CICOMP(:,inde(i))
         old_ENERGX(i) = ENERGX(inde(i))
      END DO
      old_de = dee

C  If mdstep .LE. 1 (2 for numerical CC) then store H in Hold and 
C  rd in rdold and exit

      IF (.NOT. tully_rest) THEN
         IF ((.NOT. an_CC .AND. mdstep-startStep .LE. 2) .OR. 
     &             (an_CC .AND. mdstep-startStep .LE. 1)) THEN
            Hold = H
            rdold = rd
            C = DCMPLX(0.D0,0.D0)
            C(istat) = DCMPLX(1.D0,0.D0)
            ehrgr(:,:) = CGX(:,:,istatt)
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
      ELSE
         WRITE(*,*)Einteg
         WRITE(*,*)"Bad definition of electronic integrator!"
         WRITE(*,*)
         STOP 'Abnormal termination in SURHOP'
      END IF


      DO i = 1, IN2(159)
         IF (g(i) .GT. 1.D0) g(i) = 1.D0
      END DO
      
      CALL pop_norm(iout,IN2(159),C)



      IF (iout .GE. 2) THEN
         WRITE(*,*)
         WRITE(*,*)"    Population of excited states"
         WRITE(*,'("State",2X,"population")')
         DO i = 1, IN2(159)
            WRITE(*,'(2X,I1,5X,F7.5)')i,REAL(C(i)*CONJG(C(i)))
         END DO
         WRITE(*,*)
         WRITE(*,*)"    Transition probabilities"
         WRITE(*,'("State",2X,"probability")')
         DO i = 1, IN2(159)
            WRITE(*,'(2X,I1,5X,F8.6)')i,g(i)
         END DO
      END IF



C  update Hold and rdold

      Hold = H
      rdold = rd


C  (5)

C Build Ehrenfest effective gradient

      DO i = 1, 3
         DO j = 1, N
            eh1 = 0.D0
            eh2 = 0.D0
            DO k = 1, IN2(159)
               eh1 = eh1 + REAL(C(k)*CONJG(C(k)))*CGX(i,j,inde(k))
               DO l = 1, IN2(159)
                  eh2 = eh2 + REAL(C(l)*CONJG(C(k)))*
     &                 (ENERGX(inde(k))-ENERGX(inde(l)))*CCX(i,j,k,l)
               END DO
            END DO
               ehrgr(i,j) = eh1 + eh2 !CGX(i,j,istat) 
         END DO
      END DO




C   Write data if requested

      IF (write_hop) THEN
         OPEN(77,STATUS='REPLACE',FILE='hop.tmp')
         WRITE(77,'("mdstep: ",I4)')mdstep
         WRITE(77,'("energy",3F15.5)') (ENERGX(inde(j)),
     &                                       j=1,IN2(159))
         WRITE(77,'("energy differences")')
         DO i = 2, IN2(159)
            WRITE(77,'(4F15.5)') (ENERGX(inde(i))-ENERGX(inde(j)),
     &                                 j=1,i-1)
         END DO
         WRITE(77,'("coupling")')
         DO i = 2, IN2(159)
            WRITE(77,'(4F15.5)') (rd(i,j),j=1,i-1)
         END DO
         WRITE(77,'("population")')
         DO j = 1, IN2(159)
            write(77,'(1X,F8.5)')
     &                REAL(C(j)*DCONJG(C(j)))
         END DO
         WRITE(77,*)"---------------------------------"
         CLOSE(77)
      END IF


      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"EHRENFEST<<<"
         WRITE(*,*)
      END IF


      END SUBROUTINE

