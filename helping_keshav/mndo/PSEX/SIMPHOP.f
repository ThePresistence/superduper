      SUBROUTINE SIMPHOP(iout,ARRAY,LM5,ICALL,SCFCAL,mdstep,startStep,
     1                   dt,N,istat,ENERGX,CGX,old_CICOMP,
     2                   istatt,rnd_gen,mass,vx,vy,vz,ihop,opt)

      USE LIMIT, ONLY: LM1, LMGRD, LM1M, LMV, LMCONF


      IMPLICIT NONE

C  COMMON

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
      INTEGER :: IN2
      INTEGER :: NCICONF
      REAL*8 :: CICOMP


      COMMON
     ./ctimemd/ ctime_v,ctime_g,ctime_sf,ctime_acc
     ./wtimemd/ wtime_v,wtime_g,wtime_sf,wtime_acc
     ./ctimemd2/ ctime_eh,ctime_eint,ctime_ncc
     ./wtimemd2/ wtime_eh,wtime_eint,wtime_ncc
     ./INOPT2/ IN2(300)
     ./DYNVA/ CICOMP(LMCONF,6),NCICONF


C  INPUT VARIABLES

      INTEGER :: iout
      INTEGER :: LM5
      REAL*8, DIMENSION(LM5) :: ARRAY
      INTEGER :: ICALL
      INTEGER :: N
      INTEGER :: mdstep,startStep
      REAL*8 :: dt
      CHARACTER(8) :: rnd_gen
      REAL*8, DIMENSION(N) :: mass
      INTEGER :: opt

C  INPUT/OUTPUT VARIABLES

      REAL*8, DIMENSION(LMGRD) :: ENERGX
      REAL*8, DIMENSION(3,LM1+LM1M,LMGRD) :: CGX
      REAL*8, DIMENSION(LMCONF,6) :: old_CICOMP
      INTEGER :: istat
      REAL*8, DIMENSION(N) :: vx,vy,vz


C  OUTPUT VARIABLES

      INTEGER :: istatt
      LOGICAL :: ihop


C  LOCAL VARIABLES 

      REAL*8 :: ctime1, ctime2
      INTEGER :: wtime1, wtime2
      INTEGER :: check_ir
      REAL*8, DIMENSION(LMGRD) :: CNORMX
      INTEGER :: ICICAL
      INTEGER, DIMENSION(IN2(159)) :: inde
      REAL*8 :: DE
      INTEGER :: i,l
      REAL*8 :: a1
      REAL*8 :: a2,ia2
      REAL*8 :: a3
      REAL*8 :: csi
      LOGICAL :: shop
      REAL*8, DIMENSION(3,N) :: F2
      REAL*8 :: GAMMA
      REAL*8 :: delta
      REAL*8 :: b


C  EXTERNAL ROUTINES

      EXTERNAL SCFCAL


C  PARAMETERS

      REAL*8, PARAMETER :: Etrh = 30.D0
      REAL*8, PARAMETER :: a1tr = 0.5D0
      REAL*8, PARAMETER :: a2tr = 0.5D0
      REAL*8, PARAMETER :: sca = 418.4 ! from kcal/(mol*A*AMU) to A/ps^2

C--------------------------------------------------------------------------C

C      opt = 6


      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)">>>SIMPHOP"
         WRITE(*,*)
      END IF


      shop = .FALSE.
      ihop = .FALSE.

      DO i = 1, IN2(159)
         inde(i) = i
      END DO



      IF (mdstep-startStep .LT. 2) THEN
         DO i = 1, IN2(159)
            old_CICOMP(:,i) = CICOMP(:,inde(i))
         END DO
         RETURN
      END IF






C------------------------C
C       OPTION 1         C
C........................C
C original formulation   C
C------------------------C

      IF (opt .EQ. 1) THEN


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE
            
            a1 = 0.D0
            a2 = 0.D0
            ia2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
               ia2 = ia2 + CICOMP(l,istat-1)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0*dt)


             a1 = ABS(a1)
             a2 = ABS(a2)
             ia2 = ABS(ia2)
            
            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2
            WRITE(*,*)"ia2 = ",ia2
            write(*,*)"a3 = ",a3


            IF (a2.GE.a2tr .AND. ia2.GE.a2tr) THEN
               istat = istat - 1
               WRITE(*,*)"down hopping"
               shop = .TRUE.
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF

C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - up", DE
            
            a1 = 0.D0
            a2 = 0.D0
            ia2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
               ia2 = ia2 + CICOMP(l,istat+1)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0*dt)
 
            a1 = ABS(a1)
            a2 = ABS(a2)
            ia2 = ABS(ia2)
            

            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2
            WRITE(*,*)"ia2 = ",ia2
            write(*,*)"a3 = ",a3


            IF (a2.GE.a2tr .AND. ia2.GE.a2tr) THEN
               istat = istat + 1
               WRITE(*,*)"up hopping"
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C-----------------------------C
C          OPTION 2           C
C.............................C
C difference of cross overlap C
C-----------------------------C

      ELSE IF (opt .EQ. 2) THEN


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE            

            a1 = 0.D0
            a2 = 0.D0
            ia2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
               ia2 = ia2 + CICOMP(l,istat-1)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0*dt)



             a1 = ABS(a1)

            
            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2
            WRITE(*,*)"ia2 = ",ia2
            write(*,*)"a2-ia2 = ", ABS(a2-ia2)*0.5D0
            write(*,*)"a3 = ",a3


            IF (ABS(a2-ia2)*0.5D0 .GE. a2tr) THEN
               istat = istat - 1
               WRITE(*,*)"down hopping"
               shop = .TRUE.
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - up", DE            

            a1 = 0.D0
            a2 = 0.D0
            ia2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
               ia2 = ia2 + CICOMP(l,istat+1)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0*dt)
 
            a1 = ABS(a1)
            

            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2
            WRITE(*,*)"ia2 = ",ia2
            WRITE(*,*)"a2-ia2 = ", ABS(a2-ia2)*0.5D0
            write(*,*)"a3 = ",a3


            IF (ABS(a2-ia2)*0.5D0 .GE. a2tr) THEN
               istat = istat + 1
               WRITE(*,*)"up hopping"
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF




C-----------------------------C
C          OPTION 3           C
C.............................C
C random number compared with C
C cross overlap               C
C-----------------------------C

      ELSE IF (opt .EQ. 3) THEN


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE            

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0)
            a1 = ABS(a1)

            
            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi




            IF (csi .LE. ABS(a3)) THEN
               istat = istat - 1
               WRITE(*,*)"down hopping"
               shop = .TRUE.
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN
            
            WRITE(*,*)"try - up", DE

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0) 
            a1 = ABS(a1)
            

            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi


            IF (csi .LE. ABS(a3)) THEN
               istat = istat + 1
               WRITE(*,*)"up hopping"
            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF




C-----------------------------C
C          OPTION 4           C
C.............................C
C random number compared with C
C cross overlap               C
C with vel adjustment         C
C-----------------------------C

      ELSE IF (opt .EQ. 4) THEN

         GAMMA = 0.D0


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE            

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0)
            a1 = ABS(a1)

            
            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi




            IF (csi .LE. ABS(a3)) THEN

               WRITE(*,*)"down hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat-1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat-1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat-1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,DE,istat-1,istat,
     &                     GAMMA,mdstep)

               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

               shop = .TRUE.

            END IF


            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN
            
            WRITE(*,*)"try - up", DE

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0) 
            a1 = ABS(a1)
            

            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi


            IF (csi .LE. ABS(a3)) THEN
               WRITE(*,*)"up hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat+1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat+1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat+1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,-DE,istat+1,
     &                     istat,GAMMA,mdstep)


               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF





C-----------------------------C
C          OPTION 5           C
C.............................C
C random number compared with C
C cross overlap only negative C
C with vel adjustment         C
C-----------------------------C

      ELSE IF (opt .EQ. 5) THEN

         GAMMA = 0.D0


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE            

            a1 = 0.D0
            a2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat-1)*old_CICOMP(l,istat-1)
            END DO
            IF (a1 .LT. 0.D0) THEN
               DO l = 1, NCICONF
                  CICOMP(l,istat) = -CICOMP(l,istat)
               END DO
            END IF
            IF (a2 .LT. 0.D0) THEN
               DO l = 1, NCICONF
                  CICOMP(l,istat-1) = -CICOMP(l,istat-1)
               END DO
            END IF

            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2

            DO l = 1, NCICONF
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0)
        
            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi


            IF (csi .LE. -a3) THEN
               WRITE(*,*)"down hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat-1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat-1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat-1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,DE,istat-1,istat,
     &                     GAMMA,mdstep)

               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

               shop = .TRUE.

            END IF


            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN
            
            WRITE(*,*)"try - up", DE

            a1 = 0.D0
            a2 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a2 = a2 + CICOMP(l,istat+1)*old_CICOMP(l,istat+1)
            END DO
            IF (a1 .LT. 0.D0) THEN
               DO l = 1, NCICONF
                  CICOMP(l,istat) = -CICOMP(l,istat)
               END DO
            END IF
            IF (a2 .LT. 0.D0) THEN
               DO l = 1, NCICONF
                  CICOMP(l,istat+1) = -CICOMP(l,istat+1)
               END DO
            END IF

            WRITE(*,*)"a1 = ",a1
            WRITE(*,*)"a2 = ",a2

            DO l = 1, NCICONF
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0) 

            write(*,*)"a3 = ",a3



            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi


            IF (csi .LE. -a3) THEN
               WRITE(*,*)"up hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat+1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat+1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat+1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,-DE,istat+1,istat,
     &                     GAMMA,mdstep)

               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF






C-----------------------------C
C          OPTION 6           C
C.............................C
C random number compared with C
C cross overlap               C
C with vel adjustment         C
C random number shift         C
C-----------------------------C

      ELSE IF (opt .EQ. 6) THEN

         GAMMA = 0.D0
         b = 0.5


C     DOWN HOP
      IF (istat .GT. 1) THEN

         DE = ENERGX(istat) - ENERGX(istat-1)

         IF (DE .LE. Etrh) THEN

            WRITE(*,*)"try - down", DE            

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat-1) 
     &                 - CICOMP(l,istat-1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0)
            a1 = ABS(a1)

            
            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3


            a3 = ABS(a3)
            a3 = 1 - (1-a3)**b
            write(*,*)"a3 = ",a3


            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi
C            write(*,*)"1-b=",1.D0-b
C            write(*,*)"csi**(1.D0-b)=",csi**(1.D0-b)
C            write(*,*)"(csi+1.D0)/2.D0=",(csi+1.D0)/2.D0
C            csi = b*csi*csi +(1.D0-b)*csi
C            csi = (csi**(1.D0-b))*(csi+1.D0)/2.D0
C            csi = csi**b
C            write(*,*)"csi = ",csi



            IF (csi .LE. ABS(a3)) THEN

               WRITE(*,*)"down hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat-1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat-1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat-1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,DE,istat-1,istat,
     &                     GAMMA,mdstep)

               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

               shop = .TRUE.

            END IF


            WRITE(*,*)"istat = ",istat

         END IF

      END IF


C     UP HOP
      IF (istat .LT. IN2(159) .AND. .NOT.shop) THEN
            
         DE = ENERGX(istat+1) - ENERGX(istat)
         
         IF (DE .LE. Etrh) THEN
            
            WRITE(*,*)"try - up", DE

            a1 = 0.D0
            a3 = 0.D0
            DO l = 1, NCICONF
               a1 = a1 + CICOMP(l,istat)*old_CICOMP(l,istat)
               a3 = a3 + CICOMP(l,istat)*old_CICOMP(l,istat+1) 
     &                 - CICOMP(l,istat+1)*old_CICOMP(l,istat)
            END DO

            a3 = a3/(2.D0) 
            a1 = ABS(a1)
            

            WRITE(*,*)"a1 = ",a1
            write(*,*)"a3 = ",a3

            a3 = ABS(a3)
            a3 = 1 - (1-a3)**b
            write(*,*)"a3 = ",a3


            IF (rnd_gen .EQ. "standard") THEN
               CALL RANDOM_NUMBER(csi)
            ELSE IF (rnd_gen .EQ. "PM_BD") THEN
               CALL nran2(csi,1,.FALSE.)
            ELSE IF (rnd_gen .EQ. "knuth") THEN
               CALL nran3(csi,1,.FALSE.)
            END IF

            write(*,*)"csi = ",csi
C            write(*,*)"1-b=",1.D0-b
C            write(*,*)"csi**(1.D0-b)=",csi**(1.D0-b)
C            write(*,*)"(csi+1.D0)/2.D0=",(csi+1.D0)/2.D0
C            csi = b*csi*csi +(1.D0-b)*csi
C            csi = (csi**(1.D0-b))*(csi+1.D0)/2.D0
C            csi = csi**b
C            write(*,*)"csi = ",csi

            IF (csi .LE. ABS(a3)) THEN
               WRITE(*,*)"up hopping"
               DO l = 1, N
                  F2(1,l) = CGX(1,l,istat)-CGX(1,l,istat+1)
                  F2(2,l) = CGX(2,l,istat)-CGX(2,l,istat+1)
                  F2(3,l) = CGX(3,l,istat)-CGX(3,l,istat+1)
               END DO

               ihop = .TRUE.

               CALL ADJEK2(iout,N,F2,mass,vx,vy,vz,dt,-DE,istat+1,
     &                     istat,GAMMA,mdstep)


               DO i = 1, N
                  vx(i) =  vx(i) - GAMMA*F2(1,i)/mass(i)
                  vy(i) =  vy(i) - GAMMA*F2(2,i)/mass(i)
                  vz(i) =  vz(i) - GAMMA*F2(3,i)/mass(i)
               END DO

            END IF

            WRITE(*,*)"istat = ",istat

         END IF

      END IF












      END IF


      



      write(200,*)mdstep,istat


      istatt = inde(istat)


      DO i = 1, IN2(159)
         old_CICOMP(:,i) = CICOMP(:,inde(i))
      END DO


      IF (iout .GE. 1) THEN
         WRITE(*,*)
         WRITE(*,*)"SIMPHOP<<<"
         WRITE(*,*)
      END IF

      END SUBROUTINE

