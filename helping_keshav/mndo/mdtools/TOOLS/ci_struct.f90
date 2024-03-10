PROGRAM MAIN
  IMPLICIT NONE
  INTEGER(KIND=4)    STATE,STATE_LAST 
  INTEGER(KIND=4)    I,K,TIME_STATE,I_CI
  INTEGER(KIND=4)    NSTEP_STATE,N_ATOM,FSAV
  REAL(KIND=8)       X,Y,Z
  CHARACTER (LEN=10) ATOM_LABEL
  CHARACTER (LEN=80) FILE_NAME
  LOGICAL            FILE_EXIST
  !     --BEGIN--
  call system ("sed -n 's/.*NSTEP *=* *\([0-9][0-9]*\).*/\1/p' dynvar.in >> tmp.txt")
  call system ("sed -n 's/.*FSAV *=* *\([0-9][0-9]*\).*/\1/p' dynvar.in >> tmp.txt")
  call system ("sed -n 's/.*N_atoms *=* *\([0-9][0-9]*\).*/\1/p' dynam.restart >> tmp.txt")
  OPEN(UNIT=11, FILE="tmp.txt",STATUS="OLD")
  READ (11,*) NSTEP_STATE
  READ (11,*) FSAV
  READ (11,*) N_ATOM
  CLOSE(11) 
  call system ("rm -f tmp.txt")
  call system ("rm -f ci_struct*")
  NSTEP_STATE = NSTEP_STATE - 2

  OPEN(UNIT=12, FILE="./MD_data/state.dat",STATUS="OLD")
  STATE_LAST=0
  DO I=1,NSTEP_STATE,1 
    READ (12,*) TIME_STATE, STATE
    IF (STATE<STATE_LAST) THEN
      WRITE(FILE_NAME, '(A,I1,A,I1,A)') "ci_struct_",STATE_LAST,"_",STATE,".xyz"
      OPEN(UNIT=13, FILE="traj.out" ,STATUS="OLD")
      INQUIRE(FILE=FILE_NAME, EXIST=FILE_EXIST)
      IF (.NOT. FILE_EXIST) THEN
         OPEN(UNIT=14, FILE=FILE_NAME ,STATUS="NEW")

         I_CI=INT(I/FSAV)

         !            MOLDEN FORMAT HEADLINE
         DO K=1,N_ATOM+3,1
         READ (13,*)
         ENDDO

         !            FILTER THE PREVIOUS STRUCTURES
         DO K=1,(I_CI-1)*(N_ATOM+2),1
         READ (13,*) 
         ENDDO

         !            READ CI STRUCTURE
         READ  (13,*) 
         READ  (13,*)
         WRITE (14,*) N_ATOM
         WRITE (14,*) I 

         DO K=1,N_ATOM,1
         READ  (13,*)    ATOM_LABEL,X,Y,Z
         WRITE (14,'(A, 3(4X,F14.7))') ATOM_LABEL,X,Y,Z
         ENDDO
         CLOSE(14) 
      END IF
      CLOSE(13) 
    ENDIF
    IF (STATE .EQ. 1) THEN
      EXIT
    ENDIF
    STATE_LAST=STATE
  ENDDO
  CLOSE(12)

  STOP
END PROGRAM

