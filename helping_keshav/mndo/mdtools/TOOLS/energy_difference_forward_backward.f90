      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER(KIND=4)                           :: I,J,NSTATES,STATE,STATE_LAST,TIME_STATE,NSTEP
      logical                                   :: HOP
      REAL(KIND=8)                              :: E_DIFF,CUTOFF ,E_TOT,E_TOT_OLD
      REAL(KIND=8), dimension( : ), allocatable :: E_STATE,E_STATE_OLD
      CHARACTER                                 :: WORD*6
      character(len=120)                        :: temp,str

!     --BEGIN--
      call system ("sed -n 's/.*NSTEP *=* *\([0-9][0-9]*\).*/\1/p' dynvar.in > tmp.txt")
      call system ("sed -n 's/.*iroot *=* *\([0-9][0-9]*\).*/\1/p' mndo.inp >> tmp.txt")
      OPEN(UNIT=11, FILE="tmp.txt",STATUS="OLD")
      READ (11,*) NSTEP
      READ (11,*) NSTATES
      CLOSE(11) 
      call system ("rm -f tmp.txt")
      NSTEP = NSTEP - 2
      allocate(E_STATE(NSTATES),E_STATE_OLD(NSTATES))
      CUTOFF=10

      OPEN(UNIT=12, FILE="./MD_data/state.dat",STATUS="OLD")
      OPEN(UNIT=13, FILE="hopping.out"        ,STATUS="OLD")
      OPEN(UNIT=15, FILE="stat.out",STATUS="OLD")
      READ (15,*)

!     JUDGE THE ENERGY DIFFERENCE BETWEEN THE SEQUENTIAL TWO STEPS
      STATE_LAST=0
      DO I=1,NSTEP,1
        ! FIND THE HOPPING POINT
        HOP=.false.
        READ (12,*) TIME_STATE, STATE
        IF (STATE/=STATE_LAST) THEN
          HOP=.true.
        ENDIF

        READ (15,*) TIME_STATE, E_TOT

        READ(13,*)
        READ(13,'(A120)') str
        str=adjustl(str)
        str=str(index(str," "):)
        DO J=1,NSTATES,1
          str=adjustl(str)
          temp=str(1:index(str," "))
          read(temp, '(F)') E_STATE(J)
          str=str(index(str," "):)
        ENDDO
        DO J=1,3+NSTATES*3,1
          READ(13,*)  
        ENDDO
        IF (I .GT. 1) THEN
          IF (HOP) THEN
            E_DIFF=0.0D0
          ELSE
            E_DIFF=E_STATE(STATE)-E_STATE_OLD(STATE)
          ENDIF
          IF (ABS(E_DIFF) .GE. CUTOFF) THEN
            OPEN(UNIT=14,FILE="bad_run_pot",STATUS="NEW")
!           INDEX STARTER IS A LITTLE DIFFERENT
            WRITE(14,"(2I,F10.5,I)") I+1,I,ABS(E_DIFF),STATE
            CLOSE(14)
            EXIT
          ENDIF

          E_DIFF=E_TOT-E_TOT_OLD
          IF (ABS(E_DIFF) .GE. CUTOFF) THEN
            OPEN(UNIT=14,FILE="bad_run_tot",STATUS="NEW")
!           INDEX STARTER IS A LITTLE DIFFERENT
            WRITE(14,"(2I,F10.5,I)") I+1,I,ABS(E_DIFF),STATE
            CLOSE(14)
            EXIT
          ENDIF

        ENDIF

        DO J=1,NSTATES,1
          E_STATE_OLD(J)=E_STATE(J)
        ENDDO

        E_TOT_OLD=E_TOT
        STATE_LAST=STATE
      ENDDO

      CLOSE(12)
      CLOSE(13)

      STOP
      END PROGRAM

