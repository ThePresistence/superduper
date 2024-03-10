C     ******************************************************************
C
C     Memory management package in (almost) standard F77
C     with memory compaction.
C
C     ******************************************************************
C
C   Philosophy:
C
C      This is special-purpose memory allocator with
C      lifetime of memory blocks explicitly divided
C      between computation stages. At the end of each
C      stage, memory allocator would discard all segments
C      which are not active for the extent of next stage
C      (calling routine is responsible either for swapping
C      them to secondary storage or is prepared to recompute
C      values if they are needed on further stages). Then,
C      if necessary to accomodate new segments, remaining 
C      memory blocks will be compacted, and space for new
C      blocks will be allocated. Maximum amount of memory
C      used could be slightly (by "lax" amount) more then
C      strictly necessary to accomodate all memory blocks
C      during most-demanding stage - this should minimize
C      copying expense between stages.
C
C      As the number of blocks which are expected to be
C      allocated is rather small (10 or so would be a 
C      reasonable estimation), as well as the total
C      number of stages, we will not be as careful with
C      allocator overheads as is common in general-purpose
C      memory allocation world.
C
C      We would not permit allocation of new blocks once
C      first stage was initiated, because for the type
C      of problem we are solving, memory requirements 
C      should be known in advance. Also, there is no way
C      to deallocate block, bar the complete reset of the
C      allocator routines.
C
C      To guard against code running over it's bid of
C      dynamic memory, fenceposts are placed before the
C      first word and after the last word of each block
C      of usable memory. This also should serve as the
C      additional check of the correctness of the 
C      memory compaction part of PSMSTG.
C
C   Use:
C
C      PSMINI(MEMLAX,IPRT) will initialize memory allocator
C      PSMALC(ISIZE,ISTCNT,ISTLST,ID) will allocate memory
C         block of size ISIZE. There is no error indication,
C         call IPSMUS() to get current top of memory arena.
C      IPSMOF(ID) return current position of the memory
C         block ID. It will abort caller if memory block is 
C         not active on the current stage.
C      PSMSTG(BASE,ISTAGE) set execution stage. BASE is the
C         base address of the memory arena, in case we need
C         to shuffle blocks to make everything to fit in.
C      IPSMUS() integer function, returns current top of the
C         memory arena.
C      PSMSTS will generate status report on memory allocation.
C      PSMCHK(BASE) will check integrity of memory arena by 
C         verifying contents of guard words at the beginning 
C         and end of each memory block
C
C      All other subroutines defined in this file should be
C      considered private to the memory allocator and should
C      never be called from outside.
C
C   Data structures: 
C
C      Please do not make external modules to rely on
C      these data - it is bad programming practice, and
C      should not be really necessary - all access should
C      be performed through IDs returned by allocator
C      routines.
C
C      Constants:
C
C         MAXBLK - Maximum number of blocks which could be
C                  allocated
C         MAXSTG - Maximum number of stages which could be
C                  present in program
C         MAXFRE - Largest number of free memory blocks
C                  (Usually 2*MAXBLK+1)
C         IDBASE - Arbitrary constant to shift ID numbers
C                  by, should prevent callers from passing
C                  bogus block IDs.
C
C      Variables (in COMMON /PSMMGR/):
C
C         IPRINT - Printing level
C                  -1 = No output, except on fatal errors
C                   0 = No output, except on fatal errors
C                   1 = Report compactions only
C                   5 = Report everything
C         ISTAGE - Currently active stage, zero if none
C         MAXUSE - Largest amount of memory used during
C                  all phases
C         MEMLAX - Extra amount of memory available for
C                  fragmentation
C         MEMTOP - Top of required memory (should be
C                  MAXUSE+MEMLAX). 
C         NBLKS  - Number of allocated blocks
C         NFREE  - Number of free memory blocks
C         MEMUSE - Amount of memory used during each
C                  of the stages, MEMUSE(MAXSTG)
C         IBLLEN - Length of each allocated memory block,
C                  IBLLEN(MAXBLK). Memory actually used
C                  by calling routine is 2 words less,
C                  because two words (at the beginning
C                  and end of block) are used to guard
C                  agains memory boundary problems.
C         IBLOFF - Position of each allocated block on
C                  memory arena, or -1 if not active on
C                  this stage. User data will begin one
C                  word after this point
C         IFRLEN - Length of each free block, IFRLEN(MAXFRE)
C         IFROFF - Position of each free block on memory arena,
C                  IFROFF(MAXFRE)
C         LBLACT - Active stage list for each of the
C                  memory blocks (LOGICAL LBLAST(MAXBLK,MAXSTG)
C         LMOFL  - Memory count overflow had occured
C
C   Bugs:
C
C      Common block /PSMMGR/ have larger scope then it should.
C      This could be fixed by removing SAVE attribute from it
C      and adding definition of PSMMGR to PSDRV, but as it is
C      relatively small and such change will hurt pleasant 
C      modularity which exist now, I won't do it.
C
C      Parameter IHUGE is only suitable for 32-bit machines.
C      I see no portable way of determining it in F77. For
C      64-bit machines, set it to 9223372036854775807, and if
C      you are lucky enough to have Fortran compiler at the 
C      level of ISO/IEC 1539(1991), to HUGE(1)
C
C      Then sum of the sizes of memory blocks exceeds IHUGE,
C      it is silently truncated to IHUGE. This should cause 
C      no problems in a "real" world, because this allocation
C      pattern will be refused on the grounds of insufficient
C      memory space anyway.
C
C  Speedups possible:
C
C     Inlining IPSMOF might be a good idea.
C
      SUBROUTINE PSMINI(LAXMEM,IPRT)
C
C   Initialize memory allocator.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      LAXMEM - Lax memory for memory fragmentation. It's OK
C               to set it to zero.
C      IPRT   - Printing level
C
C   Accessed common blocks:
C
C   Modified common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      IPRINT = IPRT
      MEMLAX = LAXMEM
      MEMTOP = LAXMEM
      ISTAGE = 0
      NBLKS  = 0
      NFREE  = 0
      MAXUSE = 0
      LMOFL  = .FALSE.
      DO 10 I=1,MAXSTG
          MEMUSE(I) = 0
   10 CONTINUE
C
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6, 10000) IPRINT, MEMLAX
      ENDIF
C
      RETURN
10000 FORMAT(' MEMORY ALLOCATOR INITIALIZED WITH IPRINT = ', I4,
     .       ' MEMLAX = ', I10 )
      END
C
      SUBROUTINE PSMALC(IUSIZE,ISTCNT,ISTLST,ID)
C
C   Allocate memory block ID
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      IUSIZE - Size of the memory block requested by
C               caller.
C      ISTCNT - Number of active stages, 0 if always
C               active
C      ISTLST - List of stages, contaning 1 for the
C               stages on which block will be active and
C               0 for inactive stages. All other values
C               are reserved for future use.
C      ID     - Output memory block ID.
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      PARAMETER (IDBASE=5639)
C-AK  PARAMETER (IHUGE=2147483647)
C-AK  PARAMETER (IHUGE=9223372036854775807)
      PARAMETER (IHUGE=HUGE(1))
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      DIMENSION ISTLST(ISTCNT)
C
      ISIZE = IUSIZE + 2
      IF( IPRINT.GE.5 ) THEN
          WRITE( 6, 11000 ) ISIZE, (ISTLST(I).EQ.1,I=1,ISTCNT)
      ENDIF
      IF( ISTAGE.NE.0 ) THEN
          WRITE(NB6,10000) ISTAGE
          STOP 'PSMALC'
      ENDIF
C
C    Register memory block
C
      NBLKS = NBLKS + 1
      IF( NBLKS.GT.MAXBLK ) THEN
          WRITE(NB6,10100) MAXBLK
          STOP 'PSMALC'
      ENDIF
      IBLLEN(NBLKS) = ISIZE
      IBLOFF(NBLKS) = -1
C
C    Update stage info
C
      LBLACT(NBLKS,0) = .FALSE.
      DO 200 I=1,MAXSTG
          LBLACT(NBLKS,I) = .FALSE.
          IF( I.LE.ISTCNT ) THEN
              IF( ISTLST(I).EQ.1 ) GOTO 110
          ENDIF
          IF( ISTCNT.NE.0 ) GOTO 200
  110     CONTINUE
C
C         Memory block will be active on stage I
C
          LBLACT(NBLKS,I) = .TRUE.
          IF( IHUGE - MEMUSE(I) .LT. ISIZE ) THEN
              MEMUSE(I) = IHUGE
              LMOFL     = .TRUE.
          ELSE
              MEMUSE(I) = MEMUSE(I) + ISIZE
          ENDIF
          IF( MEMUSE(I).GT.MAXUSE ) MAXUSE = MEMUSE(I)
  200 CONTINUE
      MEMTOP = MAXUSE + MEMLAX
      ID = NBLKS + IDBASE
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,11100) NBLKS, ID
          WRITE(NB6,11200) MEMTOP
      ENDIF
C
      RETURN
10000 FORMAT(' MEMORY ALLOCATION ATTEMPTED DURING STAGE ', I3,
     .       '. PROGRAM WILL STOP' )
10100 FORMAT(' ALLOWED NUMBER OF DYNAMIC MEMORY BLOCKS (',I3, 
     .       ') EXCEEDED. PROGRAM WILL STOP' )
11000 FORMAT(' ALLOCATING ', I10, ' WORDS, ACTIVE ON STAGE(S):', 
     .       20(1X,L2))
11100 FORMAT(' ALLOCATED SLOT ', I3, ' WITH ID ', I10 )
11200 FORMAT(' MEMORY ARENA TOP SO FAR: ', I10 )
      END
C
      FUNCTION IPSMOF(ID)
C
C   Return memory arena offset of the active memory block.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      ID     - Memory block ID returned by PSMALC
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      PARAMETER (IDBASE=5639)
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      IF( ID.LE.IDBASE .OR. ID.GT.IDBASE+NBLKS ) THEN
          WRITE(NB6,10000) ID
          STOP 'IPSMOF'
      ENDIF
      NID = ID - IDBASE
      IF( IBLOFF(NID).EQ.-1 ) THEN
          WRITE(NB6,10100) ID, ISTAGE
          STOP 'IPSMOF'
      ENDIF
      IPSMOF = IBLOFF(NID) + 1
C
      RETURN
10000 FORMAT(' BOGUS MEMORY BLOCK ID ', I10, ' PASSED TO IPSMOF.',
     .       ' PROGRAM WILL STOP' )
10100 FORMAT(' ADDRESS OF MEMORY BLOCK ', I10, ', UNALLOCATED ON STAGE '
     .       ,I3,', WAS REQUESTED FROM IPSMOF. PROGRAM WILL STOP' )
      END
C
      SUBROUTINE PSMSTG(BASE,INEWST)
C
C   Setup memory arena for the new execution stage,
C   compacting allocated blocks if necessary.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      BASE   - Base address of the memory arena
C      INEWST - New execution stage. If INEWST is equal
C               to the current execution stage, call
C               is ignored.
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C
C   Local storage:
C
C      2*MAXBLK + 2*MAXFRE INTEGER cells used for uncommited
C      copies of memory control blocks and sorting of blocks
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      PARAMETER (GUARDL=19681002D0)
      PARAMETER (GUARDH=10021968D0)
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      DIMENSION BASE(*)
C
      DIMENSION NIBLOF(MAXBLK), NIFROF(MAXFRE), NIFRLE(MAXFRE)
      DIMENSION IBLIND(MAXBLK)
C
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,10000) INEWST, ISTAGE
      ENDIF
      IF( INEWST.LT.0 .OR. INEWST.GT.MAXSTG ) THEN
          WRITE(NB6,11000) INEWST
          STOP 'PSMSTG'
      ENDIF
      IF( INEWST.EQ.ISTAGE ) RETURN
      IF( ISTAGE.EQ.0 ) THEN
C
C         First stage special: initialize free list
C
          NFREE     = 1
          IFRLEN(1) = MEMTOP
          IFROFF(1) = 1
      ELSE
C
C         First, check memory blocks integrity and die if they
C         were stepped over
C
          CALL PSMCHK(BASE)
          IF( INEWST.EQ.0 ) THEN
C
C             Shut down memory allocator
C
              ISTAGE = 0
              RETURN
          ENDIF
C
C         Now, put blocks which will not be active
C         for the next phase on a free list
C
          DO 100 I=1,NBLKS
              IF( LBLACT(I,ISTAGE) .AND. .NOT. LBLACT(I,INEWST) ) THEN
                  IF( IPRINT.GE.5 ) THEN
                      WRITE(NB6,10100) I, IBLLEN(I), IBLOFF(I)
                  ENDIF
                  NFREE = NFREE + 1
                  IF( NFREE.GT.MAXFRE ) THEN
                      WRITE(NB6,11100)
                      STOP 'PSMSTG'
                  ENDIF
                  IFRLEN(NFREE) = IBLLEN(I)
                  IFROFF(NFREE) = IBLOFF(I)
                  IBLOFF(I)     = -1
              ENDIF
  100     CONTINUE
C
C         Sort free list. Use straight insertion, as it 
C         is going to be rather small sort...
C
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,10150) ( IFROFF(I), IFRLEN(I), I=1,NFREE )
          ENDIF
          DO 200 I=2,NFREE
              IFLTMP = IFRLEN(I)
              IFOTMP = IFROFF(I)
              DO 190 J=I-1,1,-1
                  IF( IFOTMP.GT.IFROFF(J) ) THEN
                      IFRLEN(J+1) = IFLTMP
                      IFROFF(J+1) = IFOTMP
                      GOTO 200
                  ENDIF
                  IFRLEN(J+1) = IFRLEN(J)
                  IFROFF(J+1) = IFROFF(J)
  190         CONTINUE
              IFRLEN(1) = IFLTMP
              IFROFF(1) = IFOTMP
  200     CONTINUE
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,10200) ( IFROFF(I), IFRLEN(I), I=1,NFREE )
          ENDIF
C
C         Coalesce adjasent free memory blocks
C
          NEWFRE=1
          DO 300 I=2,NFREE
              IPRTOP = IFROFF(NEWFRE)+IFRLEN(NEWFRE)
              IF( IFRLEN(I).NE.0 ) THEN
                  IF( IFROFF(I).GT.IPRTOP ) THEN
                      NEWFRE = NEWFRE + 1
                      IFROFF(NEWFRE) = IFROFF(I)
                      IFRLEN(NEWFRE) = IFRLEN(I)
                  ELSE IF( IFROFF(I).EQ.IPRTOP ) THEN
                      IFRLEN(NEWFRE) = IFRLEN(NEWFRE)+IFRLEN(I)
                  ELSE
                      WRITE(NB6,11200)
                      STOP 'PSMSTG'
                  ENDIF
              ENDIF
  300     CONTINUE
          NFREE = NEWFRE
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,10300) ( IFROFF(I), IFRLEN(I), I=1,NFREE )
          ENDIF
      ENDIF
C
C     Try to allocate new blocks into free areas. Do not make
C     allocation changes permanent yet.
C
      DO 1000 IPASS=1,2
          DO 500 I=1, NFREE
              NIFROF(I) = IFROFF(I)
              NIFRLE(I) = IFRLEN(I)
  500     CONTINUE
          DO 600 I=1,NBLKS
              NIBLOF(I) = IBLOFF(I)
              IF( LBLACT(I,INEWST) .AND. .NOT. LBLACT(I,ISTAGE) ) THEN
                  IF( IPRINT.GE.5 ) THEN
                      WRITE(NB6,10400) I, IBLLEN(I)
                  ENDIF
                  DO 550 J=1,NFREE
                      IF( NIFRLE(J).GE.IBLLEN(I) ) THEN
                          IF( IPRINT.GE.5 ) THEN
                              WRITE(NB6,10500) J, NIFRLE(J), NIFROF(J)
                          ENDIF
                          NIFRLE(J) = NIFRLE(J) - IBLLEN(I)
                          NIBLOF(I) = NIFROF(J)
                          NIFROF(J) = NIFROF(J) + IBLLEN(I)
                          BASE(NIBLOF(I))=GUARDL+I
                          BASE(NIBLOF(I)+IBLLEN(I)-1)=GUARDH+I
                          GOTO 560
                      ENDIF
  550             CONTINUE
C
C                 There is no free block which could satisfy our request,
C                 do compaction now and try to fit blocks into memory again
C
                  IF( IPRINT.GE.1 ) THEN
                      WRITE(NB6,11400) INEWST
                  ENDIF
C
C                 Sort table of indices to allocated memory blocks to
C                 simplify compaction pass. (BTW, who told you Fortran
C                 does not support pointers to pointers and pointers
C                 comparison? Here they are :-)
C
                  IBLC = 0
                  DO 552 IBL=1,NBLKS
                      IF( IBLOFF(IBL).NE.-1 ) THEN
                          IBLC = IBLC + 1
                          IBLIND(IBLC) = IBL
                      ENDIF
  552             CONTINUE
C
                  IF( IPRINT.GE.5 ) THEN
                      WRITE(NB6,11600) ( IBLIND(IT), IBLOFF(IBLIND(IT)), 
     .                                 IBLLEN(IBLIND(IT)), IT=1,IBLC)
C
                  ENDIF
                  DO 554 IT=2,IBLC
                      IBLTMP = IBLIND(IT)
                      KEY    = IBLOFF(IBLTMP)
                      DO 553 JT=IT-1,1,-1
                          IF( KEY.GT.IBLOFF(IBLIND(JT)) ) THEN
                              IBLIND(JT+1) = IBLTMP
                              GOTO 554
                          ENDIF
                          IBLIND(JT+1) = IBLIND(JT)
  553                 CONTINUE
                      IBLIND(1) = IBLTMP
  554             CONTINUE
C
                  IF( IPRINT.GE.5 ) THEN
                      WRITE(NB6,11700) ( IBLIND(IT), IBLOFF(IBLIND(IT)), 
     .                                 IBLLEN(IBLIND(IT)), IT=1,IBLC)
                  ENDIF
C
                  IBASE = 1
                  DO 559 IBL = 1, IBLC
                      IBLOCK = IBLIND(IBL)
                      IF( IBLOFF(IBLOCK).GT.IBASE ) THEN
C
C                         There is a gap between top of the previous block
C                         and bottom of this one, move the block down
C
                          NBASE = IBLOFF(IBLOCK)
                          DO 555 J=0,IBLLEN(IBLOCK)-1
                              BASE(IBASE+J) = BASE(NBASE+J)
  555                     CONTINUE
                          IBLOFF(IBLOCK) = IBASE
                      ELSE IF( IBLOFF(IBLOCK).LT.IBASE ) THEN
                          WRITE(NB6,11500)
                          STOP 'PSMSTG'
                      ENDIF
                      IBASE = IBASE + IBLLEN(IBLOCK)
  559             CONTINUE
                  NFREE=1
                  IFROFF(1) = IBASE
                  IFRLEN(1) = MEMTOP - IBASE + 1
                  GOTO 1000
  560             CONTINUE
              ENDIF
  600     CONTINUE
C
C         All requests for memory were satisfied, commit changes
C
          DO 700 I=1,NBLKS
              IBLOFF(I) = NIBLOF(I)
  700     CONTINUE
          DO 710 I=1,NFREE
              IFROFF(I) = NIFROF(I)
              IFRLEN(I) = NIFRLE(I)
  710     CONTINUE
          ISTAGE = INEWST
C
C         Make sure we had not damaged blocks in the process
C
          CALL PSMCHK(BASE)
          RETURN
 1000 CONTINUE
      WRITE(NB6,11300)
      STOP 'PSMSTG'
C
10000 FORMAT(' INITIALIZING MEMORY ARENA FOR THE EXECUTION STAGE ', I3,
     .       ', CURRENT STAGE ', I3 )
10100 FORMAT(' RELEASING MEMORY OCCUPIED BY SLOT ', I3, ' LEN = ',
     .       I10, ' BASE = ', I10 )
10150 FORMAT(' LIST OF FREE BLOCKS BEFORE SORT:',
     .       10(/5X,3(1X,'(',I10,',',I10,')')))
10200 FORMAT(' LIST OF SORTED FREE BLOCKS:',
     .       10(/5X,3(1X,'(',I10,',',I10,')')))
10300 FORMAT(' LIST OF COALESCED FREE BLOCKS:',
     .       10(/5X,3(1X,'(',I10,',',I10,')')))
10400 FORMAT(' LOOKING FOR PLACE FOR MEMORY BLOCK ', I3, 
     .       ' LEN = ', I10 )
10500 FORMAT(' FIT INTO FREE BLOCK ', I3, ' LEN = ', I10, ' AT ', I10 )
11000 FORMAT(' NEW MEMORY STAGE ', I10, ' IS OUTSIDE VALID RANGE.',
     .       ' PROGRAM WILL STOP.' )
11100 FORMAT(' FREE BLOCK LIST OVERFLOWED. PROGRAM WILL STOP.' )
11200 FORMAT(' OVERLAPPING ENTRIES ON FREE LIST DETECTED, MEMORY ',
     .       'CORRUPTION IS PROBABLE TO OCCURE.'/
     .       ' PROGRAM WILL STOP.' )
11300 FORMAT(' CANNOT FIT REQUESTED AREAS INTO MEMORY ARENA AFTER ',
     .       'COMPACTION PASS. PROGRAM WILL STOP.' )
11400 FORMAT(' RUNNING MEMORY COMPACTION PASS AT BEGINNING OF STAGE ',
     .       I3 )
11500 FORMAT(' OVERLAPPING ENTRIES ON OCCUPIED LIST DETECTED, MEMORY ',
     .       'CORRUPTION DID OCCUR.'/
     .       ' PROGRAM WILL STOP' )
11600 FORMAT(' LIST OF OCCUPIED BLOCKS BEFORE SORT: '/
     .         (' IND = ', I3, ' START = ', I10, ' LEN = ', I10))
11700 FORMAT(' LIST OF SORTED OCCUPIED BLOCKS: '/
     .         (' IND = ', I3, ' START = ', I10, ' LEN = ', I10))
      END
C
      FUNCTION IPSMUS()
C
C   Return current top of memory arena required to satisfy
C   all memory requests up to now.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      None.
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
      SAVE /PSMMGR/
C
      IPSMUS = MEMTOP
C
      RETURN
      END
C
      SUBROUTINE PSMSTS
C
C   Report current state of memory allocation system
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      None.
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above.
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      PARAMETER (IDBASE=5639)
      LOGICAL LBLACT, LMOFL
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      WRITE(NB6,10000) ISTAGE
      IF( LMOFL ) WRITE(NB6,10005)
      DO 100 I=1,MAXSTG
          IF( MEMUSE(I).GT.0 ) THEN
              WRITE(NB6,10010) I, MEMUSE(I)
              DO 90 J=1,NBLKS
                  IF( LBLACT(J,I) ) THEN
                      IF( ISTAGE.EQ.I ) THEN
                          WRITE(NB6,10020) J, J+IDBASE, IBLLEN(J), 
     .                                   IBLOFF(J)
                      ELSE
                          WRITE(NB6,10030) J, J+IDBASE, IBLLEN(J)
                      ENDIF
                  ENDIF
   90         CONTINUE
          ENDIF
  100 CONTINUE
      WRITE(NB6,10040) NFREE
      DO 200 I=1,NFREE
          WRITE(NB6,10050) I, IFRLEN(I), IFROFF(I)
  200 CONTINUE
C
      RETURN
10000 FORMAT(/' MEMORY ALLOCATION STATUS'
     .       /' CURRENT STAGE = ', I3)
10005 FORMAT( ' TOTAL MEMORY REQUEST COUNT OVERFLOWED' )
10010 FORMAT( ' FOR STAGE ', I3, ' MEMORY USE IS ', I10
     .       /' ACTIVE BLOCKS ARE:' )
10020 FORMAT( '      ', I2, ': ID = ', I10, ' LEN = ', I10, 
     .        ' AT ', I10 )
10030 FORMAT( '      ', I2, ': ID = ', I10, ' LEN = ', I10 )
10040 FORMAT( ' FREE LIST CONTAINS ', I3, ' BLOCKS' )
10050 FORMAT( '      ', I2, ': LEN = ', I10, ' AT ', I10 )
      END
C
      SUBROUTINE PSMCHK(BASE)
C
C   Check integrity of the memory guard words.
C
C   Coded by: Serge Pachkovsky
C
C   Parameters:
C
C      BASE   - Base address of the memory arena
C
C   Accessed common blocks:
C
C      PSMMGR - Data structures used by memory allocator.
C      PSPRT  - Printing unit.
C
C   Modified common blocks:
C
C      None.
C
C   Local storage:
C
C      None.
C
C   Module logic:
C
C      See discussion above
C
C   Bugs:
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXBLK=60)
      PARAMETER (MAXSTG=9)
      PARAMETER (MAXFRE=2*MAXBLK+1)
      PARAMETER (IDBASE=5639)
      PARAMETER (GUARDL=19681002D0)
      PARAMETER (GUARDH=10021968D0)
      LOGICAL LBLACT, LMOFL, HALT
      COMMON
     ./PSMMGR/ IPRINT, ISTAGE, MAXUSE, MEMLAX, MEMTOP, NBLKS, NFREE,
     .         MEMUSE(MAXSTG), IBLLEN(MAXBLK), IBLOFF(MAXBLK),
     .         IFRLEN(MAXFRE), IFROFF(MAXFRE),
     .         LBLACT(MAXBLK,0:MAXSTG), LMOFL
     ./PSPRT / NB6
      SAVE /PSMMGR/, /PSPRT /
C
      DIMENSION BASE(*)
C
      IF( IPRINT.GE.5 ) THEN
          WRITE(NB6,11000)
      ENDIF
      IF( ISTAGE.LT.1 .OR. ISTAGE.GT.MAXSTG ) THEN
          IF( IPRINT.GE.5 ) THEN
              WRITE(NB6,11100) ISTAGE
          ENDIF
          RETURN
      ENDIF
      HALT = .FALSE.
      DO 1000 I=1,NBLKS
          IF( ((IBLOFF(I).EQ.-1) .AND.       LBLACT(I,ISTAGE)) .OR.
     .        ((IBLOFF(I).NE.-1) .AND. .NOT. LBLACT(I,ISTAGE)) ) THEN
              WRITE(NB6,10000) I, I+IDBASE
              HALT = .TRUE.
          ENDIF
          IF( IBLOFF(I).NE.-1 ) THEN
              IF( BASE(IBLOFF(I)).NE.GUARDL+I ) THEN
                 WRITE(NB6,10100) I, I+IDBASE, GUARDL+I, BASE(IBLOFF(I))
                 HALT = .TRUE.
              ENDIF
              IF( BASE(IBLOFF(I)+IBLLEN(I)-1).NE.GUARDH+I ) THEN
                  WRITE(NB6,10200) I, I+IDBASE, GUARDH+I, 
     .                           BASE(IBLOFF(I)+IBLLEN(I)-1)
                  HALT = .TRUE.
              ENDIF
          ENDIF
 1000 CONTINUE
      IF(HALT) THEN
          CALL PSMSTS
          STOP 'PSMCHK'
      ENDIF
C
      RETURN
10000 FORMAT(' MEMORY CONTROL STRUCTURE DESTROYED FOR BLOCK ',
     .       I4, ' (ID = ', I7, ')' )
10100 FORMAT(' BLOCK ', I4, ' (ID = ', I7, ') OVERFLOWED TO ',
     .       'LOW MEMORY.'
     .      /' GUARD WORD EXPECTED: ', G16.9, 
     .       ' GUARD WORD FOUND: ', G16.9 )
10200 FORMAT(' BLOCK ', I4, ' (ID = ', I7, ') OVERFLOWED TO ',
     .       'HIGH MEMORY.'
     .      /' GUARD WORD EXPECTED: ', G16.9, 
     .       ' GUARD WORD FOUND: ', G16.9 )
C
11000 FORMAT(' CHECKING INTEGRITY OF MEMORY BLOCKS' )
11100 FORMAT(' STAGE NUMBER ', I4, ' NOT IN VALID RANGE, NO BLOCKS')
      END
C 
