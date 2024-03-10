C     ******************************************************************
C
C     Public-domain GEMM-based BLAS3 routines.
C     Author: Per Ling, University of Umea, Sweden.
C     Written in December 1993 on the basis of the
C     original BLAS3 code from Jack Dongarra (February 1989).
C     Obtained from
C
C     ******************************************************************
C
C     LIST OF AVAILABLE BLAS3 ROUTINES AND TYPICAL CALLS.
C
C     LOGIC1 = DBIGP(IP,DIM1,DIM2)
C     LOGIC2 = DCLD(LD)
C     CALL DSYMM ( SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC )
C     CALL DSYR2K( UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC )
C     CALL DTRSM ( SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,C,LDC)
C     CALL DTRMM ( SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,C,LDC )
C     CALL DSYRK ( UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC )
C
C     ******************************************************************
      LOGICAL FUNCTION DBIGP( IP, DIM1, DIM2 )
*     .. Scalar Arguments ..
      INTEGER                IP, DIM1, DIM2
*     ..
*
*  Purpose
*  =======
*
*  DBIGP determines which of two alternative code sections in a GEMM-
*  Based Level 3 BLAS routine that will be the fastest for a particular
*  problem. If the problem is considered large enough DBIGP returns
*  .TRUE., otherwise .FALSE. is returned. The input parameter IP
*  specifies the calling routine and a break point for alternative code
*  sections. The input parameters DIM1 and DIM2 are matrix dimensions.
*  The returned value is a function of the input parameters and the
*  performance characteristics of the two alternative code sections.
*
*  In this simple implementation, the returned values are determined by
*  looking at only one of the two dimensions DIM1 and DIM2. It may be
*  rewarding to rewrite the logical expressions in DBIGP so that both
*  dimensions are involved. The returned values should effectively
*  reflect the performance characteristics of the underlying BLAS
*  routines.
*
*
*  Input
*  =====
*
*  IP     - INTEGER
*           On entry, IP specifies which routine and which alternative
*           code sections that the decision is intended for.
*           Unchanged on exit.
*
*  DIM1   - INTEGER.
*           On entry, DIM1 specifies the first dimension in the calling
*           sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*  DIM2   - INTEGER.
*           On entry, DIM2 specifies the second dimension in the
*           calling sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*
*  -- Written in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. User specified parameters for DBIGP ..
      INTEGER            DIP41, DIP42, DIP81, DIP82, DIP83,
     $                   DIP91, DIP92, DIP93
      PARAMETER        ( DIP41 = 4, DIP42 = 3,
     $                   DIP81 = 4, DIP82 = 3, DIP83 = 4,
     $                   DIP91 = 4, DIP92 = 3, DIP93 = 4 )
*     ..
*     .. Executable Statements ..
      IF( IP.EQ.41 )THEN
         DBIGP = DIM1.GE.DIP41
      ELSE IF( IP.EQ.42 )THEN
         DBIGP = DIM2.GE.DIP42
      ELSE IF( IP.EQ.81 )THEN
         DBIGP = DIM2.GE.DIP81
      ELSE IF( IP.EQ.82 )THEN
         DBIGP = DIM2.GE.DIP82
      ELSE IF( IP.EQ.83 )THEN
         DBIGP = DIM1.GE.DIP83
      ELSE IF( IP.EQ.91 )THEN
         DBIGP = DIM2.GE.DIP91
      ELSE IF( IP.EQ.92 )THEN
         DBIGP = DIM2.GE.DIP92
      ELSE IF( IP.EQ.93 )THEN
         DBIGP = DIM1.GE.DIP93
      ELSE
         DBIGP = .FALSE.
      END IF
*
      RETURN
*
*     End of DBIGP.
*
      END
C     ******************************************************************
      LOGICAL FUNCTION DCLD( LD )
*     .. Scalar Arguments ..
      INTEGER                LD
*     ..
*
*  Purpose
*  =======
*
*  The size of the leading dimension of a two-dimensional array may
*  cause severe problems. Often when an array with a 'critical' leading
*  dimension is referenced, the execution time becomes significantly
*  longer than expected. This is caused by shortcomings of the memory
*  system.
*
*  The function DCLD returns .TRUE. if the leading dimension LD is
*  critical and .FALSE. if it is not critical. In this implementation
*  DCLD is designed to detect critical leading dimensions in an
*  environment with a multi-way associative cache. Parameters defining
*  cache characteristics are adjustable to match different machines.
*  It may be rewarding to rewrite DCLD for a machine with a different
*  cache policy.
*
*  The cache lines in a multi-way associative cache are divided among a
*  number of partitions, each containing the same number of lines. Each
*  address of main memory is mapped into a particular partition. The
*  number of lines in a partition equals the associativity. For example,
*  in a four way associative cache, each partition contain four cache
*  lines.
*
*  Data are transferred between the cache and main memory according to
*  an associative mapping scheme. A transfer of a data word from main
*  memory to cache is accomplished as follows. A unit of data
*  (data line) in main memory, with the size of a cache line, and
*  containing several contiguous data words including the referenced
*  one, is mapped (copied) to a certain partition in the cache memory.
*  The partition is determined by the location of the element in the
*  main memory and the associative mapping scheme. A replacement
*  algorithm makes room for the data line in one of the cache lines in
*  the selected partition. For example, an LRU-based (Least Recently
*  Used) replacement algorithm places the data line in the least
*  recently 'touched' cache line in the selected partition.
*
*
*  Input
*  =====
*
*  LD     - On entry, LD specifies the leading dimension of a
*           2-dimensional array. Unchanged on exit.
*
*
*  User specified parameters for DCLD
*  ================================
*
*  LNSZ   - Size of a cache line in number of bytes.
*
*  NPRT   - Number of partitions in the cache memory.
*
*  PRTSZ  - The number of cache lines in a partition that can be used
*           exclusively to hold a local array containing a matrix block
*           during the execution of a GEMM-Based Level 3 BLAS routine.
*           The remaining cache lines may be occupied by scalars,
*           vectors and possibly program code depending on the system.
*
*  LOLIM  - Leading dimensions smaller than or equal to LOLIM are not
*           considered critical.
*
*  DP     - Number of bytes in a double-precision word.
*
*
*  Local Variables and Parameters
*  ==============================
*
*  ONEWAY - The maximum number of double precision words that can be
*           stored in the cache memory if only a single cache line in
*           each partition may be used.
*
*  UPDIF  - The difference between the multiple of LD that is nearest
*           ONEWAY, or nearest a multiple of ONEWAY, and the nearest
*           multiple of ONEWAY that is larger than LD. In number of
*           double precision words.
*
*  MXDIF  - If both UPDIF and LD - UPDIF are less than MXDIF, and LD
*           is greater than LOLIM, then the leading dimension is
*           considered critical. Otherwise, the leading dimension is
*           considered not critical.
*
*
*  -- Written in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Variables ..
      INTEGER            UPDIF
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     .. User specified parameters for DCLD ..
      INTEGER            LOLIM, LNSZ, NPRT, PRTSZ, DP
      PARAMETER        ( LNSZ = 64, NPRT = 128, PRTSZ = 3,
     $                   LOLIM = 64, DP = 8 )
*     .. Parameters ..
      INTEGER            ONEWAY, MXDIF
      PARAMETER        ( ONEWAY = ( LNSZ*NPRT )/DP,
     $                   MXDIF = LNSZ/( DP*PRTSZ ) )
*     ..
*     .. Executable Statements ..
*
      IF( LD.LE.LOLIM )THEN
         DCLD = .FALSE.
      ELSE
         UPDIF = MOD( ( LD/ONEWAY )*ONEWAY+ONEWAY, LD )
         DCLD = MIN( UPDIF, LD-UPDIF ).LE.MXDIF
      END IF
*
      RETURN
*
*     End of DCLD.
*
      END
C     ******************************************************************
      SUBROUTINE DSYMM( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars,  A is a symmetric matrix and  B and
*  C are  m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  symmetric matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  symmetric  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  symmetric matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  symmetric matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA
      LOGICAL            LSIDE, UPPER
      INTEGER            I, II, IX, ISEC, J, JJ, JX, JSEC
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           DGEMM, DCOPY
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     .. User specified parameters for DSYMM ..
      INTEGER            RCB, CB
      PARAMETER        ( RCB = 96, CB = 48 )
*     .. Local Arrays ..
      DOUBLE PRECISION   T1( RCB, RCB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE = SIDE .EQ. 'L'
      UPPER = UPLO .EQ. 'U'
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.( SIDE .EQ. 'R' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER ).AND.( .NOT.( UPLO .EQ. 'L' ) ) )THEN
         INFO = 2
      ELSE IF( M.LT.0 )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         CALL DGEMM ( 'N', 'N', M, N, 0, ZERO, A, MAX( LDA, LDB ),
     $                                B, MAX( LDA, LDB ), BETA, C, LDC )
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( UPPER )THEN
*
*           Form  C := alpha*A*B + beta*C. Left, Upper.
*
            DO 40, II = 1, M, RCB
               ISEC = MIN( RCB, M-II+1 )
*
*              T1 := A, a upper triangular diagonal block of A is copied
*              to the upper triangular part of T1.
*
               DO 10, I = II, II+ISEC-1
                  CALL DCOPY ( I-II+1, A( II, I ), 1, T1( 1, I-II+1 ),
     $                                                               1 )
   10          CONTINUE
*
*              T1 :=  A', a strictly upper triangular diagonal block of
*              A is copied to the strictly lower triangular part of T1.
*              Notice that T1 is referenced by row and that the maximum
*              length of a vector referenced by DCOPY is CB.
*
               DO 30, JJ = II, II+ISEC-1, CB
                  JSEC = MIN( CB, II+ISEC-JJ )
                  DO 20, J = JJ+1, II+ISEC-1
                     CALL DCOPY ( MIN( JSEC, J-JJ ), A( JJ, J ), 1,
     $                                      T1( J-II+1, JJ-II+1 ), RCB )
   20             CONTINUE
   30          CONTINUE
*
*              C := alpha*T1*B + beta*C, a horizontal block of C is
*              updated using the general matrix multiply, DGEMM. T1
*              corresponds to a full diagonal block of the matrix A.
*
               CALL DGEMM ( 'N', 'N', ISEC, N, ISEC, ALPHA, T1( 1, 1 ),
     $                     RCB, B( II, 1 ), LDB, BETA, C( II, 1 ), LDC )
*
*              C := alpha*A'*B + C and C := alpha*A*B + C, general
*              matrix multiply operations involving rectangular blocks
*              of A.
*
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'T', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( 1, II ), LDA, B( 1, 1 ), LDB,
     $                                            ONE, C( II, 1 ), LDC )
               END IF
               IF( II+ISEC.LE.M )THEN
                  CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1, ALPHA,
     $                           A( II, II+ISEC ), LDA, B( II+ISEC, 1 ),
     $                                       LDB, ONE, C( II, 1 ), LDC )
               END IF
   40       CONTINUE
         ELSE
*
*           Form  C := alpha*A*B + beta*C. Left, Lower.
*
            DO 80, IX = M, 1, -RCB
               II = MAX( 1, IX-RCB+1 )
               ISEC = IX-II+1
*
*              T1 := A, a lower triangular diagonal block of A is copied
*              to the lower triangular part of T1.
*
               DO 50, I = II, II+ISEC-1
                  CALL DCOPY ( II+ISEC-I, A( I, I ), 1,
     $                                         T1( I-II+1, I-II+1 ), 1 )
   50          CONTINUE
*
*              T1 :=  A', a strictly lower triangular diagonal block of
*              A is copied to the strictly upper triangular part of T1.
*              Notice that T1 is referenced by row and that the maximum
*              length of a vector referenced by DCOPY is CB.
*
               DO 70, JX = II+ISEC-1, II, -CB
                  JJ = MAX( II, JX-CB+1 )
                  JSEC = JX-JJ+1
                  DO 60, J = II, JJ+JSEC-2
                     CALL DCOPY ( MIN( JSEC, JJ+JSEC-1-J ),
     $                                        A( MAX( JJ, J+1 ), J ), 1,
     $                       T1( J-II+1, MAX( JJ-II+1, J-II+2 ) ), RCB )
   60             CONTINUE
   70          CONTINUE
*
*              C := alpha*T1*B + beta*C, a horizontal block of C is
*              updated using the general matrix multiply, DGEMM. T1
*              corresponds to a full diagonal block of the matrix A.
*
               CALL DGEMM ( 'N', 'N', ISEC, N, ISEC, ALPHA, T1( 1, 1 ),
     $                     RCB, B( II, 1 ), LDB, BETA, C( II, 1 ), LDC )
*
*              C := alpha*A'*B + C and C := alpha*A*B + C, general
*              matrix multiply operations involving rectangular blocks
*              of A.
*
               IF( II+ISEC.LE.M )THEN
                  CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1, ALPHA,
     $                           A( II+ISEC, II ), LDA, B( II+ISEC, 1 ),
     $                                       LDB, ONE, C( II, 1 ), LDC )
               END IF
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'N', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( II, 1 ), LDA, B( 1, 1 ), LDB,
     $                                            ONE, C( II, 1 ), LDC )
               END IF
   80       CONTINUE
         END IF
      ELSE
         IF( UPPER )THEN
*
*           Form  C := alpha*B*A + beta*C. Right, Upper.
*
            DO 120, JJ = 1, N, RCB
               JSEC = MIN( RCB, N-JJ+1 )
*
*              T1 := A, a upper triangular diagonal block of A is copied
*              to the upper triangular part of T1.
*
               DO 90, J = JJ, JJ+JSEC-1
                  CALL DCOPY ( J-JJ+1, A( JJ, J ), 1, T1( 1, J-JJ+1 ),
     $                                                               1 )
   90          CONTINUE
*
*              T1 :=  A', a strictly upper triangular diagonal block of
*              A is copied to the strictly lower triangular part of T1.
*              Notice that T1 is referenced by row and that the maximum
*              length of a vector referenced by DCOPY is CB.
*
               DO 110, II = JJ, JJ+JSEC-1, CB
                  ISEC = MIN( CB, JJ+JSEC-II )
                  DO 100, I = II+1, JJ+JSEC-1
                     CALL DCOPY ( MIN( ISEC, I-II ), A( II, I ), 1,
     $                                      T1( I-JJ+1, II-JJ+1 ), RCB )
  100             CONTINUE
  110          CONTINUE
*
*              C := alpha*B*T1 + beta*C, a vertical block of C is updated
*              using the general matrix multiply, DGEMM. T1 corresponds
*              to a full diagonal block of the matrix A.
*
               CALL DGEMM ( 'N', 'N', M, JSEC, JSEC, ALPHA, B( 1, JJ ),
     $                     LDB, T1( 1, 1 ), RCB, BETA, C( 1, JJ ), LDC )
*
*              C := alpha*B*A + C and C := alpha*B*A' + C, general
*              matrix multiply operations involving rectangular blocks
*              of A.
*
               IF( JJ.GT.1 )THEN
                  CALL DGEMM ( 'N', 'N', M, JSEC, JJ-1, ALPHA,
     $                                  B( 1, 1 ), LDB, A( 1, JJ ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
               END IF
               IF( JJ+JSEC.LE.N )THEN
                  CALL DGEMM ( 'N', 'T', M, JSEC, N-JJ-JSEC+1, ALPHA,
     $                           B( 1, JJ+JSEC ), LDB, A( JJ, JJ+JSEC ),
     $                                       LDA, ONE, C( 1, JJ ), LDC )
               END IF
  120       CONTINUE
         ELSE
*
*           Form  C := alpha*B*A + beta*C. Right, Lower.
*
            DO 160, JX = N, 1, -RCB
               JJ = MAX( 1, JX-RCB+1 )
               JSEC = JX-JJ+1
*
*              T1 := A, a lower triangular diagonal block of A is copied
*              to the lower triangular part of T1.
*
               DO 130, J = JJ, JJ+JSEC-1
                  CALL DCOPY ( JJ+JSEC-J, A( J, J ), 1,
     $                                         T1( J-JJ+1, J-JJ+1 ), 1 )
  130          CONTINUE
*
*              T1 :=  A', a strictly lower triangular diagonal block of
*              A is copied to the strictly upper triangular part of T1.
*              Notice that T1 is referenced by row and that the maximum
*              length of a vector referenced by DCOPY is CB.
*
               DO 150, IX = JJ+JSEC-1, JJ, -CB
                  II = MAX( JJ, IX-CB+1 )
                  ISEC = IX-II+1
                  DO 140, I = JJ, II+ISEC-2
                     CALL DCOPY ( MIN( ISEC, II+ISEC-1-I ),
     $                                        A( MAX( II, I+1 ), I ), 1,
     $                       T1( I-JJ+1, MAX( II-JJ+1, I-JJ+2 ) ), RCB )
  140             CONTINUE
  150          CONTINUE
*
*              C := alpha*B*T1 + beta*C, a vertical block of C is
*              updated using the general matrix multiply, DGEMM. T1
*              corresponds to a full diagonal block of the matrix A.
*
               CALL DGEMM ( 'N', 'N', M, JSEC, JSEC, ALPHA,
     $                                 B( 1, JJ ), LDB, T1( 1, 1 ), RCB,
     $                                           BETA, C( 1, JJ ), LDC )
*
*              C := alpha*B*A + C and C := alpha*B*A' + C, general
*              matrix multiply operations involving rectangular blocks
*              of A.
*
               IF( JJ+JSEC.LE.N )THEN
                  CALL DGEMM ( 'N', 'N', M, JSEC, N-JJ-JSEC+1, ALPHA,
     $                           B( 1, JJ+JSEC ), LDB, A( JJ+JSEC, JJ ),
     $                                       LDA, ONE, C( 1, JJ ), LDC )
               END IF
               IF( JJ.GT.1 )THEN
                  CALL DGEMM ( 'N', 'T', M, JSEC, JJ-1, ALPHA,
     $                                  B( 1, 1 ), LDB, A( JJ, 1 ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
               END IF
  160       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYMM.
*
      END
C     ******************************************************************
      SUBROUTINE DSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR2K  performs one of the symmetric rank 2k operations
*
*     C := alpha*A*B' + alpha*B*A' + beta*C,
*
*  or
*
*     C := alpha*A'*B + alpha*B'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*  matrices in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
*                                        beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrices  A and B.  K must be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA
      INTEGER            I, II, IX, ISEC, JJ, JX, JSEC
      LOGICAL            UPPER, NOTR
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           DGEMM, DAXPY, DSCAL
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. User specified parameters for DSYR2K ..
      INTEGER            RCB, CB
      PARAMETER        ( RCB = 96, CB = 48 )
*     .. Local Arrays ..
      DOUBLE PRECISION   T1( RCB, RCB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      UPPER = UPLO .EQ. 'U'
      NOTR = TRANS .EQ. 'N'
      IF( NOTR )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      INFO = 0
      IF( ( .NOT.UPPER ).AND.( .NOT.( UPLO .EQ. 'L' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.( TRANS .EQ. 'T' ) ).AND.
     $                                ( .NOT.( TRANS .EQ. 'C' ) ) )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( K.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2K', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when alpha.eq.zero or k.eq.0.
*
      IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
         IF( UPPER )THEN
            DO 10, I = 1, N
               CALL DSCAL ( I, BETA, C( 1, I ), 1 )
   10       CONTINUE
         ELSE
            DO 20, I = 1, N
               CALL DSCAL ( N-I+1, BETA, C( I, I ), 1 )
   20       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( UPPER )THEN
         IF( NOTR )THEN
*
*           Form  C := alpha*A*B' + alpha*B*A' + beta*C. Upper, Notr.
*
            DO 70, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
*
*              T1 := alpha*A*B', general matrix multiply on rectangular
*              blocks of A and B. T1 is square.
*
               CALL DGEMM ( 'N', 'T', ISEC, ISEC, K, ALPHA, A( II, 1 ),
     $                     LDA, B( II, 1 ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C :=  beta*C, a upper triangular diagonal block of C is
*              updated with beta.
*
               IF( BETA.NE.ONE )THEN
                  DO 30, I = II, II+ISEC-1
                     CALL DSCAL ( I-II+1, BETA, C( II, I ), 1 )
   30             CONTINUE
               END IF
*
*              C := T1 + C, the upper triangular part of T1 is added to
*              the upper triangular diagonal block of C.
*
               DO 40, I = II, II+ISEC-1
                  CALL DAXPY ( I-II+1, ONE, T1( 1, I-II+1 ), 1,
     $                                                   C( II, I ), 1 )
   40          CONTINUE
*
*              C := T1' + C, the transpose of the lower triangular part
*              of T1 is added to the upper triangular diagonal block
*              of C. Notice that T1 is referenced by row and that the
*              maximum length of a vector referenced by DAXPY is CB.
*
               DO 60, JJ = II, II+ISEC-1, CB
                  JSEC = MIN( CB, II+ISEC-JJ )
                  DO 50, I = JJ, II+ISEC-1
                     CALL DAXPY ( MIN( JSEC, I-JJ+1 ), ONE,
     $                       T1( I-II+1, JJ-II+1 ), RCB, C( JJ, I ), 1 )
   50             CONTINUE
   60          CONTINUE
*
*              C := alpha*A*B' + beta*C  and  C := alpha*B*A' + C,
*              general matrix multiply on upper vertical blocks of C.
*
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                            A( 1, 1 ), LDA, B( II, 1 ), LDB, BETA,
     $                                                 C( 1, II ), LDC )
                  CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                             B( 1, 1 ), LDB, A( II, 1 ), LDA, ONE,
     $                                                 C( 1, II ), LDC )
               END IF
   70       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + alpha*B'*A + beta*C. Upper, Trans.
*
            DO 120, II = 1, N, RCB
               ISEC = MIN( RCB, N-II+1 )
*
*              T1 := alpha*A'*B, general matrix multiply on rectangular
*              blocks of A and B. T1 is square.
*
               CALL DGEMM ( 'T', 'N', ISEC, ISEC, K, ALPHA, A( 1, II ),
     $                     LDA, B( 1, II ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C :=  beta*C, a upper triangular diagonal block of C is
*              updated with beta.
*
               IF( BETA.NE.ONE )THEN
                  DO 80, I = II, II+ISEC-1
                     CALL DSCAL ( I-II+1, BETA, C( II, I ), 1 )
   80             CONTINUE
               END IF
*
*              C := T1 + C, the upper triangular part of T1 is added to
*              the upper triangular diagonal block of C.
*
               DO 90, I = II, II+ISEC-1
                  CALL DAXPY ( I-II+1, ONE, T1( 1, I-II+1 ), 1,
     $                                                   C( II, I ), 1 )
   90          CONTINUE
*
*              C := T1' + C, the transpose of the lower triangular part
*              of T1 is added to the upper triangular diagonal block
*              of C. Notice that T1 is referenced by row and that the
*              maximum length of a vector referenced by DAXPY is CB.
*
               DO 110, JJ = II, II+ISEC-1, CB
                  JSEC = MIN( CB, II+ISEC-JJ )
                  DO 100, I = JJ, II+ISEC-1
                     CALL DAXPY ( MIN( JSEC, I-JJ+1 ), ONE,
     $                       T1( I-II+1, JJ-II+1 ), RCB, C( JJ, I ), 1 )
  100             CONTINUE
  110          CONTINUE
*
*              C := alpha*A'*B + beta*C  and  C := alpha*B'*A + C,
*              general matrix multiply on upper vertical blocks of C.
*
               IF( II.GT.1 )THEN
                  CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                            A( 1, 1 ), LDA, B( 1, II ), LDB, BETA,
     $                                                 C( 1, II ), LDC )
                  CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                             B( 1, 1 ), LDB, A( 1, II ), LDA, ONE,
     $                                                 C( 1, II ), LDC )
               END IF
  120       CONTINUE
         END IF
      ELSE
         IF( NOTR )THEN
*
*           Form  C := alpha*A*B' + alpha*B*A' + beta*C. Lower, Notr.
*
            DO 170, IX = N, 1, -RCB
               II = MAX( 1, IX-RCB+1 )
               ISEC = IX-II+1
*
*              T1 := alpha*A*B', general matrix multiply on rectangular
*              blocks of A and B. T1 is square.
*
               CALL DGEMM ( 'N', 'T', ISEC, ISEC, K, ALPHA, A( II, 1 ),
     $                     LDA, B( II, 1 ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C :=  beta*C, a lower triangular diagonal block of C is
*              updated with beta.
*
               IF( BETA.NE.ONE )THEN
                  DO 130, I = II, II+ISEC-1
                     CALL DSCAL ( II+ISEC-I, BETA, C( I, I ), 1 )
  130             CONTINUE
               END IF
*
*              C := T1 + C, the lower triangular part of T1 is added to
*              the lower triangular diagonal block of C.
*
               DO 140, I = II, II+ISEC-1
                  CALL DAXPY ( II+ISEC-I, ONE, T1( I-II+1, I-II+1 ), 1,
     $                                                    C( I, I ), 1 )
  140          CONTINUE
*
*              C := T1' + C, the transpose of the upper triangular part
*              of T1 is added to the lower triangular diagonal block
*              of C. Notice that T1 is referenced by row and that the
*              maximum length of a vector referenced by DAXPY is CB.
*
               DO 160, JX = II+ISEC-1, II, -CB
                  JJ = MAX( II, JX-CB+1 )
                  JSEC = JX-JJ+1
                  DO 150, I = II, JJ+JSEC-1
                     CALL DAXPY ( MIN( JSEC, JJ+JSEC-I ), ONE,
     $                        T1( I-II+1, MAX( JJ-II+1, I-II+1 ) ), RCB,
     $                                         C( MAX( JJ, I ), I ), 1 )
  150             CONTINUE
  160          CONTINUE
*
*              C := alpha*A*B' + beta*C  and  C := alpha*B*A' + C,
*              general matrix multiply on lower vertical blocks of C.
*
               IF( II+ISEC.LE.N )THEN
                  CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K, ALPHA,
     $                            A( II+ISEC, 1 ), LDA, B( II, 1 ), LDB,
     $                                     BETA, C( II+ISEC, II ), LDC )
                  CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K, ALPHA,
     $                            B( II+ISEC, 1 ), LDB, A( II, 1 ), LDA,
     $                                      ONE, C( II+ISEC, II ), LDC )
               END IF
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + alpha*B'*A + beta*C. Lower, Trans.
*
            DO 220, IX = N, 1, -RCB
               II = MAX( 1, IX-RCB+1 )
               ISEC = IX-II+1
*
*              T1 := alpha*A*B', general matrix multiply on rectangular
*              blocks of A and B. T1 is square.
*
               CALL DGEMM ( 'T', 'N', ISEC, ISEC, K, ALPHA, A( 1, II ),
     $                     LDA, B( 1, II ), LDB, ZERO, T1( 1, 1 ), RCB )
*
*              C :=  beta*C, a lower triangular diagonal block of C is
*              updated with beta.
*
               IF( BETA.NE.ONE )THEN
                  DO 180, I = II, II+ISEC-1
                     CALL DSCAL ( II+ISEC-I, BETA, C( I, I ), 1 )
  180             CONTINUE
               END IF
*
*              C := T1 + C, the lower triangular part of T1 is added to
*              the lower triangular diagonal block of C.
*
               DO 190, I = II, II+ISEC-1
                  CALL DAXPY ( II+ISEC-I, ONE, T1( I-II+1, I-II+1 ), 1,
     $                                                    C( I, I ), 1 )
  190          CONTINUE
*
*              C := T1' + C, the transpose of the upper triangular part
*              of T1 is added to the lower triangular diagonal block
*              of C. Notice that T1 is referenced by row and that the
*              maximum length of a vector referenced by DAXPY is CB.
*
               DO 210, JX = II+ISEC-1, II, -CB
                  JJ = MAX( II, JX-CB+1 )
                  JSEC = JX-JJ+1
                  DO 200, I = II, JJ+JSEC-1
                     CALL DAXPY ( MIN( JSEC, JJ+JSEC-I ), ONE,
     $                        T1( I-II+1, MAX( JJ-II+1, I-II+1 ) ), RCB,
     $                                         C( MAX( JJ, I ), I ), 1 )
  200             CONTINUE
  210          CONTINUE
*
*              C := alpha*A'*B + beta*C  and  C := alpha*B'*A + C,
*              general matrix multiply on lower vertical blocks of C.
*
               IF( II+ISEC.LE.N )THEN
                  CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K, ALPHA,
     $                            A( 1, II+ISEC ), LDA, B( 1, II ), LDB,
     $                                     BETA, C( II+ISEC, II ), LDC )
                  CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K, ALPHA,
     $                            B( 1, II+ISEC ), LDB, A( 1, II ), LDA,
     $                                      ONE, C( II+ISEC, II ), LDC )
               END IF
  220       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYR2K.
*
      END
C     ******************************************************************
      SUBROUTINE DTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDC
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*C,   or   X*op( A ) = alpha*C,
*
*  where alpha is a scalar, X and C are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on C.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*C.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*C.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of C. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  C need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry,  the leading  m by n part of the array  C must
*           contain  the  right-hand  side  matrix  C,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA, OFFD
      LOGICAL            LSIDE, UPPER, NOTR, NOUNIT, CLDC, SMALLN,
     $                   TINYN, TINYM
      INTEGER            I, II, IX, ISEC, J, JJ, JX, JSEC, TSEC
      DOUBLE PRECISION   GAMMA, DELTA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     .. External Functions ..
      LOGICAL            LSAME, DBIGP, DCLD
      EXTERNAL           LSAME, DBIGP, DCLD
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           DGEMM, DGEMV, DTRSV, DCOPY
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      INTEGER            DIP91, DIP92, DIP93
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   DIP91 = 91, DIP92 = 92, DIP93 = 93 )
*     .. User specified parameters for DTRSM ..
      INTEGER            RB, CB, RCB
      PARAMETER        ( RCB = 48, RB = 48, CB = 48 )
      DOUBLE PRECISION   T1( RB, CB ), T2( CB, CB ), T3( RCB, RCB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE = SIDE  .EQ. 'L'
      UPPER = UPLO  .EQ. 'U'
      NOTR = TRANSA .EQ. 'N'
      NOUNIT = DIAG .EQ. 'N'
      IF( NOUNIT )THEN
         OFFD = 0
      ELSE
         OFFD = 1
      END IF
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.( SIDE .EQ. 'R' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER ).AND.( .NOT.( UPLO .EQ. 'L' ) ) )THEN
         INFO = 2
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.( TRANSA .EQ. 'T' ) ).AND.
     $                               ( .NOT.( TRANSA .EQ. 'C' ) ) )THEN
         INFO = 3
      ELSE IF( ( .NOT.NOUNIT ).AND.( .NOT.( DIAG .EQ. 'U' ) ) )THEN
         INFO = 4
      ELSE IF( M.LT.0 )THEN
         INFO = 5
      ELSE IF( N.LT.0 )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ) )
     $   RETURN
*
*     And when alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         CALL DGEMM ( 'N', 'N', M, N, 0, ZERO, C, MAX( LDA, LDC ), C,
     $                                   MAX( LDA, LDC ), ZERO, C, LDC )
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( UPPER )THEN
            IF( NOTR )THEN
*
*              Solve  A*X = alpha*C. Left, Upper, No transpose.
*
               SMALLN = .NOT.DBIGP( DIP91, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP92, M, N )
                  DO 40, II = M-MOD( M-1, RCB ), 1, -RCB
                     ISEC = MIN( RCB, M-II+1 )
*
*                    C := -1*A*C + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1, -ONE,
     $                      A( II, II+ISEC ), LDA, C( II+ISEC, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       Solve A*X = C, triangular system solve involving
*                       a upper triangular diagonal block of A. The
*                       block of X is overwritten on C.
*
                        DO 10, J = 1, N
                           CALL DTRSV ( 'U', 'N', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
   10                   CONTINUE
                     ELSE
*
*                       T3 := A, a upper unit or non-unit triangular
*                       diagonal block of A is copied to the upper
*                       triangular part of T3.
*
                        DO 20, I = II+OFFD, II+ISEC-1
                           CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                              T3( 1, I-II+1 ), 1 )
   20                   CONTINUE
*
*                       Solve T3*X = C, triangular system solve
*                       involving a upper triangular diagonal block of A
*                       stored in T3. The block of X is overwritten
*                       on C.
*
                        DO 30, J = 1, N
                           CALL DTRSV ( 'U', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
   30                   CONTINUE
                     END IF
   40             CONTINUE
               ELSE
                  DELTA = ONE
                  CLDC = DCLD( LDC )
                  DO 110, II = M-MOD( M-1, CB ), 1, -CB
                     ISEC = MIN( CB, M-II+1 )
*
*                    C := -1*A*C + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1, -ONE,
     $                      A( II, II+ISEC ), LDA, C( II+ISEC, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
*
*                    T2 := A', the transpose of a upper unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    lower triangular part of T2.
*
                     DO 50, I = II+OFFD, II+ISEC-1
                        CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                             T2( I-II+1, 1 ), CB )
   50                CONTINUE
                     DO 100, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 60, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
   60                      CONTINUE
                        ELSE
                           DO 70, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
   70                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*T2 + delta*T1, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether T2 stores a unit or non-unit
*                       triangular block. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 80, I = II+ISEC-1, II, -1
                           IF( NOUNIT )THEN
                              DELTA = ONE/T2( I-II+1, I-II+1 )
                           END IF
                           GAMMA = -DELTA
                           TSEC = II+ISEC-1-I
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                        T1( 1, I-II+2 ), RB, T2( I-II+2, I-II+1 ),
     $                                    1, DELTA, T1( 1, I-II+1 ), 1 )
   80                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 90, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
   90                   CONTINUE
  100                CONTINUE
  110             CONTINUE
               END IF
            ELSE
*
*              Solve  A'*X = alpha*C. Left, Upper, Transpose.
*
               SMALLN = .NOT.DBIGP( DIP91, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP92, M, N )
                  DO 150, II = 1, M, RCB
                     ISEC = MIN( RCB, M-II+1 )
*
*                    C := -1*A'*C + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'T', 'N', ISEC, N, II-1, -ONE,
     $                                  A( 1, II ), LDA, C( 1, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       Solve A'*X = C, triangular system solve
*                       involving the transpose of a upper triangular
*                       diagonal block of A. The block of X is
*                       overwritten on C.
*
                        DO 120, J = 1, N
                           CALL DTRSV ( 'U', 'T', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  120                   CONTINUE
                     ELSE
*
*                       T3 :=  A', the transpose of a upper unit or
*                       non-unit triangular diagonal block of A is
*                       copied to the lower triangular part of T3.
*
                        DO 130, I = II+OFFD, II+ISEC-1
                           CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                            T3( I-II+1, 1 ), RCB )
  130                   CONTINUE
*
*                       Solve T3*X = C, triangular system solve
*                       involving the transpose of a upper triangular
*                       diagonal block of A stored in T3. The block of X
*                       is overwritten on C.
*
                        DO 140, J = 1, N
                           CALL DTRSV ( 'L', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  140                   CONTINUE
                     END IF
  150             CONTINUE
               ELSE
                  DELTA = ONE
                  CLDC = DCLD( LDC )
                  DO 210, II = 1, M, CB
                     ISEC = MIN( CB, M-II+1 )
*
*                    C := -1*A'*C + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'T', 'N', ISEC, N, II-1, -ONE,
     $                                  A( 1, II ), LDA, C( 1, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     DO 200, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 160, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  160                      CONTINUE
                        ELSE
                           DO 170, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  170                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*A + delta*T1, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether A is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 180, I = II, II+ISEC-1
                           IF( NOUNIT )THEN
                              DELTA = ONE/A( I, I )
                           END IF
                           GAMMA = -DELTA
                           TSEC = I-II
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                                    T1( 1, 1 ), RB, A( II, I ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  180                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 190, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  190                   CONTINUE
  200                CONTINUE
  210             CONTINUE
               END IF
            END IF
         ELSE
            IF( NOTR )THEN
*
*              Solve  A*X = alpha*C. Left, Lower, No transpose.
*
               SMALLN = .NOT.DBIGP( DIP91, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP92, M, N )
                  DO 250, IX = MOD( M-1, RCB )+1, M, RCB
                     II = MAX( 1, IX-RCB+1 )
                     ISEC = IX-II+1
*
*                    C := -1*A*C + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', ISEC, N, II-1, -ONE,
     $                                  A( II, 1 ), LDA, C( 1, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       Solve A*X = C, triangular system solve involving
*                       a lower triangular diagonal block of A. The
*                       block of X is overwritten on C.
*
                        DO 220, J = 1, N
                           CALL DTRSV ( 'L', 'N', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  220                   CONTINUE
                     ELSE
*
*                       T3 := A, a lower unit or non-unit triangular
*                       diagonal block of A is copied to the lower
*                       triangular part of T3. The block of X is
*                       overwritten on C.
*
                        DO 230, I = II, II+ISEC-1-OFFD
                           CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                                 1, T3( I-II+1+OFFD, I-II+1 ), 1 )
  230                   CONTINUE
*
*                       Solve T3*X = C, triangular system solve
*                       involving a lower triangular diagonal block of A
*                       stored in T3. The block of X is overwritten
*                       on C.
*
                        DO 240, J = 1, N
                           CALL DTRSV ( 'L', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  240                   CONTINUE
                     END IF
  250             CONTINUE
               ELSE
                  DELTA = ONE
                  CLDC = DCLD( LDC )
                  DO 320, IX = MOD( M-1, CB )+1, M, CB
                     II = MAX( 1, IX-CB+1 )
                     ISEC = IX-II+1
*
*                    C := -1*A*C + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', ISEC, N, II-1, -ONE,
     $                                  A( II, 1 ), LDA, C( 1, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
*
*                    T2 := A', the transpose of a lower unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    upper triangular part of T2.
*
                     DO 260, I = II, II+ISEC-1-OFFD
                        CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                                1, T2( I-II+1, I-II+1+OFFD ), CB )
  260                CONTINUE
                     DO 310, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 270, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  270                      CONTINUE
                        ELSE
                           DO 280, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  280                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*T2 + delta*T1, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether T2 stores a unit or non-unit
*                       triangular block. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 290, I = II, II+ISEC-1
                           IF( NOUNIT )THEN
                              DELTA = ONE/T2( I-II+1, I-II+1 )
                           END IF
                           GAMMA = -DELTA
                           TSEC = I-II
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                               T1( 1, 1 ), RB, T2( 1, I-II+1 ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  290                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 300, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  300                   CONTINUE
  310                CONTINUE
  320             CONTINUE
               END IF
            ELSE
*
*              Solve  A'*X = alpha*C. Left, Lower, Transpose.
*
               SMALLN = .NOT.DBIGP( DIP91, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP92, M, N )
                  DO 360, IX = M, 1, -RCB
                     II = MAX( 1, IX-RCB+1 )
                     ISEC = IX-II+1
*
*                    C := -1*A'*C + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1, -ONE,
     $                      A( II+ISEC, II ), LDA, C( II+ISEC, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       Solve A'*X = C, triangular system solve
*                       involving the transpose of a lower triangular
*                       diagonal block of A. The block of X is
*                       overwritten on C.
*
                        DO 330, J = 1, N
                           CALL DTRSV ( 'L', 'T', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  330                   CONTINUE
                     ELSE
*
*                       T3 :=  A', the transpose of a lower unit or
*                       non-unit triangular diagonal block of A is
*                       copied to the upper triangular part of T3.
*
                        DO 340, I = II, II+ISEC-1-OFFD
                           CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                               1, T3( I-II+1, I-II+1+OFFD ), RCB )
  340                   CONTINUE
*
*                       Solve T3*X = C, triangular system solve
*                       involving the transpose of a lower triangular
*                       diagonal block of A stored in T3. The block of X
*                       is overwritten on C.
*
                        DO 350, J = 1, N
                           CALL DTRSV ( 'U', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  350                   CONTINUE
                     END IF
  360             CONTINUE
               ELSE
                  DELTA = ONE
                  CLDC = DCLD( LDC )
                  DO 420, IX = M, 1, -CB
                     II = MAX( 1, IX-CB+1 )
                     ISEC = IX-II+1
*
*                    C := -1*A'*C + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1, -ONE,
     $                      A( II+ISEC, II ), LDA, C( II+ISEC, 1 ), LDC,
     $                                          ALPHA, C( II, 1 ), LDC )
                     DO 410, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 370, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  370                      CONTINUE
                        ELSE
                           DO 380, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  380                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*A + delta*T1, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether A is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 390, I = II+ISEC-1, II, -1
                           IF( NOUNIT )THEN
                              DELTA = ONE/A( I, I )
                           END IF
                           GAMMA = -DELTA
                           TSEC = II+ISEC-1-I
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                              T1( 1, I-II+2 ), RB, A( I+1, I ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  390                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 400, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  400                   CONTINUE
  410                CONTINUE
  420             CONTINUE
               END IF
            END IF
         END IF
      ELSE
         IF( UPPER )THEN
            IF( NOTR )THEN
*
*              Solve  X*A = alpha*C. Right, Upper, No transpose.
*
               TINYM = .NOT.DBIGP( DIP93, M, N )
               IF( TINYM )THEN
                  DO 440, JJ = 1, N, RCB
                     JSEC = MIN( RCB, N-JJ+1 )
*
*                    C := -1*C*A + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', M, JSEC, JJ-1, -ONE,
     $                                  C( 1, 1 ), LDC, A( 1, JJ ), LDA,
     $                                          ALPHA, C( 1, JJ ), LDC )
*
*                    Solve X*A = C, triangular system solve involving
*                    a upper triangular diagonal block of A. The block
*                    of X is overwritten on C.
*
                     DO 430, I = 1, M
                        CALL DTRSV ( 'U', 'T', DIAG, JSEC, A( JJ, JJ ),
     $                                            LDA, C( I, JJ ), LDC )
  430                CONTINUE
  440             CONTINUE
               ELSE
                  DELTA = ONE
                  DO 490, JJ = 1, N, CB
                     JSEC = MIN( CB, N-JJ+1 )
*
*                    C := -1*C*A + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', M, JSEC, JJ-1, -ONE,
     $                                  C( 1, 1 ), LDC, A( 1, JJ ), LDA,
     $                                          ALPHA, C( 1, JJ ), LDC )
                     DO 480, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 450, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  450                   CONTINUE
*
*                       C := gamma*T1*A + delta*C, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether A is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 460, J =  JJ, JJ+JSEC-1
                           IF( NOUNIT )THEN
                              DELTA = ONE/A( J, J )
                           END IF
                           GAMMA = -DELTA
                           TSEC = J-JJ
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                                    T1( 1, 1 ), RB, A( JJ, J ), 1,
     $                                       DELTA, T1( 1, J-JJ+1 ), 1 )
  460                   CONTINUE
*
*                       C := T1, T1 is copied back to C.
*
                        DO 470, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( 1, J-JJ+1 ), 1,
     $                                                   C( II, J ), 1 )
  470                   CONTINUE
  480                CONTINUE
  490             CONTINUE
               END IF
            ELSE
*
*              Solve  X*A' = alpha*C. Right, Upper, Transpose.
*
               TINYM = .NOT.DBIGP( DIP93, M, N )
               IF( TINYM )THEN
                  DO 510, JJ = N-MOD( N-1, RCB ), 1, -RCB
                     JSEC = MIN( RCB, N-JJ+1 )
*
*                    C := -1*C*A' + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'N', 'T', M, JSEC, N-JJ-JSEC+1, -ONE,
     $                           C( 1, JJ+JSEC ), LDC, A( JJ, JJ+JSEC ),
     $                                     LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    Solve X*A' = C, triangular system solve involving
*                    the transpose of a upper triangular diagonal block
*                    of A. The block of X is overwritten on C.
*
                     DO 500, I = 1, M
                        CALL DTRSV ( 'U', 'N', DIAG, JSEC, A( JJ, JJ ),
     $                                            LDA, C( I, JJ ), LDC )
  500                CONTINUE
  510             CONTINUE
               ELSE
                  DELTA = ONE
                  DO 570, JJ = N-MOD( N-1, CB ), 1, -CB
                     JSEC = MIN( CB, N-JJ+1 )
*
*                    C := -1*C*A' + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'N', 'T', M, JSEC, N-JJ-JSEC+1, -ONE,
     $                           C( 1, JJ+JSEC ), LDC, A( JJ, JJ+JSEC ),
     $                                     LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    T2 := A', the transpose of a upper unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    lower triangular part of T2.
*
                     DO 520, J = JJ+OFFD, JJ+JSEC-1
                        CALL DCOPY ( J-JJ+1-OFFD, A( JJ, J ), 1,
     $                                             T2( J-JJ+1, 1 ), CB )
  520                CONTINUE
                     DO 560, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 530, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  530                   CONTINUE
*
*                       C := gamma*T1*T2 + delta*C, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether T2 is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 540, J = JJ+JSEC-1, JJ, -1
                           IF( NOUNIT )THEN
                              DELTA = ONE/T2( J-JJ+1, J-JJ+1 )
                           END IF
                           GAMMA = -DELTA
                           TSEC = JJ+JSEC-1-J
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                        T1( 1, J-JJ+2 ), RB, T2( J-JJ+2, J-JJ+1 ),
     $                                    1, DELTA, T1( 1, J-JJ+1 ), 1 )
  540                   CONTINUE
*
*                       C := T1, T1 is copied back to C.
*
                        DO 550, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( 1, J-JJ+1 ), 1,
     $                                                   C( II, J ), 1 )
  550                   CONTINUE
  560                CONTINUE
  570             CONTINUE
               END IF
            END IF
         ELSE
            IF( NOTR )THEN
*
*              Solve  X*A = alpha*C. Right, Lower, No transpose.
*
               TINYM = .NOT.DBIGP( DIP93, M, N )
               IF( TINYM )THEN
                  DO 590, JX = N, 1, -RCB
                     JJ = MAX( 1, JX-RCB+1 )
                     JSEC = JX-JJ+1
*
*                    C := -1*C*A + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', M, JSEC, N-JJ-JSEC+1, -ONE,
     $                           C( 1, JJ+JSEC ), LDC, A( JJ+JSEC, JJ ),
     $                                     LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    Solve X*A = C, triangular system solve involving
*                    a lower triangular diagonal block of A. The block
*                    of X is overwritten on C.
*
                     DO 580, I = 1, M
                        CALL DTRSV ( 'L', 'T', DIAG, JSEC, A( JJ, JJ ),
     $                                            LDA, C( I, JJ ), LDC )
  580                CONTINUE
  590             CONTINUE
               ELSE
                  DELTA = ONE
                  DO 640, JX = N, 1, -CB
                     JJ = MAX( 1, JX-CB+1 )
                     JSEC = JX-JJ+1
*
*                    C := -1*C*A + alpha*C, general matrix multiply
*                    involving a rectangular block of A.
*
                     CALL DGEMM ( 'N', 'N', M, JSEC, N-JJ-JSEC+1, -ONE,
     $                           C( 1, JJ+JSEC ), LDC, A( JJ+JSEC, JJ ),
     $                                     LDA, ALPHA, C( 1, JJ ), LDC )
                     DO 630, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 600, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  600                   CONTINUE
*
*                       C := gamma*T1*A + delta*C, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether A is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 610, J = JJ+JSEC-1, JJ, -1
                           IF( NOUNIT )THEN
                              DELTA = ONE/A( J, J )
                           END IF
                           GAMMA = -DELTA
                           TSEC = JJ+JSEC-1-J
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                              T1( 1, J-JJ+2 ), RB, A( J+1, J ), 1,
     $                                       DELTA, T1( 1, J-JJ+1 ), 1 )
  610                   CONTINUE
*
*                       C := T1, T1 is copied back to C.
*
                        DO 620, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( 1, J-JJ+1 ), 1,
     $                                                   C( II, J ), 1 )
  620                   CONTINUE
  630                CONTINUE
  640             CONTINUE
               END IF
            ELSE
*
*              Solve  X*A' = alpha*C. Right, Lower, Transpose.
*
               TINYM = .NOT.DBIGP( DIP93, M, N )
               IF( TINYM )THEN
                  DO 660, JX = MOD( N-1, RCB )+1, N, RCB
                     JJ = MAX( 1, JX-RCB+1 )
                     JSEC = JX-JJ+1
*
*                    C := -1*C*A' + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'N', 'T', M, JSEC, JJ-1, -ONE,
     $                                  C( 1, 1 ), LDC, A( JJ, 1 ), LDA,
     $                                          ALPHA, C( 1, JJ ), LDC )
*
*                    Solve X*A' = C, triangular system solve involving
*                    the transpose of a lower triangular diagonal block
*                    of A. The block of X is overwritten on C.
*
                     DO 650, I = 1, M
                        CALL DTRSV ( 'L', 'N', DIAG, JSEC, A( JJ, JJ ),
     $                                            LDA, C( I, JJ ), LDC )
  650                CONTINUE
  660             CONTINUE
               ELSE
                  DELTA = ONE
                  DO 720, JX = MOD( N-1, CB )+1, N, CB
                     JJ = MAX( 1, JX-CB+1 )
                     JSEC = JX-JJ+1
*
*                    C := -1*C*A' + alpha*C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     CALL DGEMM ( 'N', 'T', M, JSEC, JJ-1, -ONE,
     $                                  C( 1, 1 ), LDC, A( JJ, 1 ), LDA,
     $                                          ALPHA, C( 1, JJ ), LDC )
*
*                    T2 := A', the transpose of a lower unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    upper triangular part of T2.
*
                     DO 670, J = JJ, JJ+JSEC-1-OFFD
                        CALL DCOPY ( JJ+JSEC-J-OFFD, A( J+OFFD, J ),
     $                                1, T2( J-JJ+1, J-JJ+1+OFFD ), CB )
  670                CONTINUE
                     DO 710, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 680, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  680                   CONTINUE
*
*                       C := gamma*T1*T2 + delta*C, triangular matrix
*                       multiply where the values of gamma and delta
*                       depend on whether T2 is a unit or non-unit
*                       triangular matrix. Gamma and tsec are also used
*                       to compensate for a deficiency in DGEMV that
*                       appears if the second dimension (tsec) is zero.
*
                        DO 690, J = JJ, JJ+JSEC-1
                           IF( NOUNIT )THEN
                              DELTA = ONE/T2( J-JJ+1, J-JJ+1 )
                           END IF
                           GAMMA = -DELTA
                           TSEC = J-JJ
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                               T1( 1, 1 ), RB, T2( 1, J-JJ+1 ), 1,
     $                                       DELTA, T1( 1, J-JJ+1 ), 1 )
  690                   CONTINUE
*
*                       C := T1, T1 is copied back to C.
*
                        DO 700, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( 1, J-JJ+1 ), 1,
     $                                                   C( II, J ), 1 )
  700                   CONTINUE
  710                CONTINUE
  720             CONTINUE
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM.
*
      END
C     ******************************************************************
      SUBROUTINE DTRMM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDC
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*C,   or   C := alpha*C*op( A ),
*
*  where  alpha  is a scalar,  C  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies C from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*op( A )*C.
*
*              SIDE = 'R' or 'r'   C := alpha*C*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of C. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  C need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry,  the leading  m by n part of the array  C must
*           contain the matrix  C,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA, OFFD
      LOGICAL            LSIDE, UPPER, NOTR, NOUNIT, CLDC, SMALLN,
     $                   TINYN, TINYM
      INTEGER            I, II, IX, ISEC, J, JJ, JX, JSEC, TSEC
      DOUBLE PRECISION   GAMMA, DELTA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     .. External Functions ..
      LOGICAL            LSAME, DBIGP, DCLD
      EXTERNAL           LSAME, DBIGP, DCLD
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           DGEMM, DGEMV, DTRMV, DCOPY
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      INTEGER            DIP81, DIP82, DIP83
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   DIP81 = 81, DIP82 = 82, DIP83 = 83 )
*     .. User specified parameters for DTRMM ..
      INTEGER            RB, CB, RCB
      PARAMETER        ( RCB = 64, RB = 64, CB = 64 )
      DOUBLE PRECISION   T1( RB, CB ), T2( CB, CB ), T3( RCB, RCB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      NOTR = LSAME( TRANSA, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( NOUNIT )THEN
         OFFD = 0
      ELSE
         OFFD = 1
      END IF
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      INFO = 0
      IF( ( .NOT.LSIDE ).AND.( .NOT.LSAME( SIDE, 'R' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER ).AND.( .NOT.LSAME( UPLO, 'L' ) ) )THEN
         INFO = 2
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $                               ( .NOT.LSAME( TRANSA, 'C' ) ) )THEN
         INFO = 3
      ELSE IF( ( .NOT.NOUNIT ).AND.( .NOT.LSAME( DIAG, 'U' ) ) )THEN
         INFO = 4
      ELSE IF( M.LT.0 )THEN
         INFO = 5
      ELSE IF( N.LT.0 )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ) )
     $   RETURN
*
*     And when alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         CALL DGEMM ( 'N', 'N', M, N, 0, ZERO, C, MAX( LDA, LDC ), C,
     $                                   MAX( LDA, LDC ), ZERO, C, LDC )
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( UPPER )THEN
            IF( NOTR )THEN
*
*              Form  C := alpha*A*C. Left, Upper, No transpose.
*
               SMALLN = .NOT.DBIGP( DIP81, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP82, M, N )
                  DO 40, II = 1, M, RCB
                     ISEC = MIN( RCB, M-II+1 )
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'N', ISEC, N, 0, ZERO, A, LDA,
     $                                  C, LDC, ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       C := A*C, triangular matrix multiply involving
*                       a upper triangular diagonal block of A.
*
                        DO 10, J = 1, N
                           CALL DTRMV ( 'U', 'N', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
   10                   CONTINUE
                     ELSE
*
*                       T3 := A, a upper unit or non-unit triangular
*                       diagonal block of A is copied to the upper
*                       triangular part of T3.
*
                        DO 20, I = II+OFFD, II+ISEC-1
                           CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                              T3( 1, I-II+1 ), 1 )
   20                   CONTINUE
*
*                       C := T3*C, triangular matrix multiply involving
*                       a upper triangular diagonal block of A stored
*                       in T3.
*
                        DO 30, J = 1, N
                           CALL DTRMV ( 'U', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
   30                   CONTINUE
                     END IF
*
*                    C := alpha*A*C + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( II+ISEC.LE.M )THEN
                        CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1,
     $                                     ALPHA, A( II, II+ISEC ), LDA,
     $                                        C( II+ISEC, 1 ), LDC, ONE,
     $                                                 C( II, 1 ), LDC )
                     END IF
   40             CONTINUE
               ELSE
                  DELTA = ALPHA
                  CLDC = DCLD( LDC )
                  DO 110, II = 1, M, CB
                     ISEC = MIN( CB, M-II+1 )
*
*                    T2 := A', the transpose of a upper unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    lower triangular part of T2.
*
                     DO 50, I = II+OFFD, II+ISEC-1
                        CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                             T2( I-II+1, 1 ), CB )
   50                CONTINUE
                     DO 100, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 60, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
   60                      CONTINUE
                        ELSE
                           DO 70, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
   70                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*T2 + delta*T1, triangular matrix
*                       multiply where the value of delta depends on
*                       whether T2 stores a unit or non-unit triangular
*                       block. Gamma and tsec are used to compensate for
*                       a deficiency in DGEMV that appears if the second
*                       dimension (tsec) is zero.
*
                        DO 80, I = II, II+ISEC-1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*T2( I-II+1, I-II+1 )
                           END IF
                           GAMMA = ALPHA
                           TSEC = II+ISEC-1-I
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                        T1( 1, I-II+2 ), RB, T2( I-II+2, I-II+1 ),
     $                                    1, DELTA, T1( 1, I-II+1 ), 1 )
   80                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 90, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
   90                   CONTINUE
  100                CONTINUE
*
*                    C := alpha*A*C + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( II+ISEC.LE.M )THEN
                        CALL DGEMM ( 'N', 'N', ISEC, N, M-II-ISEC+1,
     $                                     ALPHA, A( II, II+ISEC ), LDA,
     $                                        C( II+ISEC, 1 ), LDC, ONE,
     $                                                 C( II, 1 ), LDC )
                     END IF
  110             CONTINUE
               END IF
            ELSE
*
*              Form  C := alpha*A'*C. Left, Upper, Transpose.
*
               SMALLN = .NOT.DBIGP( DIP81, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP82, M, N )
                  DO 150, II = M-MOD( M-1, RCB ), 1, -RCB
                     ISEC = MIN( RCB, M-II+1 )
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'T', 'N', ISEC, N, 0, ZERO, A, LDA,
     $                                  C, LDC, ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       C := A'*C, triangular matrix multiply involving
*                       a upper triangular diagonal block of A.
*
                        DO 120, J = 1, N
                           CALL DTRMV ( 'U', 'T', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  120                   CONTINUE
                     ELSE
*
*                       T3 :=  A', the transpose of a upper unit or
*                       non-unit triangular diagonal block of A is
*                       copied to the lower triangular part of T3.
*
                        DO 130, I = II+OFFD, II+ISEC-1
                           CALL DCOPY ( I-II+1-OFFD, A( II, I ), 1,
     $                                            T3( I-II+1, 1 ), RCB )
  130                   CONTINUE
*
*                       C := T3*C, triangular matrix multiply involving
*                       the transpose of a upper triangular diagonal
*                       block of A stored in T3.
*
                        DO 140, J = 1, N
                           CALL DTRMV ( 'L', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  140                   CONTINUE
                     END IF
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( II.GT.1 )THEN
                        CALL DGEMM ( 'T', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( 1, II ), LDA, C( 1, 1 ), LDC,
     $                                            ONE, C( II, 1 ), LDC )
                     END IF
  150             CONTINUE
               ELSE
                  DELTA = ALPHA
                  CLDC = DCLD( LDC )
                  DO 210, II = M-MOD( M-1, CB ), 1, -CB
                     ISEC = MIN( CB, M-II+1 )
                     DO 200, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 160, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  160                      CONTINUE
                        ELSE
                           DO 170, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  170                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*A + delta*T1, triangular matrix
*                       multiply where the value of delta depends on
*                       whether A is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 180, I = II+ISEC-1, II, -1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*A( I, I )
                           END IF
                           GAMMA = ALPHA
                           TSEC = I-II
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                                    T1( 1, 1 ), RB, A( II, I ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  180                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 190, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  190                   CONTINUE
  200                CONTINUE
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( II.GT.1 )THEN
                        CALL DGEMM ( 'T', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( 1, II ), LDA, C( 1, 1 ), LDC,
     $                                            ONE, C( II, 1 ), LDC )
                     END IF
  210             CONTINUE
               END IF
            END IF
         ELSE
            IF( NOTR )THEN
*
*              Form  C := alpha*A*C. Left, Lower, No transpose.
*
               SMALLN = .NOT.DBIGP( DIP81, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP82, M, N )
                  DO 250, IX = M, 1, -RCB
                     II = MAX( 1, IX-RCB+1 )
                     ISEC = IX-II+1
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'N', ISEC, N, 0, ZERO, A, LDA,
     $                                  C, LDC, ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       C := A*C, triangular matrix multiply involving
*                       a lower triangular diagonal block of A.
*
                        DO 220, J = 1, N
                           CALL DTRMV ( 'L', 'N', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  220                   CONTINUE
                     ELSE
*
*                       T3 := A, a lower unit or non-unit triangular
*                       diagonal block of A is copied to the lower
*                       triangular part of T3.
*
                        DO 230, I = II, II+ISEC-1-OFFD
                           CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                                 1, T3( I-II+1+OFFD, I-II+1 ), 1 )
  230                   CONTINUE
*
*                       C := T3*C, triangular matrix multiply involving
*                       a lower triangular diagonal block of A stored
*                       in T3.
*
                        DO 240, J = 1, N
                           CALL DTRMV ( 'L', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  240                   CONTINUE
                     END IF
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( II.GT.1 )THEN
                        CALL DGEMM ( 'N', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( II, 1 ), LDA, C( 1, 1 ), LDC,
     $                                            ONE, C( II, 1 ), LDC )
                     END IF
  250             CONTINUE
               ELSE
                  DELTA = ALPHA
                  CLDC = DCLD( LDC )
                  DO 320, IX = M, 1, -CB
                     II = MAX( 1, IX-CB+1 )
                     ISEC = IX-II+1
*
*                    T2 := A', the transpose of a lower unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    upper triangular part of T2.
*
                     DO 260, I = II, II+ISEC-1-OFFD
                        CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                                1, T2( I-II+1, I-II+1+OFFD ), CB )
  260                CONTINUE
                     DO 310, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 270, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  270                      CONTINUE
                        ELSE
                           DO 280, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  280                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*T2 + delta*T1, triangular matrix
*                       multiply where the value of delta depends on
*                       whether T2 stores a unit or non-unit triangular
*                       block. Gamma and tsec are used to compensate for
*                       a deficiency in DGEMV that appears if the second
*                       dimension (tsec) is zero.
*
                        DO 290, I = II+ISEC-1, II, -1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*T2( I-II+1, I-II+1 )
                           END IF
                           GAMMA = ALPHA
                           TSEC = I-II
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                               T1( 1, 1 ), RB, T2( 1, I-II+1 ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  290                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 300, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  300                   CONTINUE
  310                CONTINUE
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( II.GT.1 )THEN
                        CALL DGEMM ( 'N', 'N', ISEC, N, II-1, ALPHA,
     $                                  A( II, 1 ), LDA, C( 1, 1 ), LDC,
     $                                            ONE, C( II, 1 ), LDC )
                     END IF
  320             CONTINUE
               END IF
            ELSE
*
*              Form  C := alpha*A'*C. Left, Lower, Transpose.
*
               SMALLN = .NOT.DBIGP( DIP81, M, N )
               IF( SMALLN )THEN
                  TINYN = .NOT.DBIGP( DIP82, M, N )
                  DO 360, IX = MOD( M-1, RCB )+1, M, RCB
                     II = MAX( 1, IX-RCB+1 )
                     ISEC = IX-II+1
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'T', 'N', ISEC, N, 0, ZERO, A, LDA,
     $                                  C, LDC, ALPHA, C( II, 1 ), LDC )
                     IF( TINYN )THEN
*
*                       C := A'*C, triangular matrix multiply involving
*                       a lower triangular diagonal block of A.
*
                        DO 330, J = 1, N
                           CALL DTRMV ( 'L', 'T', DIAG, ISEC,
     $                                 A( II, II ), LDA, C( II, J ), 1 )
  330                   CONTINUE
                     ELSE
*
*                       T3 :=  A', the transpose of a lower unit or
*                       non-unit triangular diagonal block of A is
*                       copied to the upper triangular part of T3.
*
                        DO 340, I = II, II+ISEC-1-OFFD
                           CALL DCOPY ( II+ISEC-I-OFFD, A( I+OFFD, I ),
     $                               1, T3( I-II+1, I-II+1+OFFD ), RCB )
  340                   CONTINUE
*
*                       C := alpha*T3*C, triangular matrix multiply
*                       involving the transpose of a lower triangular
*                       diagonal block of A stored in T3.
*
                        DO 350, J = 1, N
                           CALL DTRMV ( 'U', 'N', DIAG, ISEC,
     $                                  T3( 1, 1 ), RCB, C( II, J ), 1 )
  350                   CONTINUE
                     END IF
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( II+ISEC.LE.M )THEN
                        CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1,
     $                                     ALPHA, A( II+ISEC, II ), LDA,
     $                                        C( II+ISEC, 1 ), LDC, ONE,
     $                                                 C( II, 1 ), LDC )
                     END IF
  360             CONTINUE
               ELSE
                  DELTA = ALPHA
                  CLDC = DCLD( LDC )
                  DO 420, IX = MOD( M-1, CB )+1, M, CB
                     II = MAX( 1, IX-CB+1 )
                     ISEC = IX-II+1
                     DO 410, JJ = 1, N, RB
                        JSEC = MIN( RB, N-JJ+1 )
*
*                       T1 := C', the transpose of a rectangular block
*                       of C is copied to T1.
*
                        IF( CLDC )THEN
                           DO 370, J = JJ, JJ+JSEC-1
                              CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                             T1( J-JJ+1, 1 ), RB )
  370                      CONTINUE
                        ELSE
                           DO 380, I = II, II+ISEC-1
                              CALL DCOPY ( JSEC, C( I, JJ ), LDC,
     $                                              T1( 1, I-II+1 ), 1 )
  380                      CONTINUE
                        END IF
*
*                       T1 := gamma*T1*A + delta*T1, triangular matrix
*                       multiply where the value of delta depends on
*                       whether A is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 390, I = II, II+ISEC-1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*A( I, I )
                           END IF
                           GAMMA = ALPHA
                           TSEC = II+ISEC-1-I
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', JSEC, TSEC, GAMMA,
     $                              T1( 1, I-II+2 ), RB, A( I+1, I ), 1,
     $                                       DELTA, T1( 1, I-II+1 ), 1 )
  390                   CONTINUE
*
*                       C := T1', the transpose of T1 is copied back
*                       to C.
*
                        DO 400, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, T1( J-JJ+1, 1 ), RB,
     $                                                   C( II, J ), 1 )
  400                   CONTINUE
  410                CONTINUE
*
*                    C := alpha*A'*C + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( II+ISEC.LE.M )THEN
                        CALL DGEMM ( 'T', 'N', ISEC, N, M-II-ISEC+1,
     $                                     ALPHA, A( II+ISEC, II ), LDA,
     $                                        C( II+ISEC, 1 ), LDC, ONE,
     $                                                 C( II, 1 ), LDC )
                     END IF
  420             CONTINUE
               END IF
            END IF
         END IF
      ELSE
         IF( UPPER )THEN
            IF( NOTR )THEN
*
*              Form  C := alpha*C*A. Right, Upper, No transpose.
*
               TINYM = .NOT.DBIGP( DIP83, M, N )
               IF( TINYM )THEN
                  DO 440, JJ = N-MOD( N-1, RCB ), 1, -RCB
                     JSEC = MIN( RCB, N-JJ+1 )
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'N', M, JSEC, 0, ZERO, C, LDC,
     $                                  A, LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    C := C*A, triangular matrix multiply involving a
*                    upper triangular diagonal block of A.
*
                     DO 430, I = 1, M
                        CALL DTRMV ( 'U', 'T', DIAG, JSEC,
     $                               A( JJ, JJ ), LDA, C( I, JJ ), LDC )
  430                CONTINUE
*
*                    C := alpha*C*A + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( JJ.GT.1 )THEN
                        CALL DGEMM ( 'N', 'N', M, JSEC, JJ-1, ALPHA,
     $                                  C( 1, 1 ), LDC, A( 1, JJ ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
                     END IF
  440             CONTINUE
               ELSE
                  DELTA = ALPHA
                  DO 480, JJ = N-MOD( N-1, CB ), 1, -CB
                     JSEC = MIN( CB, N-JJ+1 )
                     DO 470, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 450, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  450                   CONTINUE
*
*                       C := gamma*T1*A + delta*C, triangular matrix
*                       multiply where the value of delta depends on
*                       whether A is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 460, J = JJ+JSEC-1, JJ, -1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*A( J, J )
                           END IF
                           GAMMA = ALPHA
                           TSEC = J-JJ
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                                    T1( 1, 1 ), RB, A( JJ, J ), 1,
     $                                            DELTA, C( II, J ), 1 )
  460                   CONTINUE
  470                CONTINUE
*
*                    C := alpha*C*A + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( JJ.GT.1 )THEN
                        CALL DGEMM ( 'N', 'N', M, JSEC, JJ-1, ALPHA,
     $                                  C( 1, 1 ), LDC, A( 1, JJ ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
                     END IF
  480             CONTINUE
               END IF
            ELSE
*
*              Form  C := alpha*C*A'. Right, Upper, Transpose.
*
               TINYM = .NOT.DBIGP( DIP83, M, N )
               IF( TINYM )THEN
                  DO 500, JJ = 1, N, RCB
                     JSEC = MIN( RCB, N-JJ+1 )
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'T', M, JSEC, 0, ZERO, C, LDC,
     $                                  A, LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    C := C*A', triangular matrix multiply involving a
*                    upper triangular diagonal block of A.
*
                     DO 490, I = 1, M
                        CALL DTRMV ( 'U', 'N', DIAG, JSEC,
     $                               A( JJ, JJ ), LDA, C( I, JJ ), LDC )
  490                CONTINUE
*
*                    C := alpha*C*A' + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( JJ+JSEC.LE.N )THEN
                        CALL DGEMM ( 'N', 'T', M, JSEC, N-JJ-JSEC+1,
     $                                      ALPHA, C( 1, JJ+JSEC ), LDC,
     $                                       A( JJ, JJ+JSEC ), LDA, ONE,
     $                                                 C( 1, JJ ), LDC )
                     END IF
  500             CONTINUE
               ELSE
                  DELTA = ALPHA
                  DO 550, JJ = 1, N, CB
                     JSEC = MIN( CB, N-JJ+1 )
*
*                    T2 := A', the transpose of a upper unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    lower triangular part of T2.
*
                     DO 510, J = JJ+OFFD, JJ+JSEC-1
                        CALL DCOPY ( J-JJ+1-OFFD, A( JJ, J ), 1,
     $                                             T2( J-JJ+1, 1 ), CB )
  510                CONTINUE
                     DO 540, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 520, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  520                   CONTINUE
*
*                       C := gamma*T1*T2 + delta*C, triangular matrix
*                       multiply where the value of delta depends on
*                       whether T2 is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 530, J = JJ, JJ+JSEC-1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*T2( J-JJ+1, J-JJ+1 )
                           END IF
                           GAMMA = ALPHA
                           TSEC = JJ+JSEC-1-J
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                        T1( 1, J-JJ+2 ), RB, T2( J-JJ+2, J-JJ+1 ),
     $                                         1, DELTA, C( II, J ), 1 )
  530                   CONTINUE
  540                CONTINUE
*
*                    C := alpha*C*A' + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( JJ+JSEC.LE.N )THEN
                        CALL DGEMM ( 'N', 'T', M, JSEC, N-JJ-JSEC+1,
     $                                      ALPHA, C( 1, JJ+JSEC ), LDC,
     $                                       A( JJ, JJ+JSEC ), LDA, ONE,
     $                                                 C( 1, JJ ), LDC )
                     END IF
  550             CONTINUE
               END IF
            END IF
         ELSE
            IF( NOTR )THEN
*
*              Form  C := alpha*C*A. Right, Lower, No transpose.
*
               TINYM = .NOT.DBIGP( DIP83, M, N )
               IF( TINYM )THEN
                  DO 570, JX = MOD( N-1, RCB )+1, N, RCB
                     JJ = MAX( 1, JX-RCB+1 )
                     JSEC = JX-JJ+1
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'N', M, JSEC, 0, ZERO, C, LDC,
     $                                  A, LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    C := C*A, triangular matrix multiply involving a
*                    lower triangular diagonal block of A.
*
                     DO 560, I = 1, M
                        CALL DTRMV ( 'L', 'T', DIAG, JSEC,
     $                               A( JJ, JJ ), LDA, C( I, JJ ), LDC )
  560                CONTINUE
*
*                    C := alpha*C*A + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( JJ+JSEC.LE.N )THEN
                        CALL DGEMM ( 'N', 'N', M, JSEC, N-JJ-JSEC+1,
     $                                      ALPHA, C( 1, JJ+JSEC ), LDC,
     $                                       A( JJ+JSEC, JJ ), LDA, ONE,
     $                                                 C( 1, JJ ), LDC )
                     END IF
  570             CONTINUE
               ELSE
                  DELTA = ALPHA
                  DO 610, JX = MOD( N-1, CB )+1, N, CB
                     JJ = MAX( 1, JX-CB+1 )
                     JSEC = JX-JJ+1
                     DO 600, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 580, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  580                   CONTINUE
*
*                       C := gamma*T1*A + delta*C, triangular matrix
*                       multiply where the value of delta depends on
*                       whether A is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 590, J = JJ, JJ+JSEC-1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*A( J, J )
                           END IF
                           GAMMA = ALPHA
                           TSEC = JJ+JSEC-1-J
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                              T1( 1, J-JJ+2 ), RB, A( J+1, J ), 1,
     $                                            DELTA, C( II, J ), 1 )
  590                   CONTINUE
  600                CONTINUE
*
*                    C := alpha*C*A + C, general matrix multiply
*                    involving a rectangular block of A.
*
                     IF( JJ+JSEC.LE.N )THEN
                        CALL DGEMM ( 'N', 'N', M, JSEC, N-JJ-JSEC+1,
     $                                      ALPHA, C( 1, JJ+JSEC ), LDC,
     $                                       A( JJ+JSEC, JJ ), LDA, ONE,
     $                                                 C( 1, JJ ), LDC )
                     END IF
  610             CONTINUE
               END IF
            ELSE
*
*              Form  C := alpha*C*A'. Right, Lower, Transpose.
*
               TINYM = .NOT.DBIGP( DIP83, M, N )
               IF( TINYM )THEN
                  DO 630, JX = N, 1, -RCB
                     JJ = MAX( 1, JX-RCB+1 )
                     JSEC = JX-JJ+1
*
*                    C := alpha*C, scale the rectangular block of C
*                    with alpha.
*
                     IF( ALPHA.NE.ONE )
     $                  CALL DGEMM ( 'N', 'T', M, JSEC, 0, ZERO, C, LDC,
     $                                  A, LDA, ALPHA, C( 1, JJ ), LDC )
*
*                    C := C*A', triangular matrix multiply involving
*                    a lower triangular diagonal block of A.
*
                     DO 620, I = 1, M
                        CALL DTRMV ( 'L', 'N', DIAG, JSEC,
     $                               A( JJ, JJ ), LDA, C( I, JJ ), LDC )
  620                CONTINUE
*
*                    C := alpha*C*A' + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( JJ.GT.1 )THEN
                        CALL DGEMM ( 'N', 'T', M, JSEC, JJ-1, ALPHA,
     $                                  C( 1, 1 ), LDC, A( JJ, 1 ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
                     END IF
  630             CONTINUE
               ELSE
                  DELTA = ALPHA
                  DO 680, JX = N, 1, -CB
                     JJ = MAX( 1, JX-CB+1 )
                     JSEC = JX-JJ+1
*
*                    T2 := A', the transpose of a lower unit or non-unit
*                    triangular diagonal block of A is copied to the
*                    upper triangular part of T2.
*
                     DO 640, J = JJ, JJ+JSEC-1-OFFD
                        CALL DCOPY ( JJ+JSEC-J-OFFD, A( J+OFFD, J ),
     $                                1, T2( J-JJ+1, J-JJ+1+OFFD ), CB )
  640                CONTINUE
                     DO 670, II = 1, M, RB
                        ISEC = MIN( RB, M-II+1 )
*
*                       T1 := C, a rectangular block of C is copied
*                       to T1.
*
                        DO 650, J = JJ, JJ+JSEC-1
                           CALL DCOPY ( ISEC, C( II, J ), 1,
     $                                              T1( 1, J-JJ+1 ), 1 )
  650                   CONTINUE
*
*                       C := gamma*T1*T2 + delta*C, triangular matrix
*                       multiply where the value of delta depends on
*                       whether T2 is a unit or non-unit triangular
*                       matrix. Gamma and tsec are used to compensate
*                       for a deficiency in DGEMV that appears if the
*                       second dimension (tsec) is zero.
*
                        DO 660, J = JJ+JSEC-1, JJ, -1
                           IF( NOUNIT )THEN
                              DELTA = ALPHA*T2( J-JJ+1, J-JJ+1 )
                           END IF
                           GAMMA = ALPHA
                           TSEC = J-JJ
                           IF( TSEC.EQ.0 )THEN
                              TSEC = 1
                              GAMMA = ZERO
                           END IF
                           CALL DGEMV ( 'N', ISEC, TSEC, GAMMA,
     $                               T1( 1, 1 ), RB, T2( 1, J-JJ+1 ), 1,
     $                                            DELTA, C( II, J ), 1 )
  660                   CONTINUE
  670                CONTINUE
*
*                    C := alpha*C*A' + C, general matrix multiply
*                    involving the transpose of a rectangular block
*                    of A.
*
                     IF( JJ.GT.1 )THEN
                        CALL DGEMM ( 'N', 'T', M, JSEC, JJ-1, ALPHA,
     $                                  C( 1, 1 ), LDC, A( JJ, 1 ), LDA,
     $                                            ONE, C( 1, JJ ), LDC )
                     END IF
  680             CONTINUE
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM.
*
      END
C     ******************************************************************
      SUBROUTINE DSYRK( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYRK  performs one of the symmetric rank k operations
*
*     C := alpha*A*A' + beta*C,
*
*  or
*
*     C := alpha*A'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*  in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrix  A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*  -- Rewritten in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. Local Scalars ..
      INTEGER            INFO, NROWA
      INTEGER            I, II, IX, ISEC, L, LL, LSEC
      LOGICAL            UPPER, NOTR, CLDA, SMALLN, TINYK
      DOUBLE PRECISION   DELTA
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     .. External Functions ..
      LOGICAL            LSAME, DBIGP, DCLD
      EXTERNAL           LSAME, DBIGP, DCLD
*     .. External Subroutines ..
      EXTERNAL           XERBLA
      EXTERNAL           DGEMM, DGEMV, DSYR, DCOPY, DSCAL
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      INTEGER            DIP41, DIP42
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0,
     $                   DIP41 = 41, DIP42 = 42 )
*     .. User specified parameters for DSYRK ..
      INTEGER            RB, CB, RCB
      PARAMETER        ( RCB = 64, RB = 64, CB = 64 )
*     .. Local Arrays ..
      DOUBLE PRECISION   T1( RB, CB ), T2( RCB, RCB ), T3( RCB, RCB )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      UPPER = LSAME( UPLO, 'U' )
      NOTR = LSAME( TRANS, 'N' )
      IF( NOTR )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      INFO = 0
      IF( ( .NOT.UPPER ).AND.( .NOT.LSAME( UPLO, 'L' ) ) )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTR ).AND.( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $                               ( .NOT.LSAME( TRANS, 'C' ) ) )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( K.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYRK ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when alpha.eq.zero or k.eq.0.
*
      IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
         IF( UPPER )THEN
            DO 10, I = 1, N
               CALL DSCAL ( I, BETA, C( 1, I ), 1 )
   10       CONTINUE
         ELSE
            DO 20, I = 1, N
               CALL DSCAL ( N-I+1, BETA, C( I, I ), 1 )
   20       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( UPPER )THEN
         IF( NOTR )THEN
*
*           Form  C := alpha*A*A' + beta*C. Upper, Notr.
*
            SMALLN = .NOT.DBIGP( DIP41 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.DBIGP( DIP42 , N, K )
               DO 90, II = 1, N, RCB
                  ISEC = MIN( RCB, N-II+1 )
*
*                 C := alpha*A*A' + beta*C, general matrix multiply on
*                 upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                                 A( 1, 1 ), LDA, A( II, 1 ), LDA,
     $                                          BETA, C( 1, II ), LDC )
                  END IF
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a upper triangular diagonal block
*                    of C is updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 30, I = II, II+ISEC-1
                           CALL DSCAL ( I-II+1, BETA, C( II, I ), 1 )
   30                   CONTINUE
                     END IF
*
*                    C := alpha*A*A' + C, symmetric matrix multiply.
*                    C is a symmetric diagonal block having upper
*                    triangular storage format.
*
                     DO 40, L = 1, K
                        CALL DSYR  ( 'U', ISEC, ALPHA, A( II, L ), 1,
     $                                               C( II, II ), LDC )
   40                CONTINUE
                  ELSE
*
*                    T2 := C, a upper triangular diagonal block of the
*                    symmetric matrix C is copied to the upper
*                    triangular part of T2.
*
                     DO 50, I = II, II+ISEC-1
                        CALL DCOPY ( I-II+1, C( II, I ), 1,
     $                                             T2( 1, I-II+1 ), 1 )
   50                CONTINUE
*
*                    T2 :=  beta*T2, the upper triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 60, I = II, II+ISEC-1
                           CALL DSCAL ( I-II+1, BETA,
     $                                             T2( 1, I-II+1 ), 1 )
   60                   CONTINUE
                     END IF
*
*                    T2 := alpha*A*A' + T2, symmetric matrix multiply.
*                    T2 contains a symmetric block having upper
*                    triangular storage format.
*
                     DO 70, L = 1, K
                        CALL DSYR  ( 'U', ISEC, ALPHA, A( II, L ), 1,
     $                                                T2( 1, 1 ), RCB )
   70                CONTINUE
*
*                    C := T2, the upper triangular part of T2 is copied
*                    back to C.
*
                     DO 80, I = II, II+ISEC-1
                        CALL DCOPY ( I-II+1, T2( 1, I-II+1 ), 1,
     $                                                  C( II, I ), 1 )
   80                CONTINUE
                  END IF
   90          CONTINUE
            ELSE
               DO 130, II = 1, N, RB
                  ISEC = MIN( RB, N-II+1 )
*
*                 C := alpha*A*A' + beta*C, general matrix multiply on
*                 upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL DGEMM ( 'N', 'T', II-1, ISEC, K, ALPHA,
     $                                 A( 1, 1 ), LDA, A( II, 1 ), LDA,
     $                                          BETA, C( 1, II ), LDC )
                  END IF
                  DELTA = BETA
                  DO 120, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A, a rectangular block of A is copied to T1.
*
                     DO 100, L = LL, LL+LSEC-1
                        CALL DCOPY ( ISEC, A( II, L ), 1,
     $                                             T1( 1, L-LL+1 ), 1 )
  100                CONTINUE
*
*                    C := alpha*T1*T1' + delta*C, C is symmetric having
*                    triangular storage format. Delta is used instead
*                    of beta to avoid updating the block of C with beta
*                    multiple times.
*
                     DO 110, I = II, II+ISEC-1
                        CALL DGEMV ( 'N', I-II+1, LSEC, ALPHA,
     $                             T1( 1, 1 ), RB, T1( I-II+1, 1 ), RB,
     $                                           DELTA, C( II, I ), 1 )
  110                CONTINUE
                     DELTA = ONE
  120             CONTINUE
  130          CONTINUE
            END IF
         ELSE
*
*           Form  C := alpha*A'*A + beta*C. Upper, Trans.
*
            SMALLN = .NOT.DBIGP( DIP41 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.DBIGP( DIP42 , N, K )
               DO 220, II = 1, N, RCB
                  ISEC = MIN( RCB, N-II+1 )
*
*                 C := alpha*A'*A + beta*C, general matrix multiply on
*                 upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                                 A( 1, 1 ), LDA, A( 1, II ), LDA,
     $                                          BETA, C( 1, II ), LDC )
                  END IF
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a upper triangular diagonal block
*                    of C is updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 140, I = II, II+ISEC-1
                           CALL DSCAL ( I-II+1, BETA, C( II, I ), 1 )
  140                   CONTINUE
                     END IF
*
*                    C := alpha*A*A' + C, symmetric matrix multiply.
*                    C is a symmetric diagonal block having upper
*                    triangular storage format.
*
                     DO 150, L = 1, K
                        CALL DSYR  ( 'U', ISEC, ALPHA, A( L, II ), LDA,
     $                                               C( II, II ), LDC )
  150                CONTINUE
                  ELSE
*
*                    T2 := C, a upper triangular diagonal block of the
*                    symmetric matrix C is copied to the upper
*                    triangular part of T2.
*
                     DO 160, I = II, II+ISEC-1
                        CALL DCOPY ( I-II+1, C( II, I ), 1,
     $                                             T2( 1, I-II+1 ), 1 )
  160                CONTINUE
*
*                    T2 :=  beta*T2, the upper triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 170, I = II, II+ISEC-1
                           CALL DSCAL ( I-II+1, BETA,
     $                                             T2( 1, I-II+1 ), 1 )
  170                   CONTINUE
                     END IF
                     DO 200, LL = 1, K, RCB
                        LSEC = MIN( RCB, K-LL+1 )
*
*                       T3 :=  A', the transpose of a square block of A
*                       is copied to T3.
*
                        DO 180, I = II, II+ISEC-1
                           CALL DCOPY ( LSEC, A( LL, I ), 1,
     $                                           T3( I-II+1, 1 ), RCB )
  180                   CONTINUE
*
*                       T2 := alpha*T3*T3' + T2, symmetric matrix
*                       multiply. T2 contains a symmetric block having
*                       upper triangular storage format.
*
                        DO 190, L = LL, LL+LSEC-1
                           CALL DSYR  ( 'U', ISEC, ALPHA,
     $                            T3( 1, L-LL+1 ), 1, T2( 1, 1 ), RCB )
  190                   CONTINUE
  200                CONTINUE
*
*                    C := T2, the upper triangular part of T2 is copied
*                    back to C.
*
                     DO 210, I = II, II+ISEC-1
                        CALL DCOPY ( I-II+1, T2( 1, I-II+1 ), 1,
     $                                                  C( II, I ), 1 )
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE
               CLDA = DCLD( LDA )
               DO 270, II = 1, N, RB
                  ISEC = MIN( RB, N-II+1 )
*
*                 C := alpha*A'*A + beta*C, general matrix multiply on
*                 upper vertical blocks of C.
*
                  IF( II.GT.1 )THEN
                     CALL DGEMM ( 'T', 'N', II-1, ISEC, K, ALPHA,
     $                                 A( 1, 1 ), LDA, A( 1, II ), LDA,
     $                                          BETA, C( 1, II ), LDC )
                  END IF
                  DELTA = BETA
                  DO 260, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A', the transpose of a rectangular block
*                    of A is copied to T1.
*
                     IF( CLDA )THEN
                        DO 230, I = II, II+ISEC-1
                           CALL DCOPY ( LSEC, A( LL, I ), 1,
     $                                            T1( I-II+1, 1 ), RB )
  230                   CONTINUE
                     ELSE
                        DO 240, L = LL, LL+LSEC-1
                           CALL DCOPY ( ISEC, A( L, II ), LDA,
     $                                             T1( 1, L-LL+1 ), 1 )
  240                   CONTINUE
                     END IF
*
*                    C := alpha*T1*T1' + delta*C, C is symmetric having
*                    triangular storage format. Delta is used instead
*                    of beta to avoid updating the block of C with beta
*                    multiple times.
*
                     DO 250, I = II, II+ISEC-1
                        CALL DGEMV ( 'N', I-II+1, LSEC, ALPHA,
     $                             T1( 1, 1 ), RB, T1( I-II+1, 1 ), RB,
     $                                           DELTA, C( II, I ), 1 )
  250                CONTINUE
                     DELTA = ONE
  260             CONTINUE
  270          CONTINUE
            END IF
         END IF
      ELSE
         IF( NOTR )THEN
*
*           Form  C := alpha*A*A' + beta*C. Lower, Notr.
*
            SMALLN = .NOT.DBIGP( DIP41 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.DBIGP( DIP42 , N, K )
               DO 340, IX = N, 1, -RCB
                  II = MAX( 1, IX-RCB+1 )
                  ISEC = IX-II+1
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a lower triangular diagonal block
*                    of C is updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 280, I = II, II+ISEC-1
                           CALL DSCAL ( II+ISEC-I, BETA, C( I, I ), 1 )
  280                   CONTINUE
                     END IF
*
*                    C := alpha*A*A' + C, symmetric matrix multiply.
*                    C is a symmetric diagonal block having lower
*                    triangular storage format.
*
                     DO 290, L = 1, K
                        CALL DSYR  ( 'L', ISEC, ALPHA, A( II, L ), 1,
     $                                               C( II, II ), LDC )
  290                CONTINUE
                  ELSE
*
*                    T2 := C, a lower triangular diagonal block of the
*                    symmetric matrix C is copied to the lower
*                    triangular part of T2.
*
                     DO 300, I = II, II+ISEC-1
                        CALL DCOPY ( II+ISEC-I, C( I, I ), 1,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  300                CONTINUE
*
*                    T2 :=  beta*T2, the lower triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 310, I = II, II+ISEC-1
                           CALL DSCAL ( II+ISEC-I, BETA,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  310                   CONTINUE
                     END IF
*
*                    T2 := alpha*A*A' + T2, symmetric matrix multiply.
*                    T2 contains a symmetric block having lower
*                    triangular storage format.
*
                     DO 320, L = 1, K
                        CALL DSYR  ( 'L', ISEC, ALPHA, A( II, L ), 1,
     $                                                T2( 1, 1 ), RCB )
  320                CONTINUE
*
*                    C := T2, the lower triangular part of T2 is copied
*                    back to C.
*
                     DO 330, I = II, II+ISEC-1
                        CALL DCOPY ( II+ISEC-I, T2( I-II+1, I-II+1 ),
     $                                                1, C( I, I ), 1 )
  330                CONTINUE
                  END IF
*
*                 C := alpha*A*A' + beta*C, general matrix multiply on
*                 lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K,
     $                         ALPHA, A( II+ISEC, 1 ), LDA, A( II, 1 ),
     $                               LDA, BETA, C( II+ISEC, II ), LDC )
                  END IF
  340          CONTINUE
            ELSE
               DO 380, IX = N, 1, -RB
                  II = MAX( 1, IX-RB+1 )
                  ISEC = IX-II+1
                  DELTA = BETA
                  DO 370, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A, a rectangular block of A is copied to T1.
*
                     DO 350, L = LL, LL+LSEC-1
                        CALL DCOPY ( ISEC, A( II, L ), 1,
     $                                             T1( 1, L-LL+1 ), 1 )
  350                CONTINUE
*
*                    C := alpha*T1*T1' + delta*C, C is symmetric having
*                    triangular storage format. Delta is used instead
*                    of beta to avoid updating the block of C with beta
*                    multiple times.
*
                     DO 360, I = II, II+ISEC-1
                        CALL DGEMV ( 'N', II+ISEC-I, LSEC, ALPHA,
     $                        T1( I-II+1, 1 ), RB, T1( I-II+1, 1 ), RB,
     $                                            DELTA, C( I, I ), 1 )
  360                CONTINUE
                     DELTA = ONE
  370             CONTINUE
*
*                 C := alpha*A*A' + beta*C, general matrix multiply on
*                 lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL DGEMM ( 'N', 'T', N-II-ISEC+1, ISEC, K,
     $                         ALPHA, A( II+ISEC, 1 ), LDA, A( II, 1 ),
     $                               LDA, BETA, C( II+ISEC, II ), LDC )
                  END IF
  380          CONTINUE
            END IF
         ELSE
*
*           Form  C := alpha*A'*A + beta*C. Lower, Trans.
*
            SMALLN = .NOT.DBIGP( DIP41 , N, K )
            IF( SMALLN )THEN
               TINYK = .NOT.DBIGP( DIP42 , N, K )
               DO 470, IX = N, 1, -RCB
                  II = MAX( 1, IX-RCB+1 )
                  ISEC = IX-II+1
                  IF( TINYK )THEN
*
*                    C :=  beta*C, a lower triangular diagonal block
*                    of C is updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 390, I = II, II+ISEC-1
                           CALL DSCAL ( II+ISEC-I, BETA, C( I, I ), 1 )
  390                   CONTINUE
                     END IF
*
*                    C := alpha*A*A' + C, symmetric matrix multiply.
*                    C is a symmetric diagonal block having lower
*                    triangular storage format.
*
                     DO 400, L = 1, K
                        CALL DSYR  ( 'L', ISEC, ALPHA, A( L, II ), LDA,
     $                                               C( II, II ), LDC )
  400                CONTINUE
                  ELSE
*
*                    T2 := C, a lower triangular diagonal block of the
*                    symmetric matrix C is copied to the lower
*                    triangular part of T2.
*
                     DO 410, I = II, II+ISEC-1
                        CALL DCOPY ( II+ISEC-I, C( I, I ), 1,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  410                CONTINUE
*
*                    T2 :=  beta*T2, the lower triangular part of T2 is
*                    updated with beta.
*
                     IF( BETA.NE.ONE )THEN
                        DO 420, I = II, II+ISEC-1
                           CALL DSCAL ( II+ISEC-I, BETA,
     $                                        T2( I-II+1, I-II+1 ), 1 )
  420                   CONTINUE
                     END IF
                     DO 450, LL = 1, K, RCB
                        LSEC = MIN( RCB, K-LL+1 )
*
*                       T3 :=  A', the transpose of a square block of A
*                       is copied to T3.
*
                        DO 430, I = II, II+ISEC-1
                           CALL DCOPY ( LSEC, A( LL, I ), 1,
     $                                           T3( I-II+1, 1 ), RCB )
  430                   CONTINUE
*
*                       T2 := alpha*T3*T3' + T2, symmetric matrix
*                       multiply. T2 contains a symmetric block having
*                       lower triangular storage format.
*
                        DO 440, L = LL, LL+LSEC-1
                           CALL DSYR  ( 'L', ISEC, ALPHA,
     $                            T3( 1, L-LL+1 ), 1, T2( 1, 1 ), RCB )
  440                   CONTINUE
  450                CONTINUE
*
*                    C := T2, the lower triangular part of T2 is copied
*                    back to C.
*
                     DO 460, I = II, II+ISEC-1
                        CALL DCOPY ( II+ISEC-I, T2( I-II+1, I-II+1 ),
     $                                                1, C( I, I ), 1 )
  460                CONTINUE
                  END IF
*
*                 C := alpha*A'*A + beta*C, general matrix multiply on
*                 lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K,
     $                         ALPHA, A( 1, II+ISEC ), LDA, A( 1, II ),
     $                               LDA, BETA, C( II+ISEC, II ), LDC )
                  END IF
  470          CONTINUE
            ELSE
               CLDA = DCLD( LDA )
               DO 520, IX = N, 1, -RB
                  II = MAX( 1, IX-RB+1 )
                  ISEC = IX-II+1
                  DELTA = BETA
                  DO 510, LL = 1, K, CB
                     LSEC = MIN( CB, K-LL+1 )
*
*                    T1 := A', the transpose of a rectangular block
*                    of A is copied to T1.
*
                     IF( CLDA )THEN
                        DO 480, I = II, II+ISEC-1
                           CALL DCOPY ( LSEC, A( LL, I ), 1,
     $                                            T1( I-II+1, 1 ), RB )
  480                   CONTINUE
                     ELSE
                        DO 490, L = LL, LL+LSEC-1
                           CALL DCOPY ( ISEC, A( L, II ), LDA,
     $                                             T1( 1, L-LL+1 ), 1 )
  490                   CONTINUE
                     END IF
*
*                    C := alpha*T1*T1' + delta*C, C is symmetric having
*                    triangular storage format. Delta is used instead
*                    of beta to avoid updating the block of C with beta
*                    multiple times.
*
                     DO 500, I = II, II+ISEC-1
                        CALL DGEMV ( 'N', II+ISEC-I, LSEC, ALPHA,
     $                        T1( I-II+1, 1 ), RB, T1( I-II+1, 1 ), RB,
     $                                            DELTA, C( I, I ), 1 )
  500                CONTINUE
                     DELTA = ONE
  510             CONTINUE
*
*                 C := alpha*A'*A + beta*C, general matrix multiply on
*                 lower vertical blocks of C.
*
                  IF( II+ISEC.LE.N )THEN
                     CALL DGEMM ( 'T', 'N', N-II-ISEC+1, ISEC, K,
     $                         ALPHA, A( 1, II+ISEC ), LDA, A( 1, II ),
     $                               LDA, BETA, C( II+ISEC, II ), LDC )
                  END IF
  520          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DSYRK.
*
      END
