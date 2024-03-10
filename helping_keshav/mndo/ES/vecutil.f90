      ! C = A + B
      SUBROUTINE VecAdd(C,A,B,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),B(*),C(*)
      M = MOD(N,4)
      IF(M.NE.0) THEN
         DO I = 1,M
            C(I) = A(I) + B(I)
         END DO
      ENDIF
      IF(N.LT.4) RETURN
      DO I = M+1,N,4
         C(I)   = A(I)   + B(I)
         C(I+1) = A(I+1) + B(I+1)
         C(I+2) = A(I+2) + B(I+2)
         C(I+3) = A(I+3) + B(I+3)
      ENDDO
      RETURN
      END

      ! C = A - B
      SUBROUTINE Vecsub(C,A,B,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),B(*),C(*)
      M = MOD(N,4)
      IF(M.NE.0) THEN
         DO I = 1,M
            C(I) = A(I) - B(I)
         END DO
      ENDIF
      IF(N.LT.4) RETURN
      DO I = M+1,N,4
         C(I)   = A(I)   - B(I)
         C(I+1) = A(I+1) - B(I+1)
         C(I+2) = A(I+2) - B(I+2)
         C(I+3) = A(I+3) - B(I+3)
      ENDDO
      RETURN
      END

      ! C = x*A + B
      SUBROUTINE Vecadd2(C,x,A,B,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),B(*),C(*)
      M = MOD(N,4)
      IF(M.NE.0) THEN
         DO I = 1,M
            C(I) = x*A(I) + B(I)
         END DO
      ENDIF
      IF(N.LT.4) RETURN
      DO I = M+1,N,4
         C(I)   = x*A(I)   + B(I)
         C(I+1) = x*A(I+1) + B(I+1)
         C(I+2) = x*A(I+2) + B(I+2)
         C(I+3) = x*A(I+3) + B(I+3)
      ENDDO
      RETURN
      END

      ! C = x*A + C
      SUBROUTINE Vecadd3(C,x,A,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*), C(*)
      M = MOD(N,4)
      IF(M.NE.0) THEN
         DO I = 1,M
            C(I) = x*A(I) + C(I)
         END DO
      ENDIF
      IF(N.LT.4) RETURN
      DO I = M+1,N,4
         C(I)   = x*A(I)   + C(I)
         C(I+1) = x*A(I+1) + C(I+1)
         C(I+2) = x*A(I+2) + C(I+2)
         C(I+3) = x*A(I+3) + C(I+3)
      ENDDO
      RETURN
      END


      SUBROUTINE VecCopy(B,A,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),B(*)
      M = MOD(N,7)
      IF(M.NE.0) THEN
         DO I = 1,M
            B(I) = A(I) 
         END DO
      ENDIF
      IF(N.LT.7) RETURN
      DO I = M+1,N,7
         B(I)   = A(I)   
         B(I+1) = A(I+1) 
         B(I+2) = A(I+2) 
         B(I+3) = A(I+3) 
         B(I+4) = A(I+4) 
         B(I+5) = A(I+5) 
         B(I+6) = A(I+6) 
      ENDDO
      RETURN
      END


      ! C = A \dot B
      SUBROUTINE Vecdot(Q,A,B,N)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),B(*),Q
      Q = 0d0
      M = MOD(N,5)
      IF(M.NE.0) THEN
         DO I = 1,M
            Q = Q + A(I)*B(I)
         END DO
         IF (N.LT.5) RETURN
      END IF
      DO I=M+1,N,5
         Q=Q + A(I)   * B(I)   & 
             + A(I+1) * B(I+1) &
             + A(I+2) * B(I+2) &
             + A(I+3) * B(I+3) &
             + A(I+4) * B(I+4) 
      END DO
      RETURN
      END
 
      ! A = Value
      SUBROUTINE VecInit(A,N,Value)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),Value
      M = MOD(N,4)
      IF(M.NE.0) THEN
         DO I = 1,M
            A(I) = Value
         END DO
      ENDIF
      IF(N.LT.4) RETURN
      DO I = M+1,N,4
         A(I)   = Value 
         A(I+1) = Value
         A(I+2) = Value
         A(I+3) = Value
      ENDDO
      RETURN
      END

      ! A = Scale*A
      SUBROUTINE VecScale(A,N,x)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 A(*),x
      M = MOD(N,5)
      IF(M.NE.0) THEN
         DO I = 1,M
            A(I) = x*A(I) 
         END DO
      ENDIF
      IF(N.LT.5) RETURN
      DO I = M+1,N,5
         A(I)   = x*A(I)   
         A(I+1) = x*A(I+1) 
         A(I+2) = x*A(I+2) 
         A(I+3) = x*A(I+3) 
         A(I+4) = x*A(I+4) 
      ENDDO
      RETURN
      END


