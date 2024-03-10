      SUBROUTINE SORTIX (IND,N)
C     *
C     SORT NUMBERS IN IND(N) SUCH THAT THEY ARE IN ASCENDING ORDER.
C     *
      DIMENSION IND(N)
      DO I=1,N-1
         K = I
         DO J=I+1,N
            IF(IND(J).LT.IND(K))  K = J
         ENDDO
         IF(K.NE.I) THEN
            J      = IND(I)
            IND(I) = IND(K)
            IND(K) = J
         ENDIF
      ENDDO
      RETURN
      END
