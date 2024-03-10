      CHARACTER*80 FUNCTION GETNAM(NAMEIN)
C     *
C     ON A UNIX SYSTEM, GETENV WILL CONSULT THE ENVIRONMENT FOR THE
C     CURRENT ALIAS OF THE CHARACTER STRING CONTAINED IN 'NAMEIN'.
C     THE ALIAS, IF IT EXISTS, OR THE ORIGINAL NAME IN 'NAMEIN'
C     WILL BE RETURNED.
C     ON NON-UNIX SYSTEMS, PLEASE COMMENT OUT THE CALL TO GETENV
C     SO THAT THE ORIGINAL NAME IN 'NAMEIN' IS RETURNED.
C     ADAPTED FROM MOPAC(6.0) WRITTEN BY J.J.P.STEWART.
C     *
      CHARACTER*(*) NAMEIN
      CHARACTER*(80) NAMOUT
      NAMOUT=' '
      CALL GETENV(NAMEIN,NAMOUT)
      IF (NAMOUT.EQ.'  ') NAMOUT=NAMEIN
      GETNAM = NAMOUT
      RETURN
      END