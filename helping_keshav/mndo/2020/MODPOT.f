      SUBROUTINE MODPOT
C     *
C     MODIFY PARAMETERS FOR THE ELECTROSTATIC POTENTIAL.
C     *
C     THESE PARAMETERS MAY BE DEFINED SEPARATELY FOR EACH MOLECULE.
C     EXPLICIT INPUT FOR THE PARAMETERS (OPTION IPAROK) THUS REQUIRES
C     A REDEFINITION OF THE DEFAULT VALUES (OPTIONS IPOT,MMPOT),
C     WHICH IS IMPLEMENTED BY CALLING MODPOT AFTER CALLING INIPOT.
C     SEE SUBROUTINES INIPOT AND MODPAR FOR FURTHER INFORMATION.
C     *
      USE LIMIT, ONLY: LMZ, LMP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON
     ./COSMO6/ SOLV1(LMZ,6)
     ./PAR20 / XXPAR(LMP),NFPAR(LMP),NSPAR(LMP),NNPAR
     ./QMMM3 / DELTAM(LMZ,2)
      IF(NNPAR.LE.0) RETURN
C *** LOOP OVER INPUT PARAMETERS AND REPLACE PARAMETER VALUES.
      DO 10 I=1,NNPAR
      NF     = NFPAR(I)
      NS     = NSPAR(I)
      A      = XXPAR(I)
C     PARAMETERS FOR SOLVATION MODELS.
      IF(NF.GE.47 .AND. NF.LE.48) THEN
         SOLV1(NS,NF-42) = A
C     PARAMETERS FOR ELECTROSTATIC POTENTIALS.
      ELSE IF(NF.GE.53 .AND. NF.LE.54) THEN
         DELTAM(NS,NF-52) = A
      ENDIF
   10 CONTINUE
      RETURN
      END
