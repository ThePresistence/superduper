SUBROUTINE ScosmoAddFck(F,PA,PB,COORD,LM4,NUMAT,&
     NAT,NFIRST,NLAST,TORE,IP,IP1,IP2,DD,NITER,Epol)
!!!*** 
!!! NAME
!!!     ScosmoAddFck - A f95 interface routine to obtain COSMO modified Fock matrix 
!!! CALLED BY
!!!     ITER in MNDO
!!! SCRATCH
!!!     IP(LMI) -  INDICES OF UNIQUE ONE-CENTER AO PAIRS = I*(I-1)/2+J
!!!     LM4 - size of PA,PB
!!!     
!!!  
!!!***

  USE DataTypes
  USE ConstantsMod
  USE SCosmoMNDOMod
  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: LM1=600
  INTEGER(I4B), PARAMETER :: LMI=45*LM1
  INTEGER(I4B), PARAMETER :: LMZ=86
  INTEGER(I4B), INTENT(IN) :: LM4, NUMAT, NITER

!!! **** Cant declare these as assumed-shape or assumed-size arrays; hence explicit dimension.
  INTEGER(I4B), INTENT(IN):: NAT(LM1), NFIRST(LM1),NLAST(LM1), IP(LMI),IP1(LMI),IP2(LMI)
  REAL(SP), INTENT(IN) :: PA(LM4), PB(LM4), COORD(3,LM1), TORE(LMZ), DD(6,0:LMZ)
  REAL(SP), INTENT(INOUT):: F(LM4)
  REAL(SP), INTENT(OUT) :: Epol
  REAL(SP) :: Ecdr
  REAL(SP) :: P(LM4)

!  WRITE(6,*)'### ScosmoAddFck',LM4

  P(:) = PA(:) + PB(:)

  CALL AddFckNew(F,P,COORD(1:3,1:NUMAT),NUMAT,NAT, &
       NFIRST,NLAST,TORE,IP,IP1,IP2,DD,NITER,Epol, &
       Ecdr)

  Epol = Epol/ELECTRON_VOLT

END SUBROUTINE SCosmoAddFck

SUBROUTINE InitSC(NQ,NUMAT,NAT, ELEMNT, SCRQ0, SCRQ1, SCRQ2, SCNET, SCNEB, JPRINT)
  USE DataTypes
  USE SCosmoMNDOMod, only :InitSCosmo
  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: LM1=600
  INTEGER(I4B), PARAMETER :: LMZ=86
  INTEGER(I4B), INTENT(IN) :: NQ, NUMAT, NAT(LM1)
  REAL(SP), INTENT(IN) :: SCRQ0(LMZ),SCRQ1(LMZ),SCRQ2(LMZ),SCNET(LMZ),SCNEB(LMZ)
  CHARACTER(LEN=2), INTENT(IN) :: ELEMNT(107) 
  INTEGER(I4B), INTENT(IN) :: JPRINT
  
  CALL InitSCosmo(NQ, NUMAT,NAT, ELEMNT, SCRQ0, SCRQ1, SCRQ2, SCNET, SCNEB)
  
END SUBROUTINE InitSC

! Allows MNDO-SCOSMO Code to use surface construced externally
! (i,e. from an MD code)
subroutine mndo_scosmo_from_external(cdata,surf,Qtot,prs,grad)
  Use Datatypes
  USE SurfaceMod
  USE CosmoMod
  USE SCosmoMNDOMod
  IMPLICIT NONE
  
  TYPE(cosmo_t),TARGET :: cdata
  TYPE(surface_t), TARGET :: surf
  REAL(SP), INTENT(IN) :: Qtot
  REAL(SP), pointer :: prs(:)
  REAL(SP), pointer :: grad(:,:)

  QEXTERN = Qtot 
  cosmo_mndo => cdata
  esurf_mndo => surf

  phi_rho_sigma_cons => prs
  grad_surf => grad
  iradii=1
      
  EXTSURF=.TRUE.

end subroutine mndo_scosmo_from_external

SUBROUTINE SCosmoDieleCav(Epol, Ecdr, PA, &
     PB,COORD,LM4,NUMAT,NAT,NFIRST,NLAST,TORE,IP,IP1,IP2,DD)
!!!***
!!! NAME
!!!     SCosmoDielen - A f95 interface routine to obtain Smooth COSMO
!!!                    dielectric screening energy Epol
!!! CALLED BY
!!!     SCFCAL in MNDO
!!!***

  USE DataTypes
  USE SCosmoMNDOMod
  IMPLICIT NONE

!!! **** assumed-shape or assumed-size arrays don't work; hence explicit dimension.
  INTEGER(I4B), PARAMETER :: LM1=600
  INTEGER(I4B), PARAMETER :: LMI=45*LM1
  INTEGER(I4B), PARAMETER :: LMZ=86
  INTEGER(I4B), INTENT(IN) :: LM4, NUMAT
  INTEGER(I4B), INTENT(IN):: NAT(LM1), NFIRST(LM1),NLAST(LM1), IP(LMI),IP1(LMI),IP2(LMI)
  REAL(SP), INTENT(IN) :: PA(LM4), PB(LM4), COORD(3,LM1), TORE(LMZ), DD(6,0:LMZ)
  REAL(SP), INTENT(OUT) :: Epol,Ecdr
  REAL(SP) :: P(LM4)

  P(:) = PA(:) + PB(:)
  CALL CptSCosmoE(P, COORD(:,1:NUMAT), NUMAT, NAT, NFIRST, NLAST, TORE, IP, IP1, IP2, DD, &
       .TRUE., -1 , Epol, Ecdr)


END SUBROUTINE SCosmoDieleCav


SUBROUTINE SCosmoPrtcsm(ESolut,NAT,ELEMNT)
!!!***
!!! NAME
!!!     ScosmoPrtcsm - A f95 routine to print smooth COSMO results
!!! CALLED BY
!!!     SCFCAL in MNDO
!!!***

  USE DataTypes
  USE ConstantsMod, ONLY : AU=>AU_DISTANCE_IN_ANGSTROM, KCAL_PER_MOL
  USE SCosmoMNDOMod
  USE SurfaceMod
  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: LM1=600
  REAL(SP), INTENT(IN) :: ESolut
  INTEGER(I4B), INTENT(IN) :: NAT(LM1)
  CHARACTER*2, INTENT(IN) :: ELEMNT(107)
  INTEGER(I4B) :: i

  WRITE(6,500)
  IF(iRadii >= 2) THEN
    WRITE(6,510)
    WRITE(6,520)
    DO i=1, SIZE(radii_main)
       WRITE(6,'(3X,I6,9X,A,4x,F11.4,3X,F11.4)')i,ELEMNT(NAT(i)),radii_main(i)*AU,AtomCharges(i)
    ENDDO
    WRITE(6,*)
  ENDIF
!!!**** PRINT SUMMARY OF SCOSMO RESULTS.

  WRITE(6,540)
  WRITE(6,560) ESOLUT
  WRITE(6,600) cosmo_mndo%Ecdr/KCAL_PER_MOL
  WRITE(6,561) ESOLUT + cosmo_mndo%Ecdr/KCAL_PER_MOL
  WRITE(6,610) cosmo_mndo%Epol/KCAL_PER_MOL
  WRITE(6,615) (cosmo_mndo%Epol+cosmo_mndo%Ecdr)/KCAL_PER_MOL
  WRITE(6,620) AU*AU*esurf_mndo%SA
  WRITE(6,630) SIZE(esurf_mndo%SECRD,2)

500 FORMAT(///5X,'*** SUMMARY OF COSMO RESULTS ***')
510 FORMAT(/  5X,'SUMMARY OF VARIABLE RADII USED:')
520 FORMAT(   5X,'ATOM      ELEMENT      RADIUS(A)      CHARGE')
540 FORMAT(/  5X,'SMOOTH COSMO CALCULATION')
560 FORMAT(/  5X,'HEAT OF FORMATION IN SOLUTION ',ES15.5,' KCAL/MOL')
561 FORMAT(/  5X,'HEAT OF FORM IN SOL /W NONELEC',ES15.5,' KCAL/MOL')
590 FORMAT(   5X,'ELECTROSTATIC SCOSMO TERM     ',ES15.5,' KCAL/MOL')
600 FORMAT(   5X,'NONELECTROSTATIC TERM         ',F15.5,' KCAL/MOL')
610 FORMAT(   5X,'DIELECTRIC SCREENING ENERGY   ',ES15.5,' KCAL/MOL')
615 FORMAT(   5X,'DIELECTRIC SCREENING + NONELEC',ES15.5,' KCAL/MOL')
620 FORMAT(   5X,'SURFACE OF THE MOLECULE       ',F15.5,' A**2')
630 FORMAT(   5X,'TOTAL NUMBER OF SEGMENTS      ',I15)

END SUBROUTINE SCosmoPrtcsm


SUBROUTINE SCosmoDGrad(DXYZ, COORD, NUMAT, DD, NAT, NFIRST, NLAST)
!!!***
!!! DESCRIPTION
!!!     an interface routine to compute COSMO gradients by calling 
!!!     AnalytGradient in MinimizeMod
!!! AUTHOR
!!!     Jana Khandogin
!!! CALLED BY
!!!     SCF in MNDO
!!!***

  USE DataTypes
  USE ConstantsMod, AU=>AU_DISTANCE_IN_ANGSTROM
  USE UtilitiesMod
  USE SurfaceMod
  USE MinimizeMod
  USE SCosmoMNDOMod
  Use cosmoMod
  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: LM1=600
  INTEGER(I4B), PARAMETER :: LMZ=86
  INTEGER(I4B), INTENT(IN) :: NUMAT, NAT(LM1), NFIRST(LM1),NLAST(LM1)
  INTEGER(I4B) :: n, Nse, i

  REAL(SP), INTENT(IN) :: COORD(3,*), DD(6,0:LMZ)
  REAL(SP), INTENT(INOUT) :: DXYZ(3,*)
  REAL(SP), allocatable :: Grad_rho_sigma(:,:), Grad_sigma_rho(:,:)
  REAL(SP), allocatable :: CavGrad(:,:),Grad(:,:)
  REAL(SP) :: cosmo_f,  Epol, Ecdr, iCDR
  LOGICAL(LGD) :: success
  REAL(SP), parameter :: conv = ANGSTROM/KCAL_PER_MOL !au to kcal/mol/Ang


  Epol = Cosmo_mndo%epol
  iCDR = esurf_mndo%icdr
  Ecdr = Cosmo_mndo%ecdr
  cosmo_f = Cosmo_mndo%f

  IF (esurf_mndo%VERBOSE) THEN
     WRITE(6,*) 'compute gradient with Gaussians ', Ifgaus
     WRITE(6,*) 'Gradient Without Correction'
     DO i = 1, NUMAT
        WRITE(6,'(a5,i5, 3f12.5)') 'Ggas', i, DXYZ(:,i)
     END DO
  END IF

  IF(cosmo_mndo%igrad==0) RETURN
  ALLOCATE(CavGrad(3,NUMAT))
  ALLOCATE(Grad(3,NUMAT))

  IF (cosmo_mndo%igrad==2 .OR. cosmo_mndo%igrad ==3) THEN ! numerical gradients
     CALL MNDONumGradientFast(DD, NAT, NFIRST, NLAST, cosmo_f, &
          ANGSTROM*COORD(1:3,1:NUMAT), Cosmo_mndo%surface_charge, &
          cosmo_mndo%minpara%FDstep, IfGaus, iCDR, Epol, Ecdr, Grad, CavGrad, success)

     IF (.NOT. success) THEN
        
        CALL MNDONumGradientSlow(DD, NAT, NFIRST, NLAST, &
             ANGSTROM*COORD(1:3,1:NUMAT), IfGaus, &
             Epol, Ecdr, Grad, CavGrad)
     END IF

     IF (esurf_mndo%VERBOSE) THEN
        WRITE(6,*) 'Numerical Gradient due to solvent'
        DO i = 1, NUMAT
           WRITE(6,'(a8,i5, 3f12.5)') 'EpolGrad', i, Grad(:,i) * conv
        END DO

        IF (iCDR /= 0 .AND. esurf_mndo%gammaswitch > 0.0) THEN
           WRITE(6,*) 'Numerical Gradient due to cavitation dispersion repulsion term'
           DO i = 1, NUMAT
              WRITE(6,'(a8,i5, 3f12.5)') 'ECavGrad', i, CavGrad(:,i) * conv
           END DO
        END IF
     END IF
  END IF

  IF (cosmo_mndo%igrad == 1 .OR. cosmo_mndo%igrad ==3) THEN ! analytic gradients

     NSE = SIZE(esurf_mndo%SECRD,2)

     ALLOCATE(Grad_rho_sigma(3,Nse))
     ALLOCATE(Grad_sigma_rho(3,NUMAT))
     IF (IfGaus) THEN
        CALL CptGradMNDO(COORD(1:3,1:NUMAT), AU * esurf_mndo%SECRD, NUMAT, &
             DD, NAT, NFIRST, NLAST, GRAD_RHO_SIGMA=Grad_rho_sigma, &
             GRAD_SIGMA_RHO=Grad_sigma_rho, RHO=rho, SIGMA=Cosmo_mndo%surface_charge, & 
             ZETA=ANGSTROM*esurf_mndo%COSMOZETA)
     ELSE
        CALL CptGradMNDO(COORD(1:3,1:NUMAT), AU * esurf_mndo%SECRD, NUMAT, &
             DD, NAT, NFIRST, NLAST, GRAD_RHO_SIGMA=Grad_rho_sigma, &
             GRAD_SIGMA_RHO=Grad_sigma_rho, RHO=rho, SIGMA=Cosmo_mndo%surface_charge)
     END IF
     
!!! **** convert Angstrom-2 to AU-2
     Grad_rho_sigma = Grad_rho_sigma * AU*AU
     Grad_sigma_rho = Grad_sigma_rho * AU*AU
 
!!!**** convert Angstrom to Bohr first
     n = assert_eq(SIZE(Cosmo_mndo%surface_charge),SIZE(esurf_mndo%SECRD,2), SIZE(Grad_rho_sigma,2), 'SCosmoDgrad')

     IF(EXTSURF)THEN 
        grad_surf = Grad_rho_sigma
     ELSE
        CALL AnalytGradient(esurf_mndo, Cosmo_mndo%surface_charge, cosmo_f, &
             cosmo_mndo%minpara, Grad_sigma_rho, Grad_rho_sigma)
     ENDIF
     DEALLOCATE(Grad_rho_sigma)

     IF (iCDR /= 0 .AND. esurf_mndo%gammaswitch > 0.0) THEN
        CALL AnalytMolNonelecGradient(esurf_mndo, CavGrad)
     END IF

     IF (esurf_mndo%VERBOSE) THEN
        WRITE(6,*) 'Analytic Gradient due to solvent'
        DO i = 1, NUMAT
           WRITE(6,'(a8,i5, 3f12.5)') 'EpolGrad', i, Grad_sigma_rho(:,i) * conv
        END DO

        IF (iCDR /= 0 .AND. esurf_mndo%gammaswitch > 0.0) THEN
           WRITE(6,*) 'Analytic Gradient due to cavitation dispersion repulsion term'
           DO i = 1, NUMAT
              WRITE(6,'(a8,i5, 3f12.5)') 'ECavGrad', i, CavGrad(:,i) * conv
           END DO
        END IF
     END IF

  END IF

!!! **** convert from AU-2 to kcal/(mol*Ang) and add to DXYZ
  DXYZ(:,1:NUMAT) = DXYZ(:,1:NUMAT) + Grad_sigma_rho(:,1:NUMAT) * conv

  IF (icdr /= 0 .AND. esurf_mndo%gammaswitch > 0.0) THEN 
     DXYZ(:,1:NUMAT) = DXYZ(:,1:NUMAT) + CavGrad(:,1:NUMAT) * conv
  END IF

  DEALLOCATE(Grad_sigma_rho,CavGrad)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  SUBROUTINE MNDONumGradientFast(DD, NAT, NFIRST, NLAST, cosmof, &
       crd, surface_charge, FDstep, IfGaus, iCDR, &
       Epol, Ecdr, EGrad, CavGrad, success)
!!!****f* MNDONumGradientFast
!!!
!!! Description
!!!     fast way to compute the numerical gradients of Epol and Ecav
!!!     using frozen density approximation (surface_charge remains the same)
!!! Scratch
!!!     SigmaFldSigma - potential and field at srf points due to srf charges
!!!     SigmaFldRho   - potential and field at srf points due to solute charges
!!!     RhoFld        - potential and field at solute atom positions due to srf charges
!!!
!!!***

    USE CosmoMod
    USE SurfaceMod
    USE RBM
    USE cosmo_util
    USE SCosmoMNDOMod
    USE ConstantsMod
    IMPLICIT NONE

    TYPE(surface_t) :: Surface_local
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    INTEGER(I4B) :: nse, natm, iatm, k, NQ

    REAL(SP), INTENT(IN) :: crd(:,:),surface_charge(:)
    REAL(SP), INTENT(IN) :: cosmof, FDstep, Epol, Ecdr, iCDR
    REAL(SP), INTENT(IN) :: DD(:,:)

    REAL(SP), INTENT(OUT) :: EGrad(:,:)
    REAL(SP), INTENT(OUT) :: Cavgrad(:,:)
    REAL(SP) :: CrdNew(3,SIZE(crd,2)), EpolNew, EcdrNew, AU, f
    REAL(SP), ALLOCATABLE :: SigmaCoord(:,:), phi_sigma_sigma(:), phi_sigma_rho(:)

    LOGICAL, INTENT(IN) :: IfGAUS
    LOGICAL, INTENT(OUT) :: success

    AU = AU_DISTANCE_IN_ANGSTROM
    nse = SIZE(surface_charge)
    Natm = SIZE(crd,2)
    NQ = SIZE(rho)

    EGrad = 0.0_SP
    CavGrad = 0.0_SP

    IF (cosmof/=0.0_SP) THEN
       f = 1.0_SP/cosmof
    ELSE
       RETURN
    END IF

    CALL copy_surface_scalar_inputs(esurf_mndo, Surface_local)

    IF (iCDR /= 0) THEN
       WRITE(6,*) 'Cavitation Dispersion Repulsion ON'
    END IF

    ALLOCATE(SigmaCoord(1:5,1:nse), phi_sigma_sigma(1:nse))
    allocate(phi_sigma_rho(1:NQ))

    SigmaCoord(1,:) =  surface_charge(:)

    DO iatm = 1, Natm
       DO k = 1, 3
          CrdNew(:,:) = crd(:,:)
          CrdNew(k,iatm) = CrdNew(k,iatm) + FDstep
          CALL SASurface(CrdNew, radii_main, Surface_local)

          IF (SIZE(Surface_local%cosmozeta) /= nse) THEN
             IF (surface_local%verbose) &
                  WRITE(6,*) 'surface charge vector length changes. do slow gradient'

             CALL DeallocateSurfaceArrays(Surface_local)
             DEALLOCATE(SigmaCoord, phi_sigma_sigma,phi_sigma_rho)
             Success = .FALSE.
             RETURN
          END IF

          SigmaCoord(2:4,:) = Surface_local%SECRD(1:3,:)
          SigmaCoord(5,:) = Surface_local%cosmozeta(:)

          CALL DirectCoulomb(SigmaCoord,phi_sigma_sigma)

          IF (IfGaus) THEN  
             CALL Cpt_bmat_Phi_MNDO(AU*CrdNew, AU*Surface_local%SECRD, &
                  DD, NAT, Natm, NFIRST, NLAST, PHI_SIGMA_RHO=phi_sigma_rho, &
                  SIGMA=surface_charge, ZETA=ANGSTROM *Surface_local%cosmozeta)
          ELSE
             CALL Cpt_bmat_Phi_MNDO(AU*CrdNew, AU*Surface_local%SECRD, &
                  DD, NAT, Natm, NFIRST, NLAST, &
                  PHI_SIGMA_RHO=phi_sigma_rho, SIGMA=surface_charge)
          END IF

          phi_sigma_rho = phi_sigma_rho * AU !convert Angstrom-1 to Bohr-1

          IF (surface_local%GammaSwitch > 0) THEN
             CALL switch_pot(phi_sigma_sigma, surface_charge, Surface_local%Sescale, &
                  Surface_local%cosmozeta)
          END IF

!!!**** gradient due to 0.5*sigma(T)*A*sigma + sigma(T)*B*rho
          EpolNew = 0.5_SP*f*DOT_PRODUCT(surface_charge,phi_sigma_sigma) + &
               DOT_PRODUCT(rho, phi_sigma_rho)

          EGrad(k,iatm) = (EpolNew - Epol)/FDStep

!!!**** cavitation energy
          IF (iCDR /= 0  .AND. Surface_local%GammaSwitch > 0) THEN
             CALL MolecularNonelecEnergy(Surface_local, EcdrNew)
             CavGrad(k,iatm) = (EcdrNew - Ecdr)/FDStep 
          END IF

          CALL DeallocateSurfaceArrays(Surface_local)

       END DO ! k
    END DO ! iatm

    Success = .TRUE.

    DEALLOCATE(SigmaCoord, phi_sigma_sigma,phi_sigma_rho)

  END SUBROUTINE MNDONumGradientFast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MNDONumGradientSlow(DD, NAT, NFIRST, NLAST, &
       crd, IfGaus, Epol, Ecdr, &
       EGrad, CavGrad)
!!!****f* MNDONumGradientSlow
!!!
!!! Description
!!!     Slow way to compute the numerical gradients of Epol and Ecav
!!!     by recomputing Epol and Ecav
!!!     This necessary when surface_charge changes dimension when small FDstep is applied
!!!***
    USE ConstantsMod
    USE SurfaceMod
    USE CosmoMod
    USE cosmo_util
    USE SCosmoMNDOMod
    USE DATATYPES
    IMPLICIT NONE
 
    TYPE(cosmo_t) :: cosmo_data_local
    TYPE(surface_t) :: surface_local

    INTEGER(I4B) :: natm, iatm, k, NQ
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)

    REAL(SP), INTENT(IN) :: crd(:,:), DD(:,:)
    REAL(SP), INTENT(IN) :: Epol, Ecdr
    REAL(SP), INTENT(OUT) :: EGrad(:,:)
    REAL(SP), INTENT(OUT) :: Cavgrad(:,:)
    LOGICAL(LGD), intent(in) :: IFGaus
    REAL(SP),allocatable :: phi(:)
    REAL(SP) :: CrdNew(3,SIZE(crd,2)), EpolNew, EcdrNew, AU


    AU = AU_DISTANCE_IN_ANGSTROM

    Natm = SIZE(crd,2)
    NQ = SIZE(rho)

    ! Retrieve variables used for cosmo calcs.
    call copy_cosmo_scalar_inputs(cosmo_mndo, cosmo_data_local)
    CALL copy_surface_scalar_inputs(esurf_mndo, surface_local)

    EGrad = 0.0_SP
    CavGrad = 0.0_SP

    DO iatm = 1, Natm
       DO k = 1, 3

!!!**** assign CrdNew to  be the same as Crd, then move a FDstep in k direction
          CrdNew(:,:) = crd(:,:)
          CrdNew(k,iatm) = CrdNew(k,iatm) + cosmo_data_local%minpara%FDstep

          CALL SASurface(CrdNew, Radii_main, surface_local)
          allocate(phi(SIZE(surface_local%SECRD,2)))

          IF (ifGaus) THEN
             CALL Cpt_bmat_Phi_MNDO(AU*CrdNew, AU*surface_local%SECRD, &
                  DD, NAT, NUMAT, NFIRST, NLAST, &
                  PHI_RHO_SIGMA=phi, RHO=rho, &
                  ZETA=ANGSTROM * surface_local%COSMOZETA)
          ELSE
             CALL Cpt_bmat_Phi_MNDO(AU*CrdNew, AU*surface_local%SECRD, &
                  DD, NAT, NUMAT, NFIRST, NLAST, &
                  PHI_RHO_SIGMA=phi, RHO=rho)
          END IF
          phi =  AU * phi ! convert to Bohr-1

!!!**** run SCOSMO calculation
          CALL do_cosmo(cosmo_data_local, surface_local, SUM(rho),phi)
          EpolNew = cosmo_data_local%Epol

          EGrad(k,iatm) = (EpolNew - Epol)/cosmo_data_local%minpara%FDStep

!!!**** nonelectrostatic energy
          IF (surface_local%iCDR /= 0  .AND. surface_local%GammaSwitch > 0) THEN
             CALL MolecularNonelecEnergy(Surface_local, EcdrNew)
             CavGrad(k,iatm) = (EcdrNew - ECDR)/cosmo_data_local%minpara%FDStep 
         END IF
          CALL DeallocateSurfaceArrays(surface_local)
          CALL DeallocateCosmoArrays(cosmo_data_local)
          Deallocate(phi)

       END DO ! k
    END DO ! iatm

  END SUBROUTINE MNDONumGradientSlow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE SCosmoDGrad
