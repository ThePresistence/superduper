MODULE MinimizeMod
!!!****h* MinimizeMod
!!!
!!! NAME
!!!     MinimizeMod -- an energy minimization module for COSMO
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Jana Khandogin and Darrin York
!!! CREATION DATE
!!!     2002
!!! DESCRIPTION
!!!     Solves for the surface polarization charge distribution
!!!     by the preconditioned conjugate gradient technique
!!!***
  USE DataTypes
  USE CosmoMod
  IMPLICIT NONE  

  PRIVATE
  PUBLIC :: AnalytGradient, NumGradientFast, NumGradientSlow      

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AnalytGradient(surface,surface_charge,&
       cosmof, minpara, GRAD_SIGMA_RHO, GRAD_RHO_SIGMA)
!!!****f* MinimizeMod/CptAnalytGradient
!!!
!!! DESCRIPTION
!!!     computes the analytic gradients of Energy (negative of force)
!!! USES
!!!     recursive bisection method to evaluate electric fields
!!! INPUT
!!!     cosmof        - dielectric scale factor
!!!     q             - solute charges
!!!     surface       - surface type data
!!!     surface_charge- charges on the COSMO surface
!!!     GRAD_SIGMA_RHO - Energy Gradients (negative of force) at rho due to sigma points
!!!     GRAD_RHO_SIGMA - Energy Gradients (negative of force) at sigma points due to rho
!!! OUTPUT
!!!     GRAD_SIGMA_RHO - Gradient of Epol (negative of force) at rho positions
!!! SCRATCH
!!!     SigmaCoord    - charges and coordinates of srf points
!!!     SigmaFldSigma - potential and electric field at srf points due to srf charges 
!!!     A0(ik)        - A_0_ik,ik as in Eq 74, where ik is the double index
!!!     Diag          - diag part of gradient
!!!     dr_hat        - grad r (Eq 76 without Delta functions)
!!!     dSwitch       - dS_wf/dr (Eq 75)
!!!     GradSwf       - dr_hat * dSwitch/Swf
!!!***

    USE UtilitiesMod
    USE ConstantsMod
    USE SurfaceMod
    USE RBM 
    USE ErrorMod
    USE Cosmo_util

    TYPE(surface_t) :: surface
    TYPE(minpara_t) :: minpara

    INTEGER(I4B) :: nse, NseAtm, nkatm, first, last, ikatm, k, ik, stat_alloc,index
    REAL(SP), INTENT(IN)  :: cosmof
    REAL(SP), INTENT(IN)  :: surface_charge(:)
    REAL(SP), INTENT(INOUT)  :: GRAD_SIGMA_RHO(:,:)
    REAL(SP), INTENT(IN)  ::  GRAD_RHO_SIGMA(:,:)
    REAL(SP), POINTER :: SigmaCoord(:,:) => NULL()
    REAL(SP), POINTER :: SigmaPotSigma(:) => NULL(),SigmaFldSigma(:,:) => NULL()
    REAL(SP) :: f, Diag, A0, dS, zeta, sigma2
    REAL(SP), PARAMETER :: SQRT_TWO_over_PI = SQRT2/SQRT_PI

    nse = assert_eq(SIZE(surface%SECRD,2), SIZE(Surface_charge), 'AnalytGradient')
    nkatm = SIZE(Surface%ASAINDEX)    

    IF (cosmof/=0.0_SP) THEN
       f = 1/cosmof
    ELSE
       RETURN
    END IF

    ALLOCATE(SigmaCoord(5,nse),SigmaPotSigma(nse),SigmaFldSigma(3,nse), STAT=stat_alloc)

    SigmaCoord(1,:) =  surface_charge(:)
    SigmaCoord(2:4,:) = Surface%SECRD(1:3,:)
    SigmaCoord(5,:) = Surface%cosmozeta(:)
    SigmaPotSigma = 0.0_SP 
    SigmaFldSigma = 0.0_SP 

!!!**** calculate electric field at srf points due to grad(A)*Sigma and grad(B)*Rho 
    SELECT CASE(minpara%multipole)
    CASE(0)
       CALL DirectCoulomb(SigmaCoord,SigmaPotSigma,SigmaFldSigma)
    CASE(1)
       CALL RecursiveBisectionMethod(SigmaCoord, SigmaPotSigma,SigmaFldSigma, &
            WS=minpara%WS,L=minpara%L)
    CASE DEFAULT
       CALL error('illegal Ifmultipole in CptAnalytGradient')
    END SELECT

!!!**** sum up Gradients (-forces) on atoms with exposed surface elements 
    first=1
    AtomLoop: DO ikatm = 1, Nkatm
       index = surface%keptIndex(ikatm)
       last = surface%asaindex(ikatm)
       NseAtm = last - first + 1
       DO k = 1, 3
          Diag = 0.0_SP
!!!**** non-switching gradients due to sigma*A*sigma
          IF (NseAtm /=0) THEN
             GRAD_SIGMA_RHO(k,index) = GRAD_SIGMA_RHO(k,index) + &
                  f*DOT_PRODUCT(SigmaCoord(1,first:last), SigmaFldSigma(k,first:last)) + & 
                  SUM(GRAD_RHO_SIGMA(k,first:last))
          END IF
!!!**** gradient due to switching: since gradient of A_diag is zero when there is
!!!     no switching. See Eq 74 for grad(A_diag). Refer to my notes for the 
!!!     corrected Eq. 74
          IF (surface%GammaSwitch > 0) THEN 
!!!**** first term in my notes
             IF (NseAtm /=0) THEN
                DO ik = first, last
                   zeta = surface%cosmozeta(ik)
                   IF(org_Sik_in_A0_diag)THEN
                      A0 = SQRT_TWO_over_PI*zeta/surface%sescale(ik)
                   ELSE
                      A0 = 0.5_SP*SQRT_TWO_over_PI*zeta/SQRT(surface%sescale(ik))
                   ENDIF
                   dS = SUM(Surface%GradSwf(ik,:,k))
                   sigma2 = SigmaCoord(1,ik)*SigmaCoord(1,ik)
                   sigma2 = sigma2 * A0 * dS 
                   Diag = Diag - sigma2
                END DO
             END IF
!!!**** second term in my notes
             DO ik = 1, nse
                zeta = surface%cosmozeta(ik)
                IF(org_Sik_in_A0_diag)THEN
                   A0 = SQRT_TWO_over_PI*zeta/surface%sescale(ik)
                ELSE
                   A0 = 0.5_SP*SQRT_TWO_over_PI*zeta/SQRT(surface%sescale(ik))
                ENDIF
                dS =Surface%GradSwf(ik,ikatm,k)
                sigma2 = SigmaCoord(1,ik)*SigmaCoord(1,ik)
                sigma2 = sigma2 * A0 * dS
                Diag = Diag + sigma2 
             END DO
             Diag = Diag * Surface%PEXP
          END IF ! ifswitch

!!!**** add the grad(A_diag) to the non-switching gradient 
          GRAD_SIGMA_RHO(k,index) = GRAD_SIGMA_RHO(k,index) + 0.5_SP*f*Diag 
       END DO !k

       first = last + 1

    END DO AtomLoop

    DEALLOCATE(SigmaCoord,SigmaPotSigma,SigmaFldSigma)

  END SUBROUTINE AnalytGradient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE NumGradientSlow(cosmo_data, surface, Crd,radii,rho)
!!!****f* MinimizeMod/NumGradientSlow
!!!
!!! Description
!!!     computes the numerical gradients of Epol and Ecav
!!!     using frozen surface charge approach
!!! Uses
!!!     recursive bisection method to evaluate electric fields
!!! Scratch
!!!     SigmaFldSigma - potential at srf points due to srf charges 
!!!     SigmaFldRho   - potential at srf points due to solute charges
!!!     RhoFld        - potential at solute atom positions due to srf charges
!!!
!!!***

    USE CosmoMod
    USE SurfaceMod
    USE cosmo_helper

    ! this should be intent(in)
    ! but we then cant point to subcomponents
    TYPE(cosmo_t), intent(inOUT) :: cosmo_data
    TYPE(Surface_t), intent(in) :: surface

    TYPE(cosmo_t) :: cosmo_data_local
    TYPE(Surface_t) :: surface_local

    REAL(SP), INTENT(IN) :: Crd(:,:), radii(:),rho(:)
    REAL(SP) ::  EpolNew, EcdrNew
    REAL(SP) :: CrdNew(3,SIZE(crd,2))
    REAL(SP), pointer :: phi(:)
    INTEGER(I4B) :: nse, natm, iatm, k

    nse = SIZE(surface%cosmozeta)
    Natm = SIZE(crd,2)

    call copy_cosmo_scalar_inputs(cosmo_data,cosmo_data_local)
    call copy_surface_scalar_inputs(surface,surface_local)

    IF (ASSOCIATED(cosmo_data%EGrad)) DEALLOCATE(cosmo_data%EGrad)
    ALLOCATE(cosmo_data%EGrad(1:3,1:natm))
    cosmo_data%EGrad = 0.0_SP

    IF (ASSOCIATED(cosmo_data%CDRGrad)) DEALLOCATE(cosmo_data%CDRGrad)
    ALLOCATE(cosmo_data%CDRGrad(1:3,1:natm))
    cosmo_data%CDRGrad = 0.0_SP

    DO iatm = 1, Natm
       DO k = 1, 3
          CrdNew(:,:) = crd(:,:)
          CrdNew(k,iatm) = CrdNew(k,iatm) + cosmo_data%minpara%FDstep
          CALL SASurface(CrdNew, radii, Surface_local)

!!!**** run COSMO calculation
          ! phi allocated in cpt_phirhosigma_for_pt
          CALL cpt_phirhosigma_for_pt(CrdNew,rho, &
               surface_local%secrd, surface_local%cosmozeta, phi)

          CALL do_cosmo(cosmo_data_local, surface_local, SUM(rho), phi)

          EpolNew = cosmo_data_local%Epol
          cosmo_data%EGrad(k,iatm) = (EpolNew - cosmo_data%Epol)/cosmo_data%minpara%FDStep

!!!**** cavitation energy
          IF (surface%iCDR == 1 .AND. surface%GammaSwitch > 0.0_SP) THEN
             ! CALL SASurface(CrdNew, radii,cdrsurf_local)
             CALL MolecularNonelecEnergy(Surface_local, EcdrNew)
             cosmo_data%CDRGrad(k,iatm) = (EcdrNew - cosmo_data%ECDR)/cosmo_data%minpara%FDStep 
             ! call deallocateSurfaceArrays(cdrsurf_local)
          ELSE IF (surface%iCDR >= 2 .AND. surface%GammaSwitch > 0.0_SP) THEN
             WRITE(6,*)'Not implemented for SCOSMO alone'
          END IF
          CALL deallocateSurfaceArrays(Surface_local)
          CALL DeallocateCosmoArrays(cosmo_data_local)
          deallocate(phi)
       END DO ! k
    END DO ! iatm

  END SUBROUTINE NumGradientSlow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE NumGradientFast(cosmo_data, surface, Crd, radii, rho,Success)
!!!****f* MinimizeMod/NumGradientFast
!!!
!!! Description
!!!     computes the numerical gradients of Epol and Ecav
!!! Uses
!!!     recursive bisection method to evaluate electric fields
!!! Scratch
!!!     SigmaFldSigma - potential at srf points due to srf charges
!!!     SigmaFldRho   - potential at srf points due to solute charges
!!!     RhoFld        - potential at solute atom positions due to srf charges
!!!
!!!***
    USE ErrorMod
    USE CosmoMod
    USE SurfaceMod
    USE RBM
    USE cosmo_util

    TYPE(cosmo_t), INTENT(INOUT) :: cosmo_data
    TYPE(surface_t), INTENT(IN) :: surface  
    TYPE(surface_t) :: TmpSurface
 
    INTEGER(I4B) :: nse, natm, iatm, k

    REAL(SP), INTENT(IN) :: Crd(:,:),radii(:),rho(:)
    REAL(SP) :: CrdNew(3,SIZE(crd,2)), f, EpolNew, EcdrNew
    REAL(SP), POINTER :: SigmaCoord(:,:), SigmaPotSigma(:), &
         SigmaPotRho(:), RhoCoord(:,:), RhoPot(:)

    LOGICAL(LGD), INTENT(OUT) :: Success

    nse = SIZE(surface%cosmozeta)
    Natm = SIZE(crd,2)

    IF (ASSOCIATED(cosmo_data%EGrad)) DEALLOCATE(cosmo_data%EGrad)
    ALLOCATE(cosmo_data%EGrad(1:3,1:natm))
    cosmo_data%EGrad = 0.0_SP

    IF (ASSOCIATED(cosmo_data%CDRGrad)) DEALLOCATE(cosmo_data%CDRGrad)
    ALLOCATE(cosmo_data%CDRGrad(1:3,1:natm))
    cosmo_data%CDRGrad = 0.0_SP

    IF (cosmo_data%f /= 0.0_SP) THEN
       f = 1/cosmo_data%f
    ELSE
       RETURN
    END IF

    call copy_surface_scalar_inputs(surface,TmpSurface)

    ALLOCATE(SigmaCoord(1:5,1:nse), SigmaPotSigma(1:nse), &
         SigmaPotRho(1:nse), RhoCoord(1:4,1:natm), RhoPot(1:natm))

    SigmaCoord(1,:) =  cosmo_data%surface_charge(:)
    RhoCoord(1,:) =  rho(:)

    DO iatm = 1, Natm
       DO k = 1, 3
          CrdNew(:,:) = crd(:,:)
          CrdNew(k,iatm) = CrdNew(k,iatm) + cosmo_data%minpara%FDstep

          CALL SASurface(CrdNew, radii, TmpSurface)

          IF (SIZE(TmpSurface%cosmozeta) /= nse) THEN
             IF (cosmo_data%verbose) &
                  WRITE(6,*) 'surface sigma vector length changes. do slow gradient' 
             CALL deallocateSurfaceArrays(TmpSurface)
             DEALLOCATE(SigmaCoord, SigmaPotSigma, SigmaPotRho, RhoCoord, RhoPot)
             Success = .FALSE.
             RETURN 
          END IF

          SigmaCoord(2:4,:) = TmpSurface%SECRD(1:3,:)
          SigmaCoord(5,:) = TmpSurface%cosmozeta(:)
          RhoCoord(2:4,:) = CrdNew(1:3,:)

          SigmaPotSigma = 0.0_SP
          SigmaPotRho = 0.0_SP
          RhoPot=0.0_SP
          SELECT CASE(cosmo_data%minpara%Multipole)
          CASE(0)
             write(*,*)'numerical gradients with direct Coulomb...'
             !CALL DirectCoulomb(RhoCoord,RhoFldRho)
             CALL DirectCoulomb(SigmaCoord,SigmaPotSigma)
             CALL DirectCoulomb(SigmaCoord,SigmaPotRho,RhoCoord,RhoPot,.false.)

          CASE(1)
             write(*,*)'numerical gradients with recursive bisection...'
             !CALL RecursiveBisectionMethod(RhoCoord,RhoFldRho,minpara%WS,minpara%L)
             CALL RecursiveBisectionMethod(SigmaCoord,SigmaPotSigma,cosmo_data%minpara%WS,cosmo_data%minpara%L)
             CALL RecursiveBisectionMethod(SigmaCoord,SigmaPotRho,RhoCoord,&
                  RhoPot,cosmo_data%minpara%WS,cosmo_data%minpara%L)
          CASE DEFAULT
             CALL ERROR('illegal IFmultipole')
          END SELECT

          IF (surface%GammaSwitch > 0) THEN
             IF (Surface%verbose) WRITE(6,*)'Correcting for switching'
             CALL switch_pot(SigmaPotSigma(:), SigmaCoord(1,:), TmpSurface%Sescale, &
                  TmpSurface%cosmozeta)
          END IF

!!!**** gradient due to 0.5*sigma(T)*A*sigma + sigma(T)*B*rho
          EpolNew = 0.5_SP*f*DOT_PRODUCT(SigmaCoord(1,:),SigmaPotSigma(:)) + &
               DOT_PRODUCT(RhoCoord(1,:), RhoPot(:))

          cosmo_data%EGrad(k,iatm) = (EpolNew - cosmo_data%Epol)/cosmo_data%minpara%FDStep

!!!**** cavitation energy
          IF (ASSOCIATED(cosmo_data%CDRGrad) .AND. TmpSurface%iCDR == 1) THEN
             ! Call  SASurface(CrdNew, radii, cdrSurf_new)
             CALL MolecularNonelecEnergy(TmpSurface, EcdrNew)
             cosmo_data%CDRGrad(k,iatm) = (EcdrNew - cosmo_data%ECDR)/cosmo_data%minpara%FDStep
             !CALL deallocateSurfaceArrays(cdrSurf_new)
          ELSE IF (ASSOCIATED(cosmo_data%CDRGrad) .AND. TmpSurface%iCDR >= 2) THEN
             WRITE(6,*) 'Non implemented for Atomic Nonelectrostatics Yet'
          END IF

          CALL deallocateSurfaceArrays(TmpSurface)

       END DO ! k
    END DO ! iatm

    Success = .TRUE.
    

    DEALLOCATE(SigmaCoord, SigmaPotSigma, SigmaPotRho, RhoCoord, RhoPot)

  END SUBROUTINE NumGradientFast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE MinimizeMod


