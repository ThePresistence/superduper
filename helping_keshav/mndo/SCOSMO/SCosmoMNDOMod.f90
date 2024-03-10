MODULE SCosmoMNDOMod
!!!****h* SCosmoMNDOMod
!!!
!!! NAME
!!!    CosmoMNDOMod -- a module for computing COSMO corrections for
!!! MNDO 
!!!    or any other semiempirical program
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Jana Khandogin, Brent Gregersen
!!! CREATION DATE
!!!     2003
!!! DESCRIPTION
!!!     computes dielectric correction to the Fock matrix due to
!!! solvent 
!!!     reaction field in the smooth COSMO formulation. 
!!!***

  USE DataTypes
  USE CosmoMod
  USE SurfaceMod
  USE MinimizeMod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: InitSCOSMO, AddFckNew, CptSCosmoE, Cpt_bmat_Phi_MNDO,&
       & CptGradMNDO


!!!****s* module variables
!!! These could be put into a common block if ported to f77
  ! Variable storage for cosmo calculation
  TYPE(cosmo_t),pointer, SAVE, public :: cosmo_mndo =>NULL()
  TYPE(surface_t),pointer, SAVE,public :: esurf_mndo =>NULL()

  LOGICAL(LGD), SAVE, public :: ifGaus=.TRUE.
  LOGICAL(LGD), SAVE, public :: EXTSURF = .FALSE.
  INTEGER(I4B), SAVE, public :: iRadii=0
  REAL(SP),pointer,save, public :: radii_main(:) => NULL()

  REAL(SP), pointer,save :: vRadii(:,:) => NULL()
  REAL(SP), pointer,save :: ASAP(:,:) => NULL()
  REAL(SP), pointer,save, public :: rho(:) =>NULL()
  REAL(SP), pointer,save :: Bmat(:,:)=>NULL()
  REAL(SP), pointer,save :: phi_rho_sigma(:)=>NULL()
  REAL(SP), public, save :: QEXTERN = 0.0_SP
  REAL(SP), pointer, save, public :: phi_rho_sigma_cons(:)=>NULL()
  REAL(SP), pointer, save, public :: grad_surf(:,:)=>NULL()
  REAL(SP), pointer, save, public :: AtomCharges(:) =>NULL()
  INTEGER(I4B), save :: NQbasis = 0  ! Number of charge basis
  ! functions for current molecule
!!!***

CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AddFckNew(F,P,COORD,NUMAT,&
       NAT,NFIRST,NLAST,TORE,IP,IP1,IP2,DD,NITER,Epol,Ecdr)
!!!****f* SCosmoMNDOMod/AddFCKNew
!!!
!!! DESCRIPTION
!!!     computes the dielectric correction to the atom block diagonal 
!!!     elements of the Fock matrix due to the reaction field of the 
!!!     electronic charge
!!! CALLED BY
!!!     ScosmoAddFck - an interface f90 subroutine called by ITER
!!!                    in MNDO
!!!***
    USE ConstantsMod
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN):: NUMAT, NITER
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    INTEGER, INTENT(IN):: IP(:),IP1(:),IP2(:)
    INTEGER(I4B):: iQ, JM
    REAL(SP), INTENT(IN):: TORE(:), DD(:,:)
    REAL(SP), INTENT(IN) :: P(:), COORD(:,:)
    REAL(SP), INTENT(INOUT) :: F(:)
    REAL(SP), INTENT(OUT) :: Epol,Ecdr
    REAL(SP) :: FACTOR, AU

    AU = AU_DISTANCE_IN_ANGSTROM

!!!**** remain in Angstrom and compute smooth COSMO solvation energy Epol
    CALL CptSCosmoE(P, COORD, NUMAT, NAT, NFIRST, NLAST, TORE, IP, IP1, &
         IP2, DD, .FALSE., NITER, Epol, Ecdr)

!!!**** convert to Angstrom and compute reaction field potential phi_sigma_rho
    IF (.NOT. ASSOCIATED(cosmo_mndo%phi_sigma_rho)) THEN
       ALLOCATE(cosmo_mndo%phi_sigma_rho(NQbasis))
    END IF
    IF (ASSOCIATED(Bmat)) THEN
       cosmo_mndo%phi_sigma_rho = MATMUL(cosmo_mndo%surface_charge, &
            Bmat)
    ELSE
       IF (ifGaus) THEN
          CALL Cpt_bmat_Phi_MNDO(COORD, AU * esurf_mndo%SECRD, &
               DD, NAT, NUMAT, NFIRST, NLAST, &
               PHI_SIGMA_RHO=cosmo_mndo%phi_sigma_rho, &
               SIGMA=cosmo_mndo%surface_charge, &
               ZETA=ANGSTROM * esurf_mndo%COSMOZETA)
       ELSE
          CALL Cpt_bmat_Phi_MNDO(COORD, AU * esurf_mndo%SECRD, &
               DD, NAT, NUMAT, NFIRST, NLAST, &
               PHI_SIGMA_RHO=cosmo_mndo%phi_sigma_rho, &
               SIGMA=cosmo_mndo%surface_charge)
       END IF
       ! convert to Bohr-1
       cosmo_mndo%phi_sigma_rho = AU * cosmo_mndo%phi_sigma_rho
    END IF

!!!**** compute dielectric corrected Fock matrix
!!!**** convert Angstrom-1 to EV
    FACTOR = 1/ELECTRON_VOLT

!!!**** add reaction field potential to the atom-block diagonal 
!!!     part of Fock matrix
    DO iQ=1,NQbasis
       JM     = IP(iQ)
       F(JM) = F(JM) - FACTOR * cosmo_mndo%phi_sigma_rho(iQ)
    END DO

  END SUBROUTINE AddFckNew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CptSCosmoE(P, COORD, NUMAT, NAT, NFIRST, NLAST, &
       TORE, IP, IP1, IP2, DD, IFLAST,SCFITER, EPOL, ECDR)
!!!***
!!! NAME
!!!     SCosmoE - A f95 routine to obtain smooth COSMO solvation Energy
!!! CALLED BY
!!!     AddFckNew
!!!     externally: SCFCAL in MNDO
!!! INPUT (from MNDO)
!!!     IP(LMI)  - INDICES OF UNIQUE ONE-CENTER AO PAIRS = I*(I-1)/2+J
!!!     IP1(LMI) - INDEX OF FIRST  AO IN THE ONE-CENTER PAIR = I (COUL)
!!!     IP2(LMI) - INDEX OF SECOND AO IN THE ONE-CENTER PAIR = J (COUL)
!!!     DD       - 
!!!     COORD, RADII - In angstrom
!!! OUTPUT (to MNDO)
!!!     Epol, Ecav - optional, in hartree
!!!***

    USE ConstantsMod, ONLY: ANGSTROM, AU_DISTANCE_IN_ANGSTROM
    USE MolecularPropertiesMod
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN):: NUMAT
    INTEGER(I4B), INTENT(IN):: SCFITER ! SCF iter #
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    INTEGER(I4B), INTENT(IN):: IP(:),IP1(:),IP2(:)
    LOGICAL, INTENT(IN) :: IfLast

    INTEGER(I4B):: Nse, iQ, JM, i

    REAL(SP), INTENT(IN):: COORD(:,:), TORE(:), DD(:,:), P(:)
    REAL(SP) :: PCosmo(NQbasis)
    REAL(SP), INTENT(OUT) :: EPOL, ECDR
    REAL(SP),parameter :: AU = AU_DISTANCE_IN_ANGSTROM

!!!**** extract atom-block diagonal elements from the density matrix 
!!!      to form the electronic part of the COSMO vector rho
    DO iQ=1,NQbasis
       JM     = IP(iQ)
       PCosmo(iQ)   = P(JM)
!!!**** copied from MNDO/ADDFCK
       IF (IP1(iQ) /= IP2(iQ)) THEN
          PCosmo(iQ) = Pcosmo(iQ)*2.0_SP
       END IF
    END DO

!!!**** modify radii if asked
    IF (iRadii >= 2 ) THEN ! varying radii
       !###########################################################
       ! Experimental Code Segment for charge dependent radii
       !###########################################################
       CALL DeallocateSurfaceArrays(esurf_mndo)
       CALL DeallocateCosmoArrays(cosmo_mndo)
       IF(associated(Bmat))deallocate(Bmat)
       IF(associated(phi_rho_sigma))deallocate(phi_rho_sigma)

!!!**** Add the core charges to (S,S) and obtain rho
       CALL CptRho(PCosmo, TORE, NUMAT, NAT, NFIRST, NLAST, &
            rho, AtomCharges)
       ! variable radii

       DO i=1,NUMAT

          SELECT CASE(iRadii)
          CASE (2) !direct radii
             radii_main(i) = vRadii(1,i) & 
                  + vRadii(2,i)*AtomCharges(i) &
                  + vRadii(3,i)*AtomCharges(i)*AtomCharges(i)
             IF(radii_main(i) <= 0)THEN
                WRITE(6,*)'Variable Radii Negative for atom ',i,&
                     ' setting radii to 0.00:SCosmoMNDOMod'
                radii_main(i) = 0.0_SP
             END IF
          CASE(3) !inverse radii
             radii_main(i) = ( 1.0_SP/vRadii(1,i) + vRadii(2,i)*AtomCharges(i) &
                  + vRadii(3,i)*AtomCharges(i)*AtomCharges(i) ) ** (-1.0_SP)
             IF(radii_main(i) <= 0)THEN
                WRITE(6,*)'Variable Radii Negative for atom ',i,&
                     ' setting radii to 0.00:SCosmoMNDOMod'
                radii_main(i) = 0.0_SP
             END IF
          CASE(4) !exponential radii
             radii_main(i) = vRadii(1,i) * EXP(-1.0_SP * &
                  (vRadii(2,i)*AtomCharges(i) &
                  + vRadii(3,i)*AtomCharges(i)*AtomCharges(i)))
          END SELECT
       END DO

       CALL SASurface(ANGSTROM*COORD, radii_main, esurf_mndo)
       Nse = SIZE(esurf_mndo%SECRD,2)
       ALLOCATE(phi_rho_sigma(Nse))

       IF (ifGaus) THEN
          CALL Cpt_bmat_Phi_MNDO(COORD, AU*esurf_mndo%SECRD, &
               DD, NAT, NUMAT, NFIRST, NLAST, &
               PHI_RHO_SIGMA=phi_rho_sigma, RHO=rho, &
               ZETA=ANGSTROM * esurf_mndo%COSMOZETA)
       ELSE
          CALL Cpt_bmat_Phi_MNDO(COORD, AU*esurf_mndo%SECRD, &
               DD, NAT, NUMAT, NFIRST, NLAST, &
               PHI_RHO_SIGMA=phi_rho_sigma, RHO=rho)
       END IF
       ! convert to Bohr-1
       phi_rho_sigma =  AU * phi_rho_sigma 
    ELSE ! fixed radii
       !###########################################################
       ! Radii are NOT charge dependant, we can store surface during SCF
       !###########################################################
       IF (SCFITER==0) THEN
          IF (cosmo_mndo%verbose) THEN
             WRITE(6,*) 'SCFITER= ',SCFITER, &
                  ': generate surface and compute Bmat'
          END IF
          IF(associated(Bmat))deallocate(Bmat)
          IF(associated(phi_rho_sigma))deallocate(phi_rho_sigma)
          IF(.NOT.EXTSURF)THEN
             ! We dont want to regenerate the surface if it was 
             ! created externally (i.e. for RxnFld calculations in QM/MM)
             CALL DeallocateSurfaceArrays(esurf_mndo)
             CALL DeallocateCosmoArrays(cosmo_mndo)
             CALL SASurface(ANGSTROM*COORD, radii_main, esurf_mndo)
          ENDIF

          Nse = SIZE(esurf_mndo%SECRD,2)
          ALLOCATE(Bmat(Nse,NQbasis))
          ALLOCATE(phi_rho_sigma(Nse))

          ! Compute surface-density interaction matrix
          IF (ifGaus) THEN
             CALL Cpt_bmat_Phi_MNDO(COORD, AU*esurf_mndo%SECRD, &
                  DD, NAT, NUMAT, NFIRST, NLAST, BMAT=Bmat, &
                  ZETA=ANGSTROM * esurf_mndo%COSMOZETA)
          ELSE
             CALL Cpt_bmat_Phi_MNDO(COORD, AU*esurf_mndo%SECRD, &
                  DD, NAT, NUMAT, NFIRST, NLAST, BMAT=Bmat)
          END IF
          Bmat = AU * Bmat ! convert to Bohr-1
       ELSE
          IF (cosmo_mndo%verbose) THEN
             WRITE(6,*) 'SCFITER', SCFITER, &
                  'use saved surface and Bmat'
          END IF
       END IF

!!!**** Add the core charges to (S,S) and obtain rho
       CALL CptRho(PCosmo, TORE, NUMAT, NAT, NFIRST, NLAST, &
            rho, AtomCharges)
       phi_rho_sigma = MATMUL(Bmat, rho)
    END IF

!!! Add in phi at sigma from EXTERNAL calculations with current surface
!!! This is VERY useful for QM/MM type calculations where the surface
!!! responds to
    IF(EXTSURF .AND. associated(phi_rho_sigma_cons))THEN
       phi_rho_sigma = phi_rho_sigma + phi_rho_sigma_cons   
    ENDIF

!!!**** run SCOSMO calculation
    CALL do_cosmo(cosmo_mndo, esurf_mndo, SUM(AtomCharges)+QEXTERN, phi_rho_sigma)   

    Epol =  cosmo_mndo%Epol

    ! cavitation dispersion repulsion terms
    if(esurf_mndo%iCDR == 1)THEN  
       CALL MolecularNonelecEnergy(esurf_mndo, cosmo_mndo%Ecdr)
    else if(esurf_mndo%iCDR == 2)THEN
       CALL AtomicNonelecEnergy(esurf_mndo,ASAP(1,:),ASAP(2,:),cosmo_mndo%Ecdr)
    else
       cosmo_mndo%Ecdr = 0.0_SP
    ENDIF

    Ecdr = cosmo_mndo%Ecdr

    IF (IFLAST) THEN
       IF(LEN_TRIM(esurf_mndo%surf_file) /= 0) & 
            & call visualize_surface(esurf_mndo,&
            trim(esurf_mndo%surf_file) // '.xyz')
    ENDIF

  END SUBROUTINE CptSCosmoE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE CptRho(PCosmo, TORE, NUMAT, NAT, NFIRST, NLAST, Rho, QI)
!!!****f* SCosmoMNDOMod/CptRho
!!!
!!! DESCRIPTION
!!!     This subroutine fills COSMO rho vector with the
!!!     total charges (reduced core+electronic) 
!!!     where reducred core charges are added to the "ss" matrix elements
!!! INPUT
!!!     PCosmo - density matrix in COSMO vector form
!!! OUTPUT
!!!     Rho    - COSMO charge vector (total charge)
!!!     QI     - Mulliken charge vector  
!!!
!!!***  
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN):: NUMAT
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    INTEGER(I4B) :: iQ, i, ii, jj
    REAL(SP), INTENT(IN):: TORE(:)
    REAL(SP), INTENT(IN):: PCosmo(:)
    REAL(SP), INTENT(OUT) :: Rho(:)
    REAL(SP), INTENT(OUT) :: QI(:)

    Rho = 0.0_SP
    QI = 0.0_SP

    iQ = 0
    DO  i=1,NUMAT
       DO ii = NFIRST(I), NLAST(I)
          DO jj = NFIRST(i),ii
             iQ = iQ + 1
             Rho(iQ) = -PCosmo(iQ)
             IF ( (ii==NFIRST(i)) .AND. (jj==ii) ) THEN
                Rho(iQ) = Rho(iQ) + TORE(NAT(i))  ! s,s element
                QI(i) = QI(i) + TORE(NAT(i)) 
             END IF
             IF ( (ii==jj) ) THEN
                QI(i) = QI(i) - PCosmo(iQ)
             END IF
          END DO
       END DO
    END DO ! i

  END SUBROUTINE CptRho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE cpt_Bmat_phi_MNDO(COORD, SECRD, DD, NAT, NUMAT, NFIRST, NLAST, &
       Bmat, PHI_RHO_SIGMA, RHO, PHI_SIGMA_RHO, SIGMA, ZETA)
!!!****f* SCosmoMNDOMod/CptBMNDO
!!!
!!! DESCRIPTION
!!!     compute Bmat in units of Angtrom-1
!!!**
    USE ConstantsMod, ONLY: AU_DISTANCE_IN_ANGSTROM
    USE ErrorMod
    USE UtilitiesMod
    IMPLICIT NONE

    REAL(SP), PARAMETER :: ONE=1.0_SP, TWO=2.0_SP, THREE=3.0_SP, &
         ZERO=0.0_SP, PT5=0.5_SP, SQ3=1.7320508075689_SP

    INTEGER(I4B), INTENT(IN) :: NUMAT
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    REAL(SP), INTENT(IN) :: SECRD(:,:), COORD(:,:), DD(:,:)
    REAL(SP), OPTIONAL, INTENT(IN) :: rho(:), sigma(:),  zeta(:)
    REAL(SP), optional, INTENT(OUT) :: Bmat(:,:),phi_rho_sigma(:),phi_sigma_rho(:)
    REAL(SP) :: Brow(1:45)
    REAL(SP) :: x, y, z, r, r1, r3, r5, ddi, qqi2, AU, zetai, QC, QB,DB, &
         dipolX, dipolY, dipolZ, quadru, Phi0, Phi1, Phi1x, Phi1y, Phi1z, &
         Phi2x, Phi2y, Phi2z, Phi2xy, Phi2xz, Phi2yz
    INTEGER(I4B) :: Nse, ise, iQ, i, nati, nbase, npairs
    LOGICAL(LGD) :: cpt_phi_rho_sigma,cpt_phi_sigma_rho,cpt_Bmat

    Nse = SIZE(SECRD,2)
    AU = AU_DISTANCE_IN_ANGSTROM

    cpt_phi_rho_sigma=.FALSE.
!!!**** compute B*rho ( phi_rho_sigma)
    IF (PRESENT(PHI_RHO_SIGMA) .AND. PRESENT(RHO)) THEN
       ! should be size NSE
       IF(SIZE(phi_rho_sigma) /= Nse) THEN
          WRITE(6,*)'Fatal Error In CptPhiMNDO: phi_rho_sigma was not correct length'
          STOP
       ENDIF
       cpt_phi_rho_sigma = .TRUE.
       phi_rho_sigma = 0.0_sp
    ENDIF

    cpt_phi_sigma_rho=.FALSE.
!!!**** compute B(T)*sigma (phi_sigma_rho)
!!!     atom-block elements of B is multiplied by surface charges from
!!!     the same atom 
    IF (PRESENT(phi_sigma_rho) .AND. PRESENT(sigma)) THEN
       ! should be size NQbasis
       IF(SIZE(phi_sigma_rho) /= NQbasis) THEN
          WRITE(6,*)'Fatal Error In CptPhiMNDO: phi_sigma_rho was not correct length'
          STOP
       ENDIF
       cpt_phi_sigma_rho = .TRUE.
       phi_sigma_rho = 0.0_sp
    ENDIF

    cpt_Bmat=.FALSE.
    IF (PRESENT(Bmat)) THEN
       cpt_Bmat = .TRUE.
       Bmat = 0.0_SP
    ENDIF

    iQ = 0
!!!**** compute Bmat
    IF (PRESENT(ZETA)) THEN ! Gaussian-atomic multipole
       DO i=1,NUMAT
          nati=NAT(i)
          ddi=DD(2,nati)*AU
          qqi2=(DD(3,nati)*AU)**2
          QB  = (DD(4,nati)*AU)**2*PT5
          DB  = DD(5,nati)*AU
          QC  = (DD(6,nati)*AU)**2*PT5             
          nbase = NLAST(i) - NFIRST(i) + 1
          npairs = (nbase*(nbase+1))/2
          IF (nbase == 1) THEN
             DO ise = 1, NSe
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)
                Brow(1) = GaussianMultipole(zetai,r,0)

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) +  &
                        DOT_PRODUCT(Brow(1:1),rho(iQ+1:iQ+1))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+1) = phi_sigma_rho(iQ+1:iQ+1) &
                        + Brow(1:1) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+1) = Brow(1:1)
                ENDIF

             END DO
          ELSE IF (nbase == 4) THEN ! SP basis
             DO ise = 1, Nse
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)

                r1=ONE/r
                r3=r1**3
                r5=r1**5

                Phi0 = GaussianMultipole(zetai,r,0)
                Phi1 = GaussianMultipole(zetai,r,1)
                Phi1x = GaussianMultipole(zetai,r,1,RX=x)
                Phi1y = GaussianMultipole(zetai,r,1,RY=y)
                Phi1z = GaussianMultipole(zetai,r,1,RZ=z)
                Phi2x = GaussianMultipole(zetai,r,2,RX=x)
                Phi2y = GaussianMultipole(zetai,r,2,RY=y)
                Phi2z = GaussianMultipole(zetai,r,2,RZ=z)
                Phi2xy = GaussianMultipole(zetai,r,2,RX=x,RY=y)
                Phi2xz = GaussianMultipole(zetai,r,2,RX=x,RZ=z)
                Phi2yz = GaussianMultipole(zetai,r,2,RY=y,RZ=z)

                Brow(1) = Phi0 

                Brow(2) = ddi*Phi1x
                Brow(4) = ddi*Phi1y
                Brow(7) = ddi*Phi1z

                Brow(3) = Phi0 + THREE*qqi2*Phi2x - qqi2*Phi1
                Brow(6) = Phi0 + THREE*qqi2*Phi2y - qqi2*Phi1
                Brow(10)= Phi0 + THREE*qqi2*Phi2z - qqi2*Phi1

                Brow(5) = THREE*qqi2*Phi2xy
                Brow(8) = THREE*qqi2*Phi2xz 
                Brow(9) = THREE*qqi2*Phi2yz 

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) + &
                        DOT_PRODUCT(Brow(1:10),rho(iQ+1:iQ+10))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+10) = phi_sigma_rho(iQ+1:iQ+10) &
                        + Brow(1:10) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+10) = Brow(1:10)
                ENDIF

             END DO ! Nse
          ELSE ! SPD basis
             DO ise = 1, Nse
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)

                r1=ONE/r
                r3=r1**3
                r5=r1**5

                Phi0 = GaussianMultipole(zetai,r,0)
                Phi1 = GaussianMultipole(zetai,r,1)
                Phi1x = GaussianMultipole(zetai,r,1,RX=x)
                Phi1y = GaussianMultipole(zetai,r,1,RY=y)
                Phi1z = GaussianMultipole(zetai,r,1,RZ=z)
                Phi2x = GaussianMultipole(zetai,r,2,RX=x)
                Phi2y = GaussianMultipole(zetai,r,2,RY=y)
                Phi2z = GaussianMultipole(zetai,r,2,RZ=z)
                Phi2xy = GaussianMultipole(zetai,r,2,RX=x,RY=y)
                Phi2xz = GaussianMultipole(zetai,r,2,RX=x,RZ=z)
                Phi2yz = GaussianMultipole(zetai,r,2,RY=y,RZ=z)

                Brow(1) = Phi0 

                Brow(2) = ddi*Phi1x
                Brow(4) = ddi*Phi1y
                Brow(7) = ddi*Phi1z

                Brow(3) = Phi0 + THREE*qqi2*Phi2x - qqi2*Phi1
                Brow(6) = Phi0 + THREE*qqi2*Phi2y - qqi2*Phi1
                Brow(10)= Phi0 + THREE*qqi2*Phi2z - qqi2*Phi1

                Brow(5) = THREE*qqi2*Phi2xy
                Brow(8) = THREE*qqi2*Phi2xz 
                Brow(9) = THREE*qqi2*Phi2yz 

                Brow(11) = 1.5_SP*QB*(Phi2x-Phi2y)
                Brow(16) = THREE*QB*Phi2xz
                Brow(22) = PT5*SQ3*QB*(TWO*Phi2z -Phi2x -Phi2y)
                Brow(29) = THREE*QB*Phi2yz
                Brow(37) = THREE*QB*Phi2xy

                Brow(12) = DB*Phi1x
                Brow(13) =-DB*Phi1y
                Brow(17) = DB*Phi1z
                Brow(19) = Brow(12)
                Brow(23) =-Brow(12)/SQ3
                Brow(24) = Brow(13)/SQ3
                Brow(25) = TWO*Brow(17)/SQ3
                Brow(31) = Brow(17)
                Brow(32) = -Brow(13)
                Brow(38) = -Brow(13)
                Brow(39) = Brow(12)

                Brow(15) = Phi0 - THREE*QC*Phi2z + QC*Phi1
                Brow(21) = Phi0 - THREE*QC*Phi2y + QC*Phi1
                Brow(28) = Phi0 + THREE*QC*Phi2z - QC*Phi1
                Brow(36) = Phi0 - THREE*QC*Phi2x + QC*Phi1
                Brow(45) = Brow(15)

                Brow(20) = THREE*QC*Phi2xz 
                Brow(26) =-SQ3*QC*(Phi2x-Phi2y)
                Brow(27) = SQ3*QC*Phi2xz
                Brow(33) =-THREE*QC*Phi2yz
                Brow(34) = THREE*QC*Phi2xy
                Brow(35) = SQ3*QC*Phi2yz
                Brow(42) = -Brow(33)
                Brow(43) =-TWO*SQ3*QC*Phi2xy
                Brow(44) = Brow(20)

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) + &
                        DOT_PRODUCT(Brow(1:45),rho(iQ+1:iQ+45))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+45) = phi_sigma_rho(iQ+1:iQ+45) &
                        + Brow(1:45) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+45) = Brow(1:45)
                ENDIF

             END DO ! Nse
          END IF
          iQ = iQ + npairs
       END DO ! NUMAT
    ELSE ! point charge-atomic multipole
       DO i=1,NUMAT
          nati=NAT(i)
          ddi  = DD(2,nati)*AU
          qqi2 = (DD(3,nati)*AU)**2
          QB  = (DD(4,nati)*AU)**2*PT5
          DB  = DD(5,nati)*AU
          QC  = (DD(6,nati)*AU)**2*PT5
          nbase = NLAST(i) - NFIRST(i) + 1
          npairs = (nbase*(nbase+1))/2

          IF (nbase == 1) THEN !s basis
             DO ise = 1, NSe
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)
                r1=ONE/r
                Brow(1) =  r1

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) + &
                        DOT_PRODUCT(Brow(1:1),rho(iQ+1:iQ+1))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+1) = phi_sigma_rho(iQ+1:iQ+1) &
                        + Brow(1:1) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+1) = Brow(1:1)
                ENDIF

             END DO
          ELSE IF (nbase == 4) THEN ! sp basis
             DO ise = 1, Nse
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)
                r1=ONE/r
                r3=r1**3
                r5=r1**5
                Brow(1) = r1

                Brow(2) = x*ddi*r3
                Brow(4) = y*ddi*r3
                Brow(7) = z*ddi*r3

                Brow(3) = r1 + THREE*x**2*qqi2*r5 - qqi2*r3
                Brow(6) = r1 + THREE*y**2*qqi2*r5 - qqi2*r3
                Brow(10)= r1 + THREE*z**2*qqi2*r5 - qqi2*r3

                Brow(5) = THREE*x*y*qqi2*r5 
                Brow(8) = THREE*x*z*qqi2*r5 
                Brow(9) = THREE*z*y*qqi2*r5 

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) + &
                        DOT_PRODUCT(Brow(1:10),rho(iQ+1:iQ+10))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+10) = phi_sigma_rho(iQ+1:iQ+10) &
                        + Brow(1:10) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+10) = Brow(1:10)
                ENDIF

             END DO ! Nse
          ELSE ! spd basis
             DO ise = 1, Nse
                x = SECRD(1,ise)-COORD(1,i)
                y = SECRD(2,ise)-COORD(2,i)
                z = SECRD(3,ise)-COORD(3,i)
                r = SQRT(x*x + y*y + z*z)
                r1=ONE/r
                r3=r1**3
                r5=r1**5

                Brow(1) = r1

                Brow(2) = x*ddi*r3
                Brow(4) = y*ddi*r3
                Brow(7) = z*ddi*r3

                Brow(3) = r1+THREE*x*x*qqi2*r5-qqi2*r3
                Brow(6) = r1+THREE*y*y*qqi2*r5-qqi2*r3
                Brow(10)= r1+THREE*z*z*qqi2*r5-qqi2*r3

                Brow(5) = THREE*x*y*qqi2*r5
                Brow(8) = THREE*x*z*qqi2*r5
                Brow(9) = THREE*y*z*qqi2*r5

                dipolX = x*DB*r3
                dipolY = y*DB*r3
                dipolZ = z*DB*r3
                quadru = QC*r5

                Brow(11) = 1.5_SP*(x*x-y*y)*QB*r5
                Brow(16) = THREE*x*z*QB*r5
                Brow(22) = PT5*SQ3*(TWO*z*z-x*x-y*y)*QB*r5
                Brow(29) = THREE*y*z*QB*r5
                Brow(37) = THREE*x*y*QB*r5
                Brow(12) = dipoLX
                Brow(13) =-dipolY
                Brow(17) = dipolZ
                Brow(19) = dipolX
                Brow(23) =-dipolX/SQ3
                Brow(24) =-dipoly/SQ3
                Brow(25) = TWO*dipolZ/SQ3
                Brow(31) = dipolZ
                Brow(32) = dipolY
                Brow(38) = dipolY
                Brow(39) = dipolX
                Brow(15) = r1- THREE*z*z*quadru + QC*r3
                Brow(21) = r1- THREE*y*y*quadru + QC*r3
                Brow(28) = r1+ THREE*z*z*quadru - QC*r3
                Brow(36) = r1- THREE*x*x*quadru + QC*r3
                Brow(45) = Brow(15)

                Brow(20) = THREE*x*z*quadru
                Brow(26) =-SQ3*(x*x-y*y)*quadru
                Brow(27) = SQ3*x*z*quadru
                Brow(33) =-THREE*y*z*quadru
                Brow(34) = THREE*x*y*quadru
                Brow(35) = SQ3*y*z*quadru
                Brow(42) = -Brow(33)
                Brow(43) =-TWO*SQ3*x*y*quadru
                Brow(44) = Brow(20)

                IF(cpt_phi_rho_sigma)THEN
                   phi_rho_sigma(ise) = phi_rho_sigma(ise) + &
                        DOT_PRODUCT(Brow(1:45),rho(iQ+1:iQ+45))
                ENDIF
                IF(cpt_phi_sigma_rho)THEN
                   phi_sigma_rho(iQ+1:iQ+45) = phi_sigma_rho(iQ+1:iQ+45) &
                        + Brow(1:45) * sigma(ise)
                ENDIF
                IF(cpt_bmat)THEN
                   Bmat(ise,iQ+1:IQ+45) = Brow(1:45)
                ENDIF

             END DO ! Nse
          END IF
          iQ = iQ + npairs
       END DO ! NUMAT
    END IF ! Zeta

  END SUBROUTINE Cpt_Bmat_phi_MNDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CptGradMNDO(COORD, SECRD, NUMAT, DD, NAT, NFIRST, NLAST, &
       GRAD_RHO_SIGMA, GRAD_SIGMA_RHO, RHO, SIGMA, ZETA)
!!!****f* SCosmoMNDOMod/CpFldMNDO
!!!
!!! DESCRIPTION
!!!     compute fld_rho_sigma or fld_sigma_rho
!!!     Point charge-point charge part are taken from Klamt COSMO (in DGRAD)
!!! INPUT
!!!     COORD                    - coordinates for vector rho in ANGSTROM
!!!     SECRD                    - coordinates for surface charges in ANGSTROM
!!!     RHO                      - solute charge vector
!!!     SIGMA                    - surface (screening) charge vector
!!!     ZETA (OPTIONAL)          - Gaussian zeta in ANGSTROM-1
!!! OUTPUT
!!!     Grad_rho_sigma           - sigma*Grad(B)*rho    (gradient at crds of sigma)
!!!     Grad_sigma_rho           - rho*Grad(B(T))*sigma (gradient at crds of atoms)
!!! SCRATCH
!!!     DB(i,j)                  - DB(k,0): only s orbital, k=1,2,3
!!!                                DB(k,1-10): for s,p orbital, k=1,2,3  
!!!***

    USE UtilitiesMod
    USE ConstantsMod, ONLY: AU_DISTANCE_IN_ANGSTROM
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: NUMAT
    INTEGER(I4B), INTENT(IN) :: NAT(:), NFIRST(:), NLAST(:)
    INTEGER(I4B) :: Nse, ise, iatm, nati, nbase, npairs, j, iQ
    REAL(SP), INTENT(IN) :: SECRD(:,:), COORD(:,:), DD(:,:)
    REAL(SP), INTENT(OUT) :: Grad_rho_sigma(:,:), Grad_sigma_rho(:,:)
    REAL(SP), INTENT(IN) :: rho(:), sigma(:)
    REAL(SP), OPTIONAL, INTENT(IN) :: zeta(:)
    REAL(SP) :: AU, x, y, z, r, r1, r2, r3, r4, r5, ddi, qqi2, qqi3r2, tmp, &
         zetai, DX, DY, DZ, DB(0:3,0:45), DXtot, DYtot, DZtot, &
         DXrho, DYrho, DZrho, d(1:3), QC, QB, DDB

!!!*v** dPhi(l,lm,k)= l,lm notation see dGaussianMultipole, k=x,y,z
    REAL(SP) :: dPhi(0:2,0:6,1:3)

    REAL(SP), PARAMETER :: ONE=1.0_SP, THREE=3.0_SP, TWO=2.0_SP, FIVE=5.0_SP, &
         PT5=0.5_SP, SQ3=1.7320508075689_SP

    Nse = SIZE(SECRD,2)
    AU = AU_DISTANCE_IN_ANGSTROM

!!!**** compute sigma*grad(B*rho) (Grad_rho_sigma) and rho*grad(B(T)*sigma)
!!!**** initialize
    Grad_rho_sigma = 0.0_sp; Grad_sigma_rho = 0.0_SP
    DB(0:3,1:45)=0.0_SP
    DB(0,1)=ONE

!!!**** compute gradients of Gaussian-multipole integrals and multiply by -charge
!!!      to obtain energy gradient
    IF (PRESENT(ZETA)) THEN
       iQ = 0
       DO iatm=1,NUMAT
          nati=NAT(iatm)
          ddi=DD(2,nati)*AU
          qqi2=(DD(3,nati)*AU)**2
          nbase = NLAST(iatm) - NFIRST(iatm) + 1
          npairs = (nbase*(nbase+1))/2
          IF ( nbase == 1) THEN ! s basis
             DO ise = 1, NSe
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,iatm)
                y = SECRD(2,ise)-COORD(2,iatm)
                z = SECRD(3,ise)-COORD(3,iatm)
                r = SQRT(x*x + y*y + z*z)

                CALL dGaussianMultipole(zetai,r,0,x,y,z,dPhi(0,0,1), &
                     dPhi(0,0,2),dPhi(0,0,3))
                DB(1:3,0) = dPhi(0,0,1:3)

                Grad_rho_sigma(1:3,ise) = Grad_rho_sigma(1:3,ise) &
                     + DB(1:3,0)*rho(iQ+1)

                Grad_sigma_rho(1:3,iatm) = Grad_sigma_rho(1:3,iatm) &
                     - DB(1:3,0)*sigma(ise)*rho(iQ+1)

             END DO
          ELSE IF (nbase == 4) THEN ! s,p basis
             DO ise = 1, Nse
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,iatm)
                y = SECRD(2,ise)-COORD(2,iatm)
                z = SECRD(3,ise)-COORD(3,iatm)
                r = SQRT(x*x + y*y + z*z)
!!! dPhi(l,lm,k) 
                CALL dGaussianMultipole(zetai,r,0,x,y,z,dPhi(0,0,1),&
                     dPhi(0,0,2),dPhi(0,0,3))     
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,1,1),&
                     dPhi(1,1,2),dPhi(1,1,3),LM=1)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,2,1),&
                     dPhi(1,2,2),dPhi(1,2,3),LM=2)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,3,1),&
                     dPhi(1,3,2),dPhi(1,3,3),LM=3)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,0,1),&
                     dPhi(1,0,2),dPhi(1,0,3))
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,1,1),&
                     dPhi(2,1,2),dPhi(2,1,3),LM=1)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,2,1),&
                     dPhi(2,2,2),dPhi(2,2,3),LM=2)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,3,1),&
                     dPhi(2,3,2),dPhi(2,3,3),LM=3)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,4,1),&
                     dPhi(2,4,2),dPhi(2,4,3),LM=4)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,5,1),&
                     dPhi(2,5,2),dPhi(2,5,3),LM=5)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,6,1),&
                     dPhi(2,6,2),dPhi(2,6,3),LM=6)

                DB(1:3,1) = dPhi(0,0,1:3)

                DB(1:3,2) = dPhi(1,1,1:3)*ddi
                DB(1:3,4) = dPhi(1,2,1:3)*ddi
                DB(1:3,7) = dPhi(1,3,1:3)*ddi

                d(1:3) = dPhi(0,0,1:3) - qqi2*dPhi(1,0,1:3)
                DB(1:3,3) = d(1:3) + THREE*qqi2*dPhi(2,1,1:3)
                DB(1:3,6) = d(1:3) + THREE*qqi2*dPhi(2,2,1:3)
                DB(1:3,10) = d(1:3) + THREE*qqi2*dPhi(2,3,1:3)

                DB(1:3,5) = dPhi(2,4,1:3)*THREE*qqi2
                DB(1:3,8) = dPhi(2,5,1:3)*THREE*qqi2
                DB(1:3,9) = dPhi(2,6,1:3)*THREE*qqi2

                dx = DOT_PRODUCT(DB(1,1:10),rho(iQ+1:iQ+10))
                dy = DOT_PRODUCT(DB(2,1:10),rho(iQ+1:iQ+10))
                dz = DOT_PRODUCT(DB(3,1:10),rho(iQ+1:iQ+10))

                Grad_rho_sigma(1,ise) = Grad_rho_sigma(1,ise) + dx
                Grad_rho_sigma(2,ise) = Grad_rho_sigma(2,ise) + dy
                Grad_rho_sigma(3,ise) = Grad_rho_sigma(3,ise) + dz

                Grad_sigma_rho(1,iatm) = Grad_sigma_rho(1,iatm) - dx * sigma(ise)
                Grad_sigma_rho(2,iatm) = Grad_sigma_rho(2,iatm) - dy * sigma(ise)
                Grad_sigma_rho(3,iatm) = Grad_sigma_rho(3,iatm) - dz * sigma(ise)
             END DO ! Nse
          ELSE ! SPD basis
             QB  = (DD(4,nati)*AU)**2*PT5
             DDB  = DD(5,nati)*AU
             QC  = (DD(6,nati)*AU)**2*PT5

             DO ise = 1, Nse
                zetai = zeta(ise)
                x = SECRD(1,ise)-COORD(1,iatm)
                y = SECRD(2,ise)-COORD(2,iatm)
                z = SECRD(3,ise)-COORD(3,iatm)
                r = SQRT(x*x + y*y + z*z)

                ! dPhi(l,lm,k) 
                CALL dGaussianMultipole(zetai,r,0,x,y,z,dPhi(0,0,1),&
                     dPhi(0,0,2),dPhi(0,0,3))     
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,1,1),&
                     dPhi(1,1,2),dPhi(1,1,3),LM=1)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,2,1),&
                     dPhi(1,2,2),dPhi(1,2,3),LM=2)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,3,1),&
                     dPhi(1,3,2),dPhi(1,3,3),LM=3)
                CALL dGaussianMultipole(zetai,r,1,x,y,z,dPhi(1,0,1),&
                     dPhi(1,0,2),dPhi(1,0,3))
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,1,1),&
                     dPhi(2,1,2),dPhi(2,1,3),LM=1)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,2,1),&
                     dPhi(2,2,2),dPhi(2,2,3),LM=2)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,3,1),&
                     dPhi(2,3,2),dPhi(2,3,3),LM=3)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,4,1),&
                     dPhi(2,4,2),dPhi(2,4,3),LM=4)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,5,1),&
                     dPhi(2,5,2),dPhi(2,5,3),LM=5)
                CALL dGaussianMultipole(zetai,r,2,x,y,z,dPhi(2,6,1),&
                     dPhi(2,6,2),dPhi(2,6,3),LM=6)

                DB(1:3,1) = dPhi(0,0,1:3)

                DB(1:3,2) = dPhi(1,1,1:3)*ddi
                DB(1:3,4) = dPhi(1,2,1:3)*ddi
                DB(1:3,7) = dPhi(1,3,1:3)*ddi

                d(1:3) = dPhi(0,0,1:3) - qqi2*dPhi(1,0,1:3)
                DB(1:3,3) = d(1:3) + THREE*qqi2*dPhi(2,1,1:3)
                DB(1:3,6) = d(1:3) + THREE*qqi2*dPhi(2,2,1:3)
                DB(1:3,10) = d(1:3) + THREE*qqi2*dPhi(2,3,1:3)

                DB(1:3,5) = dPhi(2,4,1:3)*THREE*qqi2
                DB(1:3,8) = dPhi(2,5,1:3)*THREE*qqi2
                DB(1:3,9) = dPhi(2,6,1:3)*THREE*qqi2


                DB(1:3,11) = 1.5_SP*QB*(dPhi(2,1,1:3) - dPhi(2,2,1:3))
                DB(1:3,16) = THREE*QB*dPhi(2,5,1:3)
                DB(1:3,22) = PT5*SQ3*QB*(TWO*dPhi(2,3,1:3) &
                     - dPhi(2,1,1:3) - dPhi(2,2,1:3))
                DB(1:3,29) = THREE*QB*dPhi(2,6,1:3)
                DB(1:3,37) = THREE*QB*dPhi(2,4,1:3)

                DB(1:3,12) = DDB*dPhi(1,1,1:3)
                DB(1:3,13) =-DDB*dPhi(1,2,1:3)
                DB(1:3,17) = DDB*dPhi(1,3,1:3)
                DB(1:3,19) = DB(1:3,12)

                DB(1:3,23) =-DB(1:3,12)/SQ3
                DB(1:3,24) = DB(1:3,13)/SQ3
                DB(1:3,25) = TWO*DB(1:3,17)/SQ3

                DB(1:3,31) = DB(1:3,17)
                DB(1:3,32) =-DB(1:3,13)
                DB(1:3,38) =-DB(1:3,13)
                DB(1:3,39) = DB(1:3,12)

                DB(1:3,15) = dPhi(0,0,1:3) - THREE*QC*dPhi(2,3,1:3) &
                     + QC*dPhi(1,0,1:3)

                DB(1:3,21) = dPhi(0,0,1:3) - THREE*QC*dPhi(2,2,1:3) &
                     + QC*dPhi(1,0,1:3)

                DB(1:3,28) = dPhi(0,0,1:3) + THREE*QC*dPhi(2,3,1:3) &
                     - QC*dPhi(1,0,1:3)

                DB(1:3,36) = dPhi(0,0,1:3) - THREE*QC*dPhi(2,1,1:3) &
                     + QC*dPhi(1,0,1:3)
                DB(1:3,45) = DB(1:3,15)

                DB(1:3,20) = THREE*QC*dPhi(2,5,1:3)
                DB(1:3,26) = -SQ3*QC*(dPhi(2,1,1:3) - dPhi(2,2,1:3))
                DB(1:3,27) =  SQ3*QC*dPhi(2,5,1:3)

                DB(1:3,33) = -THREE*QC*dPhi(2,6,1:3)
                DB(1:3,34) = THREE*QC*dPhi(2,4,1:3)
                DB(1:3,35) = SQ3*QC*dPhi(2,6,1:3)

                DB(1:3,42) = -DB(1:3,33)
                DB(1:3,43) = -TWO*SQ3*QC*dPhi(2,4,1:3)
                DB(1:3,44) = DB(1:3,20)

                dx = DOT_PRODUCT(DB(1,1:45),rho(iQ+1:iQ+45))
                dy = DOT_PRODUCT(DB(2,1:45),rho(iQ+1:iQ+45))
                dz = DOT_PRODUCT(DB(3,1:45),rho(iQ+1:iQ+45))

                Grad_rho_sigma(1,ise) = Grad_rho_sigma(1,ise) + dx
                Grad_rho_sigma(2,ise) = Grad_rho_sigma(2,ise) + dy
                Grad_rho_sigma(3,ise) = Grad_rho_sigma(3,ise) + dz

                Grad_sigma_rho(1,iatm) = Grad_sigma_rho(1,iatm) &
                     - dx * sigma(ise)
                Grad_sigma_rho(2,iatm) = Grad_sigma_rho(2,iatm) &
                     - dy * sigma(ise)
                Grad_sigma_rho(3,iatm) = Grad_sigma_rho(3,iatm) &
                     - dz * sigma(ise)

             END DO ! Nse
          END IF
          iQ = iQ + npairs
       END DO ! NUMAT

!!!**** point-charge multipole interaction
    ELSE 

       do ise = 1,nse
          iQ=0
          do iatm = 1,NUMAT
             nbase = NLAST(iatm)-NFIRST(iatm) + 1
             npairs = (nbase*(nbase+1))/2
             nati = NAT(iatm)

             x = SECRD(1,ise)-COORD(1,iatm)
             y = SECRD(2,ise)-COORD(2,iatm)
             z = SECRD(3,ise)-COORD(3,iatm)
             r = SQRT(x*x + y*y + z*z)
             r1 = ONE/r
             r2 = r1**2
             r3 = r1**3
             r5 = r1**5

             IF (nbase > 1) THEN ! sp basis

                ddi = DD(2,nati)* AU
                qqi2 = (DD(3,nati)*AU)**2

                tmp = ddi * THREE * r2

                DB(0,2) = x * tmp
                DB(0,4) = y * tmp
                DB(0,7) = z * tmp

                qqi3r2 = qqi2 * THREE * r2

                DB(0,3)= ONE +(FIVE* x**2 * r2 -ONE)*qqi3r2
                DB(0,6)= ONE +(FIVE* y**2 * r2 -ONE)*qqi3r2
                DB(0,10)=ONE +(FIVE* z**2 * r2 -ONE)*qqi3r2

                tmp = qqi3r2 * FIVE * r2

                DB(0,5) = tmp * x * y
                DB(0,8) = tmp * x * z
                DB(0,9) = tmp * z * y

                DB(1,2) = ddi
                DB(2,4) = DB(1,2)
                DB(3,7) = DB(1,2)

                tmp = TWO * qqi3r2

                DB(1,3)  = x * tmp
                DB(2,6)  = y * tmp
                DB(3,10) = z * tmp

                DB(1,5) = PT5 * DB(2,6)
                DB(2,5) = PT5 * DB(1,3)
                DB(1,8) = PT5 * DB(3,10)
                DB(3,8) = PT5 * DB(1,3)
                DB(2,9) = PT5 * DB(3,10)
                DB(3,9) = PT5 * DB(2,6)
             END IF

             IF (nbase > 4) THEN ! spd basis

                QB  = (DD(4,nati)*AU)**2*PT5
                DDB  = DD(5,nati)*AU
                QC  = (DD(6,nati)*AU)**2*PT5

                r4 = r2**2

                DB(0,11) = 7.5_SP*QB*(x-y)*(x+y)*r4
                DB(0,16) = 15.0_SP*QB*x*z*r4
                DB(0,22) = SQ3*QB*r2+2.5_SP*SQ3*QB*(TWO*z*z-x*x-y*y)*r4
                DB(0,29) = 15.0_SP*QB*y*z*r4
                DB(0,37) = 15.0_SP*QB*x*y*r4
                DB(0,12) = THREE*DDB*x*r2
                DB(0,13) =-THREE*DDB*y*r2
                DB(0,17) = THREE*DDB*z*r2
                DB(0,19) = DB(0,12)
                DB(0,23) =-DB(0,12)/SQ3
                DB(0,24) = DB(0,13)/SQ3
                DB(0,25) = TWO*DB(0,17)/SQ3
                DB(0,31) = DB(0,17)
                DB(0,32) =-DB(0,13)
                DB(0,38) =-DB(0,13)
                DB(0,39) = DB(0,12)
                DB(0,15) = ONE + QC*(THREE-15.0_SP*z*z*r2)*r2
                DB(0,21) = ONE + QC*(THREE-15.0_SP*y*y*r2)*r2
!!!**** corrected sign
                DB(0,28) = ONE - QC*(THREE-15.0_SP*z*z*r2)*r2

                DB(0,36) = ONE + QC*(THREE-15.0_SP*x*x*r2)*r2
                DB(0,45) = DB(0,15)
                DB(0,20) = 15.0_SP*QC*x*z*r4
                DB(0,26) = -FIVE*SQ3*QC*(x-y)*(x+y)*r4
                DB(0,27) =  FIVE*SQ3*QC*x*z*r4
                DB(0,33) =-15.0_SP*QC*y*z*r4
                DB(0,34) = 15.0_SP*QC*x*y*r4

!!!**** added QC
                DB(0,35) =  FIVE*SQ3*QC*y*z*r4

                DB(0,42) =-DB(0,33)
                DB(0,43) =-10.0_SP*SQ3*QC*x*y*r4
                DB(0,44) = DB(0,20)

                DB(1,11) = THREE*QB*x*r2
                DB(2,11) =-THREE*QB*y*r2
                DB(1,16) = THREE*QB*z*r2
                DB(3,16) = DB(1,11)
                DB(3,22) = THREE*SQ3*QB*z*r2
                DB(2,29) = THREE*QB*z*r2
                DB(3,29) = THREE*QB*y*r2
                DB(1,37) = DB(3,29)
                DB(2,37) = THREE*QB*x*r2
                DB(1,12) = DDB
                DB(2,13) =-DB(1,12)
                DB(3,17) = DB(1,12)
                DB(1,19) = DB(1,12)
                DB(1,23) =-DB(1,12)/SQ3
                DB(2,24) = DB(1,23)
                DB(3,25) = TWO*DDB/SQ3
                DB(3,31) = DB(1,12)
                DB(2,32) = DB(1,12)
                DB(2,38) = DB(1,12)
                DB(1,39) = DB(1,12)

                DB(3,15) =-6.0_SP*QC*z*r2
                DB(2,21) =-6.0_SP*QC*y*r2

!!!**** sign verified
                DB(3,28) = -DB(3,15)

                DB(1,36) =-6.0_SP*QC*x*r2
                DB(3,45) = DB(3,15)
                DB(1,20) = THREE*QC*z*r2
                DB(3,20) = THREE*QC*x*r2
                DB(1,26) =-TWO*SQ3*QC*x*r2
                DB(2,26) = TWO*SQ3*QC*y*r2
                DB(1,27) = SQ3*QC*z*r2
                DB(3,27) = SQ3*QC*x*r2
                DB(2,33) =-DB(1,20)
                DB(3,33) =-THREE*QC*y*r2
                DB(1,34) =-DB(3,33)
                DB(2,34) = THREE*QC*x*r2

!!!**** verified
                DB(2,35) = SQ3*QC*z*r2
                DB(3,35) = SQ3*QC*y*r2

                DB(2,42) = THREE*QC*z*r2
                DB(3,42) = THREE*QC*y*r2
                DB(1,43) =-TWO*SQ3*QC*y*r2
                DB(2,43) =-TWO*SQ3*QC*x*r2
                DB(1,44) = DB(1,20)
                DB(3,44) = DB(2,34)
             END IF

             DXtot =0.0_SP; DYtot=0.0_SP; DZtot=0.0_SP

             DO j = 1, npairs
                DX = (x * DB(0,j)-DB(1,j))* r3
                DY = (y * DB(0,j)-DB(2,j))* r3
                DZ = (z * DB(0,j)-DB(3,j))* r3

                DXrho = DX * rho(iQ+j)
                DYrho = DY * rho(iQ+j)
                DZrho = DZ * rho(iQ+j)

                DXtot = DXtot + DXrho
                DYtot = DYtot + DYrho
                DZtot = DZtot + DZrho

                Grad_rho_sigma(1,ise) = Grad_rho_sigma(1,ise) - DXrho
                Grad_rho_sigma(2,ise) = Grad_rho_sigma(2,ise) - DYrho
                Grad_rho_sigma(3,ise) = Grad_rho_sigma(3,ise) - DZrho
             END DO

             Grad_sigma_rho(1,iatm) = Grad_sigma_rho(1,iatm) &
                  + DXtot * sigma(ise)
             Grad_sigma_rho(2,iatm) = Grad_sigma_rho(2,iatm) &
                  + DYtot * sigma(ise)
             Grad_sigma_rho(3,iatm) = Grad_sigma_rho(3,iatm) &
                  + DZtot * sigma(ise)
             iQ = iQ + npairs

          END DO ! iatm
       END DO ! ise
    END IF

    DO j = 1,3
       Grad_rho_sigma(j,:) = Grad_rho_sigma(j,:) * sigma(:)
    END DO

  END SUBROUTINE CptGradMNDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION GaussianMultipole(zeta, r, l, RX, RY, RZ) RESULT(Phi)
!!!****f* SCosmoMNDOMod/GaussianMultipole
!!!
!!! DESCRIPTION
!!!   computes Gaussian multipole integrals (see my notes)
!!!   \int 1/ABS(r-r')^{l}*gi dr'
!!!   where l is the order of multipole (l=0,1,2,3,4)
!!!   and gi is Gaussian function centered at r_i
!!! INPUT
!!!   zeta - Gaussian zeta 
!!!   r    - distance between the multipole point and the Gaussian center
!!!   l    - order of multipole: 0(point charge), 1(dipole), 2(quadrupole)
!!!   RX,RY,RZ (optional) - components of the vector r-r'. 
!!!   Needed when computing
!!!          \int r_k/ABS(r-r')^{3}*gi dr' (r_k=X,Y,or Z)
!!!          \int r_k^2/ABS(r-r')^{5}*gi dr' (r_k=X,Y,or Z)
!!!          \int r_k*r_k'/ABS(r-r')^{5}*gi dr' (r_k /= r_k')
!!! SCRATCH
!!!   phi_0 - Gaussian-monopole interaction
!!!   dPhi_0 - 1/r*d/dr(Phi_0)
!!!  
!!!***

    USE ConstantsMod
    USE MATHMOD
    USE ErrorMod
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: l 
    REAL(SP), INTENT(IN) ::  Zeta, r
    REAL(SP), OPTIONAL, INTENT(IN) ::  RX, RY, RZ
    REAL(SP) :: Phi
    REAL(SP), PARAMETER :: TWO_D_SQRTPI = 2.0_SP/SQRT_PI, &
         TOL_HIGH = 5.0_SP, TOL_LOW = 1.E-6_SP, HUGE = 75.0_SP, &
         TWO = 2.0_SP, THREE = 3.0_SP
    REAL(SP) :: x, x2, r2, rinv, r2inv, zeta2, zeta4, zeta2r2, gaus, &
         Phi_0, dPhi_0,  d2Phi_0, Phi_1, tmp, r4, r4inv

    x = zeta * r; x2 = x*x
    r2 = r*r
    rinv = 1/r; r2inv = rinv*rinv
    zeta2 = zeta*zeta; zeta4 = zeta2*zeta2; zeta2r2 = zeta2*r2

!!!**** Gaussian-monopole l=0
    IF ( x > TOL_HIGH ) THEN
       Phi_0 = rinv
       dPhi_0 = -Phi_0*r2inv
    ELSE IF ( x > TOL_LOW ) THEN
       Phi_0 = ERF(x)*rinv
       dPhi_0 = (TWO_D_SQRTPI * EXP(-x2) - Phi_0/zeta)*r2inv*zeta
    ELSE
       x2 = x*x
       Phi_0 = Zeta* TWO_D_SQRTPI * ( 1.0_SP -x2/THREE + x2*x2/10.0_SP )
       dPhi_0 = zeta2 * TWO_D_SQRTPI * ( -TWO / THREE &
            + 4.0_SP * x2 / 10.0_SP )
    END IF

    IF (l==0) THEN
       Phi = Phi_0
       RETURN
    END IF

!!!**** If Gaussian < 1.0E-32 then we set it to be zero
    IF (zeta2r2 > HUGE) THEN
       Gaus = 0.0_SP
    ELSE
       Gaus = zeta*zeta2/(SQRT_PI*PI) * EXP(-zeta2r2)
    END IF

!!!**** Gaussian-dipole l=1 or -quadrupole l=2
    SELECT CASE (l)
    CASE(1)
       IF (PRESENT(RX)) THEN
          Phi = -RX * dPhi_0
          RETURN
       ELSE IF (PRESENT(RY)) THEN
          Phi = -RY * dPhi_0
          RETURN
       ELSE IF (PRESENT(RZ)) THEN
          Phi = -RZ * dPhi_0
          RETURN
       ELSE
          Phi = -dPhi_0 - TWO_PI*Gaus
          RETURN
       END IF
    CASE(2)
       Phi_1 = -(dPhi_0 + TWO_PI*Gaus)
       d2Phi_0 = 2*Phi_1
       IF (PRESENT(RX) .AND. PRESENT(RY)) THEN
          tmp = RX*RY*r2inv
          Phi = tmp*(d2Phi_0 - dPhi_0)
          Phi = Phi/THREE
          RETURN
       ELSE IF (PRESENT(RX) .AND. PRESENT(RZ)) THEN
          tmp = RX*RZ*r2inv
          Phi = tmp*(d2Phi_0 - dPhi_0)
          Phi = Phi/THREE
          RETURN
       ELSE IF (PRESENT(RY) .AND. PRESENT(RZ)) THEN
          tmp = RY*RZ*r2inv
          Phi = tmp*(d2Phi_0 - dPhi_0)
          Phi = Phi/THREE
          RETURN
       ELSE IF (PRESENT(RX)) THEN
          tmp = RX*RX*r2inv
          Phi = (1 - tmp)*dPhi_0 + tmp*d2Phi_0 + Phi_1
          Phi = Phi/THREE
          RETURN
       ELSE IF (PRESENT(RY)) THEN
          tmp = RY*RY*r2inv
          Phi = (1 - tmp)*dPhi_0 + tmp*d2Phi_0 + Phi_1
          Phi = Phi/THREE
          RETURN
       ELSE IF (PRESENT(RZ)) THEN
          tmp = RZ*RZ*r2inv
          Phi = (1 - tmp)*dPhi_0 + tmp*d2Phi_0 + Phi_1
          Phi = Phi/THREE
          RETURN
       ELSE
          Phi = -PI*Gaus*(zeta2 + 4.0_SP*r2inv + TWO*zeta4*r4) &
               - THREE*r2inv*dPhi_0
          Phi = Phi/THREE
          RETURN
       END IF
    CASE(3)
       r4 = r2*r2; r4inv = 1/r4
       Phi = PI*Gaus*(- 4.0_SP*zeta4*zeta4*r4 + 8.0_SP*zeta4*zeta2*r2 &
            - 9.0_SP*zeta4 &
            -24.0_SP*zeta2*r2inv - 60.0_SP*r4inv) - 45.0_SP*r4inv*dPhi_0
       Phi = Phi/45.0_SP
       RETURN
    CASE DEFAULT
       CALL ERROR('I cannot do this l')
    END SELECT ! l

  END FUNCTION GaussianMultipole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE dGaussianMultipole(zeta,r,l,rx,ry,rz,dPhix,dPhiy,dPhiz,LM)
!!!****f* SCosmoMNDOMod/GaussianMultipole
!!!
!!! DESCRIPTION
!!!   l=0:
!!!     dPhik = Grad_k {\int [1/ABS(r-r')*gi] dr']}
!!!   l=1:
!!!     dPhik = Grad_k {\int [1/ABS(r-r')^{3}*gi] dr']}
!!!   l=1, lm=1,2,3:
!!!     dPhik = Grad_k {\int [Rk/ABS(r-r')^{3}*gi] dr']}
!!!     where Rk=x,y,and z
!!!   l=2, lm=1,2,3:
!!!     dPhik = Grad_k {\int [Rk^2/ABS(r-r')^{5}*gi] dr']}
!!!     1->x^2; 2->y^2; 3->z^2
!!!     where Rk=x,y,and z
!!!   l=2, lm=4,5,6:
!!!     dPhik = Gradk {\int [Rk1*Rk2/ABS(r-r')^{5}*gi] dr']}
!!!     k1/=k2  
!!!     4->x*y; 5->x*z; 6->y*z
!!! INPUT
!!!   zeta - Gaussian zeta 
!!!   r    - distance between the multipole point and the Gaussian center
!!!   l    - order of multipole: 0(point charge), 1(dipole), 2(quadrupole) ...
!!!   LM   - optional, see above
!!! 
!!!***

    USE ConstantsMod
    USE MathMod
    USE ErrorMod
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: l
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: LM 
    REAL(SP), INTENT(IN) ::  zeta, r, rx,ry,rz
    REAL(SP), INTENT(OUT) :: dPhix, dPhiy, dPhiz
    LOGICAL :: IfLM
    REAL(SP), PARAMETER :: TWO_D_SQRTPI = 2.0_SP/SQRT_PI, &
         TOL_HIGH = 5.0_SP, TOL_LOW = 1.E-6_SP, HUGE = 75.0_SP, &
         TWO = 2.0_SP, THREE = 3.0_SP, SIX =6.0_SP, SIXinv=1.0_SP/6.0_SP, &
         THREEinv=1.0_SP/3.0_SP
    REAL(SP) :: x, x2, r2, rinv, r2inv, zeta2, zeta4, gaus, &
         Phi_0, dPhi, dPhi_0, d2Phi_0, d3Phi_0, Phi1, Phi2x, Phi2y, &
         Phi2z, Phi2xy, Phi2yz, Phi2xz,tmp1, tmp2, tmp3, tmp4, &
         rx2, ry2, rz2, r4inv

    x = zeta * r; x2 = x*x
    r2 = r*r
    rinv = 1.0_SP/r; r2inv = rinv*rinv
    zeta2 = zeta*zeta; zeta4 = zeta2*zeta2

!!!**** Gaussian-monopole l=0
    IF ( x > TOL_HIGH ) THEN
       Phi_0 = rinv
       dPhi_0 = -Phi_0*r2inv
    ELSE IF ( x > TOL_LOW ) THEN
       Phi_0 = ERF(x)*rinv
       dPhi_0 = (TWO_D_SQRTPI * EXP(-x2) - Phi_0/zeta)*r2inv*zeta
    ELSE
       x2 = x*x
       Phi_0 = Zeta* TWO_D_SQRTPI * (1.0_SP -x2/THREE + x2*x2/10.0_SP)
       dPhi_0 = zeta2 * TWO_D_SQRTPI *(- TWO/THREE + 4.0_SP * x2 / 10.0_SP)
    END IF

    IF (l==0) THEN
       dPhix = dPhi_0*RX
       dPhiy = dPhi_0*RY
       dPhiz = dPhi_0*RZ
       RETURN
    END IF

!!!**** If Gaussian < 1.0E-32 then we set it to be zero
    IF (x2 > HUGE) THEN
       Gaus = 0.0_SP
    ELSE
       Gaus = zeta*zeta2/(SQRT_PI*PI) * EXP(-x2)
    END IF

    IF (PRESENT(LM)) THEN
       IfLm = .TRUE.
    ELSE
       IfLm = .FALSE.
    END IF

!!!**** Gaussian-dipole l=1 or -quadrupole l=2
    SELECT CASE (l)
    CASE(1)
       IF (IfLm) THEN
          Phi1 = GaussianMultipole(zeta,r,1)
          Phi2xy = GaussianMultipole(zeta,r,2,RX=rx,RY=ry)
          Phi2yz = GaussianMultipole(zeta,r,2,RY=ry,RZ=rz)
          Phi2xz = GaussianMultipole(zeta,r,2,RX=rx,RZ=rz)
          SELECT CASE (LM)
          CASE(1)
             Phi2x = GaussianMultipole(zeta,r,2,RX=rx)
             dPhix = -THREE*Phi2x + Phi1 
             dPhiy = -THREE*Phi2xy
             dPhiz = -THREE*Phi2xz
             RETURN
          CASE(2)
             Phi2y = GaussianMultipole(zeta,r,2,Ry=ry)
             dPhix = -THREE*Phi2xy 
             dPhiy = -THREE*Phi2y + Phi1 
             dPhiz = -THREE*Phi2yz
             RETURN
          CASE(3)
             Phi2z = GaussianMultipole(zeta,r,2,RZ=rz)
             dPhix = -THREE*Phi2xz 
             dPhiy = -THREE*Phi2yz  
             dPhiz = -THREE*Phi2z + Phi1 
             RETURN
          CASE DEFAULT
             CALL ERROR('I cannot do this LM')
          END SELECT
       ELSE
          dPhi = FOUR_PI*Gaus*(r2inv + zeta2) + THREE*r2inv*dPhi_0
          dPhix = dPhi*rx
          dPhiy = dPhi*ry
          dPhiz = dPhi*rz
          RETURN
       END IF
    CASE(2)
       IF (IfLm) THEN
          r4inv = r2inv**2
          d2Phi_0 = -TWO*dPhi_0 - FOUR_PI*Gaus
          d3Phi_0 = FOUR_PI*Gaus*(rinv + zeta2*r) + THREE*rinv*dPhi_0
          d3Phi_0 = TWO*d3Phi_0

          SELECT CASE (LM)
          CASE(1)
             rx2=rx*rx
             tmp1 =  SIX*(r2 - rx2)
             tmp2 =  r*(r2 + TWO*rx2)
             tmp3 =  TWO*(r2 - THREE*rx2) 
             tmp4 =  tmp3*(-dPhi_0 + d2Phi_0) + tmp2*d3Phi_0 
             tmp4 =  tmp4*r4inv*SIXinv

             dPhix = tmp1*(-dPhi_0 + d2Phi_0) + tmp2*d3Phi_0
             dPhix = dPhix*r4inv*rx*SIXinv
             dPhiy = tmp4*ry
             dPhiz = tmp4*rz
             RETURN
          CASE(2)
             ry2=ry*ry
             tmp1 = SIX*(r2 - ry2)
             tmp2 = r*(r2 + TWO*ry2)
             tmp3 = TWO*(r2 - THREE*ry2) 
             tmp4 = (-dPhi_0 + d2Phi_0)*tmp3 + tmp2*d3Phi_0 
             tmp4 = tmp4*r4inv*SIXinv

             dPhix = tmp4*rx
             dPhiy = tmp1*(-dPhi_0 + d2Phi_0) + tmp2*d3Phi_0
             dPhiy = dPhiy*r4inv*ry*SIXinv
             dPhiz =tmp4*rz
             RETURN
          CASE(3)
             rz2=rz*rz
             tmp1 = SIX*(r2 - rz2)
             tmp2 = r*(r2 + TWO*rz2)
             tmp3 = TWO*(r2 - THREE*rz2) 
             tmp4 = (-dPhi_0 + d2Phi_0)*tmp3 + tmp2*d3Phi_0 
             tmp4 = tmp4*r4inv*SIXinv

             dPhix = tmp4*rx
             dPhiy = tmp4*ry
             dPhiz = tmp1*(-dPhi_0 + d2Phi_0) + tmp2*d3Phi_0
             dPhiz = dPhiz*r4inv*rz*SIXinv
             RETURN
          CASE(4)
             rx2 = rx*rx; ry2 = ry*ry
             tmp1 = r2 - THREE*rx2 
             tmp2 = r2 - THREE*ry2 

             dPhix = (-dPhi_0 + d2Phi_0)*tmp1 + r*rx2*d3Phi_0
             dPhix = dPhix*r4inv*ry*THREEinv
             dPhiy = (-dPhi_0 + d2Phi_0)*tmp2 + r*ry2*d3Phi_0
             dPhiy = dPhiy*r4inv*rx*THREEinv
             dPhiz = THREE*(dPhi_0 - d2Phi_0) + r*d3Phi_0
             dPhiz = dPhiz*r4inv*rx*ry*rz*THREEinv
             RETURN
          CASE(5)
             rx2 = rx*rx; rz2 = rz*rz
             tmp1 = r2 - THREE*rx2 
             tmp3 = r2 - THREE*rz2

             dPhix = (-dPhi_0 + d2Phi_0)*tmp1 + r*rx2*d3Phi_0
             dPhix = dPhix*r4inv*rz*THREEinv
             dPhiy = THREE*(dPhi_0 - d2Phi_0)+ r*d3Phi_0
             dPhiy = dPhiy*r4inv*rx*ry*rz*THREEinv
             dPhiz = (-dPhi_0 + d2Phi_0)*tmp3 + r*rz2*d3Phi_0
             dPhiz = dPhiz*r4inv*rx*THREEinv
             RETURN
          CASE(6)
             ry2 = ry*ry; rz2 = rz*rz
             tmp2 = r2 - THREE*ry2 
             tmp3 = r2 - THREE*rz2

             dPhix = THREE*(dPhi_0 - d2Phi_0) + r*d3Phi_0
             dPhix = dPhix*r4inv*rx*ry*rz*THREEinv
             dPhiy = (-dPhi_0 + d2Phi_0)*tmp2 + r*ry2*d3Phi_0
             dPhiy = dPhiy*r4inv*rz*THREEinv
             dPhiz = (-dPhi_0 + d2Phi_0)*tmp3 + r*rz2*d3Phi_0
             dPhiz = dPhiz*r4inv*ry*THREEinv
             RETURN
          CASE DEFAULT
             CALL ERROR('I cannot do this LM')
          END SELECT
       ELSE
       END IF
    CASE DEFAULT
       CALL error('I cannot do this l')
    END SELECT

  END SUBROUTINE dGaussianMultipole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE InitSCOSMO(NQ,NUMAT, NAT, ELEMNT, SCRQ0, SCRQ1, SCRQ2, SCNET, SCNEB)
!!!****f* SCosmoMNDOMod/InitSCOSMO
!!! DESCRIPTION
!!!     read input and initializes smooth cosmo variables
!!!     initializes radii also
!!!***
    USE ErrorMod
    USE ConstantsMod!, ONLY: ANGSTROM
    USE FileMod
    IMPLICIT NONE

    INTEGER(I4B), PARAMETER :: LMZ=86
    INTEGER(I4B), INTENT(IN) :: NQ,NUMAT,NAT(:)
    CHARACTER(LEN=2), INTENT(IN) :: ELEMNT(:)
    REAL(SP),INTENT(IN) :: SCRQ0(:),SCRQ1(:),SCRQ2(:),SCNET(:),SCNEB(:)
    INTEGER(I4B) :: iatm, I
    REAL(SP) :: RI, RBONDI(LMZ),REMS(LMZ)

!!!****  VAN-DER-WAALS RADII (IN ANGSTROM) TAKEN FROM
!!!****  A. BONDI, J.PHYS.CHEM. 68, 441 (1964), TABLE I and table XIV
    DATA RBONDI(1:LMZ) / &
         1.20_SP,   1.40_SP,   1.82_SP,   1.45_SP,  -1.00_SP, &
         1.70_SP,   1.55_SP,   1.52_SP,   1.47_SP,   1.54_SP, &
         2.27_SP,   1.73_SP,  -1.00_SP,   2.10_SP,   1.80_SP, &
         1.80_SP,   1.75_SP,   1.88_SP,   2.75_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,   1.63_SP,   1.40_SP,   1.39_SP, &
         1.87_SP,   2.19_SP,   1.85_SP,   1.90_SP,   1.85_SP, &
         2.02_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         1.63_SP,   1.72_SP,   1.58_SP,   1.93_SP,   2.27_SP, &
         -1.00_SP,   2.06_SP,   1.98_SP,   2.16_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP /

!!!****  VAN-DER-WAALS RADII (IN ANGSTROM) TAKEN FROM
!!!****  J. EMSLEY, THE ELEMENTS, OXFORD UNIVERSITY PRESS, OXFORD (1997).
!!!****  http://www.shef.ac.uk/~chem/web-elements/
    DATA REMS(1:LMZ) / &
         1.20_SP,   1.22_SP,  -1.00_SP,  -1.00_SP,   2.08_SP, & 
         1.85_SP,   1.54_SP,   1.40_SP,   1.35_SP,   1.60_SP, &
         2.31_SP,  -1.00_SP,   2.05_SP,   2.00_SP,   1.90_SP, &
         1.85_SP,   1.81_SP,   1.91_SP,   2.31_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,   2.00_SP,   2.00_SP,   1.95_SP, &
         1.98_SP,   2.44_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,   2.00_SP, &
         2.20_SP,   2.20_SP,   2.15_SP,   2.16_SP,   2.62_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP,  -1.00_SP,   2.40_SP,  -1.00_SP,  -1.00_SP, &
         -1.00_SP /

    NQbasis = NQ
    IF(ASSOCIATED(AtomCharges)) DEALLOCATE(ATOMCHARGES)
    IF(ASSOCIATED(rho))         DEALLOCATE(rho)
    ALLOCATE(AtomCharges(NUMAT))
    ALLOCATE(rho(NQbasis))

    ! Dont overwrite everything if we are getting our surface/cosmo/radii
    ! information from elsewhere.
    IF(EXTSURF)RETURN

    ! sets up main radii and cosmo values
    IF (ASSOCIATED(radii_main)) DEALLOCATE(radii_main)
    IF (ASSOCIATED(vRadii))     DEALLOCATE(vRadii)
    IF (ASSOCIATED(ASAP))       DEALLOCATE(ASAP)

    ALLOCATE(radii_main(NUMAT))

    if(associated(cosmo_mndo))THEN
       call deallocate_cosmo(cosmo_mndo)
    ELSE
       allocate(cosmo_mndo)
    endif
    if(associated(esurf_mndo))THEN
       call deallocate_Surface(esurf_mndo)
    ELSE
       allocate(esurf_mndo)
    endif

    CALL ReadWriteSCOSMO(iRadii, ifGaus)
    SELECT CASE (iRadii)

    CASE (-1) ! Emsley values
       WRITE(6,*)'Using Emsley Radii for SCOSMO, natom= ',NUMAT
       forall(iatm=1:NUMAT)
          radii_main(iatm) = REMS(NAT(iatm))
       endforall
    CASE (0) ! Bondi radii (also default)
       WRITE(6,*)'Using Bondi Radii for SCOSMO, natom=',NUMAT
       forall(iatm=1:NUMAT)
          radii_main(iatm) = RBONDI(NAT(iatm))
       end forall
    CASE (1) ! Custom STATIC radii
       WRITE(6,*)'Using Custom Radii for SCOSMO, natom=',NUMAT
       DO iatm=1,NUMAT
          READ(5,*)   I, RI
          IF (I == iatm) THEN
             radii_main(iatm) = RI
          ELSE
             CALL ERROR('wrong input radii')
          END IF
       END DO
    CASE(2:)  ! Variable radii based on atomic charge
       WRITE(6,*)'Variable Radii Requested: Model ',iRadii
       radii_main = 0.0_SP
       allocate(vRadii(3,NUMAT))
       ! *** TJG 07/20/04
       DO I=1,NUMAT
          vRadii(1,I) = SCRQ0(NAT(I))
          vRadii(2,I) = SCRQ1(NAT(I))
          vRadii(3,I) = SCRQ2(NAT(I))
       END DO
       IF(ALL(vRadii(2:3,:) == 0.0_SP))THEN
          WRITE(6,*)'### Radii have no charge dependence'
          WRITE(6,*)'### Switching to constant radii scheme'
          iradii=1
          radii_main(:) = vRadii(1,:)
          deallocate(vRadii)
       ENDIF
       ! ---
    CASE DEFAULT ! default Bondi
       forall(iatm=1:NUMAT)
          radii_main(iatm) = RBONDI(NAT(iatm))
       end forall
    END SELECT

    IF(esurf_mndo%iCDR == 2)THEN
       ALLOCATE(ASAP(2,NUMAT))
       DO I=1,NUMAT
          ASAP(1,i) = SCNET(NAT(I))
          ASAP(2,i) = SCNEB(NAT(I))
       END DO
    END IF

    IF (cosmo_mndo%verbose) THEN
       WRITE(6,*)
       IF(iradii>=2 .AND. esurf_mndo%iCDR < 2) THEN
          WRITE(6,600)
          DO i = 1, NUMAT
             WRITE(6,610)  i, ELEMNT(NAT(i)), vRadii(1,i),vRadii(2,i),vRadii(3,i)
          END DO
       ELSE IF(iradii>=2 .AND. esurf_mndo%iCDR == 2)THEN
          WRITE(6,620)
          DO i = 1, NUMAT
             WRITE(6,630) i, ELEMNT(NAT(i)), vRadii(1,i),vRadii(2,i),vRadii(3,i),SCNET(NAT(i)),SCNEB(NAT(i))
          END DO
       ELSE
          WRITE(6,700)
          DO i = 1, NUMAT
             WRITE(6,710)  i, ELEMNT(NAT(i)), radii_main(i)
          END DO
       ENDIF
       WRITE(6,*)
    END IF
    ! convert to AU
    IF(iradii == 2) THEN !direct radii, ATM 1/31/05
       WRITE(6,*) 'Direct Variable Radii Model'
       vRadii(:,:) = ANGSTROM*vRadii(:,:)
    ELSEIF (iradii == 3) THEN !inverse radii, ATM 12/6/04
       WRITE(6,*) 'Inverse Variable Radii Model'
       vRadii(1,:) = ANGSTROM*vRadii(1,:)
       vRadii(2:3,:) = vRadii (2:3,:) / ANGSTROM
    ELSEIF (iradii == 4) THEN !exponential radii, ATM 1/31/05
       WRITE(6,*) 'Exponential Variable Radii Model'
       vRadii(1,:) = ANGSTROM*vRadii(1,:)
    ELSE
       radii_main(:) = ANGSTROM*radii_main(:)
    ENDIF

600 FORMAT( 4X, 'ATOM', 6X,'ELEMENT',4X,'RADIUS R0(A)',2X,'RADIUS R1(A)',2X,'RADIUS R2(A)')
610 FORMAT( 1X,I6,10X, A, 2X,3F13.5)
620 FORMAT( 4X, 'ATOM', 6X,'ELEMENT',5X,'RADIUS R0(A)',3X,'RADIUS R1(A)',3X,'RADIUS R2(A)',2X,'Surface Tension',2X,'Buried Atom')
630 FORMAT( 1X,I6,10X, A, 2X,5F15.5)
700 FORMAT( 4X, 'ATOM', 6X,'ELEMENT',4X,'RADIUS (A)')
710 FORMAT( 1X,I6,10X, A, 2X,F13.5)

  END SUBROUTINE InitSCOSMO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ReadWriteSCOSMO(iradii, lGauss)

!!!**** name list the same as in stand-alone COSMO with additional ones:
!!!      iRadii = type of Radii: -1 (Emsley), 0 (Bondi), 1 (input), >= 2 (varying)
!!!      lGauss = use gaussians for surface charge when computing B: 
!!!      F (point charge), T(Gaussian)
!!!      

    USE FileMod
    USE cosmo_io
    IMPLICIT NONE

    INTEGER(I4B),INTENT(INOUT):: iradii
    LOGICAL(LGD),INTENT(INOUT):: lGauss
    INTEGER(I4B) :: iunit,ierr

!!!**** Initialization
    iunit = AvailableUnitNumber()
    call OPENFILE(UNIT=iunit,FILE="scosmo.in",IOSTAT=ierr,STAT='old')
    IF (ierr /= 0) STOP 'Smooth cosmo input file not found. Looking for "scosmo.in"'
    CALL read_cosmo_parms(cosmo_mndo, RadType=iRadii, &
         GaussType=lGauss, unit=iunit)
    REWIND(iunit)
    CALL read_surface_parms(esurf_mndo,unit=iunit)
    CLOSE(iunit) 
    CALL print_cosmo_parms(cosmo_mndo, RadType=iRadii, &
         GaussType=lGauss)
    CALL print_surface_parms(esurf_mndo)
  END SUBROUTINE ReadWriteSCOSMO


END MODULE SCosmoMNDOMod

 
