program cosmo_test

  use DataTypes
  use PDB_io
  use ConstantsMod !only : ANGSTROM, KCAL_PER_MOL
  use MolecularPropertiesMod
  use TimerMod, only: pptime
  USE filemod
  use ioMod, only: write_array

  use surfacemod
  use cosmomod
  use cosmo_helper
  use cosmo_util
  use cosmo_io
  USE MinimizeMod
  USE RBM
  implicit none

  ! solute
  type(pdb_dataT) :: molecule

  ! cosmo data
  type(cosmo_t) :: cosmo_data
  ! Surface Construction
  type(Surface_t) :: e_surf    ! electrostatic surface
!  type(Surface_t) :: cdr_surf  ! non-electrostatic surface

  REAL(SP),pointer :: phi(:)=>NULL(),rho(:)=>NULL(),radii(:)=>NULL(),Bmat(:,:)=>NULL(),ZMAT(:,:)=>NULL()

  real(sp) :: be, surface_dipole(1:3), solute_dipole(1:3), sur_dmag, sol_dmag
  real(sp) :: solute_qm(1:6), surface_qm(1:6)
  integer(i4b) :: i
  character(len=256) :: filename='',energy_file=''


!!!**** scratch
  INTEGER(I4B) :: iatm, unit, ierr, nrho, nse, ik,ft
  REAL(SP) :: tmp
  REAL :: time1, time2
  real :: start_time, stop_time
  REAL(SP), pointer :: ChargeCoord1(:,:)=>NULL(), ChargeCoord2(:,:)=>NULL(), &
       Pot1(:)=>NULL(), Pot2(:)=>NULL(),Field1(:,:)=>NULL(), Field2(:,:)=>NULL()
  REAL(SP), allocatable :: CDRGrad(:,:),EPolGrad(:,:)
  LOGICAL(LGD) :: Success, PQR

  WRITE(6,*)'File format: 1) PQR    2) PRQ'
  READ(5,*)ft

  IF(ft==1)THEN
    PQR=.TRUE.
  ELSE
    PQR=.FALSE.
  ENDIF

  CALL read_cosmo_parms(cosmo_data,infile=filename,energyf=energy_file)
  CALL read_surface_parms(e_surf)
  CALL print_cosmo_parms(cosmo_data,infile=filename,energyf=energy_file)
  CALL print_surface_parms(e_surf)

  CALL readpdb(filename, PDB=molecule, ALLOC=.true., UNITS='AU')

  call cpu_time(start_time)

  IF(pqr)THEN
     radii => molecule%Bvalue
     rho => molecule%occupation
  ELSE
     radii => molecule%occupation
     rho => molecule%Bvalue 
  ENDIF
!!! **** convert radii from angstrom to atomic units
   radii = radii * ANGSTROM

   CALL SASurface(molecule%crd, radii, e_surf)
   DO i=1, min(5,SIZE(molecule%Crd,2))
      WRITE(6,'(I4,A,3ES17.7)')i,'  crd= ', molecule%Crd(1:3,i)
      WRITE(6,'(4X,A,2ES17.7)')'radii, charge= ',radii(i),rho(i)
   ENDDO

   nrho = SIZE(rho)
   nse = SIZE(e_surf%secrd,2)

!!!**** compute Epol
     CALL CPU_TIME(time1)
     SELECT CASE(cosmo_data%minpara%Multipole)
     CASE (0)
        IF (cosmo_data%verbose) WRITE(6,*)'phi_rho_sigma by matrix multiplication...'
        CALL cpt_phirhosigma_for_pt(molecule%crd,rho,e_surf%secrd, e_surf%cosmozeta,phi)
     CASE (1)
        IF (cosmo_data%verbose) WRITE(6,*)'phi_rho_sigma with recursive bisection...'
        allocate(phi(nse))
        ALLOCATE(ChargeCoord1(1:4,1:nrho), ChargeCoord2(1:4,1:nse))
        ALLOCATE(Pot1(nrho), Pot2(nse))

        ChargeCoord1(1,:) = rho(:)
        ChargeCoord1(2:4,:) = molecule%Crd(1:3,:)
        ChargeCoord2(1,1:nse) = 0.0_sp
        ChargeCoord2(2:4,:) = e_surf%SECRD(1:3,:)
        Pot1 = 0.0_sp;  Pot2 = 0.0_sp;
        CALL RecursiveBisectionMethod(ChargeCoord1,Pot1, ChargeCoord2,Pot2,&
             WS=cosmo_data%minpara%WS,L=cosmo_data%minpara%L)
        phi = Pot2

        DEALLOCATE(ChargeCoord1, ChargeCoord2, Pot1, Pot2)
     CASE DEFAULT
        CALL error('illegal Multipole')
     END SELECT

   WRITE(6,*)'My lval', cosmo_data%l
   IF(cosmo_data%l<0)THEN
     call do_cosmo(cosmo_data, e_surf, SUM(rho),phi)
   ELSE
     call cpt_Z_for_pt(molecule%crd, cosmo_data%l, Zmat)        
     call do_cosmo(cosmo_data, e_surf, SUM(rho),phi,matmul(Zmat,rho))
   ENDIF

  CALL CPU_TIME(time2)
  WRITE(6,*) 'Etot CPU time=', time2-time1, 'sec'

  !**** compute Ecav (Eq 11 in York_ChemPhysLett_1996_v263_p297
  IF (e_surf%iCDR == 1 ) THEN
     ! Call SASurface(molecule%crd, radii,cdr_surf)
     CALL MolecularNonelecEnergy(e_surf, cosmo_data%Ecdr)
     cosmo_data%Etotal = cosmo_data%Epol + cosmo_data%Ecdr
  ELSE IF (e_surf%iCDR >= 2 ) THEN
     WRITE(6,*) 'Atomic Nonelectrostatics not implemented yet, no nonelec energy for you'
     cosmo_data%Etotal = cosmo_data%Epol 
  ELSE
     cosmo_data%Etotal = cosmo_data%Epol 
  END IF

  WRITE(6,*)'Epol=',cosmo_data%Epol / KCAL_PER_MOL,' kcal/mol'
  WRITE(6,*)'Ecdr=',cosmo_data%Ecdr / KCAL_PER_MOL,' kcal/mol'
  WRITE(6,*)'Esolvation=',cosmo_data%Etotal / KCAL_PER_MOL,' kcal/mol'

  call cpt_born(cosmo_data%e1,cosmo_data%e2,e_surf%sa,rho, be)
  write(6,*)
  write(6,*)'Born energy=   ',be / KCAL_PER_MOL,' kcal/mol'

  call cpt_onsager(cosmo_data%e1, cosmo_data%e2, e_surf%sa, rho, molecule%crd, be)
  write(6,*)'Onsager energy=',be / KCAL_PER_MOL,' kcal/mol'

  call cpt_GBorn(cosmo_data%e1, cosmo_data%e2,rho,molecule%crd,radii,be)
  write(6,*)'GenBorn energy=',be / KCAL_PER_MOL,' kcal/mol'

  call cpt_born2(cosmo_data%e1, cosmo_data%e2,rho,molecule%crd,radii,be)
  write(6,*)'Born2 energy=   ',be / KCAL_PER_MOL, ' kcal/mol'




! Gradients
    unit = -1
    IF(LEN_TRIM(energy_file) /= 0)THEN
      CALL OpenFile(FILE=energy_file, UNIT=unit, IOSTAT=ierr)
    ELSE
      unit=6
    ENDIF
   write(UNIT,'(5a15)') 'Epol', 'Ecdr', 'Etotal' , 'kcal/mol'
   write(UNIT,'(5es15.5)') cosmo_data%Epol/KCAL_PER_MOL, &
        cosmo_data%Ecdr/KCAL_PER_MOL, cosmo_data%Etotal/KCAL_PER_MOL


 !!!**** coverting factor tmp: au to kcal/mol/Ang

   tmp = 1/ (KCAL_PER_MOL * AU_DISTANCE_IN_ANGSTROM)

 !!!**** analytical gradients or both

   IF (cosmo_data%igrad == 1 .OR. cosmo_data%igrad == 3) THEN
      WRITE(6,*)'Computing analytical gradients...'

      CALL CPU_TIME(time1)
      IF (ASSOCIATED(ChargeCoord1)) DEALLOCATE(ChargeCoord1)
      IF (ASSOCIATED(ChargeCoord2)) DEALLOCATE(ChargeCoord2)
      IF (ASSOCIATED(Pot1)) DEALLOCATE(Pot1)
      IF (ASSOCIATED(Pot2)) DEALLOCATE(Pot2)

      ALLOCATE(ChargeCoord1(1:4,1:nrho), ChargeCoord2(1:4,1:nse))
      ALLOCATE(Pot1(nrho), Pot2(1:nse))
      ALLOCATE(Field1(3,nrho), Field2(3,nse))
      ChargeCoord1(1,:) = rho(:)
      ChargeCoord1(2:4,:) = molecule%Crd(1:3,:)
      ChargeCoord2(1,1:nse) = cosmo_data%surface_charge(:)
      ChargeCoord2(2:4,1:nse) = e_surf%SECRD(1:3,1:nse)
      Pot1 = 0.0_sp; Field1 = 0.0_SP
      Pot2 = 0.0_sp; Field2 = 0.0_SP

    SELECT CASE(cosmo_data%minpara%multipole)
    CASE(0) 
       CALL DirectCoulomb(ChargeCoord1, Pot1,Field1, ChargeCoord2, Pot2, Field2,.FALSE.)
    CASE(1) 
       CALL RecursiveBisectionMethod(ChargeCoord1, Pot1,Field1, ChargeCoord2, Pot2, Field2, &
            WS=cosmo_data%minpara%WS,L=cosmo_data%minpara%L)
    CASE DEFAULT
       CALL error('illegal Ifmultipole in CptAnalytGradient')
    END SELECT



 !!!**** obtain energy gradients (negative of forces) from potential gradients
      DO ik = 1,3
         Field1(ik,:) = Field1(ik,:) * rho(:)
         Field2(ik,:) = Field2(ik,:) * cosmo_data%surface_charge(:)
      END DO

      WRITE(UNIT,*)
      WRITE(UNIT,*)'GradEpol (negative forces) due to electric field of sigma in kcal/mol/ang'
      WRITE(UNIT,*)'GradEpol(sigma)'
      DO iatm = 1, nrho
         WRITE(UNIT,'(3es15.5)') (tmp*Field1(i,iatm), i=1,3)
      END DO
      WRITE(UNIT,*)'***********************************************************'
      WRITE(UNIT,*)
      Allocate(EPolGrad(3,nrho))
      EPolGrad = Field1

      CALL AnalytGradient(e_surf, cosmo_data%surface_charge, cosmo_data%f, &
           cosmo_data%minpara, EPolGrad, Field2)

      CALL CPU_TIME(time2)
      WRITE(6,*) 'gradient CPU time=', time2-time1, 'sec'

      WRITE(UNIT,*)'Total Analytic energy gradients (negative forces) from Epol in kcal/mol/ang'
      WRITE(UNIT,*)'GradEpol'
      DO iatm = 1, nrho
         WRITE(UNIT,'(3es15.5)') (tmp*EPolGrad(i,iatm), i=1,3)
      END DO

      allocate(CDRGrad(3,nrho))
      CDRGrad = 0.0_SP
      IF (e_surf%iCDR /= 0 .AND. e_surf%GammaSwitch > 0.0_SP) THEN
         CALL AnalytMolNonelecGradient(e_surf,CDRGrad)
      ENDIF
      WRITE(UNIT,*)'GradEcavdisrep'
      DO iatm = 1, nrho
         WRITE(UNIT,'(3es15.5)') (tmp*CDRGrad(i,iatm), i=1,3)
      END DO

      write(UNIT,*)'GradTotal'
      DO iatm = 1, nrho
         EPolGrad(:,iatm) = EPolGrad(:,iatm) + CDRGrad(:,iatm)
         WRITE(UNIT,'(3es15.5)') (tmp*EPolGrad(i,iatm), i=1,3)
      END DO
      deallocate(EPolGrad,CDRGrad)
   END IF


 !!!**** numerical gradients or both
   IF (cosmo_data%igrad == 2 .OR. cosmo_data%igrad == 3) THEN

      WRITE(6,*)'Computing numerical gradients'
      CALL CPU_TIME(time1)
      CALL NumGradientFast(cosmo_data, e_surf, molecule%crd, molecule%occupation,rho,Success)

      IF (.NOT. Success) THEN
         WRITE(6,*) 'Compute slow gradients'
         CALL NumGradientSlow(cosmo_data, e_surf, molecule%crd,molecule%occupation,rho)
      END IF

      CALL CPU_TIME(time2)
      WRITE(6,*) 'gradient CPU time=', time2-time1, 'sec'

      WRITE(UNIT,*)'Numerical Gradients (negative forces) in kcal/mol/ang'
      WRITE(UNIT,*)'GradEpol'
      DO iatm = 1, nrho
         WRITE(UNIT,'(3es15.5)') (tmp*cosmo_data%Egrad(i,iatm), i=1,3)
      END DO

      WRITE(UNIT,*)'GradEcavdisrep'
      DO iatm = 1, nrho
         WRITE(UNIT,'(3es15.5)') (tmp*cosmo_data%CDRGrad(i,iatm), i=1,3)
      END DO

      WRITE(UNIT,*)'GradTotal'
      DO iatm = 1, nrho
         cosmo_data%Egrad(:,iatm) = cosmo_data%Egrad(:,iatm) + cosmo_data%CDRGrad(:,iatm)
         WRITE(UNIT,'(3es15.5)') (tmp*cosmo_data%EGrad(i,iatm), i=1,3)
      END DO
   END IF

end program cosmo_test
