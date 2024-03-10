module cosmomod
!!!****h* cosmo
!!!
!!! NAME
!!!     cosmo -- module for smooth conductor-like screening model
!!! 
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Kevin Range
!!!     Darrin York
!!! CREATION DATE
!!!     2002
!!! DESCRIPTION
!!!     Implements the smooth conductor-like screening model as described in
!!!     York_JPhysChemA_1999_v103_p11060
!!! USES
!!!     SurfaceMod
!!!     YorkLib: DataTypes, etc.
!!!
!!!***
  use DataTypes
  implicit none
  private

  public :: do_cosmo,deallocate_cosmo, write_cosmo, &
       read_cosmo, DeallocateCosmoArrays, copy_cosmo_scalar_inputs

  TYPE, PUBLIC :: minpara_t
     REAL(SP) :: tol=1.0E-7, WS=2.0_SP, FDstep=1.0e-6_SP
     INTEGER(I4B) :: itmax=100, itol=1, precond=1, multipole=0, L=8
     LOGICAL :: IfAmatrix=.FALSE.
  END TYPE minpara_t

  INTERFACE do_cosmo
     Module Procedure do_cosmo_b_rho, do_cosmo_phi
  END INTERFACE

!!!****s* cosmo/cosmo_t
!!!
!!! NAME
!!!     cosmo_t -- cosmo data type
!!! SOURCE
!!!
  type, public :: cosmo_t
     ! inputs
     type(minpara_t) :: minpara
     ! dielectrics
     ! e1 is in the cavity (default 1), e2 is outside the cavity (default 80)
     ! setting e2 to <= 0 causes scaling factor f to be 1 (perfect conductor)
     real(sp) :: e1 = 1.0_sp
     real(sp) :: e2 = 80.0_sp
     ! Type of COSMO calculation
     INTEGER(I4B) :: icosmo = 0
     ! Compute Energy Gradients
     INTEGER(I4B) :: igrad = 0
     ! constraints
     ! l is the highest multipole to constrain.
     ! If l <= -1 there are no constraints
     integer(i4b) :: l = -1
     logical :: variational = .true.
     ! if low_mem is true, the code will throw away
     ! A matrix and other "unnecessary" data to save memory
     logical :: low_memory = .false.
     ! if print is true we will print out some info
     logical :: verbose = .false.
     !
     ! outputs
     !
     ! polarization energy, surface charge self energy, 
     ! surface charge-solute energy, solute-solute energy, total electrostatic energy
     real(sp) :: Epol=0.0_SP, Eself=0.0_SP, Ess=0.0_SP,  &
              &  Ecdr=0.0_SP, Etotal=0.0_SP
     ! surface charge vector
     real(sp), pointer :: surface_charge(:) => null()
     ! electrostatic potential due to the surface charge 
     ! at the solute coordinates
     real(sp), pointer :: phi_sigma_rho(:) => null()
     ! polarization Green's function matrix
     real(sp), pointer :: Gpol(:,:) => null()
     ! surface self-interaction matrix
     ! A, said, and, Ainv are all UNSCALED (i.e., A_0 in the paper)
     real(sp), pointer :: A(:,:) => null()
     ! inverse of A
     real(sp), pointer :: Ainv(:,:) => null()     
     ! Gradient of Epol
     REAL(SP), POINTER :: EGrad(:,:) => null()
     ! Gradient of Ecdr, i.e. nonelectrostatics
     REAL(SP), POINTER :: CDRGrad(:,:) => null()
     ! Surface constraint matrix (D)
     ! This number of surface elements by number of contraints matrix
     ! defines the constraints on the surface charge.  The transpose of this
     ! matrix operating on surface_charge should give the relevant multipoles
     real(sp), pointer :: D(:,:) => null()
     ! Auxillary matricies for constrained COSMO
     ! Q = D^T A^{-1} D
     real(sp), pointer :: Qinv(:,:) => null()
     ! R = Z + D^T A^{-1} B
     real(sp), pointer :: R(:,:) => null()
     ! lambda = (Z \rho + D^T Ainv \phi)
     real(sp), pointer :: lambda(:) => null()
     ! dielectric scale function
     real(sp) :: f=1.0_SP

  end type cosmo_t
 
!!!***
contains


  ! B:  surface element - density interaction matrix
  ! rho:  electron density (charge) vector of solute
  ! phi:  electrostaic potential due to the solute at the surface elements
  ! Z: Solute constraint matrix
  !    This number of constraints by number of solute interaction sites matrix
  !    defines the constraint values for constrained COSMO.
  !    This matrix operating on rho should give the constraint 
  !    values (minus multipoles)

  subroutine do_cosmo_phi(cosmo_data, surface, QTot,phi, zrho)
    ! Main cosmo routine for phi input
    use SurfaceMod
    use UtilitiesMod
    use ErrorMod
    use cosmo_util

    type(cosmo_t) :: cosmo_data
    type(surface_t) :: surface
    REAL(SP), INTENT(IN) :: QTot
    REAL(SP), INTENT(IN) :: phi(:)
    REAL(SP), OPTIONAL, INTENT(IN) :: zrho(:)

    ! compute dielectric scale factor
    cosmo_data%f = cosmo_f(cosmo_data%e1, cosmo_data%e2)

    if(cosmo_data%iCosmo == 1) THEN ! direct energy minimization
       IF(cosmo_data%minpara%IfAmatrix .AND. .NOT. associated(cosmo_data%A))THEN
          call cpt_A(surface%secrd,surface%cosmozeta,cosmo_data%A)
       ENDIF
       ! After this call, we will have Epol, Eself and Ess
       CALL cpt_surface_charge_minimize(cosmo_data, surface,QTot,phi)
    ELSE ! A inverse
       IF(.NOT. associated(cosmo_data%Ainv))THEN
          IF(.NOT. associated(cosmo_data%A))THEN
             call cpt_A(surface%secrd,surface%cosmozeta,cosmo_data%A)
          ENDIF
          IF(associated(surface%sescale))THEN
             call cpt_Ainv(surface%gammaswitch,surface%sescale, &
                  cosmo_data%A,cosmo_data%Ainv,cosmo_data%low_memory,cosmo_data%verbose)
          ELSE
             call cpt_Ainv(cosmo_data%A,cosmo_data%Ainv,cosmo_data%low_memory,cosmo_data%verbose)
          ENDIF
       ENDIF
       if ( cosmo_data%l >= 0 .AND. PRESENT(Zrho) .AND. .NOT. associated(cosmo_data%lambda)) then
          write(6,*)'Setting up for constrained COSMO'
          call cpt_constraints(cosmo_data,surface,phi,zrho)
       ENDIF
       ! After this call, we will have Epol, Eself and Ess
       CALL cpt_surface_charge_inverse(cosmo_data,surface,phi)
    END IF

  end subroutine do_cosmo_phi

  subroutine do_cosmo_b_rho(cosmo_data, surface, B,rho,Z,CPT_GPOL, CPT_PHI, CPT_SIGMA)
    ! Main cosmo routine
    use SurfaceMod
    use UtilitiesMod
    use ErrorMod
    use cosmo_util

    type(cosmo_t) :: cosmo_data
    type(surface_t) :: surface
    REAL(SP), INTENT(IN) :: rho(:)
    REAL(SP), INTENT(IN) :: B(:,:)
    REAL(SP), OPTIONAL, INTENT(IN) :: Z(:,:)
    REAL(SP), allocatable :: phi(:)
    ! Do you want Gpol, phi_sigma_rho, sigma
    ! I will compute the surface charge and the E's by default
    logical(lgd), optional :: CPT_GPOL, CPT_PHI, CPT_SIGMA
    logical(lgd) :: do_sigma,do_gpol,do_sigma_inv,do_phi, needAinv,needA

    ! compute phi_rho_sigma if we don't have it already
    allocate(phi(size(B,1)))
    phi = matmul(B,rho)
    ! compute dielectric scale factor
    cosmo_data%f = cosmo_f(cosmo_data%e1, cosmo_data%e2)

    do_phi = optional_flag(CPT_PHI,default=.true.)
    do_sigma = optional_flag(CPT_SIGMA,default=.true.) .OR. do_phi
    do_gpol = optional_flag(CPT_GPOL,default=.false.)
    do_sigma_inv = (cosmo_data%iCosmo == 0) .AND. do_sigma

    needAinv = (cosmo_data%l >= 0 .AND. PRESENT(Z)) .OR. do_sigma_inv .OR. do_gpol
    needA = ((cosmo_data%iCosmo == 1 .AND. cosmo_data%minpara%IfAmatrix ) &
         .OR. (needAinv .and. .NOT. associated(cosmo_data%ainv)))

    if ( needA .AND. .not. associated(cosmo_data%a))THEN
       call cpt_A(surface%secrd,surface%cosmozeta,cosmo_data%A)
    end if

    if ( do_sigma .AND. cosmo_data%iCosmo == 1)THEN
       CALL cpt_surface_charge_minimize(cosmo_data, surface,SUM(rho),phi)
    ENDIF

    if(needAinv)THEN
       IF(associated(surface%sescale))THEN
          call cpt_Ainv(surface%gammaswitch,surface%sescale, &
               cosmo_data%A,cosmo_data%Ainv,cosmo_data%low_memory,&
               cosmo_data%verbose)
       ELSE
          call cpt_Ainv(cosmo_data%A,cosmo_data%Ainv,cosmo_data%low_memory,&
               cosmo_data%verbose)
       ENDIF
    endif

    ! Set up constraint stuff if we need it
    if ( cosmo_data%l >= 0 .AND. PRESENT(Z)) then
       write(6,*)'Setting up for constrained COSMO and/or Gpol'
       IF(.NOT. associated(cosmo_data%lambda))THEN
          call cpt_constraints(cosmo_data,surface,phi,matmul(Z,rho))
       ENDIF
       if (do_gpol )call cpt_R(cosmo_data,B,Z)
    end if

    if ( do_sigma .and. cosmo_data%iCosmo == 0) then
       CALL cpt_surface_charge_inverse(cosmo_data,surface,phi)
    END IF

    if ( do_gpol ) call do_cpt_Gpol(cosmo_data,B,rho)

    if ( do_phi )  call do_cpt_phi(cosmo_data,B,Z)

  deallocate(phi)

  end subroutine do_cosmo_b_rho

 subroutine do_cpt_phi(cosmo_data,B,Z)
    ! Compute phi_sigma_rho, the potential due to the surface charges at the 
    ! solute interaction points.
    ! Requires the B matrix (or else this code would have to be WAY smarter)
    Use SurfaceMod
    use errormod
    
    type(cosmo_t), intent(inout) :: cosmo_data
    REAL(SP), INTENT(IN) :: B(:,:)
    REAL(SP), OPTIONAL,INTENT(IN) :: Z(:,:)
    integer(i4b) :: nq, stat_alloc

    IF (cosmo_data%verbose)  write(6,*)'Computing phi...'

    if ( .not. associated(cosmo_data%phi_sigma_rho) ) then
       nq = size(B,2)
       if ( nq == 0 ) then
          write(6,*)'ERROR: I cannot compute phi_sigma_rho without B'
          return
       else
          allocate(cosmo_data%phi_sigma_rho(1:nq),stat=stat_alloc)
          if ( stat_alloc > 0 ) then
             write(6,*)&
                  'ERROR allocating.  phi_sigma_rho will NOT be calculated.'
             return
          end if
       end if
    end if

    cosmo_data%phi_sigma_rho = matmul(cosmo_data%surface_charge,B)

    if (cosmo_data%variational .and. associated(cosmo_data%lambda) .AND. present(Z)) then
       cosmo_data%phi_sigma_rho = cosmo_data%phi_sigma_rho + &
            matmul(transpose(Z),cosmo_data%lambda)
    end if

  end subroutine do_cpt_phi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine do_cpt_Gpol(cosmo_data,B,rho)
    ! compute Gpol matrix
    ! Ainv and possible constraints should be calculated
    Use SurfaceMod
    use errormod

    type(cosmo_t), intent(inout) :: cosmo_data
    REAL(SP),INTENT(IN) :: B(:,:),rho(:)

    integer(i4b) :: nq, stat_alloc

    IF(cosmo_data%verbose)write(6,*)'Computing Gpol...'

    if ( .not. associated(cosmo_data%Gpol) ) then
       nq = size(B,2)
       allocate(cosmo_data%Gpol(1:nq,1:nq),stat=stat_alloc)
       if ( stat_alloc > 0 ) then
          write(6,*)&
               'ERROR allocating Gpol.  Gpol will NOT be calculated.'
          return
       end if
    end if
    
    IF(cosmo_data%l < 0)THEN
       ! Eqn. 43 in York_JPhysChemA_1999_v103_p11060
       cosmo_data%Gpol = - cosmo_data%f * matmul(transpose(B), &
            matmul(cosmo_data%Ainv, B))
    ELSE
       if ( cosmo_data%variational ) then
          ! Eqn B-6
         cosmo_data%Gpol = - cosmo_data%f * &
               matmul(transpose(B), matmul(cosmo_data%Ainv, B)) + &
               matmul(transpose(cosmo_data%R),&
               matmul(cosmo_data%Qinv,cosmo_data%R))
       else
          ! Eqn B-8
          cosmo_data%Gpol = - cosmo_data%f * &
               matmul(transpose(B), matmul(cosmo_data%Ainv, &
               B - matmul(cosmo_data%D,matmul(cosmo_data%Qinv,cosmo_data%R))))
       end if
    end if

    cosmo_data%Ess = 0.0_sp
    cosmo_data%Eself = 0.0_sp
    cosmo_data%Epol = 0.5_sp * dot_product(rho,matmul(cosmo_data%Gpol, rho))

  end subroutine do_cpt_Gpol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE cpt_surface_charge_minimize(cosmo_data,surface,Qtot,phi)
!!!****f* CosmoMod/cpt_surface_charge_minimize
!!! 
!!! DESCRIPTION
!!!    a driver subroutine to call conjugate minimization 
!!!***

    USE ErrorMod
    USE ConstantsMod!, ONLY: KCAL_PER_MOL, 
    USE SurfaceMod
    use cosmo_util

    TYPE(cosmo_t), INTENT(INOUT) :: cosmo_data
    TYPE(surface_t),INTENT(IN) :: surface
    REAL(SP), INTENT(IN) :: phi(:),Qtot

    REAL(SP), PARAMETER :: SQRT_TWO_over_PI = SQRT2/SQRT_PI

    INTEGER(I4B) :: iter, nse, ii
    REAL(SP) :: error_minimize
    REAL :: time1, time2

    IF (cosmo_data%verbose) WRITE (6,*) 'Start minimizing energy'
    nse = SIZE(Surface%secrd,2)

!!!**** compute the dielectric scale factor f
    cosmo_data%f = cosmo_f(cosmo_data%e1, cosmo_data%e2)

    IF ( nse == 0 ) CALL error('MinimizeEnergy: no surface.')
    IF ( .NOT. ASSOCIATED(cosmo_data%surface_charge) ) &
         ALLOCATE(cosmo_data%surface_charge(1:nse))
    
    ! This initial guess is special.  It obeys Gauss' Law in the dielectric.
    cosmo_data%surface_charge = -Qtot * cosmo_data%f / nse

!!!**** compute A if asked for
    IF (cosmo_data%minpara%IfAmatrix) THEN
       IF (cosmo_data%verbose) WRITE(6,*) 'Use A matrix'
       IF(.not. associated(cosmo_data%A))THEN
         WRITE(6,*)'Computing A'
         CALL cpt_A(surface%secrd,surface%cosmozeta,cosmo_data%A)
       ENDIF 

       IF(org_Sik_in_A0_diag)THEN
          DO ii=1,nse
             cosmo_data%A(ii,ii) = cosmo_data%A(ii,ii)/Surface%sescale(ii)
          END DO
       ELSE
          DO ii=1,nse
             cosmo_data%A(ii,ii) = cosmo_data%A(ii,ii)/SQRT(Surface%sescale(ii))
          END DO  
       ENDIF

       CALL CPU_TIME(time1)
       CALL ConjugateGradient(surface, cosmo_data%f, &
            cosmo_data%minpara, phi, cosmo_data%surface_charge, &
            iter,error_minimize,cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess,cosmo_data%A)
       CALL CPU_TIME(time2)

       IF (cosmo_data%verbose) WRITE(6,*) 'CPU time/iter =', (time2-time1)/iter, 'sec'

!!!**** fix A matrix: may be necessary if using A.
       IF(org_Sik_in_A0_diag)THEN
          DO ii=1,nse
             cosmo_data%A(ii,ii) = cosmo_data%A(ii,ii) * Surface%sescale(ii)
          END DO
       ELSE
          DO ii=1,nse
             cosmo_data%A(ii,ii) = cosmo_data%A(ii,ii) * SQRT(Surface%sescale(ii))
          END DO
       ENDIF

    ELSE
       CALL CPU_TIME(time1)
       CALL ConjugateGradient(surface, cosmo_data%f, &
            cosmo_data%minpara, phi, cosmo_data%surface_charge, &               
            iter,error_minimize,cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess)
       CALL CPU_TIME(time2)
       IF (cosmo_data%verbose)WRITE(6,*) 'CPU time/iter =', (time2-time1)/iter, 'sec'
    END IF

    IF (cosmo_data%verbose) THEN
       WRITE(6,*)
       WRITE(6,*) 'Epol minimization converged!'
       WRITE(6,*) 'Total number of iterations and error =', iter,  error_minimize
       WRITE(6,*) 'Epol =', cosmo_data%Epol/KCAL_PER_MOL, ' Kcal'
       WRITE(6,*) 'surface charge =',SUM(cosmo_data%surface_charge)
       WRITE(6,*)
    END IF
  END SUBROUTINE cpt_surface_charge_minimize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cpt_surface_charge_inverse(cosmo_data,surface,phi)

    ! compute surface charge
    use SurfaceMod
    use errormod

    type(cosmo_t), intent(inout) :: cosmo_data
    type(surface_t),intent(in) :: surface
    REAL(SP), INTENT(IN) :: phi(:)
    REAL(SP):: tmpE
    integer(i4b) :: nse

    IF (cosmo_data%verbose) THEN
       write(6,*)'Computing surface charge...'
    END IF

    if ( .not. associated(cosmo_data%surface_charge) ) then
       nse = size(surface%secrd,2)

       if ( nse == 0 ) then
          call error('setup_sigma cannot find the surface.')
       else
          allocate(cosmo_data%surface_charge(1:nse))
       end if
    end if

    if ( cosmo_data%l < 0 ) then
       ! No constraints
       ! Eqn. 41 in York_JPhysChemA_1999_v103_p11060
       
       cosmo_data%surface_charge = &
            - cosmo_data%f * matmul(cosmo_data%Ainv,phi)

       tmpE =  dot_product(cosmo_data%surface_charge,phi)
       cosmo_data%Eself = - 0.5_sp * tmpE
       cosmo_data%Ess = tmpE
    else
       ! constrained COSMO
       ! See Appendix B of York_JPhysChemA_1999_v103_p11060 and elsewhere
       ! We already did the constraints before we get here (We Have D and lambda)

       cosmo_data%surface_charge = &
            - cosmo_data%f * matmul(cosmo_data%Ainv,&
            (phi - matmul(cosmo_data%D, &
            cosmo_data%lambda)))

       cosmo_data%Eself = - 0.5_sp * &
            dot_product(cosmo_data%surface_charge, &
            (phi - matmul(cosmo_data%D,cosmo_data%lambda)))
       IF(cosmo_data%variational)THEN
          cosmo_data%Ess = dot_product(cosmo_data%surface_charge,phi)
       ELSE
          cosmo_data%Ess = dot_product(cosmo_data%surface_charge,phi) 
       ENDIF
    end if

    cosmo_data%Epol = cosmo_data%Eself + cosmo_data%Ess

  end subroutine cpt_surface_charge_inverse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine copy_cosmo_scalar_inputs(cdata,cdata_new)
    type(cosmo_t), intent(in) :: cdata
    type(cosmo_t), intent(out) :: cdata_new
 
     cdata_new%e1                  = cdata%e1
     cdata_new%e2                  = cdata%e2
     cdata_new%icosmo              = cdata%icosmo
     cdata_new%l                   = cdata%l
     cdata_new%variational         = cdata%variational
     cdata_new%low_memory          = cdata%low_memory
     cdata_new%verbose             = cdata%verbose

     CALL copy_minpara_scalar_inputs(cdata%minpara,cdata_new%minpara)
  end subroutine copy_cosmo_scalar_inputs

  SUBROUTINE  copy_minpara_scalar_inputs(minpara, TmpMinpara)
    TYPE(minpara_t):: minpara, TmpMinpara

       TmpMinpara%tol = Minpara%tol
       TmpMinpara%WS = Minpara%WS
       TmpMinpara%itmax = Minpara%itmax
       TmpMinpara%itol = Minpara%itol
       TmpMinpara%precond = Minpara%precond
       TmpMinpara%multipole = Minpara%multipole
       TmpMinpara%L = Minpara%L
       TmpMinpara%IfAmatrix = Minpara%IfAmatrix
       TmpMinpara%FDStep = Minpara%FDstep
  END SUBROUTINE copy_minpara_scalar_inputs


  subroutine deallocate_cosmo(cosmo_data)

    ! Deallocateds a cosmo data type and 
    ! sets all scalars back to default values
    type(cosmo_t), intent(inout) :: cosmo_data
    type(cosmo_t) :: my_cdata

    ! Make sure we deallocate surface arrays
    call deallocateCosmoArrays(cosmo_data)

    ! This initalizes cosmo_data to default values stored in the datatype
    ! In fortran, this is will copy scalar values and associate pointers
    ! however, the pointers in my_cdata should ALL be null and 
    ! the pointers in cosmo_data have just been nullified, thus there should be no
    ! memory leaks, or pointers to freed memory 
    cosmo_data=my_cdata

  end subroutine deallocate_cosmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DeallocateCosmoArrays(cosmo_data)

    TYPE(cosmo_t), INTENT(INOUT) :: cosmo_data
       if ( associated(cosmo_data%surface_charge) ) then
          deallocate(cosmo_data%surface_charge)
       end if
!       if ( associated(cosmo_data%surface_charge_cons) ) then
!          deallocate(cosmo_data%surface_charge_cons)
!       end if
       if ( associated(cosmo_data%phi_sigma_rho) ) then
          deallocate(cosmo_data%phi_sigma_rho)
       end if
       if ( associated(cosmo_data%Gpol) ) then
          deallocate(cosmo_data%Gpol)
       end if
       if ( associated(cosmo_data%A) ) then
          deallocate(cosmo_data%A)
       end if
       if ( associated(cosmo_data%Ainv) ) then
          deallocate(cosmo_data%Ainv)
       end if
        IF ( ASSOCIATED(cosmo_data%EGrad) ) THEN
          DEALLOCATE(cosmo_data%EGrad)
       END IF
       IF ( ASSOCIATED(cosmo_data%CDRGrad) ) THEN
          DEALLOCATE(cosmo_data%CDRGrad)
       END IF
      if ( associated(cosmo_data%D) ) then
          deallocate(cosmo_data%D)
       end if
       IF ( ASSOCIATED(cosmo_data%Qinv)) THEN
          DEALLOCATE(cosmo_data%Qinv)
       ENDIF
       IF ( ASSOCIATED(cosmo_data%R)) THEN
          DEALLOCATE(cosmo_data%R)
       ENDIF
       IF ( ASSOCIATED(cosmo_data%lambda)) THEN
          DEALLOCATE(cosmo_data%lambda)
       ENDIF

  END SUBROUTINE DeallocateCosmoArrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_constraints(cosmo_data,surface,phi,zrho)
  ! computes surface constraint matrix (D) and the Qinv auxillary matrix

    use UtilitiesMod
    USE MultipoleMod
    USE SurfaceMod
    use cosmo_util

    type(cosmo_t), intent(inout) :: cosmo_data
    type(surface_t),intent(in) :: surface
    REAL(SP), INTENT(IN) :: phi(:),zrho(:)

    integer(i4b) :: ldim, nse
    real(sp) :: ctr(3) = (/0.0_SP,0.0_SP,0.0_SP/)
    real(sp),allocatable :: d_temp(:,:)

    if ( cosmo_data%l < 0 ) then
       write(6,*)'cpt_constraints says: Who was the idiot who called me?'
       return
    end if

    nse =  size(surface%secrd,2)
    ldim=(cosmo_data%l+1)**2
    if ( .not. associated(cosmo_data%D) ) then
       allocate(cosmo_data%D(1:nse,1:ldim))
    end if

    ! If this becomes a bottle neck we could be faster here
    allocate(d_temp(1:ldim,1:nse))
    ! NearField Expansion about origin
    CALL NearField(surface%secrd,ctr,cosmo_data%l,d_temp)
    cosmo_data%D = TRANSPOSE(d_temp)
    deallocate(d_temp)

    ! Q
    if ( associated(cosmo_data%Qinv) ) then
       deallocate(cosmo_data%Qinv)
    end if

    allocate(cosmo_data%Qinv(1:ldim,1:ldim))
    cosmo_data%Qinv =  matmul(transpose(cosmo_data%D), &
         matmul(cosmo_data%f * cosmo_data%Ainv,cosmo_data%D))
    CALL do_syminv(ldim,cosmo_data%Qinv,cosmo_data%verbose)


    ! lambda
    if ( associated(cosmo_data%lambda) ) then
       deallocate(cosmo_data%lambda)
    end if

    allocate(cosmo_data%lambda(1:ldim))

    cosmo_data%lambda = matmul(cosmo_data%Qinv,&
         zrho + matmul(transpose(cosmo_data%D),&
         matmul(cosmo_data%f * cosmo_data%Ainv, phi)))

  end subroutine cpt_constraints


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_R(cosmo_data,B,Z)

    ! compute auxillary constraint matrix R

    type(cosmo_t), intent(inout) :: cosmo_data
    REAL(SP), INTENT(IN) :: B(:,:),Z(:,:)

    if ( associated(cosmo_data%R) ) then
       return
    end if

    allocate(cosmo_data%R(1:size(Z,1),1:size(Z,2)))

    cosmo_data%R = Z + matmul(transpose(cosmo_data%D), &
         matmul(cosmo_data%f * cosmo_data%Ainv,B ))

  end subroutine cpt_R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_cosmo(cosmo_data,surface, FILE, UNIT, FORMATTED)

    use FileMod
    use StringMod
    use errormod
    use surfacemod
    use array_io
    IMPLICIT NONE
    type(cosmo_t), intent(in) :: cosmo_data
    type(surface_t),intent(in) :: surface
    character(len=*), optional :: FILE
    integer(i4b), optional :: UNIT
    logical, optional :: FORMATTED

    integer(i4b) :: unit_l
    character(len=8) :: form

    unit_l = -1

    if (present(FILE)) then
       if (optional_flag(FORMATTED,DEFAULT=.false.)) then
          call OpenFile(FILE,unit_l,STAT='replace')
       else
          call OpenFile(FILE,UNIT=unit_l,FORM='unformatted',STAT='replace')
       end if
    elseif (present(UNIT)) then
       unit_l = UNIT
    else
       unit_l = 6
    end if

    call write_surface(surface, UNIT=unit_l)

    ! scalars
    inquire(UNIT=unit_l,FORM=form)

    if ( CompareString(form,'form') ) then
       ! Formatted I/O
       ! input
       write(unit_l,*)cosmo_data%e1, cosmo_data%e2, cosmo_data%l
       write(unit_l,*)cosmo_data%variational !, cosmo_data%volume_polarization
       ! output
       write(unit_l,*)cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess
    elseif ( CompareString(form,'unform') ) then
       ! Formatted I/O
       ! input
       write(unit_l)cosmo_data%e1, cosmo_data%e2, cosmo_data%l
       write(unit_l)cosmo_data%variational!, cosmo_data%volume_polarization
       ! output
       write(unit_l)cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess
    else
       call error('unable to determine file format in write_surface')
    end if

    ! arrays
    ! input
    ! output
    call write_array(cosmo_data%surface_charge,unit_l)
    call write_array(cosmo_data%phi_sigma_rho,unit_l)
    call write_array(cosmo_data%Gpol,unit_l)
    call write_array(cosmo_data%A,unit_l)
!    call write_array(cosmo_data%said,unit_l)
    call write_array(cosmo_data%Ainv,unit_l)
    call write_array(cosmo_data%D,unit_l)
    call write_array(cosmo_data%Qinv,unit_l)
    call write_array(cosmo_data%R,unit_l)
    call write_array(cosmo_data%lambda,unit_l)

    if ( unit_l /= 6 ) then
       close(unit_l)
    end if

  end subroutine write_cosmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_cosmo(cosmo_data,surface, FILE, UNIT, FORMATTED)

    use FileMod
    use StringMod
    use errormod
    use surfacemod
    use array_io
    IMPLICIT NONE

    type(cosmo_t), intent(out) :: cosmo_data
    type(surface_t),intent(out) :: surface
    character(len=*), optional :: FILE
    integer(i4b), optional :: UNIT
    logical, optional :: FORMATTED

    integer(i4b) :: unit_l, iostat
    character(len=8) :: form

    unit_l = -1

    if (present(FILE)) then
       if (optional_flag(FORMATTED,DEFAULT=.false.)) then
          call OpenFile(FILE,UNIT=unit_l,STAT='old',IOSTAT=iostat)
       else
          call OpenFile(&
               FILE,UNIT=unit_l,FORM='unformatted',STAT='old',IOSTAT=iostat)
       end if

       if (iostat > 0) then
          call error('Error opening file in read_surface')
       end if

    elseif (present(UNIT)) then
       unit_l = UNIT
    else
       unit_l = 5
    end if

    call read_surface(surface, UNIT=unit_l)

    ! scalars
    inquire(UNIT=unit_l,FORM=form)

    if ( CompareString(form,'form') ) then
       ! Formatted I/O
       ! input
       read(unit_l,*)cosmo_data%e1, cosmo_data%e2, cosmo_data%l
       read(unit_l,*)cosmo_data%variational!, cosmo_data%volume_polarization
       ! output
       read(unit_l,*)cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess
    elseif ( CompareString(form,'unform') ) then
       ! Formatted I/O
       ! input
       read(unit_l)cosmo_data%e1, cosmo_data%e2, cosmo_data%l
       read(unit_l)cosmo_data%variational!, cosmo_data%volume_polarization
       ! output
       read(unit_l)cosmo_data%Epol,cosmo_data%Eself,cosmo_data%Ess
    else
       call error('unable to determine file format in read_surface')
    end if

    ! arrays
    ! input
!    call read_array(cosmo_data%B,unit_l)
!    write(6,*)'Read B ',size(cosmo_data%B,1),' by ',size(cosmo_data%B,2)
!    call read_array(cosmo_data%rho,unit_l)
!    write(6,*)'Read rho ',size(cosmo_data%rho)
!    call read_array(cosmo_data%phi_rho_sigma,unit_l)
!    write(6,*)'Read phi_rho_sigma ',size(cosmo_data%phi_rho_sigma)
!    call read_array(cosmo_data%Z,unit_l)
!    write(6,*)'Read Z ',size(cosmo_data%Z,1),' by ',size(cosmo_data%Z,2)
    ! output
    call read_array(cosmo_data%surface_charge,unit_l)
    write(6,*)'Read surface charge ',size(cosmo_data%surface_charge)
    call read_array(cosmo_data%phi_sigma_rho,unit_l)
    write(6,*)'Read phi_sigma_rho ',size(cosmo_data%phi_sigma_rho)
    call read_array(cosmo_data%Gpol,unit_l)
    write(6,*)'Read G ',size(cosmo_data%Gpol,1),' by ',size(cosmo_data%Gpol,2)
    call read_array(cosmo_data%A,unit_l)
    write(6,*)'Read A ',size(cosmo_data%A,1),' by ',size(cosmo_data%A,2)
!    call read_array(cosmo_data%said,unit_l)
    call read_array(cosmo_data%Ainv,unit_l)
    write(6,*)'Read Ainv ',size(cosmo_data%Ainv,1),' by ',size(cosmo_data%Ainv,2)
    call read_array(cosmo_data%D,unit_l)
    call read_array(cosmo_data%Qinv,unit_l)
    call read_array(cosmo_data%R,unit_l)
    call read_array(cosmo_data%lambda,unit_l)
    
    close(unit_l)

  end subroutine read_cosmo

!!$  SUBROUTINE MinimizeEnergy(cosmo_data)
!!$
!!$    TYPE(cosmo_t) :: cosmo_data
!!$
!!$    INTEGER(I4B) :: iter, nse
!!$    REAL(SP) :: E, Eold, alpha_k
!!$    REAL(SP) :: xk(SIZE(cosmo_data%phi_rho_sigma))
!!$    REAL(SP) :: Deltaxk(SIZE(cosmo_data%phi_rho_sigma))
!!$    REAL(SP) :: A_xk(SIZE(cosmo_data%phi_rho_sigma))
!!$    REAL(SP) :: A_Deltaxk(SIZE(cosmo_data%phi_rho_sigma))
!!$
!!$    nse = SIZE(cosmo_data%phi_rho_sigma)
!!$
!!$    ! Initial guess
!!$    xk = -SUM(cosmo_data%rho) / nse
!!$    ! Replace matrix multiply with recursive bisection method
!!$    A_xk = MATMUL( cosmo_data%A, xk )
!!$
!!$    E = 0.0_SP
!!$
!!$    DO iter=1,nse
!!$
!!$       Eold = E
!!$       E = 0.5_SP * DOT_PRODUCT( xk, A_xk ) &
!!$            & + DOT_PRODUCT( xk, cosmo_data%phi_rho_sigma )
!!$
!!$       write(*,*)'iter, E, diff, SUM = ',iter,E,E-Eold,SUM(xk)
!!$
!!$       Deltaxk = A_xk + cosmo_data%phi_rho_sigma
!!$       ! Replace matrix multiply with recursive bisection method
!!$       A_Deltaxk = MATMUL( cosmo_data%A, Deltaxk )
!!$       alpha_k = - DOT_PRODUCT( Deltaxk, A_xk + cosmo_data%phi_rho_sigma ) / &
!!$            & DOT_PRODUCT( Deltaxk, A_Deltaxk )
!!$
!!$       xk   = xk + alpha_k * Deltaxk
!!$
!!$       A_xk = A_xk + alpha_k * A_Deltaxk
!!$
!!$    END DO
!!$
!!$  END SUBROUTINE MinimizeEnergy
!!$  

  SUBROUTINE ConjugateGradient(Surface, f, minpara, &
       b, x, iter, err, Epol,Eself,Ess, Amatrix)
!!!****f* MinimizeMod/ConjugateGradient
!!! 
!!! DESCRIPTION
!!!     Solve A*x = -b for x by minimization of a quadratic function 
!!!               T        T
!!!     f(x)=1/2 x*A*x + b*x using preconditioned conjugate gradient
!!!     algorithm described in Numerical Recipes pp. 78-82,
!!!     where the A is a symmetric and positive definite matrix, and
!!!     x and b are two vectors. In our version, 
!!!     a symmetric preconditioned matrix A' is used,
!!!     A' = A_c^{-1/2} * A * A_c^{-1/2}
!!!     where A_c is the preconditioner matrix
!!! INPUT
!!!     f - dielectric scale factor=(e2-e1)/e2*e1
!!!     b - a vector of potential at surface points due to solute charges
!!!         same as phi_rho_sigma = B * rho 
!!!     itmax - maximal iteration allowed
!!!     itol - option for convergence
!!!     icond - option for choosing preconditioner matrix
!!! INPUT and OUTPUT
!!!     x - a vector of surface charges
!!! OUTPUT
!!!     err - error in minimization energy 
!!! CALLS
!!!     atimes - computes A * x 
!!!     asolve - solves A_c * z = r for z
!!!***

    USE UtilitiesMod, ONLY: assert_eq
    USE SurfaceMod

    TYPE(surface_t), INTENT(IN) :: surface
    TYPE(minpara_t), INTENT(IN) :: minpara

    INTEGER(I4b), INTENT(OUT) :: iter
    INTEGER(I4b):: n 
    REAL(SP), INTENT(IN) :: f
    REAL(SP), INTENT(IN) :: b(:)
    REAL(SP), INTENT(INOUT) :: x(:)
    REAL(SP), INTENT(OUT) :: err, Epol,Eself,Ess
    REAL(SP), PARAMETER :: EPS=1.0E-14_SP
    REAL(SP), OPTIONAL,INTENT(IN) :: Amatrix(:,:)
    REAL(SP) ::ak,akden,bk,bkden,bknum,bnrm,dxnrm,zm1nrm,xnrm,znrm,Eold
    REAL(SP) :: p(SIZE(b)), z(SIZE(b)),r(SIZE(b)) ,AtimesX(SIZE(b))       
    REAL(SP), POINTER :: A_c(:) => NULL()

!!!**** Quick check for trivial b=0 case
    n = assert_eq(SIZE(b),SIZE(x),'ConjugateGradient')

    Epol = 0.0_SP
    Eself = 0.0_SP
    Ess = 0.0_SP
    iter=0
    IF ( ALL(b == 0.0) ) THEN
       WRITE(6,*) 'b vector is zero'
       x = 0.0_sp
       err = 0.0_sp
       RETURN
    END IF

!!!**** Compute the preconditioner matrix for later use
    CALL CptPreconditioner(surface,f,minpara%precond,A_c)!,Amatrix)

!!!**** r = 1/f(e2) *  A * x   
    IF (PRESENT(Amatrix)) THEN 
       r = MATMUL(Amatrix,x)
    ELSE
       CALL atimes(surface,x,r,minpara%multipole,minpara%WS,minpara%L) 
    END IF
    r = r/f
    AtimesX = r
    Epol = 0.5_SP*DOT_PRODUCT(x,AtimesX)+DOT_PRODUCT(x, b)

!!!**** r_1=-(b_1+A*x_0)
    r=-(b+r)

    IF(minpara%itol <=1)THEN
       bnrm=VectorNorm(b,minpara%itol)
       CALL asolve(surface,A_c,r,z,minpara%precond)
    ELSEIF(minpara%itol ==2)THEN
       call asolve(surface,A_c,b,z,minpara%precond) !solves Az=b
       bnrm=VectorNorm(z,minpara%itol)
       CALL asolve(surface,A_c,r,z,minpara%precond)
    ELSE
       call asolve(surface,A_c,b,z,minpara%precond)
       bnrm=VectorNorm(z,minpara%itol)
       call asolve(surface,A_c,r,z,minpara%precond)
       znrm=VectorNorm(z,minpara%itol)
    END IF

    IF (Surface%verbose) WRITE(6,*)'          Energy,   Diff,   Sum(X),   err'
!!!**** main loop
    Iteration: DO
       IF (iter > minpara%itmax) EXIT
       iter=iter+1

!!!**** calculate coefficient bk and direction vector
       bknum = DOT_PRODUCT(z,r)
       IF (iter == 1) then
          p=z
       ELSE
          bk=bknum/bkden
          p = bk * p + z
       END IF
       bkden=bknum

!!!**** calculate coefficient ak, new direction p
       IF (PRESENT(Amatrix)) THEN
          z = MATMUL(Amatrix,p)
       ELSE
          CALL atimes(surface,p,z,minpara%multipole,minpara%WS,minpara%L)
       END IF
       z = z/f

       akden = DOT_PRODUCT(z,p)
       ak=bknum/akden

!!!**** new x, r and E 
       x = x + ak*p             ! compute new x
       r = r - ak*z             ! compute new residual
       AtimesX = AtimesX + ak*z ! compute new A*x

!!!**** check stopping criterion
       IF(minpara%itol <=1)THEN
          err = VectorNorm(r,minpara%itol)/bnrm
       ELSEIF(minpara%itol ==2)THEN
          err = VectorNorm(z,minpara%itol)/bnrm
       ELSE
          zm1nrm=znrm
          znrm = VectorNorm(z,minpara%itol)
          IF (ABS(zm1nrm-znrm) > EPS*znrm) THEN
             dxnrm = ABS(ak)*VectorNorm(p,minpara%itol)
             err = znrm/ABS(zm1nrm-znrm)*dxnrm
          ELSE
             err = znrm/bnrm
             CYCLE
          END IF
          xnrm = VectorNorm(x,minpara%itol)
          IF (err <= 0.5_SP*xnrm) THEN
             err = err/xnrm
          ELSE
             err = znrm/bnrm
             CYCLE
          END IF
       ENDIF

       IF (Surface%verbose) THEN
          Eold = Epol
          Eself = 0.5_SP*DOT_PRODUCT(x,AtimesX)
          Ess = DOT_PRODUCT(x, b)
          Epol = Eself + Ess
          IF(err <= minpara%tol)THEN
             WRITE (*,'(a11,i4,4e17.8)') 'Last iter= ',iter,Epol,Epol-Eold,SUM(x),err
             EXIT
          ELSE
             WRITE(6,'(A,I5,4ES11.3)') 'iter= ', iter,Epol,Epol-Eold,SUM(x),err 
          endif
       ELSEIF(err <= minpara%tol)THEN
          Eself = 0.5_SP*DOT_PRODUCT(x,AtimesX)
          Ess = DOT_PRODUCT(x, b)
          Epol = Eself + Ess
          EXIT
       endif

!!!**** solve A_c*z = r and 
       call asolve(surface,A_c,r,z,minpara%precond)

    END DO Iteration
    DEALLOCATE (A_c)
  CONTAINS

    REAL(SP) FUNCTION VectorNorm(sx,itol)
!!!****f* ConjugateGradient/VectorNorm
!!!
!!! Calculate vector magnitude norm or the largest component norm
!!!***
      INTEGER(I4B), INTENT(IN):: itol
      REAL(SP), INTENT(IN) :: sx(:)

      IF(itol <=3)THEN
         ! Vector Magnitude Norm
         VectorNorm = SQRT(DOT_PRODUCT(sx,sx))
      ELSE
         ! Largest Component Norm
         VectorNorm = MAXVAL(ABS(sx))
      END IF
    END FUNCTION VectorNorm

  END SUBROUTINE ConjugateGradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE atimes(surface,x,r,IfMultipole, WS, L)
!!!****f* MinimizeMod/atimes
!!!
!!! r = A * x 
!!! INPUT: x, the surface charge vector,
!!! OUTPUT: r, the vector of potential at surface points due to surface charges
!!! Scratch:
!!! ChargeCoord(1,i): charge for i-th surface elements
!!! ChargeCoord(2:4,i): x,y,z coordinates for i-th surface elements
!!! ChargeCoord(5,1:i): Gaussian exponent for i-th surface elements
!!! r(i): potential at i-th surface element i
!!!*** 
    USE ConstantsMod, ONLY: SQRT2, SQRT_PI
    USE ErrorMod
    USE UtilitiesMod, ONLY: assert_eq
    USE RBM
    USE SurfaceMod 

    TYPE(surface_t), INTENT(IN) :: surface   
    INTEGER(I4B), INTENT(IN) :: IfMultipole,L
    REAL(SP), INTENT(IN) :: WS
    REAL(SP), INTENT(IN) :: x(:)
    REAL(SP), INTENT(OUT) :: r(:)
    REAL(SP), PARAMETER :: SQRT_TWO_over_PI = SQRT2/SQRT_PI
    REAL(SP), TARGET :: ChargeCoord(1:5,1:SIZE(x))
    INTEGER(I4B) :: n

    n=assert_eq(SIZE(x),SIZE(r),SIZE(surface%SECRD,2), 'atimes')
    ChargeCoord(1,1:n) = x(1:n)
    ChargeCoord(2:4,1:n) = Surface%SECRD(1:3,1:n)
    ChargeCoord(5,1:n) = Surface%cosmozeta(1:n)

    r = 0.0_SP
    IF ( .NOT. ALL(x == 0.0_SP) ) THEN
       SELECT CASE(IfMultipole)
!!!**** Brute force. potential due to Gaussian srf charges
       CASE (0) 
          CALL DirectCoulomb(ChargeCoord,r)
!!!**** potential at surface points with recursive bisection fast multipole
       CASE (1)
          CALL RecursiveBisectionMethod(ChargeCoord,r,WS=WS,L=L)
       CASE  DEFAULT 
          CALL error('atimes: illegal ifmultipole')
       END SELECT

       if ( associated(surface%sescale) ) then
          if ( any(surface%sescale /= 1.0_sp) ) then
             IF (surface%verbose) THEN
                WRITE(6,*)'Correcting for switching'
             END IF
             call switch_pot(r, x, surface%sescale, surface%cosmozeta)
          end if
       end if
    END IF

  END SUBROUTINE atimes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CptPreconditioner(surface,f,which_conditioner,A_c)!,Amatrix)
!!!
!!! Description
!!!     used for symmetric preconditioning:
!!!     A_c^{-1/2}*A*A_c^{-1/2} \sim 1
!!!     Therefore, we need A_c^{-1/2}
!!!     where A_c is a preconditioner matrix

    USE ConstantsMod, ONLY: SQRT2, SQRT_PI
    USE SurfaceMod    
    TYPE(surface_t), INTENT(IN) :: surface
    INTEGER(I4B), INTENT(IN) :: which_Conditioner
    INTEGER(I4B) :: nse,iatm,Natm,NseAtm,NseBlock,first,last
    REAL(SP), INTENT(IN) :: f
    REAL(SP), POINTER :: A_c(:)
    REAL(SP), PARAMETER :: SQRT_TWO_over_PI = SQRT2/SQRT_PI

    nse = SIZE(Surface%cosmozeta)

    SELECT CASE (which_Conditioner)
!!!**** no preconditioner
    CASE (0)  
       IF ( .NOT. ASSOCIATED(A_c)) ALLOCATE(A_c(1:nse))  
       A_c = 1.0_sp
       RETURN

!!!**** A_diag as preconditioner. A_c stores the A_diag^{-1}
    CASE (1)  
       IF ( .NOT. ASSOCIATED(A_c)) ALLOCATE(A_c(1:nse))
       IF(org_Sik_in_A0_diag)THEN
          A_c = Surface%sescale/(Surface%cosmozeta * SQRT_TWO_over_PI)
       ELSE
          A_c = SQRT(Surface%sescale)/(Surface%cosmozeta * SQRT_TWO_over_PI)
       ENDIF
!!!**** preconditioner  =  atom block diagonal matrix of A 
!!!     lower triangle parts stored in the sparse array A_c
    CASE (2) 
!!!**** Compute the length of A_c, which is the total number of 
!!!     elements in the trianglular atom block matrices
       Natm = SIZE(Surface%asaindex)
       NseBlock = 0; first = 1
       DO iatm = 1, Natm
          last = surface%asaindex(iatm)
          NseAtm = last - first +1
          first = last + 1
          NseBlock = NseBlock + 0.5_SP*NseAtm*(NseAtm+1)
       END DO
       ALLOCATE (A_c(NseBlock))

       CALL AsolveBlockDiag(surface,A_c)
    CASE DEFAULT 
       IF ( .NOT. ASSOCIATED(A_c)) ALLOCATE(A_c(1:nse))  
       A_c = 1.0_sp
       RETURN
    END SELECT

    A_c = f* A_c

  END SUBROUTINE CptPreconditioner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE asolve(Surface,A_c,r,z,which_conditioner)
!!!****f* MinimizeMod/asolve
!!!
!!! Description
!!!     used for symmetric preconditioning:
!!!     z = A_c^{-1/2} * r
!!!     where A_c is a preconditioner matrix
!!!***

    USE ErrorMod
    USE SurfaceMod    
    TYPE(surface_t), INTENT(IN) :: surface

    INTEGER(I4B), INTENT(IN) :: which_conditioner
    INTEGER(I4B) :: NAtm, iatm, i, j, first, last, NseAtm, NseBlock, &
         indexij, MaxPts
    REAL(SP), INTENT(IN) :: r(:)
    REAL(SP), INTENT(IN) :: A_c(:)
    REAL(SP), INTENT(OUT) :: z(:)
    REAL(SP), ALLOCATABLE :: AtomBlock(:,:)

    SELECT CASE (which_Conditioner)

!!!**** no preconditioner
    CASE (0) 
       z = r

!!!**** A_diag as Preconditioner
    CASE (1)         
       z = A_c * r

!!!**** atom block diagonal matrix as Preconditioner
    CASE (2) 

       Natm =  SIZE(surface%asaindex)
       MAXPts = 1; first = 1

       DO iatm = 1, Natm
          last = surface%asaindex(iatm)
          NseAtm = last - first +1
          MaxPts = MAX(NseAtm, MaxPts)
          first = last + 1
       END DO

       ALLOCATE (AtomBlock(MaxPts, MaxPts))

!!!**** Retrieve from the sparse storage array A_c and multiply r to obtain z
       NseBlock = 0; first = 1
       DO iatm = 1, NAtm
          last = surface%asaindex(iatm)
          NseAtm = last - first +1
          DO i = 1, NseAtm
             indexij = NseBlock + i*(i-1)/2 + i
             AtomBlock(i,i) = A_c(indexij)
             DO j = i-1, 1, -1
                indexij = NseBlock + i*(i-1)/2 + j
                AtomBlock(i,j) = A_c(indexij)
                AtomBlock(j,i) = AtomBlock(i,j)
             END DO
          END DO
          z(first:last) = MATMUL(AtomBlock(1:NseAtm, 1:NseAtm),r(first:last))
          first = last + 1 ! index for next atom
          NseBlock =  NseBlock + 0.5_SP*NseAtm*(NseAtm+1)
       END DO
       DEALLOCATE (AtomBlock)
    CASE DEFAULT 
!!!**** no preconditioner
       z = r
    END SELECT

  END SUBROUTINE asolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AsolveBlockDiag(surface,A_c)!,Amatrix)
!!!****f* MinimizeMod/AsolveBlockDiag
!!!
!!! Description
!!!     compute A_c^{-1/2} where
!!!     A_c is the atom block diagonal matrix A_ik,il
!!!***
    USE ConstantsMod, ONLY: SQRT2, SQRT_PI
    USE SurfaceMod
    USE Cosmo_util
    USE CoulombMod, only : GaussianInteraction    

    TYPE(surface_t), INTENT(IN) :: surface
    INTEGER(I4B) :: iatm, i, j, k, ii, jj, indexij, MaxPts, Natm
    INTEGER(I4B) :: first, last, NseAtm, NseBlock
    REAL(SP), ALLOCATABLE :: diag(:)
    REAL(SP), pointer :: AtomBlock(:,:)
    REAL(SP), PARAMETER :: SQRT_TWO_over_PI = SQRT2/SQRT_PI
    REAL(SP), INTENT(INOUT):: A_c(:)
    REAL(SP) :: zetai, zetaj, zeta, r2, rtmp

    Natm =  SIZE(surface%asaindex)
    MAXPts = 1; first = 1

    DO iatm = 1, Natm
       last = surface%asaindex(iatm)
       NseAtm = last - first +1
       MaxPts = MAX(NseAtm, MaxPts)
       first = last + 1
    END DO

    ALLOCATE (AtomBlock(MaxPts, MaxPts), diag(MaxPts))

    first = 1; NseBlock = 0;

!!!**** big loop over all atoms
    AtomLoop: DO iatm = 1, Natm
       last = surface%asaindex(iatm)
       NseAtm = last - first +1

!!!****  construct atom block matrix AtomBlock, and store in diag 
!!!      the inverse square root of A_diag and zero out the diagonal
!!!      elements of AtomBlock to avoid infinity

       IF (NseAtm /=0) THEN
          DO i = first, last
             zetai = surface%cosmozeta(i)
             ii = i - first + 1
             IF(org_Sik_in_A0_diag)THEN
                diag(ii) =  SQRT(Surface%sescale(i)/(SQRT_TWO_over_PI *zetai))
             ELSE
                diag(ii) =  SQRT(SQRT(Surface%sescale(i))/(SQRT_TWO_over_PI *zetai))
             ENDIF
             IF (surface%gammaswitch > 0.0_sp) THEN
                AtomBlock(ii,ii) = 0.0_sp
             ELSE 
                AtomBlock(ii,ii) = SQRT_TWO_over_PI *zetai
             END IF

             DO j = i+1, last
                zetaj = surface%cosmozeta(j)
                zeta = zetai*zetaj/SQRT(zetai*zetai+zetaj*zetaj)
                r2 = 0.0_sp
                DO k = 1, 3
                   r2 = r2 + (surface%secrd(k,i) - surface%secrd(k,j))**2
                END DO
                rtmp = SQRT(r2)      
                jj = j - first + 1
                CALL GaussianInteraction(zeta,rtmp,AtomBlock(ii,jj)) 
                AtomBlock(jj,ii) = AtomBlock(ii,jj)       
             END DO ! j
          END DO ! i


          IF (surface%gammaswitch > 0.0_sp) THEN

!!!**** construct the symmetric matrix to invert: 
!!!     I + A_diag^{-1/2} * A_off *  A_diag^{-1/2}
!!!     and store in AtomBlock 
             AtomBlock(1:NseAtm,1:NseAtm) = diagmul(diag(1:NseAtm), &
                  diagmul(AtomBlock(1:NseAtm,1:NseAtm), diag(1:NseAtm)))
             FORALL(i=1:NseAtm)
                AtomBlock(i,i) = AtomBlock(i,i) + 1.0_sp
             END FORALL

!!!**** invert the above matrix
             CALL do_syminv(NseAtm,AtomBlock(1:NseAtm,1:NseAtm),.false.)

!!!**** multiply the left and right sides of the above matrix by 
!!!     A_diag^{-1/2} to obtain A^{-1} 
             AtomBlock(1:NseAtm,1:NseAtm) = diagmul(diag(1:NseAtm),&
                  diagmul(AtomBlock(1:NseAtm,1:NseAtm),diag(1:NseAtm)))

          ELSE ! direct inversion
             CALL do_syminv(NseAtm,AtomBlock(1:NseAtm,1:NseAtm),.false.)
          END IF

!!!**** store in a big sparse storage array A_c that contains all A^{-1/2}'s     
          DO i = 1, NseAtm 
             DO j = i, 1, -1
                indexij = NseBlock + i*(i-1)/2 + j
                A_c(indexij) = AtomBlock(i,j)
             END DO
          END DO

          first = last + 1 ! index for next atom
          NseBlock = NseBlock + 0.5_SP*NseAtm*(NseAtm+1)

       END IF

    END DO AtomLoop

    DEALLOCATE(AtomBlock,diag)

  END SUBROUTINE AsolveBlockdiag



end module cosmomod
     
     
