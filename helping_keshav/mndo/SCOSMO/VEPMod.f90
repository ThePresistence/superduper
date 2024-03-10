Module VEPMod

  USE datatypes
  USE constantsMod
  implicit none

  private

  LOGICAL(LGD), public,save :: crazy_scale = .FALSE.
  
  PUBLIC :: VEP_Surf_T, deallocate_VEP_surface, Init_VEP_surface, &
       vep_direct, vep_rvm, set_r_omega, &
       direct_potential_model_surf, direct_potential_model_pts

  TYPE::VEP_Surf_T
     ! We dont need a whole lot of fancy stuff for VEP
     REAL(SP),pointer :: Crd(:,:)=>NULL()
     REAL(SP),pointer :: zeta(:)=>NULL()
     REAL(SP),pointer :: wts(:)=>NULL()
  END TYPE VEP_Surf_T

  INTERFACE vep_rvm
     Module Procedure  vep_rvm_pts,vep_rvm_surf,vep_rvm_surfpt
  END INTERFACE
contains

  subroutine Init_VEP_surface(surface,aqtype,nse,radius,ctr,debug)
    USE surfaceMod
    TYPE(VEP_Surf_T), INTENT(INOUT) :: surface 
    INTEGER(I4B),INTENT(IN) :: aqtype,nse
    REAL(SP), INTENT(IN) :: radius(:),ctr(:,:)
    LOGICAL(LGD),INTENT(IN) :: debug

    call deallocate_VEP_surface(surface)
    ! Generate Surface
    CALL SASurface(ctr, radius, AQTYPE=aqtype, DT=0, PTS_PER_SPHERE=nse, &
         PROBERADIUS=0.0_SP, GAMMASWITCH=0.0_SP, VERBOSE=debug, &
         SECRD=surface%crd,COSMOZETA=surface%zeta,WTS=surface%wts)
  end subroutine Init_VEP_surface

  subroutine deallocate_VEP_surface(surface)
    TYPE(VEP_Surf_T), INTENT(INOUT) :: surface 
    IF(Associated(surface%crd))Deallocate(surface%crd)
    IF(Associated(surface%zeta))Deallocate(surface%zeta)
    if(associated(surface%wts))deallocate(surface%wts)
  end subroutine deallocate_VEP_surface

 
  ! VEP method for surface charges
  subroutine vep_direct(crd,q,VEP_surf,vep_q,l_vep,verbose,Ainv_keep,Qinv_keep,D_keep)
    use SurfaceMod, only : cpt_A, cpt_Ainv
    use cosmo_util, only : do_syminv
    use cosmo_helper, only : cpt_phirhosigma_for_pt
    USE MultipoleMod, only : FarField,qtimesFarField
    IMPLICIT NONE
    REAL(SP),INTENT(IN) ::crd(:,:),q(:)
    TYPE(VEP_Surf_t),INTENT(IN) :: VEP_surf
    REAL(SP), INTENT(OUT) :: vep_q(:)
    INTEGER(I4B),intent(IN) :: l_vep
    LOGICAL(LGD),intent(IN) :: verbose
    REAL(SP),optional, pointer :: Ainv_keep(:,:),D_keep(:,:),Qinv_keep(:,:)
    real(sp), parameter :: ctr(3) = (/0.0_SP,0.0_SP,0.0_SP/)
    REAL(SP), pointer :: phi_tmp(:)=>NULL(),VEP_A(:,:)=>NULL(), VEP_ainv(:,:)=>NULL()
    REAL(SP), allocatable :: Cons_mat(:,:),Zrho(:),R_rho(:),lambda(:)
    REAL(SP), pointer :: D(:,:)=>NULL(),Qmat(:,:)=>NULL()
    INTEGER(I4B) :: ldim, nse

    if(verbose)WRITE(6,*)'Computing VEP Surface Charges'
    CALL cpt_phirhosigma_for_pt(Crd,q,VEP_surf%crd,VEP_surf%zeta,phi_tmp,2)

    CALL cpt_A(VEP_surf%crd,VEP_surf%zeta,VEP_A)
    ! vep_A is null on return from this routine
    CALL cpt_Ainv(VEP_A,VEP_Ainv,.TRUE., verbose)

    IF(l_vep < 0)THEN
       VEP_q = matmul(VEP_Ainv,phi_tmp)
    ELSE
       if(verbose)WRITE(6,*)'VEP Surface Charges constrained to l=',l_vep
       ! try to reduce everything to vectors so we dont store
       ! matricies
       ldim=(l_vep+1)**2

       ! 1st compute Constraints on external charges
       allocate(zrho(ldim))
       CALL qtimesFarField(q,Crd,ctr,l_vep,zrho)

       ! Now compute constraints on surface charges 
       nse=SIZE(VEP_surf%crd,2) 
       allocate(Cons_mat(1:ldim,1:nse))
       CALL FarField(VEP_surf%crd,ctr,l_vep,Cons_mat)

       ! Compute Q and R_rho
       allocate(D(1:nse,1:ldim))
       D = TRANSPOSE(Cons_mat)

       allocate(Qmat(1:ldim,1:ldim))
       Qmat = matmul(Cons_mat,matmul(VEP_Ainv,D))
       allocate(R_rho(ldim))
       
       R_rho =  ZRho - matmul(Cons_mat,matmul(VEP_Ainv, phi_tmp))
       deallocate(Cons_mat,zrho)
       allocate(lambda(1:ldim))
       ! Symmetryc inverse is more efficient/ better results than svdinverse
       CALL do_syminv(ldim,Qmat,verbose)
       lambda = matmul(Qmat,R_rho)
       deallocate(R_rho)

       IF(present(Qinv_keep))THEN
          IF(associated(Qinv_keep))deallocate(Qinv_keep)
          Qinv_keep => Qmat
          Qmat=>NULL()
       ELSE
          deallocate(qmat)
       endif

       phi_tmp = phi_tmp + matmul(D, lambda)
       deallocate(lambda)

       IF(present(d_keep))THEN
          IF(associated(d_keep)) deallocate(d_keep)
          d_keep => D
          D=>null()
       ELSE
          deallocate(D)
       ENDIF

       VEP_q = matmul(VEP_Ainv,phi_tmp)
    ENDIF
    deallocate(phi_tmp)

    IF(present(Ainv_keep))THEN
       if(associated(Ainv_keep))deallocate(Ainv_keep)
       Ainv_keep => VEP_Ainv
       VEP_Ainv => NULL()
    ELSE
       deallocate(VEP_ainv)
    endif

  end subroutine vep_direct


  ! VEP-RVM method for point charges
  subroutine vep_rvm_pts(crd,q,S_vep,rvm_crd,rvm_q,l_vep,verbose,tol)
    use cosmo_helper, only : cpt_B_for_pt
    IMPLICIT NONE
    REAL(SP),INTENT(IN) ::crd(:,:),q(:),rvm_crd(:,:)
    TYPE(VEP_Surf_t),INTENT(IN) :: S_vep
    REAL(SP), INTENT(OUT) :: rvm_q(:)
    INTEGER(I4B),INTENT(IN) :: l_vep
    LOGICAL(LGD), INTENT(IN) :: verbose
    REAL(SP),intent(in):: tol
    REAL(SP), pointer :: Maux(:,:)=>NULL()

    if(verbose)WRITE(6,*)'Computing RVM atom Charges from VEP Surface Charges'
    CALL cpt_B_for_pt(rvm_crd, S_vep%crd,S_vep%zeta, Maux, 1)
    CALL cpt_RVM_q(crd,q,S_vep, Maux, RVM_crd,RVM_q,l_vep, verbose,tol)
    deallocate(Maux)

  end subroutine vep_rvm_pts

  ! VEP-RVM method for surface charges
  subroutine vep_rvm_surf(crd,q,S_vep,S_rvm,rvm_q,l_vep,verbose,tol)
    use surfaceMod, only : cpt_A
    IMPLICIT NONE
    REAL(SP),INTENT(IN) ::crd(:,:),q(:)
    TYPE(VEP_Surf_t),INTENT(IN) :: S_vep, S_rvm
    REAL(SP), INTENT(OUT) :: rvm_q(:)
    INTEGER(I4B),INTENT(IN) :: l_vep
    LOGICAL(LGD), INTENT(IN) :: verbose
    REAL(SP),intent(in):: tol
    REAL(SP), pointer :: Maux(:,:)=>NULL()

    if(verbose)WRITE(6,*)'Computing RVM Surface Charges from VEP Surface Charges'
    CALL cpt_A(S_vep%crd,S_vep%zeta,S_rvm%crd,S_rvm%zeta,Maux)
    CALL cpt_RVM_q(crd,q,S_vep,Maux,S_rvm%crd,RVM_q,l_vep,verbose, tol)
    deallocate(Maux)
    
  end subroutine vep_rvm_surf

  subroutine vep_rvm_surfpt(org_crd,org_q,S_vep,S_rvm,rvm_surfq, &
       rvm_ptcrd,rvm_ptq,l_vep,verbose,tol)
    use cosmo_helper, only : cpt_B_for_pt
    use surfaceMod, only : cpt_A
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: org_crd(:,:),org_q(:),rvm_ptcrd(:,:)
    TYPE(VEP_Surf_t),INTENT(IN) :: S_vep, S_rvm
    REAL(SP), INTENT(OUT) :: rvm_surfq(:),rvm_ptq(:)
    INTEGER(I4B),INTENT(IN) :: l_vep
    LOGICAL(LGD), INTENT(IN) :: verbose
    REAL(SP),intent(in):: tol
    REAL(SP), pointer :: tmpM(:,:)=>NULL()
    REAL(SP), allocatable :: Maux(:,:),tmp_q(:),tmp_crd(:,:)
    INTEGER(I4B) :: nsq,npq,nqt

    nsq = SIZE(rvm_surfq)
    npq = size(rvm_ptq)
    nqt = nsq + npq

    allocate(Maux(SIZE(S_vep%zeta),nqt))

    CALL cpt_A(S_vep%crd,S_vep%zeta,S_rvm%crd,S_rvm%zeta,tmpM)
    Maux(:,1:nsq) = tmpM
    deallocate(tmpM)    

    CALL cpt_B_for_pt(rvm_ptcrd, S_vep%crd,S_vep%zeta, tmpM, 1)
    Maux(:,1+nsq:nqt) = tmpM
    deallocate(tmpM)

    allocate(tmp_q(nqt),tmp_crd(3,nqt))
    tmp_crd(1:3,1:nsq) = S_rvm%crd
    tmp_crd(1:3,nsq+1:nqt) = rvm_ptcrd

    CALL cpt_RVM_q(org_crd,org_q,S_vep,Maux,tmp_crd,tmp_q,l_vep,verbose,tol)
    deallocate(Maux,tmp_crd)

    rvm_surfq = tmp_q(1:nsq)
    rvm_ptq   = tmp_q(1+nsq:nqt)
    deallocate(tmp_q)

  end subroutine vep_rvm_surfpt

  ! subroutine to find a set of auxillary atom charges that 
  ! reproduces potential at VEP surface (in case we dont want to use the surface)
  ! We will use point charges for atoms and gaussians for surface... (AtomSurf_mode=1)
  ! We want a solution to: 
  !     sigma = A^-1 * B1 * p1 = A^-1 * B2 * p2
  ! where only p2 is unknown. Two methods we will try (sym=[true,false])
  ! 1) Symmetric matrix
  ! Multiply by transpose of B2: This has been removed. It essentially squares the
  !      Contition Number of the Matrix and these matricies are very poorly behaved
  ! 2) Assymmetric matrix
  !   p2 = (A^-1 * B2)^-1 * sigma
  !
  subroutine cpt_RVM_q(crd,q,S_vep,Maux,rvm_crd,rvm_q,l_rvm,verbose,tol)
    USE Datatypes
    USE constantsMod
    USE MultipoleMod, only : FarField
    USE Cosmo_util
    IMPLICIT NONE
    TYPE(VEP_Surf_t),INTENT(IN) :: S_vep
    REAL(SP),INTENT(IN) ::crd(:,:),q(:)
    ! this is inout so we can modify it with constraints if necessary
    REAL(SP), INTENT(INOUT) ::Maux(:,:)
    REAL(SP), INTENT(IN) :: rvm_crd(:,:)
    REAL(SP), INTENT(OUT) :: rvm_q(:)
    INTEGER(I4B),INTENT(IN) :: l_rvm
    LOGICAL(LGD), INTENT(IN) :: verbose
    REAL(SP),intent(in):: tol

    REAL(SP), allocatable :: vep_q(:)
    INTEGER(I4B) :: ldim,ncrd,nvep
    real(sp) :: ctr(3) = (/0.0_SP,0.0_SP,0.0_SP/)
    REAL(SP), pointer :: D(:,:)=>NULL(),Qinv(:,:)=>NULL(),Ainv(:,:)=>NULL()
    REAL(SP), allocatable :: Zmat(:,:)
    
    nvep = SIZE(S_vep%zeta)
    allocate(vep_q(nvep))

    IF(l_rvm<0) THEN
       CALL vep_direct(crd,q,S_vep,vep_q,l_rvm,verbose,Ainv) 
    ELSE
       ! get the constraint matricies we need from VEP_direct
       CALL vep_direct(crd,q,S_vep,vep_q,l_rvm,verbose,Ainv,Qinv,D)

       if(verbose)write(6,*)'Computing Constraints for RVM with l=',l_rvm
       ! now compute constraints on rvm points
       ldim=(l_rvm+1)**2
       ncrd = size(rvm_crd,2)
       allocate(Zmat(1:ldim,1:ncrd))
       CALL FarField(rvm_crd,ctr,l_rvm,Zmat)
       
       !   = A_io + D   *   Q^-1 *(z - D*A_ii^-1*A_io)
       Zmat=Zmat - matmul(transpose(D),matmul(Ainv,Maux))
       D = matmul(D,Qinv)
       deallocate(Qinv)
       Maux=Maux  + matmul(D,Zmat)
       deallocate(Zmat,D)
    Endif

    Maux = matmul(Ainv,Maux)
    call svdinv_times_vector(Maux,vep_q,rvm_q,verbose,tol=tol)
    deallocate(vep_q,Ainv)

  end subroutine Cpt_RVM_q

  subroutine direct_potential_model_pts(crd_in,org_crd,org_q, model_crd,model_q)
    USE cosmo_helper
    use cosmo_util

    ! solves:  model_q = (BT)^(-1) * C * org_q
    REAL(SP), intent(IN) :: crd_in(:,:),org_crd(:,:),org_q(:),model_crd(:,:)
    REAL(SP), intent(OUT) :: model_q(:)
    REAL(SP), pointer :: phi_in(:)=>NULL(), BT(:,:) => NULL()
    INTEGER(I4B) :: n_in

    
    n_in = size(Crd_in,2)

    ! Compute C*org_q
    CALL cpt_phirhosigma_for_pt(org_crd,org_q,Crd_in,spread(0.0_SP,1,n_in),phi_in,0)
    ! Compute B matrix using point charges (no gaussian zetas)
    CALL cpt_B_for_pt(model_crd,Crd_in,spread(0.0_SP,1,n_in),BT,0)
    CALL svdinv_times_vector(BT,phi_in,model_q,.true.)
    deallocate(BT,phi_in)

  end subroutine direct_potential_model_pts

  subroutine direct_potential_model_surf(crd_in,org_crd,org_q, model_crd, model_zeta, model_q)
    USE cosmo_helper
    use cosmo_util

    ! solves:  model_q = (BT)^(-1) * C * org_q
    REAL(SP), intent(IN) :: crd_in(:,:),org_crd(:,:),org_q(:),model_crd(:,:),model_zeta(:)
    REAL(SP), intent(OUT) :: model_q(:)
    REAL(SP), pointer :: phi_in(:)=>NULL(), B(:,:) => NULL()
    REAL(SP), allocatable :: BT(:,:)
    INTEGER(I4B) :: n_in

    n_in = size(Crd_in,2)

    ! Compute C*org_q
    CALL cpt_phirhosigma_for_pt(org_crd,org_q,Crd_in,spread(0.0_SP,1,n_in),phi_in,0)
    ! Compute B matrix using point charges (no gaussian zetas)
    CALL cpt_B_for_pt(Crd_in,model_crd,model_zeta,B,2)
    allocate(BT(n_in,SIZE(model_zeta)))
    BT=transpose(B)
    deallocate(B)
    CALL svdinv_times_vector(BT,phi_in,model_q,.true.)
    deallocate(BT,phi_in)

  end subroutine direct_potential_model_surf


  function set_r_omega(rgamma,nse) result (romega)
    REAL(SP), INTENT(IN) :: rgamma
    INTEGER(I4B), intent(in) :: nse
    REAL(SP) :: romega,rnse,fact
!    Original parameters... VEP Paper 1
!    REAL(SP), PARAMETER :: c0=178.459_SP, c1=17.9528_SP, c2=5.27423_SP, c3=0.011264_SP
!    New parameters... Charge Scaling (vep paper 2)
     REAL(SP), PARAMETER :: c0=187.075_SP, c1=17.33_SP, c2=4.32713_SP, c3=0.013537_SP
   
    rnse=REAL(nse,SP)
    WRITE(6,*)'Computing inner surface location for VEP-RVM'
    WRITE(6,*)'  Using nse_VEP, r_RVM:',nse,rgamma

   ! emperical formula to determine optimal location of inside surface
   ! fact= 1.0_SP - 1.0_SP/SQRT((rnse/450.0_SP) + 1.0_SP)

    fact = 1.0_SP - ( (c1-c2*exp(-c3*rnse)) / SQRT(rnse + c0) )
    romega = rgamma * fact
    WRITE(6,*)'  Location of inner surface is r=',romega

  end function set_r_omega
END Module VEPMod

