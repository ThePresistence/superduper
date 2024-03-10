module cosmo_helper

  ! Routines that users of cosmo might find useful
  use Datatypes
  use errormod
  implicit none

  private
  
  public :: cpt_B_for_pt, &
            cpt_phirhosigma_for_pt, &
            cpt_Z_for_pt, &
            cpt_born, &
            cpt_onsager, &
            cpt_gborn, &
            cpt_born2
contains
  
  subroutine cpt_B_for_pt(q_crds, se_crds,se_zeta, B, atom_surf_zeta_mode)
    USE ConstantsMod, only : SQRT2
    USE CoulombMod, only : GaussianInteraction
   ! computes B matrix for a set of point charges q_crds
    
    real(sp), intent(in) :: q_crds(:,:),se_crds(:,:),se_zeta(:)
    real(sp), pointer :: B(:,:)
    INTEGER(I4B), optional, intent(IN) :: atom_surf_zeta_mode
    INTEGER(I4B) :: pt_zeta_mode

    integer(i4b) :: nq, nse, q, se
    real(sp) :: r2, r,zeta

    if ( size(q_crds,1) /= 3 ) &
         call error('q_crds is weird in cpt_B_for_pt',i1=size(q_crds,1))

    ! Get matrix dimensions
    nq = size(q_crds,2)
    nse = size(se_zeta)

    if ( associated(B) ) deallocate(B)

    ! this is the default, treat atoms as pt charges, surface as gaussians
    pt_zeta_mode = 1
    if(present(atom_surf_zeta_mode)) pt_zeta_mode=atom_surf_zeta_mode
    IF(pt_zeta_mode < 0 .or. pt_zeta_mode > 2) pt_zeta_mode = 1

    ! allocate
    allocate(B(1:nse,1:nq))

    if (pt_zeta_mode == 1) THEN
       do q = 1, nq
          !$OMP PARALLEL DO
          do se = 1, nse
             r2 = (se_crds(1,se) - q_crds(1,q))**2 + &
                  (se_crds(2,se) - q_crds(2,q))**2 + &
                  (se_crds(3,se) - q_crds(3,q))**2
             r = sqrt(r2)
             CALL GaussianInteraction(se_zeta(se),r,B(se,q))
          end do
          !$OMP END PARALLEL DO
       end do
    ELSE if (pt_zeta_mode == 2) THEN
       WRITE(6,*)'Experimental cpt_B_for_pt: "pts have surface zetas"'
       do q = 1, nq
          !$OMP PARALLEL DO
          do se = 1, nse
             r2 = (se_crds(1,se) - q_crds(1,q))**2 + &
                  (se_crds(2,se) - q_crds(2,q))**2 + &
                  (se_crds(3,se) - q_crds(3,q))**2
             r = sqrt(r2)
             zeta = se_zeta(se)/SQRT2
             CALL GaussianInteraction(zeta,r,B(se,q))
          end do
          !$OMP END PARALLEL DO
       end do
    elseif (pt_zeta_mode == 0) THEN
    !! ** THIS IS REALLY STUPID... PEOPLE SHOULD NEVER USE THIS OPTION ***
       WRITE(6,*)'WARNING from cpt_B_for_pt: "surface dosent have zetas"'
       do q = 1, nq
          !$OMP PARALLEL DO
          do se = 1, nse
             r2 = (se_crds(1,se) - q_crds(1,q))**2 + &
                  (se_crds(2,se) - q_crds(2,q))**2 + &
                  (se_crds(3,se) - q_crds(3,q))**2
             r = sqrt(r2)
             IF(r==0)THEN
               write(6,*)'r is zero in cpt_B:', se_crds(:,se),q_crds(:,q)
             ENDIF
             B(se,q) =  1.0_SP / r
          end do
          !$OMP END PARALLEL DO
       end do
    endif

  end subroutine cpt_B_for_pt

  subroutine cpt_phirhosigma_for_pt(q_crds,q, se_crds,se_zeta, phi, atom_surf_zeta_mode)
    USE ConstantsMod, only : SQRT2
    USE CoulombMod, only : GaussianInteraction 
    real(sp), intent(in) :: q_crds(:,:),q(:),se_crds(:,:),se_zeta(:)
    real(sp), pointer :: phi(:)
    INTEGER(I4B), optional, intent(IN) :: atom_surf_zeta_mode
    INTEGER(I4B) :: pt_zeta_mode

    integer(i4b) :: nq, nse, iq, ise
    real(sp) :: r2, r,zeta,B

    if ( size(q_crds,1) /= 3 ) &
         call error('q_crds is weird in cpt_B_for_pt',i1=size(q_crds,1))

    ! Get matrix dimensions
    nq = size(q_crds,2)
    nse = size(se_zeta)

    if ( associated(phi) ) deallocate(phi)

    ! this is the default, treat atoms as pt charges, surface as gaussians
    pt_zeta_mode = 1
    if(present(atom_surf_zeta_mode)) pt_zeta_mode=atom_surf_zeta_mode
    IF(pt_zeta_mode < 0 .or. pt_zeta_mode > 2) pt_zeta_mode = 1

    ! allocate
    allocate(phi(1:nse))
    
    if (pt_zeta_mode == 1) THEN
       ! gauss-pt interaction for surface <=> atom
       ! this is the default
       do ise = 1, nse
          phi(ise)=0.0_SP
         !$OMP PARALLEL DO
          do iq = 1, nq
             r2 = (se_crds(1,ise) - q_crds(1,iq))**2 + &
                  (se_crds(2,ise) - q_crds(2,iq))**2 + &
                  (se_crds(3,ise) - q_crds(3,iq))**2
             r = sqrt(r2)
             CALL GaussianInteraction(se_zeta(ise), r, B)
             phi(ise) = phi(ise) + B * q(iq)
          end do
          !$OMP END PARALLEL DO
       end do
    ELSE if (pt_zeta_mode == 2) THEN
       ! gauss-gauss interaction for surface <=> atom
       ! atoms use same zeta as surface element its interacting with
       do ise = 1, nse
          phi(ise)=0.0_SP
          !$OMP PARALLEL DO
          do iq = 1, nq
             r2 = (se_crds(1,ise) - q_crds(1,iq))**2 + &
                  (se_crds(2,ise) - q_crds(2,iq))**2 + &
                  (se_crds(3,ise) - q_crds(3,iq))**2
             r = sqrt(r2)
             zeta = se_zeta(ise)/SQRT2
             CALL GaussianInteraction(zeta, r, B)
             phi(ise) = phi(ise) + B * q(iq)
          end do
          !$OMP END PARALLEL DO
       end do
    elseif (pt_zeta_mode == 0) THEN
       ! pt-pt interaction for surface <=> atom
       do ise = 1, nse
          phi(ise)=0.0_SP
          !$OMP PARALLEL DO
          do iq = 1, nq
             r2 = (se_crds(1,ise) - q_crds(1,iq))**2 + &
                  (se_crds(2,ise) - q_crds(2,iq))**2 + &
                  (se_crds(3,ise) - q_crds(3,iq))**2
             r = sqrt(r2)
             IF(r==0)THEN
               write(6,*)'r is zero in cpt_B:', se_crds(:,ise),q_crds(:,iq)
             ENDIF
             B =  1.0_SP / r
             phi(ise) = phi(ise) + B * q(iq)
          end do
          !$OMP END PARALLEL DO
       end do
    endif

  end subroutine cpt_phirhosigma_for_pt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_born(eps1, eps2, SA, q, energy)

    ! Compute born ion energy
    ! Needs 
    use constantsmod

    real(sp), intent(in) :: SA,eps1,eps2, q(:)
    real(sp), intent(out) :: energy

    real(sp) :: radius

    radius = sqrt(sa/FOUR_PI)

    if ( eps2 > 0 ) then
       energy =  -0.5_SP * ( (eps2 - eps1)/(eps2*eps1)) * (sum(q)**2 / radius)
    else
       energy =  -0.5_SP * (sum(q)**2 / radius)
    end if

  end subroutine cpt_born

  subroutine cpt_onsager(eps1, eps2, SA, q, coords, energy)

    ! Compute onsager dipole solvation energy
    
    use constantsmod
    use MolecularPropertiesMod

    real(sp), intent(in) :: SA, eps1, eps2, q(:), coords(:,:)
    real(sp), intent(out) :: energy

    real(sp) :: radius, dipole(1:3), d
    
    dipole = DipoleMoment(q,coords)
    d = dot_product(dipole,dipole)

    radius = sqrt(sa/FOUR_PI)

    if ( eps2 > 0 ) then
       energy = - ( ( eps2 - eps1 ) / &
            ( 2.0_sp * eps2 + eps1) ) * ( d / radius**3 )
    else
       energy = - ( 0.5_sp ) * ( d / radius**3 )
    end if
    
  end subroutine cpt_onsager

  subroutine cpt_GBorn(eps1, eps2, q, coords, radii, energy)
    real(sp), intent(in) :: eps1,eps2, q(:),coords(:,:), radii(:)
    real(sp), intent(out) :: energy
    real(sp) :: EGB,A2,r2,D,fgb,v(3)
    INTEGER(I4B) :: n,i,j

  n=size(q)
  EGB=0.0_SP
  do i=1,n
    do j=1,n
      v=(coords(1:3,i)-coords(1:3,j))
      r2=dot_product(v,v)
      A2=radii(i)*radii(j)
      IF(A2 == 0.0_SP)THEN
        fgb=SQRT(r2)
      ELSE
        D=r2/(4.0_SP*A2)
        fgb=SQRT(r2+a2*EXP(-D))
      ENDIF
      EGB=EGB+q(i)*q(j)/fgb
    enddo
  enddo
  energy=-0.5_SP*( (eps2 - eps1)/(eps2*eps1))*EGB

  end subroutine cpt_GBorn


  subroutine cpt_Born2(eps1, eps2, q, coords, radii, energy)
   IMPLICIT NONE
    real(sp), intent(in) :: eps1,eps2, q(:),coords(:,:), radii(:)
    real(sp), intent(out) :: energy
    real(sp) :: E,r,v(3)
    INTEGER(I4B) :: n,i,j

  n=size(q)
  E=0.0_SP
  do i=1,n
    E = E + 0.5_SP*q(i)*q(i)/radii(i)
    do j=i+1,n
      v = coords(:,i) - coords(:,j)
      r = sqrt(dot_product(v,v))
      E = E + q(i)*q(j)/r
    enddo
  enddo
  energy = -E * (eps2 - eps1)/(eps2*eps1)
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_Z_for_pt(coords, lmax, Zmat)
  USE MultipoleMod
    ! compute solute constraint matrix
    ! multipoles based on (irreducible) real solid harmonic expansion 
    ! matmul(Zmat,rho) should yield a vector of (lmax+1)**2 spherical-tensor multipoles
    ! NOTE: these are NOT cartesian multipoles, but are related

    real(sp), intent(in) :: coords(:,:)
    integer(i4b), intent(in) :: lmax
    real(sp), pointer :: Zmat(:,:)
    integer(i4b) :: nq, ldim
    real(sp) :: ctr(3) = (/0.0_SP,0.0_SP,0.0_SP/)

    if ( lmax < 0 ) return

    if ( associated(Zmat) ) then
       deallocate(Zmat)
    end if

    if ( size(coords,1) /= 3 ) &
         call error('coords are weird in cpt_Z_for_pt',i1=size(coords,1))

    nq = size(coords,2)
    ldim=(lmax+1)**2
    allocate(Zmat(1:ldim,1:nq))

    CALL NearField(coords,ctr,lmax,Zmat)
    ! for cosmo rxnfld, matmul(Zmat(1,:),rho(:)) = -q
    Zmat = -Zmat

  end subroutine cpt_Z_for_pt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cosmo_helper
    
