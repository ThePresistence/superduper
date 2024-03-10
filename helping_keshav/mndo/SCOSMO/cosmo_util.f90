module cosmo_util
  use DataTypes
  IMPLICIT NONE

  ! Things I wrote for cosmo that probaly should be in YorkLib.UtilitiesMod
  ! or somewhere else

  private
  ! default is machine precision
  REAL(SP),parameter :: svd_tol_default = epsilon(1.0_SP)

  public :: do_syminv, &
            diagmul, diagadd, fill_upper, &
            unpack_sym_array, pack_sym_array, &
            cosmo_f, CheckMultipoles, &
            cpt_constraints_COSMO,&
            svdinv_times_vector, svdinv_matrix


  interface diagmul
     module procedure diagmul_l, diagmul_r
  end interface

contains

  function diagadd(M, r)
    ! add r to each element of M's diagonal
    ! Only works for square matricies
    use UtilitiesMod, only : assert_eq
    IMPLICIT NONE
    real(sp) :: M(:,:), r
    real(sp) :: diagadd(1:size(M,1),1:size(M,2))
    integer(i4b) :: i, n

    n = assert_eq(size(M,1), size(M,2),&
         'Not square in diagadd.  Please fix me.')

    diagadd(1:n,1:n) = M
    
    forall(i = 1:n)
          diagadd(i,i) =  M(i,i) + r
    end forall

  end function diagadd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    
  function diagmul_l(diag, L)
    IMPLICIT NONE
    ! multiply by a diagonal matrix from the left
    real(sp) :: diag(:), L(:,:)
    real(sp) :: diagmul_l(1:size(L,1),1:size(L,2))
    integer(i4b) :: i, j, n

    n = size(diag)

    forall(i = 1:n)
       forall(j = 1:n)
          diagmul_l(j,i) = diag(j) * L(j,i)
       end forall
    end forall

  end function diagmul_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function diagmul_r(M, diag)
    IMPLICIT NONE
    ! multiply by a diagonal matrix from the right
    real(sp) :: diag(:), M(:,:)
    real(sp) :: diagmul_r(1:size(M,1),1:size(M,2))
    integer(i4b) :: i, j, n

    n = size(diag)

    forall(i = 1:n)
       forall(j = 1:n)
          diagmul_r(j,i) = diag(i) * M(j,i)
       end forall
    end forall

  end function diagmul_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine do_syminv(n,inv,verbose)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN):: n
    REAL(SP), INTENT(INOUT) :: inv(:,:)
    LOGICAL(LGD), INTENT(IN) :: verbose

    integer(i4b) :: info

    IF (verbose)write(6,*)'Cholesky decomp using LAPACK...'
    call dpotrf('l',n,inv,n,info)
    IF (verbose)write(6,*)'Constructing inverse using LAPACK...'
    call dpotri('l',n,inv,n,info)
    ! fill upper triangle of inverse
    call fill_upper(inv)
  end subroutine do_syminv

  subroutine do_syminv_packed(n,inv,verbose)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN):: n
    REAL(SP), INTENT(INOUT) :: inv(:)
    LOGICAL(LGD), INTENT(IN) :: verbose
    integer(i4b) :: info

    IF (verbose)write(6,*)'Cholesky decomp using LAPACK with packed array...'
    call dpptrf('l',n,inv,info)
    IF (verbose)write(6,*)'Constructing inverse using LAPACK with packed array...'
    call dpptri('l',n,inv,info)
  end subroutine do_syminv_packed

  subroutine make_cholinv(A, d)
    ! Take a matrix munged by choldc (and its associated vector) and creates
    ! a big old regular inverse matrix out of it.  This destroys A.
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(sp) :: A(:,:)
    real(sp) :: d(:)
    integer(i4b) :: i, j, n

    n = assert_eq(size(A,1),size(A,2),size(d),'sizes wrong in make_cholinv')

    ! NR pg. 91
    ! construct Linv lower triangle
    do i = 1, n
       A(i,i) = 1.0_sp / d(i)
       do j = i+1, n
          A(j,i) = -dot_product(A(j,i:j-1),A(i:j-1,i))/d(j)
       end do
    end do

    ! zero out junk in upper triangle
    forall (i = 1:n)
       forall (j = i+1:n)
          A(i,j) = 0.0_sp
       end forall
    end forall

    ! Ainv = transpose(Linv) times Linv
    ! A = matmul(transpose(A),A)
    A = fill_from_lower(A)

  end subroutine make_cholinv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function fill_from_lower(A) result(B)
    ! takes a lower triangular matrix and returns matmul(transpose(A),A)
    use UtilitiesMod, only : assert_eq
    IMPLICIT NONE
    real(sp) :: A(:,:)
    real(sp) :: B(1:size(A,1),1:size(A,2))
    integer(i4b) :: i, j, n

    B = 0.0_sp

    n = assert_eq(size(A,1),size(A,2),'A must be square')

    forall (i=1:n)
       forall (j=1:i)
          B(i,j) = dot_product(A(i:n,j),A(i:n,i))
          B(j,i) = B(i,j)
       end forall
    end forall

  end function fill_from_lower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_inv_bysl(A, d, Ainv)
    ! Take a matrix munged by choldc (and its associated vector) and creates
    ! a big old regular inverse matrix out of it in Ainv
    use EigenMod, only : cholsl
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(sp), pointer :: A(:,:)
    real(sp), pointer :: d(:)
    real(sp), pointer :: Ainv(:,:)
    integer(i4b) :: i, n
    real(sp), allocatable :: b(:)

    n = assert_eq(size(A,1),size(Ainv,2),size(d),'sizes wrong in make_inv_bysl')
  
    ! The is probably a smarter way to do this (like pg. 91 in NR)
    ! But I am kind slow, so:
    allocate(b(1:n))

    do i = 1, n
       ! create a column of the identity matrix
       b = 0.0_sp
       b(i) = 1.0_sp
       call cholsl(A, d, b, Ainv(:,i))
    end do

    deallocate(b)

  end subroutine make_inv_bysl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

  function VecLen_Packed_SymArray(nArray) result (nVector)
    IMPLICIT NONE
    integer(i4b), intent(in) :: nArray
    integer(i4b) :: nVector

     nVector =nArray*(nArray+1)/2
  end function VecLen_Packed_SymArray


  function SymArrayDim_unpacked_vector(nvector) result(nArray)
    IMPLICIT NONE
    integer(i4b), intent(in) :: nvector
    integer(i4b) :: nArray

  ! recover array dimensions with a little quadratic formula action...
   nArray = (-1 + sqrt(real(1+8*nvector)))/2

  end function SymArrayDim_unpacked_vector

  subroutine pack_sym_array(array, vector)
    ! Takes a symmetric array, and packs it into a vector in lower packed 
    ! format
    use ErrorMod
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(sp), INTENT(IN) :: array(:,:)
    real(sp), INTENT(OUT) :: vector(:)
    integer(i4b) :: n, i, j,v

    n = assert_eq(size(array,1),size(array,2),&
         "If it ain't square, it definitely ain't symmetric ;)")
   
    v = VecLen_Packed_SymArray(n)
    v = assert_eq(v,SIZE(vector), 'vector is not correct size') 

    ! It is faster to do something like:
    !    k = 0 
    !      do j = 1, n
    !         do i = j, n
    !            vector(k) = array(i,j)
    !         end do
    !      end do
    ! But that is not parallelizable and even the NAG compiler can optimize
    ! the loop index math quite well.

    forall(j=1:n)
       forall(i=j:n)
          vector((i + ((2*n-j)*(j-1)/2))) = array(i,j)
       end forall
    end forall
    
  end subroutine pack_sym_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine unpack_sym_array(vector, array)
    use ErrorMod
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(sp), INTENT(IN) :: vector(:)
    real(sp), INTENT(OUT) :: array(:,:)
    integer(i4b) :: n, i, j

    n = SymArrayDim_unpacked_vector(SIZE(vector))

    n = assert_eq(size(array,1),size(array,2),n,&
       "If it ain't square, it definitely ain't symmetric ;)")

    forall(j=1:n)
       forall(i=j:n)
          array(i,j) = vector((i + ((2*n-j)*(j-1)/2)))
       end forall
    end forall

    call fill_upper(array)

  end subroutine unpack_sym_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine fill_upper(array)
    ! fills the upper triangle of a so called "lower storage mode" matrix
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(sp) :: array(:,:)
    integer(i4b) :: i, j, n

    n = assert_eq(size(array,1),size(array,2),&
         "If it ain't square, it definitely ain't symmetric ;)")
    
    forall(i=1:n)
       forall(j=i:n)
          array(i,j) = array(j,i)
       end forall
    end forall
    
  end subroutine fill_upper


       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cosmo_f(e1,e2)
    IMPLICIT NONE
    real(sp), intent(in) :: e1, e2
    real(sp) :: cosmo_f

    if(e1 <= 0 .OR. e2 <=0) THEN
       ! limit as e2 => inf
       write(6,*)'setting dielectric scale factor to 1 / e1'
       cosmo_f = 1.0_sp / e1
    else
       ! factor for everything else
       cosmo_f = (e2 - e1) / (e2 * e1)
    end if

  end function cosmo_f

subroutine CheckMultipoles(ctr,q1,crd1,q2,crd2,rxn_fld)
    USE StringMod, only : RemoveSpaces 
    USE MultiPoleMod
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: ctr(3),q1(:),crd1(:,:),q2(:),crd2(:,:)
    LOGICAL(LGD) :: rxn_fld
    REAL(SP), allocatable :: multi1(:),multi2(:)
    INTEGER(I4B) :: lmax,ldim,nq1,nq2,l,ltmp1,ltmp2
    CHARACTER(LEN=20) :: fmt

    nq1=SIZE(q1)
    nq2=size(q2)
    lmax = 4
    ldim=(lmax+1)**2

    allocate(multi1(ldim),multi2(ldim))
    IF(rxn_fld)THEN
       CALL qTimesNearField(q1,Crd1,ctr,lmax,.FALSE.,multi1)
       CALL qTimesNearField(q2,Crd2,ctr,lmax,.FALSE.,multi2)
       multi2 = -1.0_SP * multi2
    ELSE
       CALL qTimesFarField(q1,Crd1,ctr,lmax,multi1)
       CALL qTimesFarField(q2,Crd2,ctr,lmax,multi2)
    ENDIF

    WRITE(6,*)'^^^^^^^^^^^^ Summary of Multipoles ^^^^^^^^^^^^'
    Do l=0,lmax
       fmt=''
       ltmp1=(l*l)+1
       ltmp2=(l+1)**2
       WRITE(fmt,'(A,I4,A)')'(A,',ldim,'ES12.4)'
       fmt = RemoveSpaces(fmt)
       WRITE(6,'(A,I4)')'^ l = ',l
       WRITE(6,fmt)'| Charge Distribution 1 ',multi1(ltmp1:ltmp2 )
       WRITE(6,fmt)'| Charge Distribution 2 ',multi2(ltmp1:ltmp2 )
       WRITE(6,fmt)'| Difference (1-2)      ',multi1(ltmp1:ltmp2 ) - multi2(ltmp1:ltmp2 )
    enddo
    WRITE(6,*)'^^^^^^^^^^^^ End Multipole Summary ^^^^^^^^^^^^'

    deallocate(multi1,multi2)
  END subroutine CheckMultipoles

  subroutine cpt_constraints_COSMO(l,crds, f,Ainv, Z_rho, B_rho, D, lambda)
    USE MultipoleMod, only : NearField
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: crds(:,:),f,Ainv(:,:),Z_rho(:), B_rho(:)
    REAL(SP), pointer :: D(:,:), lambda(:)
    INTEGER(I4B) :: l
    real(sp), parameter :: ctr(3) = (/0.0_SP,0.0_SP,0.0_SP/)
    real(sp), pointer :: Qmat(:,:)=>NULL()
    real(sp),allocatable :: dt(:,:),R_rho(:)
    integer(i4b) :: ldim,nse

    nse=SIZE(crds,2)
    ldim=(l+1)**2

    allocate(dt(1:ldim,1:nse))
    ! NearField Expansion about origin
    CALL Nearfield(crds,ctr,l,dt)
    allocate(D(1:nse,1:ldim))
    D = TRANSPOSE(dt)
    allocate(Qmat(1:ldim,1:ldim)) 
    Qmat = matmul(Dt,f*matmul(Ainv,D))
    allocate(R_rho(ldim))
    R_rho =  Z_Rho + matmul(Dt,f*matmul(Ainv, B_rho))
    deallocate(dt)

    allocate(lambda(1:ldim))
!    call svdinv_times_vector(Qmat,R_rho,lambda)
    CALL do_syminv(ldim,Qmat,.FALSE.)
    lambda = matmul(Qmat,R_rho)
    deallocate(Qmat,R_rho)

  end subroutine cpt_constraints_COSMO


  subroutine svdinv_matrix(Mat,MatInv,verbose,Tol)
    IMPLICIT NONE
    REAL(SP),INTENT(IN)   :: Mat(:,:)
    REAL(SP), INTENT(OUT) :: MatInv(:,:)
    LOGICAL(LGD),INTENT(IN) :: verbose
    REAL(SP), optional, INTENT(IN) :: Tol
    REAL(SP), allocatable :: u(:,:),S(:),V(:,:)
    REAL(SP) :: my_tol
    REAL(SP), allocatable :: work(:)
    REAL(SP) :: tmpMat(1,1),gwork(1)
    INTEGER(I4B) :: lwork,info
    INTEGER(I4B) :: i,dim1,dim2,rank,imin

    IF(present(tol))THEN
       my_tol = MAX(tol,tiny(1.0_SP))
    ELSE
       my_tol = MAX(svd_Tol_default,tiny(1.0_SP))
    ENDIF

    dim1= SIZE(MAT,1)
    dim2= SIZE(MAT,2)
    imin=min(dim1, dim2)
 
    allocate (S(dim2))
    s=0.0_SP
    allocate (V(dim2,dim2))
    v=0.0_SP
    allocate (u(dim1,dim2))
    u=Mat

    CALL dgesvd('O',   'A',  dim1,dim2, U, dim1, S, tmpMat, 1,  V, dim2, gwork, -1, INFO )
    lwork=2*CEILING(gwork(1))
    if(verbose)WRITE(6,*)'Allocating dgesvd workspace: ', lwork
    allocate(work(lwork))
    !              JOBU,  JOBVT,  M,    N,  A, LDA,  S,   U,   LDU, VT, LDVT, WORK, LWORK, INFO
    CALL dgesvd('O',   'A',  dim1,dim2, U, dim1, S, tmpMat, 1,  V, dim2, work, lwork, INFO )
    IF(INFO /= 0) WRITE(6,*)'Warning: svd failed  to  converge'
    deallocate(work)
    v=transpose(v)

    rank=1
    do i=2,dim2
       IF(s(i) > my_tol*s(1))THEN
         rank = rank + 1
       ELSE
         s(i) = 0.0_SP      
       ENDIF
    ENDDO

    if(verbose)CALL singular_value_info(s,rank,imin,my_tol)

    DO i=1,dim2
       IF(s(i) /= 0.0_SP)THEN
         u(:,i) = u(:,i)/s(i)
       ELSE
         u(:,i) = 0.0_SP
       ENDIF
    END DO
    MatInv(:,:) = matmul(V(:,:),transpose(u(:,:)))

    deallocate(S,V,u)

  end subroutine svdinv_matrix

  
  subroutine svdinv_times_vector(A,b,x,verbose,tol)
!!!****f* svdinv_times_vector
!!!
!!! NAME
!!!     svdinv_times_vector -- solves  A * x = b for vector x 
!!! USAGE
!!!     svdinv_times_vector(A,b,x)
!!! DESCRIPTION
!!!     calculates the singular value decomposition of A, forms A-1 and
!!!     returns x = A-1 * b
!!! NOTE
!!!     A is destroyed in the process
!!!***
    IMPLICIT NONE
    REAL(SP), INTENT(INOUT) :: A(:,:)
    REAL(SP), INTENT(IN) :: B(:)
    REAL(SP), INTENT(OUT) :: x(:)
    LOGICAL(LGD),INTENT(IN) :: verbose
    REAL(SP), optional, INTENT(IN) :: tol
    INTEGER(I4B) :: dim1,dim2,rank
    REAL(SP), allocatable :: S(:)
    REAL(SP) :: my_tol
    REAL(SP), allocatable :: work(:),vec_tmp(:)
    INTEGER(I4B), allocatable :: iwork(:)
    REAL(SP) :: gwork(1)
    INTEGER(I4B) :: info,lwork,LDB, imin,liwork,tiwork(1)

    ! Tol shouldnt be smaller than tiny
    IF(present(tol))THEN
       my_tol = MAX(tol,tiny(1.0_SP))
    ELSE
       my_tol = MAX(svd_Tol_default,tiny(1.0_SP))
    ENDIF

    dim1=SIZE(A,1)
    dim2=SIZE(A,2)
    imin=min(dim1, dim2)
    allocate(S(dim2))
    s(:)=0.0_SP

    LDB=max(1,max(dim1,dim2))
    allocate(vec_tmp(LDB))
    vec_tmp=0.0_SP
    vec_tmp(1:dim1) = B(1:dim1) 
    
    IF(dim1*dim2 < 100000)THEN
       if(verbose)WRITE(6,*)'Using direct svd'
       CALL dgelss(dim1, dim2, 1, A, dim1, vec_tmp, LDB, S, my_TOL, RANK, gwork, -1, INFO )
       lwork= 2*CEILING(gwork(1))
       if(verbose)WRITE(6,*)'Using workspace size:', lwork
       allocate(work(lwork))
       ! TOL=RCOND: determines effective rank of matrix, values S(i) <= RCOND*S(1) are treated as zero
       CALL dgelss(dim1, dim2, 1, A, dim1, vec_tmp, LDB, S, my_TOL, RANK, work, lwork, INFO )
    ELSE
       if(verbose)WRITE(6,*)'Using Divide and Conquer svd'
       CALL dgelsd(dim1, dim2, 1, A, dim1, vec_tmp, LDB, S, my_TOL, RANK, gwork, -1, tiwork, INFO )
       lwork= 2*CEILING(gwork(1))
       liwork = MAX( 0, INT( LOG( REAL(imin,SP)/(11.0_SP))/LOG(2.0_SP) ) + 1 )
       liwork = 2*(3 * imin * liwork + 11 * imin)
       if(verbose)WRITE(6,*)'Using workspace size:', lwork,liwork
       allocate(work(lwork))
       allocate(iwork(liwork))
       ! TOL=RCOND: determines effective rank of matrix, values S(i) <= RCOND*S(1) are treated as zero
       CALL dgelsd(dim1, dim2, 1, A, dim1, vec_tmp, LDB, S, my_TOL, RANK, work, lwork,iwork, INFO )
       deallocate(iwork)
    ENDIF
    IF(INFO /= 0) WRITE(6,*)'Warning: svd failed  to  converge'
    x(1:dim2) = vec_tmp(1:dim2)
    deallocate(vec_tmp,work)

    if(verbose)CALL singular_value_info(s,rank,imin,my_tol)

    deallocate(s)

  end subroutine svdinv_times_vector

  Subroutine singular_value_info(s,rank,mindim,tol)
    USE Datatypes
    IMPLICIT NONE
    REAL(SP),INTENT(IN) :: S(:)
    INTEGER(I4B),INTENT(IN) :: rank,mindim
    REAL(SP),INTENT(IN) :: tol
    INTEGER(I4B):: dim,i
    REAL(SP) :: rcond,esum,tot,val

    dim = SIZE(s)
    IF(mindim<dim)THEN
      WRITE(6,*)'SVD INFO, Matrix is Rank Deficient'
    ENDIF
    ! the condition number of A in the 2-norm = S(1)/S(min(m,n)).
    WRITE(6,*)'SVD INFO, Tol is : ', tol
    WRITE(6,*)'SVD INFO, Rank is: ', rank, '(',dim,mindim,')'
    esum = SUM(S(1:rank))
    open(90,FILE='SVD_EigenSpectrum.dat')
    tot=0.0_SP
    DO i=1,rank
       val=(S(i)/esum)*100_SP
       tot = tot + val
       WRITE(90,'(I5,ES12.4,2F10.3)')i, S(i), val, tot
    ENDDO
    close(90)
    WRITE(6,*)'SVD INFO, EigenSum: ', esum
    rcond=S(mindim)/S(1)
    IF(rcond == 0.0_SP) THEN
       WRITE(6,'(A,ES15.7,A)')'SVD INFO, Condition Number (2-norm): UNDEF[Singular Matrix] (',S(1)/S(rank),')'
    ELSE
       WRITE(6,'(A,ES15.7)')'SVD INFO, Condition Number (2-norm): ', 1.0_SP/rcond
       IF(rcond < 1.0E-12)THEN
          WRITE(6,'(A,ES15.7)')'SVD INFO, Warning: Matrix is Ill Conditioned. rcond= ', rcond
       ENDIF
    END IF

  end subroutine singular_value_info


end module cosmo_util
