module matrixlib
  use global
  implicit none

!------------------------------------------------------------------------------
! Constituents of a matrix object:
!   nv:    number of elements
!   np(1): number of rows
!   np(2): number of columns
!   tag:   name
!   data(np(1),np(2)): data
!------------------------------------------------------------------------------

  type matrix
     integer nv
     integer np(2)
     character (len=20) :: tag
     real(kind=8),pointer,dimension(:,:) :: data
  end type matrix

  type int_matrix
     integer np(2)
     character (len=20) :: tag
     integer, pointer, dimension(:,:) :: data
  end type int_matrix

contains

!------------------------------------------------------------------------------
! Basic matrix algebra
!------------------------------------------------------------------------------

!******************************************************************************

  integer function matrix_add (c, a, b)

! args
    type(matrix), pointer :: a,b
    type(matrix), pointer :: c

! begin
    c%data=a%data + b%data
    matrix_add = 0
  end function matrix_add

!******************************************************************************

  integer function matrix_transpose (a)

! args
    type(matrix), pointer :: a

! local vars
    integer i, idum, j, n
    real(kind=8) scr
    type(matrix), pointer :: b

! begin, row or column vector
    if ((a%np(1).eq.1) .or. (a%np(2).eq.1)) then
       idum = a%np(1)
       a%np(1) = a%np(2)
       a%np(2) = a%np(1)

! square matrix
    elseif (a%np(1).eq.a%np(2)) then
       n = a%np(1)
       do i = 1,n
          do j = i+1,n
             scr = a%data(j,i)
             a%data(j,i) = a%data(i,j)
             a%data(i,j) = scr
          end do
       end do

! generic case
    else
       b => matrix_create (a%np(2), a%np(1), a%tag)
       do i = 1,a%np(1)
          do j = 1,a%np(2)
             b%data(j,i) = a%data(i,j)
          end do
       end do
       idum = matrix_destroy (a)
       a => b
    end if
    matrix_transpose = 0
  end function matrix_transpose

!******************************************************************************

  integer function matrix_multiply (alpha, a, b, beta, c)

! args
    type(matrix), pointer :: a, b
    type(matrix), pointer :: c
    real(kind=8) alpha, beta

! local vars
    integer nlink, ixa, ixb, ixc, ii

! begin
    nlink = a%np(2)

    if (nlink .ne. b%np(1)) then
       write (stdout,'(a,i4,i4)') 'Mismatch dimensions - nlink', &
            a%np(2), b%np(1)
       call hdlc_errflag ('Matrix error', 'abort')
    end if

    if (a%np(1) .ne. c%np(1)) then
       write (stdout,'(a,i4,i4)') 'Mismatch dimensions - nrow(a) nrow(c)', &
            a%np(1), c%np(1)
       call hdlc_errflag ('Matrix error', 'abort')
    end if

    if (b%np(2) .ne. c%np(2)) then
       write (stdout,'(a,i4,i4)') 'Mismatch dimensions - ncol(b) ncol(c)', &
            b%np(2), c%np(2)
       call hdlc_errflag ('Matrix error', 'abort')
    end if

    call dgemm ('n', 'n', c%np(1), c%np(2), nlink, alpha, a%data, a%np(1), &
         b%data, b%np(1), beta, c%data, c%np(1))

    matrix_multiply = 0
  end function matrix_multiply

!******************************************************************************

!//////////////////////////////////////////////////////////////////////////////
! Matrix diagonalisation
!
! Object oriented approach: the dimensions of the passed arrays determine the
! behaviour of the routine:
!
! n     = order of the problem (dimensions of a)
! nval  = number of eigenvalues required (leading dimension of evalues)
! nvect = number of eigenvectors required (trailing dimension of evect)
! 
! Condition: n >= nval >= nvect
! If nvect or nval are less than the order of a, the lowest or the highest
! eigenvalues are taken if increasing or decreasing, respectively
!//////////////////////////////////////////////////////////////////////////////

  integer function matrix_diagonalise (a, evect, evalues, increasing)

! args
    logical increasing
    type(matrix), pointer :: a, evect, evalues

! local vars
    integer n, nval, nvect

! begin, derive sizes from the matrix objects and bomb if inconsistent
    n = a%np(1)
    if (a%np(2) .ne. n) then
       write (stdout,'(a)') 'Matrix not square in matrix_diagonalise'
       write (stdout,'(a,a)') 'Name of the matrix: ', a%tag
       matrix_diagonalise = -1
       return
    end if
    nval = evalues%np(1)
    if (nval .gt. n) then
       write (stdout,'(a)') &
            'More eigenvalues requested than the order of the matrix'
       write (stdout,'(a,a,a)') 'Names of the matrices: ', a%tag, evalues%tag
       matrix_diagonalise = -1
       return
    end if
    nvect = evect%np(2)
    if (nvect .gt. nval) then
       write (stdout,'(a)') &
            'More eigenvectors requested than eigenvalues'
       write (stdout,'(a,a,a,a)') 'Names of the matrices: ', a%tag, &
            evalues%tag, evect%tag
       matrix_diagonalise = -1
       return
    end if

! the dirty work is done in array_diagonalise
    matrix_diagonalise = array_diagonalise (a%data, evect%data, evalues%data, &
         n, nval, nvect, increasing)
    return
  end function matrix_diagonalise

!******************************************************************************

! array_diagonalise provides a non-object-oriented entry

  integer function array_diagonalise (a, evect, evalues, n, nval, nvect, &
       increasing)

! args
    logical increasing
    integer n, nval, nvect
    real(kind=8), dimension(n,n) :: a
    real(kind=8), dimension(nval) :: evalues
    real(kind=8), dimension(n,nvect) :: evect

! externals
    real(kind=8) dlamch
    external dlamch

! local vars
    character jobz
    integer flag, i, ierr, il, iu, j, k, lwork, nfound
    integer, dimension(:), allocatable :: ifail, iwork
    real(kind=8) abstol, dummy, vwork
    real(kind=8), dimension(:), allocatable :: evalwork, work
    real(kind=8), dimension(:,:), allocatable :: awork, evecwork

! begin

! copy the matrix to a work array as it would be destroyed otherwise
    allocate (awork(n,n))
    do i = 1,n
       do j = 1,n
          awork(j,i) = a(j,i)
       end do
    end do

! full diagonalisation is required
    if (nval .eq. n) then
       lwork = 3*n ! this is rather a minimum: unblocked algorithm
       allocate (work(lwork)) 
       if (nvect .eq. 0) then
          jobz = 'N'
       else
          jobz = 'V'
       end if

! diagonaliser DSYEV and error check
       call dsyev (jobz, 'L', n, awork, n, evalues, work, lwork, ierr)
       if (ierr.ne.0 .and. ctrl%printl.ge.1) then
          write (stdout,'(A,I5)') &
               'Matrix diagonaliser DSYEV failed: returned ', ierr
       end if

! clean up
       deallocate (work)
       if (increasing) then
          do i = 1,nvect
             do j = 1,n
                evect(j,i) = awork(j,i)
             end do
          end do
       else
          do i = 1,nvect
             do j = 1,n
                evect(j,i) = awork(j,n-i+1)
             end do
          end do
          k = n
          do i = 1,n/2
             vwork = evalues(i)
             evalues(i) = evalues(k)
             evalues(k) = vwork
             k = k - 1
          end do
       end if

! partial diagonalisation is required
    else
       lwork = 8*n ! this is rather a minimum: unblocked algorithm
       abstol = 2.0_8 * dlamch('S') ! this is for maximum accuracy
       allocate (work(lwork))
       allocate (iwork(5*n))
       allocate (evalwork(n))
       allocate (evecwork(n,nval)) ! note that this may be larger than nvect
       allocate (ifail(n))
       if (nvect .eq. 0) then
          jobz = 'N'
       else
          jobz = 'V'
       end if
       if (increasing) then
          il = 1
          iu = nval
       else
          il = n - nval + 1
          iu = n
       end if

! diagonaliser DSYEVX and error check
       call dsyevx (jobz, 'I', 'L', n, awork, n, dummy, dummy, il, iu, &
            abstol, nfound, evalwork, evecwork, n, work, lwork, iwork, &
            ifail, ierr)
       if (ierr.ne.0 .and. ctrl%printl.ge.1) then
          write (stdout,'(A,I5)') &
               'Matrix diagonaliser DSYEVX failed: returned ', ierr
          if (ctrl%printl.ge.3) then
             write (stdout,'(A)') 'Detailed error message (IFAIL):'
             write (stdout,'(16I5)') (ifail(i), i=1,n)
          end if
       end if

! clean up
       deallocate (iwork)
       deallocate (work)
       if (increasing) then
          do i = 1,nvect
             do j = 1,n
                evect(j,i) = evecwork(j,i)
             end do
          end do
          do i = 1,nval
             evalues(i) = evalwork(i)
          end do
       else
          do i = 1,nvect
             do j = 1,n
                evect(j,i) = evecwork(j,nval-i+1)
             end do
          end do
          do i = 1,nval
             evalues(i) = evalwork(nval-i+1)
          end do
       end if
       deallocate (ifail)
       deallocate (evecwork)
       deallocate (evalwork)
    end if

! clear working space
    deallocate (awork)
    array_diagonalise = ierr
    return
  end function array_diagonalise

!******************************************************************************

!//////////////////////////////////////////////////////////////////////////////
! Matrix inversion
!
! Object oriented entry: matrix_invert (a, det, lcdet)
! Array entry:           array_invert  (a, det, lcdet, n)
!
! Leading dimension assumed equal the order of the matrix to be inverted
! Calculate the determinant of lcdet is set
!//////////////////////////////////////////////////////////////////////////////

  integer function matrix_invert (a, det, lcdet)

! args
    logical lcdet
    real(kind=8) det
    type(matrix), pointer :: a

! local vars
    integer n

! begin
    n = a%np(1)
    if (a%np(2).ne.n) then
       call hdlc_errflag ('Matrix not square in invert', 'abort')
       matrix_invert = -1
       return
    else
       matrix_invert = array_invert (a%data, det, lcdet, n)
    end if
  end function matrix_invert

!******************************************************************************

  integer function array_invert (a, det, lcdet, n)

! args
    logical lcdet
    integer n
    real(kind=8) det
    real(kind=8), dimension(n,n) :: a
    
! local vars
    integer info, job, r
    real(kind=8), dimension(2) :: dd
    integer, dimension(:), allocatable :: ipvt
    real(kind=8), dimension(:), allocatable :: work

! begin
    r = n
    if (lcdet) then
       job = 11
    else
       job = 01
    end if

! allocate memory
    allocate (ipvt(n))
    allocate (work(n))

! get factors, replaces LINPACK: call dgefa (a, r, n, ipvt, info)
    call dgetrf (n, n, a, r, ipvt, info)

! test for singularity
    if (info.ne.0) then
       write (stdout,'(a)') &
            'Warning: attempt to invert a (probably) singular matrix'
       write (stdout,'(a)') &
            'Matrix is left unchanged and determinant is set to zero'
       if (lcdet) then
          dd(1) = 0.0_8
          dd(2) = 0.0_8
          det = dd(1) * 10.0_8**dd(2)
       end if
       array_invert = info
       return
    end if

! determinant, replaces LINPACK: call dgedi (a, r, n, ipvt, dd, work, job)
    if (lcdet) call dgedet (a, r, n, ipvt, dd, job)
    call dgetri (n, a, r, ipvt, work, n, info)

! deallocate work space
    deallocate (work)
    deallocate (ipvt)

! recompute determinant
    if (lcdet) det = dd(1) * 10.0_8**dd(2)
    array_invert = info
  end function array_invert

!******************************************************************************

  integer function matrix_scale (a, fact)

! args
    type(matrix), pointer :: a
    real(kind=8) fact

! begin
    a%data = fact * a%data
    matrix_scale = 0

  end function matrix_scale

!******************************************************************************

  function matrix_absmax (a)
    real(kind=8) matrix_absmax

! args
    type(matrix), pointer :: a

! local vars
    integer i, j
    real(kind=8) t

! begin
    t = 0.0_8
    do j = 1,a%np(2)
       do i = 1,a%np(1)
          t = max (t, abs(a%data(i,j)))
       end do
    end do
    matrix_absmax  = t

  end function matrix_absmax
    
!******************************************************************************

  function matrix_length (a)
    real(kind=8) matrix_length

! args
    type(matrix), pointer :: a

! local vars
    integer i, j
    real(kind=8) element, sum

! begin
    sum = 0.0_8
    do j = 1,a%np(2)
       do i = 1,a%np(1)
          element = a%data(i,j)
          sum = sum + element*element
       end do
    end do
!   fac = 1.0_8 / (a%np(1)*a%np(2))
    matrix_length = sqrt(sum)

  end function matrix_length

!------------------------------------------------------------------------------
! Supplements to library functions
!------------------------------------------------------------------------------

!//////////////////////////////////////////////////////////////////////////////
! calculate the determinant of a matrix using the factors computed by DGETRF
!
! this subroutine is cut out from LINPACK: SUBROUTINE DGEDI
!
! on entry:
!   a       the output from DGETRF (LAPACK ) / DGEFA (LINPACK).
!   lda     the leading dimension of the matrix a.
!   n       the order of the matrix a.
!   ipvt    the pivot vector from DGETRF / DGEFA.
!
! on return:
!   a       unchanged.
!   ipvt    unchanged.
!   det     determinant of original matrix. determinant = det(1) * 10.0**det(2)
!           with  1.0 .le. abs(det(1)) .lt. 10.0  or  det(1) .eq. 0.0 .
!//////////////////////////////////////////////////////////////////////////////

  subroutine dgedet (a, lda, n, ipvt, det, job)

! args
    integer lda, n, ipvt(*), job
    real(kind=8) a(lda,*), det(*)

! local vars
    integer i
    real(kind=8) one, ten, zero
 
! begin
    if (job .lt. 10) return
    zero = 0.0_8
    one = 1.0_8
    ten = 10.0_8
    det(1) = one
    det(2) = zero 

! loop over diagonal elements of a
    do i = 1,n
       if (ipvt(i) .ne. i) det(1) = -det(1)
       det(1) = a(i,i)*det(1)

! exit if product is zero
       if (det(1) .eq. zero) return

! cast result into the form determinant = det(1) * ten**det(2)
       do while (abs(det(1)) .lt. one)
          if (det(2).lt.-1000.0_8) then
             call hdlc_errflag ('Problems getting determinant', 'warn')
             det(1) = zero 
             det(2) = zero 
             return 
          end if
          det(1) = ten*det(1)
          det(2) = det(2) - one
       end do
       do while (abs(det(1)) .ge. ten)
          det(1) = det(1)/ten
          det(2) = det(2) + one
       end do
    end do
  end subroutine dgedet

!------------------------------------------------------------------------------
! Creation and destruction etc.
!------------------------------------------------------------------------------

  function int_matrix_create (nr, nc, name)
    type(int_matrix), pointer :: int_matrix_create
    integer, intent(in) :: nr, nc
    character*(*) name
    type(int_matrix), pointer :: temp
    allocate (temp)
    temp%np(1) = nr
    temp%np(2) = nc
    temp%tag = name
    allocate (temp%data(temp%np(1),temp%np(2)))
    int_matrix_create => temp
  end function int_matrix_create

  function int_matrix_destroy (a)
    integer int_matrix_destroy
    type(int_matrix), pointer :: a
    deallocate (a%data)
    deallocate (a)
    int_matrix_destroy = 0
  end function int_matrix_destroy

  function int_matrix_dimension (a, i)
    type(int_matrix), pointer :: a
    integer int_matrix_dimension, i
    int_matrix_dimension = a%np(i)
  end function int_matrix_dimension

  function matrix_dimension (a, i)
    type(matrix), pointer :: a
    integer matrix_dimension, i
    matrix_dimension = a%np(i)
  end function matrix_dimension

  function matrix_data_array (a)
    type(matrix), pointer :: a
    real(kind=8), dimension(:,:), pointer :: matrix_data_array
    matrix_data_array => a%data
  end function matrix_data_array

  function matrix_create(nr,nc,name)
    type(matrix), pointer :: matrix_create
    integer, intent(in) :: nr,nc
    character*(*) name
    type(matrix), pointer:: temp
    allocate(temp)
    temp%nv = nr*nc
    temp%np(1) = nr
    temp%np(2) = nc
    temp%tag = name
    allocate (temp%data(temp%np(1),temp%np(2)))
    matrix_create => temp 
  end function matrix_create

  integer function matrix_set_column (a, data, col)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer col, i
    do i = 1,a%np(1)
       a%data(i,col) = data(i)
    end do
    matrix_set_column = 0
  end function matrix_set_column

  integer function matrix_set_row (a, data, row)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer row, i
    do i = 1,a%np(2)
       a%data(row,i) = data(i)
    end do
    matrix_set_row = 0
  end function matrix_set_row

  integer function matrix_set (a, data)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer i, j, k
    k = 0
    do j = 1,a%np(2)
       do i = 1,a%np(1)
          k = k + 1
          a%data(i,j) = data(k)
       end do
    end do
    matrix_set = 0
  end function matrix_set

  integer function matrix_get (a, data)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer i, j, k
    k = 0
    do j = 1,a%np(2)
       do i = 1,a%np(1)
          k = k + 1
          data(k) = a%data(i,j)
       end do
    end do
    matrix_get = 0
  end function matrix_get

  integer function matrix_get_column (a, data, col)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer col, i
    do i = 1,a%np(1)
       data(i) = a%data(i,col)
    end do
    matrix_get_column = 0
  end function matrix_get_column

  integer function matrix_get_row (a, data, row)
    type(matrix), pointer :: a
    real(kind=8), dimension(*) :: data
    integer row, i
    do i = 1,a%np(2)
       data(i) = a%data(row,i)
    end do
    matrix_get_row = 0
  end function matrix_get_row

  integer function int_matrix_set_column (a, data, col)
    type(int_matrix), pointer :: a
    integer, dimension(*) :: data
    integer col, i
    do i = 1,a%np(1)
       a%data(i,col) = data(i)
    end do
    int_matrix_set_column = 0
  end function int_matrix_set_column

  integer function int_matrix_set_element (a, data, row, col)
    type(int_matrix), pointer :: a
    integer data, col, row
       a%data(row,col) = data
    int_matrix_set_element = 0
  end function int_matrix_set_element

  integer function matrix_assign_unit(a)
    type(matrix), pointer :: a
    integer i
    a%data=0.0_8
    do i = 1,a%np(1)
       a%data(i,i) = 1.0_8
    enddo
    matrix_assign_unit = 0
  end function matrix_assign_unit

  integer function matrix_copy (x, y)
    type(matrix), pointer :: x, y
    integer length
    length = x%np(1) * x%np(2)
    call dcopy (length, x%data, 1, y%data, 1)
    matrix_copy = 0
  end function matrix_copy

  integer function matrix_print(a)
    type(matrix), pointer  :: a
!
    logical ohi
!
    write (stdout,*) ' '
    write (stdout,'(a,a)')'Contents of matrix: ',a%tag
    write (stdout,*)   '---------------------------------------'
!
    ohi = .false.
    call prmat(a%data,a%np(1),a%np(2), a%np(1), ohi )
!
    matrix_print = 0
!
  contains
    subroutine prmat(v,m,n,ndim,ohi)
! m = number of columns
! n = number of rows
! ndim = leading dimension
! ohi print high precision
      integer n, m, ndim
      real(kind=8) v(ndim,*)
      logical ohi
      integer imin, imax, max, i, j, iwr
      ohi = .false.
!
      max = 12
      if (ohi) max=7
      imax = 0
!
100   imin = imax+1
      imax = imax+max
      if (imax .gt. n) imax = n
      write (stdout,9008)
      if (.not. ohi) write (stdout,8028) (i,i = imin,imax)
      if (ohi) write (stdout,9028) (i,i = imin,imax)
      write (stdout,9008)
      do j = 1,m
         if (ohi) write (stdout,9048) j,(v(j,i),i = imin,imax)
         if (.not.ohi) write (stdout,8048) j,(v(j,i),i = imin,imax)
      end do
      if (imax .lt. n) goto 100
      return
9008  format(1x)
9028  format(6x,7(6x,i3,6x))
9048  format(i5,1x,7f15.10)
8028  format(6x,12(3x,i3,3x))
8048  format(i5,1x,12f9.5)
    end subroutine prmat
  end function matrix_print

  integer function matrix_destroy(a)
    type(matrix), pointer :: a
    deallocate(a%data)
    deallocate(a)
    matrix_destroy = 0
  end function matrix_destroy

!------------------------------------------------------------------------------
! Checkpointing: dump to / undump from unit iunit
!------------------------------------------------------------------------------

  subroutine hdlc_wr_matrix (iunit, a, lform, lerr)
    logical lform, lerr
    integer iunit
    type(matrix), pointer :: a
    integer i,j
    lerr = .false.
    if (lform) then
       write (iunit,'(3I5)',err=98) a%nv, a%np(1), a%np(2)
       write (iunit,'(a)',err=98) a%tag
       write (iunit,'(5E20.12)',err=98) &
            ((a%data(i,j), i=1,a%np(1)), j=1,a%np(2))
    else
       write (iunit,err=98) a%nv, a%np(1), a%np(2)
       write (iunit,err=98) a%tag
       write (iunit,err=98) &
            ((a%data(i,j), i=1,a%np(1)), j=1,a%np(2))
    end if
    return
98  lerr = .true.
  end subroutine hdlc_wr_matrix

  subroutine hdlc_rd_matrix (iunit, a, lform, lerr)
    logical lform, lerr
    integer iunit
    type(matrix), pointer :: a
    integer i, j
    lerr = .false.
    if (lform) then
       read (iunit,'(3I5)',err=98) a%nv, a%np(1), a%np(2)
       read (iunit,'(a)',err=98) a%tag
    else
       read (iunit,err=98) a%nv, a%np(1), a%np(2)
       read (iunit,err=98) a%tag
    end if
    allocate (a%data (a%np(1), a%np(2)))
    if (lform) then
       read (iunit,'(5E20.12)',err=98) &
            ((a%data(i,j), i=1,a%np(1)), j=1,a%np(2))
    else
       read (iunit,err=98) &
            ((a%data(i,j), i=1,a%np(1)), j=1,a%np(2))
    end if
    return
98  lerr = .true.
  end subroutine hdlc_rd_matrix

end module matrixlib
