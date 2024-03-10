!
! Subroutines for printing matrices and vectors.
!
! Written by Axel Koslowski in 2000-2003 at MPI Muelheim.
!

Module matutils

  Private

  Public :: PrintMatrix       ! Print rectangular matrix.
  Public :: PrintUpper        ! Print upper triangle of symmetric square matrix.
  Public :: PrintPackedUpper  ! Print upper triangle of symmetric packed matrix.
  Public :: PrintPackedLower  ! Print lower triangle of symmetric packed matrix.
  Public :: PrintVector       ! Print vector.

Contains

  !
  ! This subroutine prints a rectangular matrix on unit NB6:
  !
  Subroutine PrintMatrix(A, NB6)
    Implicit None
    Double Precision, Dimension(:,:) :: A
    Integer                          :: NB6
    ! end of dummy parameters
    Integer                          :: i, j, k
    Integer                          :: ni, nj, nk, nl, nm

    ni = UBound(A,1)
    nj = UBound(A,2)
    nk = (nj+9) / 10

    Do k=1,nk
       nl = 10*k-9
       nm = Min(nj, nl+9)
       Write(NB6,'(4X,10I12)') (j,j=nl,nm)
       Do i=1,ni
          Write(NB6,'(1X,I6,10F12.6)') i, (A(i,j),j=nl,nm)
       End Do
       Write(NB6,'(//)')
    End Do

    Return
  End Subroutine PrintMatrix



  !
  ! This subroutine prints the upper triangle of a
  ! square matrix on unit NB6:
  !
  Subroutine PrintUpper(A, NB6)
    Implicit None
    Double Precision, Dimension(:,:) :: A
    Integer                          :: NB6
    ! end of dummy parameters
    Character(Len=120)               :: sp
    Integer                          :: i, j, k
    Integer                          :: ni, nk, nl, nm

    sp = ""
    ni = UBound(A,1)
    nk = (ni+9) / 10

    Do k=1,nk
       nl = 10*k-9
       nm = Min(ni, nl+9)
       Write(NB6,'(4X,10I12)') (j,j=nl,nm)
       Do i=1,nl
          Write(NB6,'(1X,I6,10F12.6)') i, (A(i,j),j=nl,nm)
       End Do
       Do i=nl+1,nm
          Write(NB6,'(1X,I6,A,10F12.6)') i, sp(1:12*(i-nl)), (A(i,j),j=i,nm)
       End Do
       Write(NB6,'(//)')
    End Do

    Return
  End Subroutine PrintUpper



  !
  ! This subroutine prints the upper triangle of a matrix
  ! in packed storage on unit NB6:
  !
  Subroutine PrintPackedUpper(A, ni, NB6)
    Implicit None
    Double Precision, Dimension(:) :: A
    Integer                        :: ni, NB6
    ! end of dummy parameters
    Character(Len=120)             :: sp
    Integer                        :: i, j, k
    Integer                        :: nk, nl, nm

    sp = ""
    nk = (ni+9) / 10

    Do k=1,nk
       nl = 10*k-9
       nm = Min(ni, nl+9)
       Write(NB6,'(4X,10I12)') (j,j=nl,nm)
       Do i=1,nl
          Write(NB6,'(1X,I6,10F12.6)') i, (A(j*(j-1)/2+i),j=nl,nm)
       End Do
       Do i=nl+1,nm
          Write(NB6,'(1X,I6,A,10F12.6)') i, sp(1:12*(i-nl)), (A(j*(j-1)/2+i),j=i,nm)
       End Do
       Write(NB6,'(//)')
    End Do

    Return
  End Subroutine PrintPackedUpper



  !
  ! This subroutine prints the lower triangle of a matrix
  ! in packed storage on unit NB6:
  !
  Subroutine PrintPackedLower(A, ni, NB6, ioffset)
    Implicit None
    Double Precision, Dimension(:) :: A
    Integer                        :: ni, NB6
    Integer, Optional              :: ioffset
    ! end of dummy parameters
    Integer                        :: i, j, k
    Integer                        :: nk, nl, nm
    Integer                        :: ioffs

    If (Present(ioffset)) Then
       ioffs = ioffset
    Else
       ioffs = 0
    End If

    nk = (ni+9) / 10

    Do k=1,nk
       nl = 10*k-9
       nm = Min(ni, nl+9)
       Write(NB6,'(4X,10I12)') (j,j=nl,nm)
       Do i=nl,nm
          Write(NB6,'(1X,I6,10F12.6)') i+ioffs, (A(i*(i-1)/2+j),j=nl,i)
       End Do
       Do i=nm+1,ni
          Write(NB6,'(1X,I6,10F12.6)') i+ioffs, (A(i*(i-1)/2+j),j=nl,nm)
       End Do
       Write(NB6,'(//)')
    End Do

    Return
  End Subroutine PrintPackedLower



  !
  ! This subroutine prints a vector on unit NB6:
  !
  Subroutine PrintVector(A, NB6)
    Implicit None
    Double Precision, Dimension(:) :: A
    Integer                        :: NB6
    ! end of dummy parameters
    Integer                        :: i

    Do i=1,UBound(A,1)
       Write(NB6,'(1X,I6,F12.6)') i, A(i)
    End Do
    Write(NB6,'(//)')

    Return
  End Subroutine PrintVector


End Module matutils
