!
! This module provides an interface to subroutine mkl_dcsrmm that is included
! in the Intel Math Kernel Library.
!
! Written by Axel Koslowski in 2012 at MPI Muelheim.
!

Module gugasparse

Implicit None

Private

Public :: dcsrsymm

Contains

  ! This subroutine computes the matrix product C = A * B (similar to DSYMM).
  ! A is a symmetric sparse matrix with the upper triangle stored in compressed
  ! sparse row (CSR) format using 1-based indexing with the following additional
  ! requirements. Every row must begin with its diagonal element, and this
  ! diagonal element must always be stored even if it is zero. B and C are
  ! dense matrices. m is the number of rows of all matrices and the number of
  ! columns of the matrix A. n is the number of columns of the matrices B and C.
  !
  Subroutine dcsrsymm(m,n,a,ia,ja,b,c)
    Integer          :: m, n, ia(:), ja(:)
    Double Precision :: a(:), b(:,:), c(:,:)
    ! End of dummy parameters.
    Integer          :: ldb, ldc, ldget

    If (n > 1) Then
       ldb = ldget(b(1,1), b(1,2))
       ldc = ldget(c(1,1), c(1,2))
    Else
       ldb = 1
       ldc = 1
    End If

    Call mkl_dcsrmm('N', m, n, m, 1.0D0, 'SUNF56', a(1), ja(1), ia(1), ia(2), b(1,1), ldb, 0.0D0, c(1,1), ldc)

  End Subroutine dcsrsymm

End Module gugasparse
