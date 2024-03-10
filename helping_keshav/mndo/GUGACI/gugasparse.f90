!
! This module provides a replacement for subroutine mkl_dcsrmm which is needed
! if the Intel Math Kernel Library is not available.
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
    Integer          :: i, j, ij

    ! Initialize C.
    c(:,:) = 0.0D0

    ! Loop over all rows of A.
    Do i=1,m
       ! Diagonal element of row i of A.
       ij     = ia(i)
       c(i,:) = c(i,:) + a(ij) * b(i,:)

       ! Loop over all off-diagonal elements of row i of A.
       Do ij=ia(i)+1, ia(i+1)-1
          j      = ja(ij)
          c(i,:) = c(i,:) + a(ij) * b(j,:)
          c(j,:) = c(j,:) + a(ij) * b(i,:)
       End Do
    End Do
  End Subroutine dcsrsymm

End Module gugasparse
