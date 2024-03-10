Subroutine priv_dcsrmm(transa,m,n,k,alpha,matdescra,a,ja,ia,ia2,b,ldb,beta,c,ldc)
  ! m: Number of rows of A.
  ! n: Number of columns of C.
  ! k: Number of columns of A.
  Character*1      :: transa
  Integer          :: m, n, k
  Double Precision :: alpha, beta
  Character*4      :: matdescra
  Integer          :: ia(m+1), ia2(m), ja(ia2(m)-1)
  Double Precision :: a(ia2(m)-1), b(ldb,n), c(ldc,n)
  ! end of dummy parameters
  Integer          :: i, j, ij

  ! Print:
  Write(*,'(1x,"m=",i10)') m
  Write(*,'(1x,"n=",i10)') n
  Write(*,'(1x,"k=",i10)') k
  Write(*,'(1x,"alpha=",f3.1)') alpha
  Write(*,'(1x,"beta =",f3.1)') beta
  Write(*,'(1x,"matdescra=",a4)') matdescra
  Write(*,'(1x,"ldb=",i10)') ldb
  Write(*,'(1x,"ldc=",i10)') ldc

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
End Subroutine priv_dcsrmm
