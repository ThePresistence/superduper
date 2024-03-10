!
! Subroutines performing a Davidson diagonalization.
!
! Written by Axel Koslowski in 2000-2003 at MPI Muelheim.
!
! Shape-driven algorithm and interface to C++ and CUDA C
! implementation added by Axel Koslowski in 2012-2013.
!

Module davidson

  Use gugaglobal, Only: ShavittControl
  Use gugaloops,  Only: CalcProdVec, CalcProdMat
  Use gugashape3, Only: ShapeDrivenProdVec
  Use gugashape4, Only: ShapeDrivenProdVecIrrep
  Use gugashape5, Only: ShapeDrivenProdMat
  Use gugashape6, Only: ShapeDrivenProdMatIrrep

  Private

  Public :: Davidson1, Davidson2, Davidson3, Davidson4

Contains

  !
  ! The following subroutine calculates the nroots lowest eigenstates of the CI Hamiltonian
  ! by the original Davidson procedure (E. R. Davidson, J. Comp. Phys. 17 (1975), 87-94;
  ! P. D. Dacre, Theor. Chim. Acta (Berl.) 43 (1976), 197). The mindav lowest-energy
  ! GUGA-CI basis functions (Gelfand states) plus the references are used as starting vectors.
  !
  ! Parameters:
  !
  ! First  letter: A = Always accessed,  I = InCore,  D = Direct
  ! Second letter: I = Input,  O = Output
  !
  ! Flags    (A,I):  ShavittControl structure.
  ! Hii      (A,I):  Diagonal elements of the CI Hamiltonian.
  ! Hij      (I,I):  Non-zero elements of the CI Hamiltonian.
  ! ifirst   (I,I):  ifirst(m) is the index of the first element of row m in Hij and icol.
  ! icol     (I,I):  Column index of the corresponding matrix element.
  ! MOInt1   (D,I):  One-electron MO integrals.
  ! MOInt2   (D,I):  Two-electron MO integrals.
  ! Occup    (D,I):  Occupation numbers of active orbitals.
  ! Irrep    (D,I):  Irreducible representation of Gelfand states.
  ! iRefConf (A,I):  Indices of CI basis functions corresponding to the reference configurations.
  ! i0       (A,I):  Number to be subtracted from the CI basis function indices to match the C block.
  ! SymLabel (A,I):  Labels of irreducible representations.
  ! E        (A,O):  Eigenvalues.
  ! C        (A,O):  Eigenvectors.
  ! nroots   (A,I):  Number of eigenstates to be calculated.
  ! mindav   (A,I):  Minimum number of columns of the B matrix.
  ! maxdav   (A,I):  Maximum number of columns of the B matrix.
  ! kitdav   (A,I):  Maximum number of Davidson iterations.
  ! qtol     (A,I):  Convergence criterion.
  ! cidir    (A,I):  Mode of calculation (in-core, semi-direct, or fully direct).
  ! stdout   (A,I):  Unit number of output file.
  ! prtlevel (A,I):  Defines output size. No output (except error messages) for prtlevel < 2.
  ! icall    (A,IO): Set to -1 to indicate error.
  !
  Subroutine Davidson1(Flags, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, iRefConf, i0, SymLabel, &
                       E, C, nroots, mindav0, maxdav0, kitdav0, qtol, cidir, stdout, prtlevel, icall)
    Implicit None
    Type(ShavittControl)                      :: Flags
    Double Precision, Dimension(:)            :: Hii, Hij
    Integer,          Dimension(:)            :: ifirst, icol
    Double Precision, Dimension(:,:)          :: MOInt1, MOInt2
    Double Precision, Dimension(:)            :: Occup
    Integer                                   :: Irrep
    Integer,          Dimension(:)            :: iRefConf
    Integer                                   :: i0
    Character(Len=4), Dimension(:)            :: SymLabel
    Double Precision, Dimension(:)            :: E
    Double Precision, Dimension(:,:)          :: C
    Integer                                   :: nroots, mindav0, maxdav0, kitdav0
    Double Precision                          :: qtol
    Integer                                   :: stdout, prtlevel, cidir, icall
    ! End of dummy parameters.
    Integer,          Dimension(:),   Pointer :: jRefConf
    Double Precision, Dimension(:,:), Pointer :: A, B, P, alpha, tmp
    Double Precision, Dimension(:),   Pointer :: lambda, q
    Integer                                   :: N, nref, ndav, ndavold, i, k
    Integer                                   :: mindav, maxdav, kitdav
    Integer                                   :: ldc, ldget
    Double Precision                          :: qnorm
    Double Precision                          :: totalcpu0, totalcpu1, totalcpu2
    Double Precision                          :: wallclock0, wallclock2
    Double Precision                          :: tuser, tsys, twall

    If (prtlevel >= 2) Then
       Write(stdout,'(//A)') " Davidson's original algorithm."
       Call GetTime(tuser, tsys, twall)
       wallclock0 = twall
       totalcpu0  = tuser + tsys
       totalcpu1  = totalcpu0
    End If

    If (mindav0 > 0) Then
       mindav = mindav0
    Else
       mindav = Max(2*nroots, 30)
    End If

    If (maxdav0 > 0) Then
       maxdav = maxdav0
    Else
       maxdav = mindav + 30
    End If

    If (kitdav0 > 0) Then
       kitdav = kitdav0
    Else
       kitdav = 40 * nroots
    End If

    If (prtlevel >= 2) Then
       Write(stdout,'(/1X,"MINDAV=",I4,4X,"MAXDAV=",I4,4X,"KITDAV=",I4)') mindav, maxdav, kitdav
    End If

    N = UBound(Hii,1)
    Allocate(A(maxdav,maxdav))
    Allocate(B(N,maxdav))
    Allocate(P(N,maxdav))
    Allocate(alpha(maxdav,maxdav))
    Allocate(lambda(maxdav))
    Allocate(q(N))
    Allocate(jRefConf(mindav + UBound(iRefConf,1)))
    Call InitRefConf(jRefConf, iRefConf, i0, Hii, mindav, nref)
    If (prtlevel >= 2)  Call PrintRefConf(Flags, Hii, jRefConf(1:nref), SymLabel, i0, stdout)
    Call InitB(jRefConf(1:nref), 0, B, ndav, 1, stdout)
    Deallocate(jRefConf)
    Call CalcP(Flags, P, B, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
    Call InitA(B, P, A, ndav)
    i = 1
    k = 1

    If (prtlevel >= 2) &
      & Write(stdout,'(/A//A)') &
      &   ' Davidson iterations:',    &
      &   ' Iter.    k        lambda(k)              Norm          CPU time      Total CPU    Wall Clock'

    iter: Do
       Call Eigen(A, lambda, alpha, ndav, nroots, stdout, icall)
       If (icall < 0) Exit iter
       ! Find next root that has not converged:
       k = 1
       Do
          Call FormQVec(B, P, alpha(:,k), lambda(k), q, qnorm, ndav)
          If (qnorm > qtol) Exit
          If (k >= nroots) Then
             Write(stdout,'(1X,"All roots converged.")')
             Exit iter
          End If
          k = k + 1
       End Do
       If (prtlevel >= 2) Then
          Call GetTime(tuser, tsys, twall)
          totalcpu2  = tuser + tsys
          wallclock2 = twall
          Write(stdout,'(I4,I7,F23.15,E16.4,3(F12.2,A2))') &
            &   i, k, lambda(k), qnorm,                    &
            &   totalcpu2  - totalcpu1,  ' s',             &
            &   totalcpu2  - totalcpu0,  ' s',             &
            &   wallclock2 - wallclock0, ' s'
          totalcpu1 = totalcpu2
       End If
       ! Check if maximum number of iterations has been reached:
       If (i == kitdav) Then
          Write(stdout,'(1X,A)') 'The Davidson procedure failed to converge.'
          icall = -1
          Exit iter
       End If
       If (ndav < maxdav) Then
          ! Calculate next b vector:
          ndavold = ndav
          Call UpdateB(Hii, q, lambda(k), qtol, B, ndav)
          If (ndav <= ndavold) Then
             Write(stdout,'(1X,A)') 'No new basis vector found. Terminating.'
             icall = -1
             Exit iter
          End If
          ! Form H * b; this is the rate determining step:
          Select Case (cidir)
          Case (-2,-1,1,2)
             Call FormPVec(Hij, ifirst, icol, N, B(:,ndav), P(:,ndav))
          Case (3)
             If (Irrep == 0) Then
                Call ShapeDrivenProdVec(Flags, P(:,ndav), B(:,ndav), Hii, MOInt1, MOInt2)
             Else
                Call ShapeDrivenProdVecIrrep(Flags, P(:,ndav), B(:,ndav), Hii, MOInt1, MOInt2, Irrep, i0)
             End If
          Case (-3)
             Call CalcProdVec(Flags, P(:,ndav), B(:,ndav), MOInt1, MOInt2, Occup, Irrep, i0)
          Case Default
             Write(stdout,'(1X,A,I2,A)') 'Unexpected value for cidir (', cidir, ').'
             Stop 'davidson::Davidson1'
          End Select
          ! Update A matrix:
          Call UpdateA(B, P, A, ndav)
          ! Next iteration:
          i = i + 1
       Else
          ! Perform a restart:
          !-------------------
          ndavold = ndav
          ndav    = Max(nroots, mindav)
          ! Calculate more roots if necessary:
          If (ndav > nroots) Then
             Call Eigen(A, lambda, alpha, ndavold, ndav, stdout, icall)
             If (icall < 0) Exit iter
          End If
          ! Calculate new B matrix in P and orthonormalize:
          ! P(:,1:ndav) = Matmul(B(:,1:ndavold), alpha(1:ndavold,1:ndav))
          Call DGEMM('N', 'N', N, ndav, ndavold, 1.D0, B(1,1), N, alpha(1,1), maxdav, 0.D0, P(1,1), N)
          Call Schmidt(P, 1, ndav, qtol)
          Call Schmidt(P, 1, ndav, qtol)
          ! Calculate new P matrix in B:
          Call CalcP(Flags, B, P, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
          ! Exchange B and P:
          tmp => P
          P   => B
          B   => tmp
          ! Calculate new A matrix:
          Call InitA(B, P, A, ndav)
          If (prtlevel >= 2) &
            & Write(stdout,'(A)') ' The Davidson procedure has been restarted.'
       End If
    End Do iter

    If (prtlevel >= 2) Write(stdout,'(//)')

    If (icall >= 0) Then
       E(1:nroots) = lambda(1:nroots)
       If (nroots > 1) Then
          ldc = ldget(C(1,1), C(1,2))
       Else
          ldc = UBound(C,1)
       End If
       ! C(:,1:nroots) = Matmul(B(:,1:ndav), alpha(1:ndav,1:nroots))
       Call DGEMM('N', 'N', N, nroots, ndav, 1.D0, B(1,1), N, alpha(1,1), maxdav, 0.D0, C(1,1), ldc)
    End If

    Deallocate(q)
    Deallocate(lambda)
    Deallocate(alpha)
    Deallocate(P)
    Deallocate(B)
    Deallocate(A)

    Return
  End Subroutine Davidson1



  !
  ! The following subroutine calculates a single eigenstate of the CI Hamiltonian
  ! by the modified Davidson procedure by Butscher and Kammer (W. Butscher,
  ! W. E. Kammer, J. Comp. Phys. 20 (1976), 313-325). All parameters except lroot
  ! are identical to subroutine Davidson1.
  !
  ! Parameters:
  !
  ! First  letter: A = Always accessed,  I = InCore,  D = Direct
  ! Second letter: I = Input,  O = Output
  !
  ! Flags    (A,I):  ShavittControl structure.
  ! Hii      (A,I):  Diagonal elements of the CI Hamiltonian.
  ! Hij      (I,I):  Off-diagonal elements of the CI Hamiltonian.
  ! ifirst   (I,I):  ifirst(m) is the index of the first element of row m in Hij and icol.
  ! icol     (I,I):  Column index of the corresponding matrix element.
  ! icol0    (I,I):  
  ! MOInt1   (D,I):  One-electron MO integrals.
  ! MOInt2   (D,I):  Two-electron MO integrals.
  ! Occup    (D,I):  Occupation numbers of active orbitals.
  ! Irrep    (D,I):  Irreducible representation of Gelfand states.
  ! iRefConf (A,I):  Indices of CI basis functions corresponding to the reference configurations.
  ! i0       (A,I):  Number to be subtracted from the CI basis function indices to match the C block.
  ! SymLabel (A,I):  Labels of irreducible representations.
  ! E        (A,O):  Eigenvalues.
  ! C        (A,O):  Eigenvectors.
  ! lroot    (A,I):  > 0: The lroot-th eigenvector when doing a CI with only the
  !                       reference Gelfand states is taken as the starting vector.
  !                  < 0: The eigenvector in which the abs(lroot)-th reference
  !                       configuration has the largest norm (sum of squares of
  !                       CI coefficients of corresponding reference Gelfand states)
  !                       is taken as the starting vector.
  ! mindav   (A,I):  Minimum number of columns of the B matrix after a restart.
  ! maxdav   (A,I):  Maximum number of columns of the B matrix.
  ! kitdav   (A,I):  Maximum number of Davidson iterations.
  ! qtol     (A,I):  Convergence criterion.
  ! cidir    (A,I):  Mode of calculation (in-core, semi-direct, or fully direct).
  ! stdout   (A,I):  Unit number of output file.
  ! prtlevel (A,I):  Defines output size. No output (except error messages) for prtlevel < 2.
  ! icall    (A,IO): Set to -1 to indicate error.
  !
  Subroutine Davidson2(Flags, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, iFirstRef, iLastRef, iRefConf, &
                       i0, SymLabel, E, C, lroot, mindav0, maxdav0, kitdav0, qtol, cidir, stdout, prtlevel, icall)
    Implicit None
    Type(ShavittControl)                      :: Flags
    Double Precision, Dimension(:)            :: Hii, Hij
    Integer,          Dimension(:)            :: ifirst, icol
    Double Precision, Dimension(:,:)          :: MOInt1, MOInt2
    Double Precision, Dimension(:)            :: Occup
    Integer                                   :: Irrep
    Integer,          Dimension(:)            :: iFirstRef
    Integer,          Dimension(:)            :: iLastRef
    Integer,          Dimension(:)            :: iRefConf
    Integer                                   :: i0
    Character(Len=4), Dimension(:)            :: SymLabel
    Double Precision                          :: E
    Double Precision, Dimension(:)            :: C
    Integer                                   :: lroot, mindav0, maxdav0, kitdav0
    Double Precision                          :: qtol
    Integer                                   :: stdout, prtlevel, cidir, icall
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Pointer :: A, B, P, alpha, tmp
    Double Precision, Dimension(:),   Pointer :: lambda, delta, q, alpha0
    Integer                                   :: N, ndav, ndavold, nref, ierr, i, k, i1, i2
    Integer                                   :: mindav, maxdav, kitdav
    Double Precision                          :: qnorm, overlap, overlap0, fk, frac
    Double Precision                          :: totalcpu0, totalcpu1, totalcpu2
    Double Precision                          :: wallclock0, wallclock2
    Double Precision                          :: tuser, tsys, twall

    If (prtlevel >= 2) Then
       Write(stdout,'(//A)') " Davidson's algorithm with modification by Butscher and Kammer."
       Call GetTime(tuser, tsys, twall)
       wallclock0 = twall
       totalcpu0  = tuser + tsys
       totalcpu1  = totalcpu0
    End If

    If (mindav0 > 0) Then
       mindav = mindav0
    Else
       mindav = UBound(iRefConf,1)
    End If

    If (maxdav0 > 0) Then
       maxdav = maxdav0
    Else
       maxdav = mindav + 30
    End If

    If (kitdav0 > 0) Then
       kitdav = kitdav0
    Else
       kitdav = 40
    End If

    If (prtlevel >= 2) Then
       Write(stdout,'(/1X,"MINDAV=",I4,4X,"MAXDAV=",I4,4X,"KITDAV=",I4)') mindav, maxdav, kitdav
    End If

    N = UBound(Hii,1)
    Allocate(A(maxdav,maxdav))
    Allocate(B(N,maxdav))
    Allocate(P(N,maxdav))
    Allocate(alpha(maxdav,maxdav))
    Allocate(lambda(maxdav))
    Allocate(delta(maxdav))
    Allocate(q(N))
    If (prtlevel >= 2)  Call PrintRefConf(Flags, Hii, iRefConf, SymLabel, 0, stdout)
    Call InitB(iRefConf, i0, B, ndav, Abs(lroot), stdout)
    Call CalcP(Flags, P, B, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
    Call InitA(B, P, A, ndav)

    ! Diagonalize submatrix:
    Call TRED2(maxdav, ndav, A, lambda, delta, alpha)
    Call TQL2(maxdav, ndav, lambda, delta, alpha, ierr)
    If (ierr /= 0) Then
       Write(stdout, '(A)') ' The diagonalization of the submatrix failed.'
       Deallocate(A)
       Deallocate(B)
       Deallocate(P)
       Deallocate(alpha)
       Deallocate(lambda)
       Deallocate(delta)
       Deallocate(q)
       icall = -1
       Return
    End If

    ! Print initial eigenvectors and headline for iterations.
    If (prtlevel >= 2) Then
       Call PrintEigenVectors(iRefConf, alpha, stdout)
       Write(stdout,'(/A//A)')      &
         ' Davidson iterations:',   &
         ' Iter.    k        lambda(k)              Overlap         Norm          CPU time      Total CPU    Wall Clock'
    End If

    ! Save the eigenvector of interest:
    Allocate(alpha0(maxdav))
    If (lroot > 0) Then
       alpha0(1:ndav) = alpha(1:ndav,lroot)
       overlap0       = 1.D0
       overlap        = 1.D0
       k              = lroot
    Else
       If (-lroot > UBound(iFirstRef,1)) Then
          Write(stdout,'(1X,"Negative lroot inconsistent with number of references.")')
          Stop 'davidson::Davidson2'
       End If
       i1 = iFirstRef(-lroot)
       i2 = iLastRef(-lroot)
       If (i1 > i2) Then
          Write(stdout,'(1X,"Negative lroot points to reference not represented by the DRT.")')
          Stop 'davidson::Davidson2'
       End If
       k  = 1
       fk = Dot_Product(alpha(i1:i2,1), alpha(i1:i2,1))
       Do i=2, ndav
          frac = Dot_Product(alpha(i1:i2,i), alpha(i1:i2,i))
          If (frac > fk) Then
             k  = i
             fk = frac
          End If
       End Do
       alpha0(1:ndav) = alpha(1:ndav,k)
       overlap0       = Sum(alpha(i1:i2,k))
       overlap        = 1.D0
    End If
    nref = ndav

    i = 1
    iter: Do
       ! Check if root k has converged:
       Call FormQVec(B, P, alpha(:,k), lambda(k), q, qnorm, ndav)
       If (prtlevel >= 2) Then
          Call GetTime(tuser, tsys, twall)
          totalcpu2  = tuser + tsys
          wallclock2 = twall
          Write(stdout,'(I4,I7,F23.15,F16.8,E16.4,3(F12.2,A2))') &
            &   i, k, lambda(k), overlap*overlap0, qnorm,        &
            &   totalcpu2  - totalcpu1,  ' s',                   &
            &   totalcpu2  - totalcpu0,  ' s',                   &
            &   wallclock2 - wallclock0, ' s'
          totalcpu1 = totalcpu2
       End If
       If (qnorm < qtol) Exit iter
       ! Check if maximum number of iterations has been reached:
       If (i == kitdav) Then
          Write(stdout,'(A)') ' The Davidson procedure failed to converge.'
          icall = -1
          Exit iter
       End If
       If (ndav < maxdav) Then
          ! Calculate next b vector:
          ndavold = ndav
          Call UpdateB(Hii, q, lambda(k), qtol, B, ndav)
          If (ndav <= ndavold) Then
             Write(stdout,'(A)') ' No new basis vector found. Terminating.'
             icall = -1
             Exit iter
          End If
          ! Form H * b; this is the rate determining step:
          Select Case (cidir)
          Case (-1,1)
             Call FormPVec(Hij, ifirst, icol, N, B(:,ndav), P(:,ndav))
          Case (3)
             If (Irrep == 0) Then
                Call ShapeDrivenProdVec(Flags, P(:,ndav), B(:,ndav), Hii, MOInt1, MOInt2)
             Else
                Call ShapeDrivenProdVecIrrep(Flags, P(:,ndav), B(:,ndav), Hii, MOInt1, MOInt2, Irrep, i0)
             End If
          Case (-3)
             Call CalcProdVec(Flags, P(:,ndav), B(:,ndav), MOInt1, MOInt2, Occup, Irrep, i0)
          Case Default
             Write(stdout,'(1X,A,I2,A)') 'Unexpected value for cidir (', cidir, ').'
             Stop 'davidson::Davidson2'
          End Select
          ! Update A matrix:
          Call UpdateA(B, P, A, ndav)
          ! Next iteration:
          i = i + 1
       Else
          ! Perform a restart:
          !-------------------
          ! There are now maxdav valid columns in B.
          Allocate(tmp(N,maxdav))
          ndav = 1
          ! Calculate new B matrix in tmp:
          ! tmp(:,1) = Matmul(B, alpha(:,k))
          Call DGEMV('N', N, maxdav, 1.D0, B(1,1), N, alpha(1,k), 1, 0.D0, tmp(1,1), 1)
          ! Calculate new P matrix in old B (which is no longer needed):
          ! B(:,1)   = Matmul(P, alpha(:,k))
          Call DGEMV('N', N, maxdav, 1.D0, P(1,1), N, alpha(1,k), 1, 0.D0,   B(1,1), 1)
          ! Free storage of old P (which is also no longer needed):
          Deallocate(P)
          ! Restore pointers:
          P => B
          B => tmp
          ! Restore reference vector:
          overlap0 = overlap0 * overlap
          alpha0(1) = 1.D0
          nref = 1
          ! Calculate new A matrix:
          Call InitA(B, P, A, ndav)
          Write(stdout,'(A)') ' The Davidson procedure has been restarted.'
       End If
       ! Diagonalize submatrix:
       Call TRED2(maxdav, ndav, A, lambda, delta, alpha)
       Call TQL2(maxdav, ndav, lambda, delta, alpha, ierr)
       If (ierr /= 0) Then
          Write(stdout, '(A)') ' The diagonalization of the submatrix failed.'
          icall = -1
          Exit iter
       End If
       Call FindK(alpha, alpha0, ndav, nref, k, overlap)
    End Do iter

    If (prtlevel >= 2) Write(stdout,'(//)')

    If (icall >= 0) Then
       E = lambda(k)
       ! C = Matmul(B(:,1:ndav), alpha(1:ndav,k))
       Call DGEMV('N', N, ndav, 1.D0, B(1,1), N, alpha(1,k), 1, 0.D0, C, 1)
    End If

    Deallocate(alpha0)
    Deallocate(q)
    Deallocate(delta)
    Deallocate(lambda)
    Deallocate(alpha)
    Deallocate(P)
    Deallocate(B)
    Deallocate(A)

    Return
  End Subroutine Davidson2



  !
  ! The following subroutine calculates the nroots lowest eigenstates of the CI Hamiltonian
  ! by the Davidson procedure with the modification by Liu (E. R. Davidson, J. Comp. Phys.
  ! 17 (1975), 87-94; P. D. Dacre, Theor. Chim. Acta (Berl.) 43 (1976), 197; B. Liu, in
  ! C. Moler, I. Shavitt, Numerical Algorithms in Chemistry: Algebraic Methods, LBL-8158
  ! Lawrence Berkeley Laboratory (1978), 49-53). The mindav lowest-energy GUGA-CI
  ! many-electron basis functions plus the references are used as starting vectors.
  !
  ! Parameters:
  !
  ! First  letter: A = Always accessed,  I = InCore,  D = Direct
  ! Second letter: I = Input,  O = Output
  !
  ! Flags    (A,I):  ShavittControl structure.
  ! Hii      (A,I):  Diagonal elements of the CI Hamiltonian.
  ! Hij      (I,I):  Off-diagonal elements of the CI Hamiltonian.
  ! ifirst   (I,I):  ifirst(m) is the index of the first element of row m in Hij and icol.
  ! icol     (I,I):  Column index of the corresponding matrix element.
  ! MOInt1   (D,I):  One-electron MO integrals.
  ! MOInt2   (D,I):  Two-electron MO integrals.
  ! Occup    (D,I):  Occupation numbers of active orbitals.
  ! Irrep    (D,I):  Irreducible representation of Gelfand states.
  ! iRefConf (A,I):  Indices of CI basis functions corresponding to the reference configurations.
  ! i0       (A,I):  Number to be subtracted from the CI basis function indices to match the C block.
  ! SymLabel (A,I):  Labels of irreducible representations.
  ! E        (A,O):  Eigenvalues.
  ! C        (A,O):  Eigenvectors.
  ! nroots   (A,I):  Number of eigenstates to be calculated.
  ! mindav   (A,I):  Minimum number of columns of the B matrix.
  ! maxdav   (A,I):  Maximum number of columns of the B matrix.
  ! kitdav   (A,I):  Maximum number of Davidson iterations.
  ! qtol     (A,I):  Convergence criterion.
  ! cidir    (A,I):  Mode of calculation (in-core, semi-direct, or fully direct).
  ! stdout   (A,I):  Unit number of output file.
  ! prtlevel (A,I):  Defines output size. No output (except error messages) for prtlevel < 2.
  ! icall    (A,IO): Set to -1 to indicate error.
  !
  Subroutine Davidson3(Flags, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, iRefConf, i0, SymLabel, &
                       E, C, nroots, mindav0, maxdav0, kitdav0, qtol, cidir, stdout, prtlevel, icall)
    Implicit None
    Type(ShavittControl)                      :: Flags
    Double Precision, Dimension(:)            :: Hii, Hij
    Integer,          Dimension(:)            :: ifirst, icol
    Double Precision, Dimension(:,:)          :: MOInt1, MOInt2
    Double Precision, Dimension(:)            :: Occup
    Integer                                   :: Irrep
    Integer,          Dimension(:)            :: iRefConf
    Integer                                   :: i0
    Character(Len=4), Dimension(:)            :: SymLabel
    Integer                                   :: nroots, mindav0, maxdav0, kitdav0
    Double Precision, Dimension(:)            :: E
    Double Precision, Dimension(:,:)          :: C
    Double Precision                          :: qtol
    Integer                                   :: stdout, prtlevel, cidir, icall
    ! End of dummy parameters.
    Integer,          Dimension(:),   Pointer :: jRefConf
    Double Precision, Dimension(:,:), Pointer :: A, B, P, alpha, Q, tmp
    Double Precision, Dimension(:),   Pointer :: lambda
    Double Precision, Dimension(nroots)       :: qnorm
    Integer                                   :: N, nref, ndav, ndavold, nconv, kbad, new, i, k
    Integer                                   :: mindav, maxdav, kitdav
    Integer                                   :: ldc, ldget
    Double Precision                          :: tuser1, tsys1, twall1, tcpu1  ! Initial timing.
    Double Precision                          :: tuser2, tsys2, twall2, tcpu2  ! Previous timing.
    Double Precision                          :: tuser3, tsys3, twall3, tcpu3  ! Current timing.

    If (prtlevel >= 2) Then
       Write(stdout,'(/A)') " Davidson's algorithm with modification by Liu."
       Call GetTime(tuser1, tsys1, twall1)
       tcpu1  = tuser1 + tsys1
       tcpu2  = tcpu1
       twall2 = twall1
    End If

    If (mindav0 > 0) Then
       mindav = mindav0
    Else
       mindav = Max(2*nroots, 30)
    End If

    If (maxdav0 > 0) Then
       maxdav = maxdav0
    Else
       maxdav = mindav + 15 * (nroots+1)
    End If

    If (kitdav0 > 0) Then
       kitdav = kitdav0
    Else
       kitdav = 40 * nroots
    End If

    If (prtlevel >= 2) Then
       Write(stdout,'(/1X,"MINDAV=",I4,4X,"MAXDAV=",I4,4X,"KITDAV=",I4)') mindav, maxdav, kitdav
    End If

    N = UBound(Hii,1)
    Allocate(A(maxdav,maxdav))
    Allocate(B(N,maxdav))
    Allocate(P(N,maxdav))
    Allocate(alpha(maxdav,maxdav))
    Allocate(lambda(maxdav))
    Allocate(Q(N,nroots))
    Allocate(jRefConf(mindav + UBound(iRefConf,1)))
    Call InitRefConf(jRefConf, iRefConf, i0, Hii, mindav, nref)
    If (prtlevel >= 2) Then
       Call PrintRefConf(Flags, Hii, jRefConf(1:nref), SymLabel, i0, stdout)
       Write(stdout,'(/A//A)')     &
         ' Davidson iterations:',  &
         ' Iter.   Basis   To go   Worst      Eigenvalue           Norm      CPU time   Wall clock   Total wall'
!            1       30      5       5     3.123456789012345    0.1234E+12    12.34 s     1.23 s      123.45 s
    End If
    Call InitB(jRefConf(1:nref), 0, B, ndav, 1, stdout)
    Deallocate(jRefConf)
    Call CalcP(Flags, P, B, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
    Call InitA(B, P, A, ndav)
    i = 1

    iter: Do
       Call Eigen(A, lambda, alpha, ndav, nroots, stdout, icall)
       If (icall < 0) Exit iter
       ! Find roots that have not converged:
       Call FormQMat(B, P, alpha(:,1:nroots), lambda(1:nroots), Q, qnorm, ndav, nroots)
       nconv = 0
       kbad  = 1
       Do k=1,nroots
          If (qnorm(k) < qtol)         nconv = nconv + 1
          If (qnorm(k) > qnorm(kbad))  kbad  = k
       End Do
       new = nroots - nconv
       If (prtlevel >= 2) Then
          Call GetTime(tuser3, tsys3, twall3)
          tcpu3 = tuser3 + tsys3
          Write(stdout,'(I4,I9,I7,I8,F22.15,E14.4,F9.2,A2,F9.2,A2,F12.2,A2)')  &
            &   i, ndav, new, kbad, lambda(kbad), qnorm(kbad),                 &
            &   tcpu3  - tcpu2,  ' s',                                         &
            &   twall3 - twall2, ' s',                                         &
            &   twall3 - twall1, ' s'
          tcpu2  = tcpu3
          twall2 = twall3
       End If
       If (new == 0) Exit iter
       ! Check if maximum number of iterations has been reached:
       If (i == kitdav) Then
          Write(stdout,'(A)') ' The Davidson procedure failed to converge.'
          icall = -1
          Exit iter
       End If
       If (ndav+new <= maxdav) Then
          ndavold = ndav
          ! Calculate new part of B matrix and adjust ndav:
          Call UpdateBMat(Hii, Q, lambda(1:nroots), qnorm, qtol, B, ndav)
          If (ndav <= ndavold) Then
             Write(stdout,'(A)') ' No new basis vectors found. Convergence achieved.'
             Exit iter
          End If
          ! Calculate corresponding part of P:
          Call CalcP(Flags, P(:,ndavold+1:ndav), B(:,ndavold+1:ndav), Hii, Hij, ifirst, icol, &
                     MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav-ndavold, stdout)
          ! Update A matrix:
          Call UpdateAMat(B, P, A, ndavold, ndav)
          ! Next iteration:
          i = i + 1
       Else
          ! Perform a restart:
          !-------------------
          ndavold = ndav
          ndav    = Max(nroots, mindav)
          ! Calculate more roots if necessary:
          If (ndav > nroots) Then
             Call Eigen(A, lambda, alpha, ndavold, ndav, stdout, icall)
             If (icall < 0) Exit iter
          End If
          ! Calculate new B matrix in P and orthonormalize:
          ! P(:,1:ndav) = Matmul(B(:,1:ndavold), alpha(1:ndavold,1:ndav))
          Call DGEMM('N', 'N', N, ndav, ndavold, 1.D0, B(1,1), N, alpha(1,1), maxdav, 0.D0, P(1,1), N)
          Call Schmidt(P, 1, ndav, qtol)
          Call Schmidt(P, 1, ndav, qtol)
          ! Calculate new P matrix in B:
          Call CalcP(Flags, B, P, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
          ! Exchange B and P:
          tmp => P
          P   => B
          B   => tmp
          ! Calculate new A matrix:
          Call InitA(B, P, A, ndav)
          If (prtlevel >= 2) &
            & Write(stdout,'(A)') ' The Davidson procedure has been restarted.'
       End If
    End Do iter

    If (prtlevel >= 2) Write(stdout,'(//)')

    If (icall >= 0) Then
       E(1:nroots)   = lambda(1:nroots)
       If (nroots > 1) Then
          ldc = ldget(C(1,1), C(1,2))
       Else
          ldc = 1
       End If
       ! C(:,1:nroots) = Matmul(B(:,1:ndav), alpha(1:ndav,1:nroots))
       Call DGEMM('N', 'N', N, nroots, ndav, 1.D0, B(1,1), N, alpha(1,1), maxdav, 0.D0, C(1,1), ldc)
    End If

    Deallocate(Q)
    Deallocate(lambda)
    Deallocate(alpha)
    Deallocate(P)
    Deallocate(B)
    Deallocate(A)

    Return
  End Subroutine Davidson3



  !
  ! The following subroutine is an interface to the routine DavLiu that calculates the
  ! nroots lowest eigenstates of the CI Hamiltonian by a GPU-adapted version of the
  ! Davidson procedure with the modification of Liu (E. R. Davidson, J. Comp. Phys.
  ! 17 (1975), 87-94; P. D. Dacre, Theor. Chim. Acta (Berl.) 43 (1976), 197; B. Liu, in
  ! C. Moler, I. Shavitt, Numerical Algorithms in Chemistry: Algebraic Methods, LBL-8158
  ! Lawrence Berkeley Laboratory (1978), 49-53). This interface routine provides the
  ! initial basis vectors for the subspace algorithm (the mindav lowest-energy
  ! GUGA-CI many-electron basis functions plus the references) and calls the routine
  ! DavLiu which may be applied in the context of in-core calculations only. Most of
  ! the linear algebra in DavLiu is performed on the GPU.
  !
  ! Parameters:
  ! I = Input,  O = Output  (In-core only.)
  !
  ! Flags    (I):  ShavittControl structure.
  ! Hii      (I):  Diagonal elements of the CI Hamiltonian.
  ! Hij      (I):  Off-diagonal elements of the CI Hamiltonian.
  ! ifirst   (I):  ifirst(m) is the index of the first element of row m in Hij and icol.
  ! icol     (I):  Column index of the corresponding matrix element.
  ! iRefConf (I):  Indices of CI basis functions corresponding to the reference configurations.
  ! i0       (I):  Number to be subtracted from the CI basis function indices to match the C block.
  ! SymLabel (I):  Labels of irreducible representations.
  ! E        (O):  Eigenvalues.
  ! C        (O):  Eigenvectors.
  ! nroots   (I):  Number of eigenstates to be calculated.
  ! mindav   (I):  Minimum number of columns of the B matrix.
  ! maxdav   (I):  Maximum number of columns of the B matrix.
  ! kitdav   (I):  Maximum number of Davidson iterations.
  ! qtol     (I):  Convergence criterion.
  ! stdout   (I):  Unit number of output file.
  ! prtlevel (I):  Defines output size. No output (except error messages) for prtlevel < 2.
  ! icall    (IO): Set to -1 to indicate error.
  !
  Subroutine Davidson4(Flags, Hii, Hij, ifirst, icol, iRefConf, i0, SymLabel, E, C, nroots, &
                       mindav0, maxdav0, kitdav0, qtol, stdout, prtlevel, icall, DavLiu)
    Implicit None
    Type(ShavittControl)                      :: Flags
    Double Precision, Dimension(:)            :: Hii, Hij
    Integer,          Dimension(:)            :: ifirst, icol
    Integer                                   :: i0
    Integer,          Dimension(:)            :: iRefConf
    Character(Len=4), Dimension(:)            :: SymLabel
    Integer                                   :: nroots, mindav0, maxdav0, kitdav0
    Double Precision, Dimension(:)            :: E
    Double Precision, Dimension(:,:)          :: C
    Double Precision                          :: qtol
    Integer                                   :: stdout, prtlevel, icall
    ! End of dummy parameters.
    Integer,          Dimension(:),   Pointer :: jRefConf
    Integer                                   :: mindav, maxdav, kitdav
    Integer                                   :: nci, inidav
    Integer                                   :: ldc, ldget
    External                                  :: DavLiu

    If (mindav0 > 0) Then
       mindav = mindav0
    Else
       mindav = Max(2*nroots, 30)
    End If

    If (maxdav0 > 0) Then
       maxdav = maxdav0
    Else
       maxdav = mindav + 15 * (nroots+1)
    End If

    If (kitdav0 > 0) Then
       kitdav = kitdav0
    Else
       kitdav = 40 * nroots
    End If

    Allocate(jRefConf(mindav + UBound(iRefConf,1)))
    Call InitRefConf(jRefConf, iRefConf, i0, Hii, mindav, inidav)
    If (prtlevel >= 2)  Call PrintRefConf(Flags, Hii, jRefConf(1:inidav), SymLabel, i0, stdout)
    If (nroots > 1) Then
       ldc = ldget(C(1,1), C(1,2))
    Else
       ldc = 1
    End If
    nci = UBound(Hii,1)
    Flush(stdout)  ! Flush output to unit 6 so Fortran and C output will be well separated.
    Call DavLiu(Hii(1), Hij(1), ifirst(1), icol(1), nci, jRefConf(1), inidav, E(1), &
                C(1,1), ldc, nroots, mindav, maxdav, kitdav, qtol, prtlevel, icall)
    Deallocate(jRefConf)
  End Subroutine Davidson4



  !
  ! Select the mindav lowest-energy Gelfand states plus the original references.
  ! While the indices in iRefConf refer to the global CI problem, jRefConf
  ! contains local indices (i.e. i0 is subtracted from the global indices).
  ! On output, inidav is the resulting number of initial basis vectors for the
  ! Davidson algorithm.
  !
  Subroutine InitRefConf(jRefConf, iRefConf, i0, Hii, mindav, inidav)
    Implicit None
    Integer,          Dimension(:)   :: jRefConf, iRefConf
    Integer                          :: i0
    Double Precision, Dimension(:)   :: Hii
    Integer                          :: mindav, inidav
    ! End of dummy parameters.
    Integer                          :: nci, mref, i, j

    nci    = UBound(Hii,1)
    inidav = mindav

    If (inidav > nci)  inidav = nci

    ! Initially, use the first inidav CI basis functions:
    Do i=1, inidav
       jRefConf(i) = i
    End Do

    ! Sort by increasing energy:
    Call SortRefConf(Hii, jRefConf, inidav)

    If (inidav == nci)  Return

    ! Replace the element with the highest energy in the set of
    ! selected CI basis functions with a CI basis function of
    ! lower energy currently not in the set:
    Do i=inidav+1, nci
       If (Hii(i) < Hii(jRefConf(inidav))) Then
          jRefConf(inidav) = i
          Call SortRefConf(Hii, jRefConf, inidav)
       End If
    End Do

    ! If appropriate, add reference configurations
    ! that are not already in the set:
    mref = inidav
    iref: Do i=1, UBound(iRefConf,1)
       If (iRefConf(i) <= i0 .Or. iRefConf(i)-i0 > nci)  Cycle iref
       Do j=1, mref
          If (iRefConf(i)-i0 == jRefConf(j))  Cycle iref
       End Do
       inidav = inidav + 1
       jRefConf(inidav) = iRefConf(i) - i0
    End Do iref

    If (inidav > mref)  Call SortRefConf(Hii, jRefConf, inidav)

    Return
  End Subroutine InitRefConf



  !
  ! Sort the reference configurations according to increasing energy.
  !
  Subroutine SortRefConf(Hii, jRefConf, nref)
    Implicit None
    Double Precision, Dimension(:)   :: Hii
    Integer,          Dimension(:)   :: jRefConf
    Integer                          :: nref
    ! End of dummy parameters.
    Integer                          :: i, j, k

    Do i=1, nref-1
       k = i
       Do j=i+1, nref
          If (Hii(jRefConf(j)) < Hii(jRefConf(k)))  k = j
       End Do
       If (k /= i) Then
          j           = jRefConf(i)
          jRefConf(i) = jRefConf(k)
          jRefConf(k) = j
       End If
    End Do

    Return
  End Subroutine SortRefConf



  !
  ! Print the reference configurations.
  !
  Subroutine PrintRefConf(Flags, Hii, jRefConf, SymLabel, i0, stdout)
    Implicit None
    Type(ShavittControl)             :: Flags
    Double Precision, Dimension(:)   :: Hii
    Integer,          Dimension(:)   :: jRefConf
    Character(Len=4), Dimension(:)   :: SymLabel
    Integer                          :: i0, stdout
    ! End of dummy parameters.
    Integer                          :: i, j, s

    Write(stdout,'(//A)') ' Starting vectors:'
    Write(stdout,'(/ A)') '   #     Index   Symmetry     Energy'
    Do i=1, UBound(jRefConf,1)
       j = jRefConf(i)
       If (j > 0 .And. j <= UBound(Hii,1)) Then
          s = Flags%SymVec(j+i0)
          Write(stdout,'(I4,I10,4X,A,A,I1,A,F12.6)')  i, j, SymLabel(s), '(', s, ')', Hii(j)
       End If
    End Do
    Write(stdout,*)

    Return
  End Subroutine PrintRefConf



  !
  ! Print the initial eigenvectors.
  !
  Subroutine PrintEigenVectors(jRefConf, alpha, stdout)
    Implicit None
    Integer,          Dimension(:)   :: jRefConf
    Double Precision, Dimension(:,:) :: alpha
    Integer                          :: stdout
    ! End of dummy parameters.
    Integer                          :: i, j, k, ndim, nblocks, maxcol

    Write(stdout,'(/A)') ' Initial eigenvectors:'
    ndim    = UBound(jRefConf, 1)
    nblocks = (ndim + 9) / 10

    Do k=1, nblocks
       maxcol = Min(10 * k, ndim)
       Write(stdout, '(/1X,"  #     Index",I8,9I11)') (j, j=10*k-9, maxcol)
       Do i=1, ndim
          Write(stdout, '(1X,I3,I10,10F11.6)') i, jRefConf(i), (alpha(i,j), j=10*k-9, maxcol)
       End Do
    End Do
    Write(stdout,*)

    Return
  End Subroutine PrintEigenVectors



  !
  ! Initialize the B matrix (i.e. the matrix containing the orthogonal set
  ! of trial vectors) with the unit vectors corresponding to the (reference)
  ! CI basis functions indexed by iRefConf. i0 is to be subtracted from the
  ! global CI basis function indices to convert them to the indices of the
  ! local CI problem (with the appropriate symmetry). On output, ndav will be
  ! the number of basis vectors for the local CI problem. The program will be
  ! aborted (and this subroutine will not return) if ndav is less than mindav.
  !
  Subroutine InitB(iRefConf, i0, B, ndav, mindav, stdout)
    Implicit None
    Integer,          Dimension(:)   :: iRefConf
    Integer                          :: i0
    Double Precision, Dimension(:,:) :: B
    Integer                          :: ndav, mindav, stdout
    ! End of dummy parameters.
    Integer                          :: i, j

    ndav = 0
    Do j=1, UBound(iRefConf,1)
       i = iRefConf(j) - i0
       If (i > 0 .And. i <= Ubound(B,1)) Then
          ndav = ndav + 1
          B(:,ndav) = 0.D0
          B(i,ndav) = 1.D0
       End If
    End Do

    If (ndav < mindav) Then
       Write(stdout,'(1X,A)') 'Too few appropriate reference configurations.'
       Stop 'davidson::InitB'
    End If

    Return
  End Subroutine InitB



  !
  ! Initialize the matrix A, which contains the projection of the
  ! CI Hamiltonian onto the trial vectors in matrix B, by forming
  ! the product A = BT * P, where P = H * B.
  !
  Subroutine InitA(B, P, A, ndav)
    Implicit None
    Double Precision, Dimension(:,:) :: B, P, A
    Integer                          :: ndav
    ! End of dummy parameters.

    ! Although only the lower triangle of A is needed for TRED2/TQL2,
    ! a single matrix multiplication is probably faster than several
    ! individual matrix-vector operations:
    !
    ! A(1:ndav,1:ndav) = Matmul(Transpose(B(:,1:ndav)), P(:,1:ndav))
    Call DGEMM('T', 'N', ndav, ndav, UBound(B,1), 1.D0, B(1,1), UBound(B,1), P(1,1), UBound(P,1), 0.D0, A(1,1), UBound(A,1))

    Return
  End Subroutine InitA



  !
  ! Diagonalizer interface.
  !
  Subroutine Eigen(A, lambda, alpha, ndav, nroots, stdout, icall)
    Implicit None
    Double Precision, Dimension(:,:)          :: A
    Double Precision, Dimension(:)            :: lambda
    Double Precision, Dimension(:,:)          :: alpha
    Integer                                   :: ndav, nroots, stdout, icall
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Pointer :: A2
    Double Precision, Dimension(:),   Pointer :: work
    Integer,          Dimension(:),   Pointer :: iwork, ifail
    Integer                                   :: lwork, nfound, info, k, NB, ILAENV
    Double Precision                          :: abstol, DLAMCH

    Allocate(A2(ndav,ndav))
    Do k=1,ndav
       A2(k:ndav,k) = A(k:ndav,k)
    End Do
    NB = ILAENV(1, 'DSYEVX', 'VIU', ndav, -1, -1, -1)
    If (NB > 5) Then
       lwork = (NB + 3) * ndav
    Else
       lwork = 8 * ndav
    End If
    Allocate(work(lwork))
    Allocate(iwork(5 * ndav))
    Allocate(ifail(ndav))
    abstol = 2.D0 * DLAMCH('S')
    Call DSYEVX('V', 'I', 'L', ndav, A2(1,1), ndav, 0.D0, 0.D0, 1, nroots, abstol, nfound, &
                lambda, alpha(1,1), UBound(alpha,1), work, lwork, iwork, ifail, info)
    Deallocate(ifail)
    Deallocate(iwork)
    Deallocate(work)
    Deallocate(A2)
    If (info /= 0  .Or.  nfound < nroots) Then
       Write(stdout, '(1X,A)') 'The diagonalization of the submatrix failed.'
       icall = -1
       Return
    End If

    Return
  End Subroutine Eigen



  !
  ! Calculate the q vector and its norm.
  !
  Subroutine FormQVec(B, P, alvec, lambdak, q, qnorm, ndav)
    Implicit None
    Double Precision, Dimension(:,:) :: B, P
    Double Precision, Dimension(:)   :: alvec
    Double Precision                 :: lambdak
    Double Precision, Dimension(:)   :: q
    Double Precision                 :: qnorm
    Integer                          :: ndav
    ! End of dummy parameters.
    Double Precision                 :: qi
    Integer                          :: i, j

    qnorm = 0.D0

    Do i=1, UBound(q,1)
       qi = 0.D0
       Do j=1, ndav
          qi = qi + (P(i,j) - lambdak * B(i,j)) * alvec(j)
       End Do
       q(i) = qi
       qnorm = qnorm + qi * qi
    End Do

    qnorm = Sqrt(qnorm)

    Return
  End Subroutine FormQVec



  !
  ! Calculate the entire Q matrix and the column norms.
  !
  Subroutine FormQMat(B, P, alpha, lambda, Q, qnorm, ndav, nroots)
    Implicit None
    Double Precision, Dimension(:,:)          :: B, P, alpha
    Double Precision, Dimension(:)            :: lambda
    Double Precision, Dimension(:,:)          :: Q
    Double Precision, Dimension(:)            :: qnorm
    Integer                                   :: ndav, nroots
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Pointer :: alpha2
    Integer                                   :: k

    Allocate(alpha2(ndav,nroots))

    Do k=1,nroots
       alpha2(:,k) = lambda(k) * alpha(1:ndav,k)
    End Do

    ! Q = P * alpha
    ! Q = Q - B * alpha2
    Call DGEMM('N', 'N', UBound(Q,1), nroots, ndav,  1.D0, &
         &     P(1,1), UBound(P,1), alpha(1,1),  UBound(alpha,1), 0.D0, Q(1,1), UBound(Q,1))
    Call DGEMM('N', 'N', UBound(Q,1), nroots, ndav, -1.D0, &
         &     B(1,1), UBound(B,1), alpha2(1,1), ndav,            1.D0, Q(1,1), UBound(Q,1))

    Deallocate(alpha2)

    Do k=1,nroots
       qnorm(k) = Sqrt(Dot_Product(Q(:,k), Q(:,k)))
    End Do

    Return
  End Subroutine FormQMat


  !
  ! Add a new trial vector to the B matrix.
  !
  Subroutine UpdateB(Hii, q, lambdak, qtol, B, ndav)
    Implicit None
    Double Precision, Dimension(:)           :: Hii, q
    Double Precision                         :: lambdak, qtol
    Double Precision, Dimension(:,:)         :: B
    Integer                                  :: ndav, new1
    ! End of dummy parameters.
    Double Precision                         :: x
    Integer                                  :: i

    ndav = ndav + 1
    new1 = ndav
    Do i=1, UBound(Hii,1)
       x = lambdak - Hii(i)
       If (Dabs(x) >= qtol) Then
          B(i,ndav) =  q(i) / x
       Else If (x >= 0.D0) Then
          B(i,ndav) =  q(i) / qtol
       Else
          B(i,ndav) = -q(i) / qtol
       End If
    End Do
    x = 1.D0 / Sqrt(Dot_Product(B(:,ndav), B(:,ndav)))
    B(:,ndav) = B(:,ndav) * x

    Call Schmidt(B, new1, ndav, qtol)
    Call Schmidt(B, new1, ndav, qtol)

    Return
  End Subroutine UpdateB


  !
  ! Add set of new trial vectors to the B matrix.
  !
  Subroutine UpdateBMat(Hii, Q, lambda, qnorm, qtol, B, ndav)
    Implicit None
    Double Precision, Dimension(:)            :: Hii
    Double Precision, Dimension(:,:)          :: Q
    Double Precision, Dimension(:)            :: lambda, qnorm
    Double Precision                          :: qtol
    Double Precision, Dimension(:,:)          :: B
    Integer                                   :: ndav
    ! End of dummy parameters.
    Integer                                   :: i, k, ndavold, nci
    Double Precision                          :: x

    nci     = UBound(B,1)
    ndavold = ndav

    Do k=1,UBound(Q,2)
       If (qnorm(k) >= qtol) Then
          ndav = ndav + 1
          Do i=1,nci
             x = lambda(k) - Hii(i)
             If (Dabs(x) >= qtol) Then
                B(i,ndav) =  Q(i,k) / x
             Else If (x >= 0.D0) Then
                B(i,ndav) =  Q(i,k) / qtol
             Else
                B(i,ndav) = -Q(i,k) / qtol
             End If
          End Do
          x = 1.D0 / Sqrt(Dot_Product(B(:,ndav), B(:,ndav)))
          B(:,ndav) = B(:,ndav) * x
       End If
    End Do

    Call Schmidt(B, ndavold+1, ndav, qtol)
    Call Schmidt(B, ndavold+1, ndav, qtol)

    Return
  End Subroutine UpdateBMat



  !
  ! Perform Schmidt orthogonalization of the new basis vectors
  ! stored in B(:,new1:ndav) against the established vectors
  ! in B(:,1:new1-1) and among each other, dropping those new
  ! vectors with a norm less than qtol after projecting out
  ! the previous vectors. On return, ndav contains the total
  ! number of established and surviving new basis vectors.
  !
  Subroutine Schmidt(B, new1, ndav, qtol)
    Implicit None
    Double Precision, Dimension(:,:)          :: B
    Integer                                   :: new1, ndav
    Double Precision                          :: qtol
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Pointer :: BTxQ
    Double Precision, Dimension(:),   Pointer :: proj
    Integer                                   :: nci, ndavold, new, k
    Double Precision                          :: qtol2, x

    If (new1 > ndav) Return

    nci     = UBound(B,1)     ! Length of basis vectors.
    ndavold = new1 - 1        ! Number of established basis vectors.
    new     = ndav - ndavold  ! Number of new basis vectors.
    qtol2   = qtol * qtol

    ! In the following comments, B refers to the matrix of the
    ! established basis vectors and Q to the new basis vectors.

    If (new1 > 1) Then
       ! Projecting out the established basis vectors.
       !   BTxQ = B**T * Q
       !   Q    = Q - B * BTxQ
       Allocate(BTxQ(ndavold, new))
       Call DGEMM('T', 'N', ndavold, new, nci,      1.D0, B(1,1), nci, B(1,new1), nci,     0.D0, BTxQ(1,1), ndavold)
       Call DGEMM('N', 'N', nci,     new, ndavold, -1.D0, B(1,1), nci, BTxQ(1,1), ndavold, 1.D0, B(1,new1), nci)
       Deallocate(BTxQ)
    End If

    Allocate(proj(new - 1))  ! instead of 'ndav - ndavold'
    k = new1  ! Index of the current new basis vector (q).
    Do While (k <= ndav)
       If (k > new1) Then
          ! Projecting out the previous new basis vectors.
          !   proj = Q**T * q
          !   q    = q - Q * proj
          Call DGEMV('T', nci, k-new1,  1.D0, B(1,new1), nci, B(1,k), 1, 0.D0, proj,   1)
          Call DGEMV('N', nci, k-new1, -1.D0, B(1,new1), nci, proj,   1, 1.D0, B(1,k), 1)
       End If
       x = Dot_Product(B(:,k), B(:,k))
       If (x < qtol2) Then
          ! The norm of the current new basis vector after
          ! projecting out all previous vectors is too small.
          ! Drop the current vector and replace it with the
          ! last new basis vector, if any.
          ndav = ndav - 1
          If (ndav < k) Exit
          B(:,k) = B(:,ndav+1)
          Cycle
       End If
       ! Normalize the current vector.
       x = 1.D0 / Sqrt(x)
       B(:,k) = B(:,k) * x
       k = k + 1
    End Do
    Deallocate(proj)

    Return
  End Subroutine Schmidt



  !
  ! Form the product vector p = H * b.
  !
  Subroutine FormPVec(Hij, ifirst, icol, n, b, p)
    Implicit None
    Double Precision, Dimension(:) :: Hij
    Integer,          Dimension(:) :: ifirst, icol
    Integer                        :: n
    Double Precision, Dimension(:) :: b, p
    ! End of dummy parameters.
    Integer                        :: i, j, ij

    p = 0.D0

    Do i=1, n
       ij   = ifirst(i)
       p(i) = p(i) + Hij(ij) * b(i)
       Do ij=ifirst(i)+1, ifirst(i+1)-1
          j = icol(ij)
          p(i) = p(i) + Hij(ij) * b(j)
          p(j) = p(j) + Hij(ij) * b(i)
       End Do
    End Do

    Return
  End Subroutine FormPVec



  !
  ! Form the matrix product P = H * B.
  ! B and P are passed in transposed arrangement.
  !
  Subroutine FormPMat(Hij, ifirst, icol, m, B, P)
    Use gugasparse, Only: dcsrsymm
    Implicit None
    Double Precision, Dimension(:)   :: Hij
    Integer,          Dimension(:)   :: ifirst, icol
    Integer                          :: m
    Double Precision, Dimension(:,:) :: B, P
    ! End of dummy parameters.
    ! Double Precision                 :: tuser1, tsys1, twall1
    ! Double Precision                 :: tuser2, tsys2, twall2

    ! Write(Stdout,'(1X,"H: nnz=",I10)') ifirst(m+1)-1
    ! Write(Stdout,'(1X,"B: ",I8," x ",I8)') m, UBound(B,2)
    ! Call GetTime(tuser1, tsys1, twall1)
    Call dcsrsymm(m, UBound(B,2), Hij, ifirst, icol, B, P)
    ! Call GetTime(tuser2, tsys2, twall2)
    ! Write(Stdout,'(1X,"twall=",F6.2," s")') twall2 - twall1

    Return
  End Subroutine FormPMat



  !
  ! Driver routine for the formation of the matrix product P = H * B.
  !
  Subroutine CalcP(Flags, P, B, Hii, Hij, ifirst, icol, MOInt1, MOInt2, Occup, Irrep, i0, cidir, ndav, stdout)
    Implicit None
    Type(ShavittControl)                      :: Flags
    Double Precision, Dimension(:,:)          :: P, B
    Double Precision, Dimension(:)            :: Hii, Hij
    Integer,          Dimension(:)            :: ifirst, icol
    Double Precision, Dimension(:,:)          :: MOInt1, MOInt2
    Double Precision, Dimension(:)            :: Occup
    Integer                                   :: Irrep, i0, cidir, ndav, stdout
    ! End of dummy parameters.
    Double Precision, Dimension(:,:), Pointer :: Btrans, Ptrans
    Integer                                   :: N, k

    N = UBound(B,1)

    Select Case (cidir)
    Case (-2,-1,1,2)
       Call FormPMat(Hij, ifirst, icol, N, B(:,1:ndav), P(:,1:ndav))
    Case (3)
       Allocate(Btrans(ndav,N))
       Allocate(Ptrans(ndav,N))
       Do k=1,ndav
          Btrans(k,:) = B(:,k)
       End Do
       If (Irrep == 0) Then
          Call ShapeDrivenProdMat(Flags, Ptrans, Btrans, Hii, MOInt1, MOInt2)
       Else
          Call ShapeDrivenProdMatIrrep(Flags, Ptrans, Btrans, Hii, MOInt1, MOInt2, Irrep, i0)
       End If
       Do k=1,ndav
          P(:,k) = Ptrans(k,:)
       End Do
       Deallocate(Ptrans)
       Deallocate(Btrans)
    Case (-3)
       Allocate(Btrans(ndav,N))
       Allocate(Ptrans(ndav,N))
       Do k=1,ndav
          Btrans(k,:) = B(:,k)
       End Do
       Call CalcProdMat(Flags, Ptrans, Btrans, MOInt1, MOInt2, Occup, Irrep, i0)
       Do k=1,ndav
          P(:,k) = Ptrans(k,:)
       End Do
       Deallocate(Ptrans)
       Deallocate(Btrans)
    Case Default
       Write(stdout,'(1X,A,I2,A)') 'Unexpected value for cidir (', cidir, ').'
       Stop 'davidson::CalcP'
    End Select

    Return
  End Subroutine CalcP



  !
  ! Add the projection of the new trial vector to the A matrix.
  !
  Subroutine UpdateA(B, P, A, ndav)
    Implicit None
    Double Precision, Dimension(:,:) :: B, P, A
    Integer                          :: ndav
    ! End of dummy parameters.

    ! New row is added to A:
    Call DGEMV('T', UBound(B,1), ndav, 1.D0, B, UBound(B,1), P(1,ndav), 1, 0.D0, A(ndav,1), UBound(A,1))

    Return
  End Subroutine UpdateA


  !
  ! Add the projection of the new trial vectors to the A matrix.
  !
  Subroutine UpdateAMat(B, P, A, ndavold, ndav)
    Implicit None
    Double Precision, Dimension(:,:) :: B, P, A
    Integer                          :: ndavold, ndav
    ! End of dummy parameters.
    Integer                          :: new, nci

    new = ndav - ndavold
    nci = UBound(B,1)
    ! New part of lower triangle of A is computed:
    Call DGEMM('T', 'N', new, ndav, nci, 1.D0, P(1,ndavold+1), nci, B(1,1), nci, 0.D0, A(ndavold+1,1), UBound(A,1))

    Return
  End Subroutine UpdateAMat


  !
  ! Find the vector in alpha which overlaps most strongly
  ! with the starting vector alpha0.
  !
  Subroutine FindK(alpha, alpha0, ndav, nref, k, overlap)
    Implicit None
    Double Precision, Dimension(:,:)   :: alpha
    Double Precision, Dimension(:)     :: alpha0
    Integer                            :: ndav, nref, k
    Double Precision                   :: overlap
    ! End of dummy parameters.
    Double Precision, Dimension(ndav)  :: overlapvec
    Integer                            :: i

    ! overlapvec = Matmul(Transpose(alpha(1:nref,1:ndav)), alpha0(1:nref))
    Call DGEMV('T', nref, ndav, 1.D0, alpha, UBound(alpha,1), alpha0, 1, 0.D0, overlapvec, 1)

    k = 1
    Do i=2, ndav
       If (Dabs(overlapvec(i)) > Dabs(overlapvec(k)))  k = i
    End Do
    overlap = overlapvec(k)

    Return
  End Subroutine FindK


End Module davidson
