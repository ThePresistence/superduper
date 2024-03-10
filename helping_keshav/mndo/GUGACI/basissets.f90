Module basissets

Implicit None

Private

Public :: Shell, BasisSet, InitSTO6G, InitECP3G, Renormalize, PrintBasisSet, FreeBasisSet


Type Shell
  Integer,          Dimension(:),   Pointer :: l, iao  ! Dimension(N)
  Double Precision, Dimension(:),   Pointer :: alpha   ! Dimension(K)
  Double Precision, Dimension(:,:), Pointer :: D       ! Dimension(K,N)
  Integer                                   :: K, N, lmin, lmax
End Type Shell


Type BasisSet
  Type(Shell), Dimension(:), Pointer :: shell
  Integer                            :: shells
End Type BasisSet


Contains

  Subroutine InitSTO6G(bs, zs, zp, zd, ns, np, nd)
    Implicit None
    Type(BasisSet),   Dimension(:), Pointer :: bs
    Double Precision, Dimension(:)          :: zs, zp, zd
    Integer,          Dimension(:)          :: ns, np, nd
    ! End of dummy parameters.
    Double Precision, Dimension(6,4)        :: EXX, CS, CP, EXXD, CD, EXXF, CF
    Integer                                 :: LMZ, Z

    ! Principle quantum number 1:
    Data EXX (:,1) /  2.310303149D+01,  4.235915534D+00,  1.185056519D+00,  4.070988982D-01,  1.580884151D-01,  6.510953954D-02 /
    Data CS  (:,1) /  9.163596280D-03,  4.936149294D-02,  1.685383049D-01,  3.705627997D-01,  4.164915298D-01,  1.303340841D-01 /
    Data CP  (:,1) /  6*0.D0 /
    Data EXXD(:,1) /  6*0.D0 /
    Data CD  (:,1) /  6*0.D0 /
    Data EXXF(:,1) /  6*0.D0 /
    Data CF  (:,1) /  6*0.D0 /

    ! Principle quantum number 2:
    Data EXX (:,2) /  1.030869372D+01,  2.040359519D+00,  6.341422177D-01,  2.439773685D-01,  1.059595374D-01,  4.856900860D-02 /
    Data CS  (:,2) / -1.325278809D-02, -4.699171014D-02, -3.378537151D-02,  2.502417861D-01,  5.951172526D-01,  2.407061763D-01 /
    Data CP  (:,2) /  3.759696623D-03,  3.767936984D-02,  1.738967435D-01,  4.180364347D-01,  4.258595477D-01,  1.017082955D-01 /
    Data EXXD(:,2) /  6*0.D0 /
    Data CD  (:,2) /  6*0.D0 /
    Data EXXF(:,2) /  6*0.D0 /
    Data CF  (:,2) /  6*0.D0 /

    ! Principle quantum number 3:
    Data EXX (:,3) /  3.080165240D+00,  8.248959202D-01,  3.093447349D-01,  1.384683897D-01,  6.852094951D-02,  3.531333690D-02 /
    Data CS  (:,3) / -7.943126362D-03, -7.100264172D-02, -1.785026925D-01,  1.510635058D-01,  7.354914767D-01,  2.760593123D-01 /
    Data CP  (:,3) / -7.139358907D-03, -1.829277070D-02,  7.621621428D-02,  4.145098597D-01,  4.889621471D-01,  1.058816521D-01 /
    Data EXXD(:,3) /  2.488296923D+00,  7.981487853D-01,  3.311327490D-01,  1.559114463D-01,  7.877734732D-02,  4.058484363D-02 /
    Data CD  (:,3) /  7.283828112D-03,  5.386799363D-02,  2.072139149D-01,  4.266269092D-01,  3.843100204D-01,  8.902827546D-02 /
    Data EXXF(:,3) /  6*0.D0 /
    Data CF  (:,3) /  6*0.D0 /

    ! Principle quantum number 4:
    Data EXX (:,4) /  3.080165240D+00,  8.248959202D-01,  3.093447349D-01,  1.384683897D-01,  6.852094951D-02,  3.531333690D-02 /
    Data CS  (:,4) / -7.943126362D-03, -7.100264172D-02, -1.785026925D-01,  1.510635058D-01,  7.354914767D-01,  2.760593123D-01 /
    Data CP  (:,4) / -7.139358907D-03, -1.829277070D-02,  7.621621428D-02,  4.145098597D-01,  4.889621471D-01,  1.058816521D-01 /
    Data EXXD(:,4) /  4.634239420D+00,  1.341648295D+00,  2.209593028D-01,  1.101467943D-01,  5.904190370D-02,  3.232628887D-02 /
    Data CD  (:,4) / -4.749842876D-04, -3.566777891D-03,  1.108670481D-01,  4.159646930D-01,  4.621672517D-01,  1.081250196D-01 /
    Data EXXF(:,4) /  1.357718039D+00,  5.004907278D-01,  2.296565064D-01,  1.173146814D-01,  6.350097171D-02,  3.474556673D-02 /
    Data CF  (:,4) /  6.930234381D-03,  5.634263745D-02,  2.217065797D-01,  4.411388883D-01,  3.688112625D-01,  7.787514504D-02 /

    LMZ = Size(zs)

    Allocate(bs(LMZ))

    Do Z=1, LMZ
      If (ns(Z) == 0) Then
         ! Element not used.
         bs(Z)%shells = 0
      Else If (np(Z) == 0) Then
         ! Create s shell.
         bs(Z)%shells = 1
         Allocate(bs(Z)%shell(1))
         Call SShell(bs(Z)%shell(1), EXX(:,ns(Z)), CS(:,ns(Z)), zs(Z))
      Else If (nd(Z) == 0) Then
         If (ns(Z) == np(Z) .And. zs(Z) == zp(Z)) Then
            ! Create sp shell.
            bs(Z)%shells = 1
            Allocate(bs(Z)%shell(1))
            Call SPShell(bs(Z)%shell(1), EXX(:,ns(Z)), CS(:,ns(Z)), CP(:,ns(Z)), zs(Z))
         Else
            ! Create s and p shell.
            bs(Z)%shells = 2
            Allocate(bs(Z)%shell(2))
            Call SShell(bs(Z)%shell(1), EXX(:,ns(Z)), CS(:,ns(Z)), zs(Z))
            Call PShell(bs(Z)%shell(2), EXX(:,np(Z)), CP(:,np(Z)), zp(Z))
         End If
      Else
         If (ns(Z) == np(Z) .And. zs(Z) == zp(Z)) Then
            ! Create sp and d shell.
            bs(Z)%shells = 2
            Allocate(bs(Z)%shell(2))
            Call SPShell(bs(Z)%shell(1), EXX (:,ns(Z)), CS(:,ns(Z)), CP(:,ns(Z)), zs(Z))
            Call DShell (bs(Z)%shell(2), EXXD(:,nd(Z)), CD(:,nd(Z)),              zd(Z))
         Else
            ! Create s, p, and d shell.
            bs(Z)%shells = 3
            Allocate(bs(Z)%shell(3))
            Call SShell(bs(Z)%shell(1), EXX (:,ns(Z)), CS(:,ns(Z)), zs(Z))
            Call PShell(bs(Z)%shell(2), EXX (:,np(Z)), CP(:,np(Z)), zp(Z))
            Call DShell(bs(Z)%shell(3), EXXD(:,nd(Z)), CD(:,nd(Z)), zd(Z))
         End If
      End If
    End Do
  End Subroutine InitSTO6G



  Subroutine InitECP3G(bs, zsp, nsp)
    Implicit None
    Type(BasisSet),   Dimension(:), Pointer :: bs
    Double Precision, Dimension(:)          :: zsp
    Integer,          Dimension(:)          :: nsp
    ! End of dummy parameters.
    Integer                                 :: LMZ, Z, i

    LMZ = Size(zsp)

    Allocate(bs(LMZ))

    ! Indicate that the basis set for
    ! element Z is not initialized:
    Do Z=1, LMZ
      bs(Z)%shells = 0
    EndDo

    ! H
    !----
    If (nsp(1) > 0) Then
      bs(1)%shells = 1
      Allocate(bs(1)%shell(1))
      !
      bs(1)%shell(1)%K    = 3
      bs(1)%shell(1)%N    = 1
      bs(1)%shell(1)%lmin = 0
      bs(1)%shell(1)%lmax = 0
      Allocate(bs(1)%shell(1)%l(1))
      Allocate(bs(1)%shell(1)%iao(1))
      Allocate(bs(1)%shell(1)%alpha(3))
      Allocate(bs(1)%shell(1)%D(3,1))
      bs(1)%shell(1)%alpha(1) = 2.227660584D+00
      bs(1)%shell(1)%alpha(2) = 4.057711562D-01
      bs(1)%shell(1)%alpha(3) = 1.098175104D-01
      bs(1)%shell(1)%l(1) = 0
      bs(1)%shell(1)%iao(1) = 0
      bs(1)%shell(1)%D(1,1) =  1.543289673D-01
      bs(1)%shell(1)%D(2,1) =  5.353281423D-01
      bs(1)%shell(1)%D(3,1) =  4.446345422D-01
    End If

    ! Li
    !----
    If (nsp(3) > 0) Then
      bs(3)%shells = 1
      Allocate(bs(3)%shell(1))
      !
      bs(3)%shell(1)%K    = 3
      bs(3)%shell(1)%N    = 2
      bs(3)%shell(1)%lmin = 0
      bs(3)%shell(1)%lmax = 1
      Allocate(bs(3)%shell(1)%l(2))
      Allocate(bs(3)%shell(1)%iao(2))
      Allocate(bs(3)%shell(1)%alpha(3))
      Allocate(bs(3)%shell(1)%D(3,2))
      bs(3)%shell(1)%alpha(1) = 4.422700000D-01
      bs(3)%shell(1)%alpha(2) = 7.847000000D-02
      bs(3)%shell(1)%alpha(3) = 2.434000000D-02
      bs(3)%shell(1)%l(1) = 0
      bs(3)%shell(1)%iao(1) = 0
      bs(3)%shell(1)%D(1,1) = -1.880500000D-01
      bs(3)%shell(1)%D(2,1) =  6.655200000D-01
      bs(3)%shell(1)%D(3,1) =  4.776200000D-01
      bs(3)%shell(1)%l(2) = 1
      bs(3)%shell(1)%iao(2) = 1
      bs(3)%shell(1)%D(1,2) =  1.167100000D-01
      bs(3)%shell(1)%D(2,2) =  5.262800000D-01
      bs(3)%shell(1)%D(3,2) =  5.228200000D-01
    End If

    ! Be
    !----
    If (nsp(4) > 0) Then
      bs(4)%shells = 1
      Allocate(bs(4)%shell(1))
      !
      bs(4)%shell(1)%K    = 3
      bs(4)%shell(1)%N    = 2
      bs(4)%shell(1)%lmin = 0
      bs(4)%shell(1)%lmax = 1
      Allocate(bs(4)%shell(1)%l(2))
      Allocate(bs(4)%shell(1)%iao(2))
      Allocate(bs(4)%shell(1)%alpha(3))
      Allocate(bs(4)%shell(1)%D(3,2))
      bs(4)%shell(1)%alpha(1) = 1.060390000D+00
      bs(4)%shell(1)%alpha(2) = 2.075800000D-01
      bs(4)%shell(1)%alpha(3) = 5.933000000D-02
      bs(4)%shell(1)%l(1) = 0
      bs(4)%shell(1)%iao(1) = 0
      bs(4)%shell(1)%D(1,1) = -1.852500000D-01
      bs(4)%shell(1)%D(2,1) =  5.241800000D-01
      bs(4)%shell(1)%D(3,1) =  6.232000000D-01
      bs(4)%shell(1)%l(2) = 1
      bs(4)%shell(1)%iao(2) = 1
      bs(4)%shell(1)%D(1,2) =  1.512100000D-01
      bs(4)%shell(1)%D(2,2) =  5.594000000D-01
      bs(4)%shell(1)%D(3,2) =  4.776500000D-01
    End If

    ! B
    !----
    If (nsp(5) > 0) Then
      bs(5)%shells = 1
      Allocate(bs(5)%shell(1))
      !
      bs(5)%shell(1)%K    = 3
      bs(5)%shell(1)%N    = 2
      bs(5)%shell(1)%lmin = 0
      bs(5)%shell(1)%lmax = 1
      Allocate(bs(5)%shell(1)%l(2))
      Allocate(bs(5)%shell(1)%iao(2))
      Allocate(bs(5)%shell(1)%alpha(3))
      Allocate(bs(5)%shell(1)%D(3,2))
      bs(5)%shell(1)%alpha(1) = 1.724270000D+00
      bs(5)%shell(1)%alpha(2) = 3.500900000D-01
      bs(5)%shell(1)%alpha(3) = 9.394000000D-02
      bs(5)%shell(1)%l(1) = 0
      bs(5)%shell(1)%iao(1) = 0
      bs(5)%shell(1)%D(1,1) = -1.935900000D-01
      bs(5)%shell(1)%D(2,1) =  6.149700000D-01
      bs(5)%shell(1)%D(3,1) =  5.497600000D-01
      bs(5)%shell(1)%l(2) = 1
      bs(5)%shell(1)%iao(2) = 1
      bs(5)%shell(1)%D(1,2) =  1.829200000D-01
      bs(5)%shell(1)%D(2,2) =  5.443700000D-01
      bs(5)%shell(1)%D(3,2) =  4.830300000D-01
    End If

    ! C
    !----
    If (nsp(6) > 0) Then
      bs(6)%shells = 1
      Allocate(bs(6)%shell(1))
      !
      bs(6)%shell(1)%K    = 3
      bs(6)%shell(1)%N    = 2
      bs(6)%shell(1)%lmin = 0
      bs(6)%shell(1)%lmax = 1
      Allocate(bs(6)%shell(1)%l(2))
      Allocate(bs(6)%shell(1)%iao(2))
      Allocate(bs(6)%shell(1)%alpha(3))
      Allocate(bs(6)%shell(1)%D(3,2))
      bs(6)%shell(1)%alpha(1) = 2.644860000D+00
      bs(6)%shell(1)%alpha(2) = 5.421500000D-01
      bs(6)%shell(1)%alpha(3) = 1.446600000D-01
      bs(6)%shell(1)%l(1) = 0
      bs(6)%shell(1)%iao(1) = 0
      bs(6)%shell(1)%D(1,1) = -1.918800000D-01
      bs(6)%shell(1)%D(2,1) =  6.162800000D-01
      bs(6)%shell(1)%D(3,1) =  5.489600000D-01
      bs(6)%shell(1)%l(2) = 1
      bs(6)%shell(1)%iao(2) = 1
      bs(6)%shell(1)%D(1,2) =  2.025900000D-01
      bs(6)%shell(1)%D(2,2) =  5.583000000D-01
      bs(6)%shell(1)%D(3,2) =  4.551400000D-01
    End If

    ! N
    !----
    If (nsp(7) > 0) Then
      bs(7)%shells = 1
      Allocate(bs(7)%shell(1))
      !
      bs(7)%shell(1)%K    = 3
      bs(7)%shell(1)%N    = 2
      bs(7)%shell(1)%lmin = 0
      bs(7)%shell(1)%lmax = 1
      Allocate(bs(7)%shell(1)%l(2))
      Allocate(bs(7)%shell(1)%iao(2))
      Allocate(bs(7)%shell(1)%alpha(3))
      Allocate(bs(7)%shell(1)%D(3,2))
      bs(7)%shell(1)%alpha(1) = 3.688490000D+00
      bs(7)%shell(1)%alpha(2) = 7.753400000D-01
      bs(7)%shell(1)%alpha(3) = 2.049800000D-01
      bs(7)%shell(1)%l(1) = 0
      bs(7)%shell(1)%iao(1) = 0
      bs(7)%shell(1)%D(1,1) = -1.926900000D-01
      bs(7)%shell(1)%D(2,1) =  6.188800000D-01
      bs(7)%shell(1)%D(3,1) =  5.492600000D-01
      bs(7)%shell(1)%l(2) = 1
      bs(7)%shell(1)%iao(2) = 1
      bs(7)%shell(1)%D(1,2) =  2.228100000D-01
      bs(7)%shell(1)%D(2,2) =  5.603200000D-01
      bs(7)%shell(1)%D(3,2) =  4.385900000D-01
    End If

    ! O
    !----
    If (nsp(8) > 0) Then
      bs(8)%shells = 1
      Allocate(bs(8)%shell(1))
      !
      bs(8)%shell(1)%K    = 3
      bs(8)%shell(1)%N    = 2
      bs(8)%shell(1)%lmin = 0
      bs(8)%shell(1)%lmax = 1
      Allocate(bs(8)%shell(1)%l(2))
      Allocate(bs(8)%shell(1)%iao(2))
      Allocate(bs(8)%shell(1)%alpha(3))
      Allocate(bs(8)%shell(1)%D(3,2))
      bs(8)%shell(1)%alpha(1) = 4.784990000D+00
      bs(8)%shell(1)%alpha(2) = 9.986000000D-01
      bs(8)%shell(1)%alpha(3) = 2.568700000D-01
      bs(8)%shell(1)%l(1) = 0
      bs(8)%shell(1)%iao(1) = 0
      bs(8)%shell(1)%D(1,1) = -1.924800000D-01
      bs(8)%shell(1)%D(2,1) =  6.695200000D-01
      bs(8)%shell(1)%D(3,1) =  5.027000000D-01
      bs(8)%shell(1)%l(2) = 1
      bs(8)%shell(1)%iao(2) = 1
      bs(8)%shell(1)%D(1,2) =  2.415800000D-01
      bs(8)%shell(1)%D(2,2) =  5.589000000D-01
      bs(8)%shell(1)%D(3,2) =  4.316000000D-01
    End If

    ! F
    !----
    If (nsp(9) > 0) Then
      bs(9)%shells = 1
      Allocate(bs(9)%shell(1))
      !
      bs(9)%shell(1)%K    = 3
      bs(9)%shell(1)%N    = 2
      bs(9)%shell(1)%lmin = 0
      bs(9)%shell(1)%lmax = 1
      Allocate(bs(9)%shell(1)%l(2))
      Allocate(bs(9)%shell(1)%iao(2))
      Allocate(bs(9)%shell(1)%alpha(3))
      Allocate(bs(9)%shell(1)%D(3,2))
      bs(9)%shell(1)%alpha(1) = 6.017830000D+00
      bs(9)%shell(1)%alpha(2) = 1.253150000D+00
      bs(9)%shell(1)%alpha(3) = 3.176000000D-01
      bs(9)%shell(1)%l(1) = 0
      bs(9)%shell(1)%iao(1) = 0
      bs(9)%shell(1)%D(1,1) = -1.885000000D-01
      bs(9)%shell(1)%D(2,1) =  6.980000000D-01
      bs(9)%shell(1)%D(3,1) =  4.742700000D-01
      bs(9)%shell(1)%l(2) = 1
      bs(9)%shell(1)%iao(2) = 1
      bs(9)%shell(1)%D(1,2) =  2.566700000D-01
      bs(9)%shell(1)%D(2,2) =  5.601300000D-01
      bs(9)%shell(1)%D(3,2) =  4.213900000D-01
    End If

    ! Ne
    !----
    If (nsp(10) > 0) Then
      bs(10)%shells = 1
      Allocate(bs(10)%shell(1))
      !
      bs(10)%shell(1)%K    = 3
      bs(10)%shell(1)%N    = 2
      bs(10)%shell(1)%lmin = 0
      bs(10)%shell(1)%lmax = 1
      Allocate(bs(10)%shell(1)%l(2))
      Allocate(bs(10)%shell(1)%iao(2))
      Allocate(bs(10)%shell(1)%alpha(3))
      Allocate(bs(10)%shell(1)%D(3,2))
      bs(10)%shell(1)%alpha(1) = 7.478310000D+00
      bs(10)%shell(1)%alpha(2) = 1.554880000D+00
      bs(10)%shell(1)%alpha(3) = 3.905700000D-01
      bs(10)%shell(1)%l(1) = 0
      bs(10)%shell(1)%iao(1) = 0
      bs(10)%shell(1)%D(1,1) = -1.869200000D-01
      bs(10)%shell(1)%D(2,1) =  7.098800000D-01
      bs(10)%shell(1)%D(3,1) =  4.624800000D-01
      bs(10)%shell(1)%l(2) = 1
      bs(10)%shell(1)%iao(2) = 1
      bs(10)%shell(1)%D(1,2) =  2.658600000D-01
      bs(10)%shell(1)%D(2,2) =  5.613400000D-01
      bs(10)%shell(1)%D(3,2) =  4.143900000D-01
    End If

    ! Original ECP-4G basis used for second-row elements.

    ! Na
    !----
    If (nsp(11) > 0) Then
      bs(11)%shells = 1
      Allocate(bs(11)%shell(1))
      !
      bs(11)%shell(1)%K    = 4
      bs(11)%shell(1)%N    = 2
      bs(11)%shell(1)%lmin = 0
      bs(11)%shell(1)%lmax = 1
      Allocate(bs(11)%shell(1)%l(2))
      Allocate(bs(11)%shell(1)%iao(2))
      Allocate(bs(11)%shell(1)%alpha(4))
      Allocate(bs(11)%shell(1)%D(4,2))
      bs(11)%shell(1)%alpha(1) = 0.42990D+00
      bs(11)%shell(1)%alpha(2) = 0.08897D+00
      bs(11)%shell(1)%alpha(3) = 0.03550D+00
      bs(11)%shell(1)%alpha(4) = 0.01455D+00
      bs(11)%shell(1)%l(1) = 0
      bs(11)%shell(1)%iao(1) = 0
      bs(11)%shell(1)%D(1,1) = -0.20874D+00
      bs(11)%shell(1)%D(2,1) =  0.31206D+00
      bs(11)%shell(1)%D(3,1) =  0.70300D+00
      bs(11)%shell(1)%D(4,1) =  0.11648D+00
      bs(11)%shell(1)%l(2) = 1
      bs(11)%shell(1)%iao(2) = 1
      bs(11)%shell(1)%D(1,2) = -0.02571D+00
      bs(11)%shell(1)%D(2,2) =  0.21608D+00
      bs(11)%shell(1)%D(3,2) =  0.54196D+00
      bs(11)%shell(1)%D(4,2) =  0.35484D+00
    End If

    ! Mg
    !----
    If (nsp(12) > 0) Then
      bs(12)%shells = 1
      Allocate(bs(12)%shell(1))
      !
      bs(12)%shell(1)%K    = 4
      bs(12)%shell(1)%N    = 2
      bs(12)%shell(1)%lmin = 0
      bs(12)%shell(1)%lmax = 1
      Allocate(bs(12)%shell(1)%l(2))
      Allocate(bs(12)%shell(1)%iao(2))
      Allocate(bs(12)%shell(1)%alpha(4))
      Allocate(bs(12)%shell(1)%D(4,2))
      bs(12)%shell(1)%alpha(1) = 0.66060D+00
      bs(12)%shell(1)%alpha(2) = 0.18450D+00
      bs(12)%shell(1)%alpha(3) = 0.06983D+00
      bs(12)%shell(1)%alpha(4) = 0.02740D+00
      bs(12)%shell(1)%l(1) = 0
      bs(12)%shell(1)%iao(1) = 0
      bs(12)%shell(1)%D(1,1) = -0.24451D+00
      bs(12)%shell(1)%D(2,1) =  0.25323D+00
      bs(12)%shell(1)%D(3,1) =  0.69720D+00
      bs(12)%shell(1)%D(4,1) =  0.21655D+00
      bs(12)%shell(1)%l(2) = 1
      bs(12)%shell(1)%iao(2) = 1
      bs(12)%shell(1)%D(1,2) = -0.04421D+00
      bs(12)%shell(1)%D(2,2) =  0.27323D+00
      bs(12)%shell(1)%D(3,2) =  0.57626D+00
      bs(12)%shell(1)%D(4,2) =  0.28152D+00
    End If

    ! Al
    !----
    If (nsp(13) > 0) Then
      bs(13)%shells = 1
      Allocate(bs(13)%shell(1))
      !
      bs(13)%shell(1)%K    = 4
      bs(13)%shell(1)%N    = 2
      bs(13)%shell(1)%lmin = 0
      bs(13)%shell(1)%lmax = 1
      Allocate(bs(13)%shell(1)%l(2))
      Allocate(bs(13)%shell(1)%iao(2))
      Allocate(bs(13)%shell(1)%alpha(4))
      Allocate(bs(13)%shell(1)%D(4,2))
      bs(13)%shell(1)%alpha(1) = 0.90110D+00
      bs(13)%shell(1)%alpha(2) = 0.44950D+00
      bs(13)%shell(1)%alpha(3) = 0.14050D+00
      bs(13)%shell(1)%alpha(4) = 0.04874D+00
      bs(13)%shell(1)%l(1) = 0
      bs(13)%shell(1)%iao(1) = 0
      bs(13)%shell(1)%D(1,1) = -0.30377D+00
      bs(13)%shell(1)%D(2,1) =  0.13382D+00
      bs(13)%shell(1)%D(3,1) =  0.76037D+00
      bs(13)%shell(1)%D(4,1) =  0.32232D+00
      bs(13)%shell(1)%l(2) = 1
      bs(13)%shell(1)%iao(2) = 1
      bs(13)%shell(1)%D(1,2) = -0.07929D+00
      bs(13)%shell(1)%D(2,2) =  0.16540D+00
      bs(13)%shell(1)%D(3,2) =  0.53015D+00
      bs(13)%shell(1)%D(4,2) =  0.47724D+00
    End If

    ! Si
    !----
    If (nsp(14) > 0) Then
      bs(14)%shells = 1
      Allocate(bs(14)%shell(1))
      !
      bs(14)%shell(1)%K    = 4
      bs(14)%shell(1)%N    = 2
      bs(14)%shell(1)%lmin = 0
      bs(14)%shell(1)%lmax = 1
      Allocate(bs(14)%shell(1)%l(2))
      Allocate(bs(14)%shell(1)%iao(2))
      Allocate(bs(14)%shell(1)%alpha(4))
      Allocate(bs(14)%shell(1)%D(4,2))
      bs(14)%shell(1)%alpha(1) = 1.16700D+00
      bs(14)%shell(1)%alpha(2) = 0.52680D+00
      bs(14)%shell(1)%alpha(3) = 0.18070D+00
      bs(14)%shell(1)%alpha(4) = 0.06480D+00
      bs(14)%shell(1)%l(1) = 0
      bs(14)%shell(1)%iao(1) = 0
      bs(14)%shell(1)%D(1,1) = -0.32403D+00
      bs(14)%shell(1)%D(2,1) =  0.18438D+00
      bs(14)%shell(1)%D(3,1) =  0.77737D+00
      bs(14)%shell(1)%D(4,1) =  0.26767D+00
      bs(14)%shell(1)%l(2) = 1
      bs(14)%shell(1)%iao(2) = 1
      bs(14)%shell(1)%D(1,2) = -0.08450D+00
      bs(14)%shell(1)%D(2,2) =  0.23786D+00
      bs(14)%shell(1)%D(3,2) =  0.56532D+00
      bs(14)%shell(1)%D(4,2) =  0.37433D+00
    End If

    ! P
    !----
    If (nsp(15) > 0) Then
      bs(15)%shells = 1
      Allocate(bs(15)%shell(1))
      !
      bs(15)%shell(1)%K    = 4
      bs(15)%shell(1)%N    = 2
      bs(15)%shell(1)%lmin = 0
      bs(15)%shell(1)%lmax = 1
      Allocate(bs(15)%shell(1)%l(2))
      Allocate(bs(15)%shell(1)%iao(2))
      Allocate(bs(15)%shell(1)%alpha(4))
      Allocate(bs(15)%shell(1)%D(4,2))
      bs(15)%shell(1)%alpha(1) = 1.45900D+00
      bs(15)%shell(1)%alpha(2) = 0.65490D+00
      bs(15)%shell(1)%alpha(3) = 0.22560D+00
      bs(15)%shell(1)%alpha(4) = 0.08115D+00
      bs(15)%shell(1)%l(1) = 0
      bs(15)%shell(1)%iao(1) = 0
      bs(15)%shell(1)%D(1,1) = -0.34091D+00
      bs(15)%shell(1)%D(2,1) =  0.21535D+00
      bs(15)%shell(1)%D(3,1) =  0.79578D+00
      bs(15)%shell(1)%D(4,1) =  0.23092D+00
      bs(15)%shell(1)%l(2) = 1
      bs(15)%shell(1)%iao(2) = 1
      bs(15)%shell(1)%D(1,2) = -0.09378D+00
      bs(15)%shell(1)%D(2,2) =  0.29205D+00
      bs(15)%shell(1)%D(3,2) =  0.58688D+00
      bs(15)%shell(1)%D(4,2) =  0.30631D+00
    End If

    ! S
    !----
    If (nsp(16) > 0) Then
      bs(16)%shells = 1
      Allocate(bs(16)%shell(1))
      !
      bs(16)%shell(1)%K    = 4
      bs(16)%shell(1)%N    = 2
      bs(16)%shell(1)%lmin = 0
      bs(16)%shell(1)%lmax = 1
      Allocate(bs(16)%shell(1)%l(2))
      Allocate(bs(16)%shell(1)%iao(2))
      Allocate(bs(16)%shell(1)%alpha(4))
      Allocate(bs(16)%shell(1)%D(4,2))
      bs(16)%shell(1)%alpha(1) = 1.81700D+00
      bs(16)%shell(1)%alpha(2) = 0.83790D+00
      bs(16)%shell(1)%alpha(3) = 0.28540D+00
      bs(16)%shell(1)%alpha(4) = 0.09939D+00
      bs(16)%shell(1)%l(1) = 0
      bs(16)%shell(1)%iao(1) = 0
      bs(16)%shell(1)%D(1,1) = -0.34015D+00
      bs(16)%shell(1)%D(2,1) =  0.19601D+00
      bs(16)%shell(1)%D(3,1) =  0.82666D+00
      bs(16)%shell(1)%D(4,1) =  0.21652D+00
      bs(16)%shell(1)%l(2) = 1
      bs(16)%shell(1)%iao(2) = 1
      bs(16)%shell(1)%D(1,2) = -0.10096D+00
      bs(16)%shell(1)%D(2,2) =  0.31244D+00
      bs(16)%shell(1)%D(3,2) =  0.57906D+00
      bs(16)%shell(1)%D(4,2) =  0.30748D+00
    End If

    ! Cl
    !----
    If (nsp(17) > 0) Then
      bs(17)%shells = 1
      Allocate(bs(17)%shell(1))
      !
      bs(17)%shell(1)%K    = 4
      bs(17)%shell(1)%N    = 2
      bs(17)%shell(1)%lmin = 0
      bs(17)%shell(1)%lmax = 1
      Allocate(bs(17)%shell(1)%l(2))
      Allocate(bs(17)%shell(1)%iao(2))
      Allocate(bs(17)%shell(1)%alpha(4))
      Allocate(bs(17)%shell(1)%D(4,2))
      bs(17)%shell(1)%alpha(1) = 2.22500D+00
      bs(17)%shell(1)%alpha(2) = 1.17300D+00
      bs(17)%shell(1)%alpha(3) = 0.38510D+00
      bs(17)%shell(1)%alpha(4) = 0.13010D+00
      bs(17)%shell(1)%l(1) = 0
      bs(17)%shell(1)%iao(1) = 0
      bs(17)%shell(1)%D(1,1) = -0.33098D+00
      bs(17)%shell(1)%D(2,1) =  0.11528D+00
      bs(17)%shell(1)%D(3,1) =  0.84717D+00
      bs(17)%shell(1)%D(4,1) =  0.26534D+00
      bs(17)%shell(1)%l(2) = 1
      bs(17)%shell(1)%iao(2) = 1
      bs(17)%shell(1)%D(1,2) = -0.12604D+00
      bs(17)%shell(1)%D(2,2) =  0.29952D+00
      bs(17)%shell(1)%D(3,2) =  0.58357D+00
      bs(17)%shell(1)%D(4,2) =  0.34097D+00
    End If

    ! Ar
    !----
    If (nsp(18) > 0) Then
      bs(18)%shells = 1
      Allocate(bs(18)%shell(1))
      !
      bs(18)%shell(1)%K    = 4
      bs(18)%shell(1)%N    = 2
      bs(18)%shell(1)%lmin = 0
      bs(18)%shell(1)%lmax = 1
      Allocate(bs(18)%shell(1)%l(2))
      Allocate(bs(18)%shell(1)%iao(2))
      Allocate(bs(18)%shell(1)%alpha(4))
      Allocate(bs(18)%shell(1)%D(4,2))
      bs(18)%shell(1)%alpha(1) = 2.70600D+00
      bs(18)%shell(1)%alpha(2) = 1.27800D+00
      bs(18)%shell(1)%alpha(3) = 0.43540D+00
      bs(18)%shell(1)%alpha(4) = 0.14760D+00
      bs(18)%shell(1)%l(1) = 0
      bs(18)%shell(1)%iao(1) = 0
      bs(18)%shell(1)%D(1,1) = -0.31286D+00
      bs(18)%shell(1)%D(2,1) =  0.11821D+00
      bs(18)%shell(1)%D(3,1) =  0.86786D+00
      bs(18)%shell(1)%D(4,1) =  0.22264D+00
      bs(18)%shell(1)%l(2) = 1
      bs(18)%shell(1)%iao(2) = 1
      bs(18)%shell(1)%D(1,2) = -0.10927D+00
      bs(18)%shell(1)%D(2,2) =  0.32601D+00
      bs(18)%shell(1)%D(3,2) =  0.57952D+00
      bs(18)%shell(1)%D(4,2) =  0.30349D+00
    End If

    ! Scaling:
    Do Z=1, LMZ
       Do i=1, bs(Z)%shells
          bs(Z)%shell(i)%alpha = bs(Z)%shell(i)%alpha * (zsp(Z)*zsp(Z))
       End Do
    End Do
  End Subroutine InitECP3G



  Subroutine SShell(sh, EXX, CS, zs)
    Implicit None
    Type(Shell)                    :: sh
    Double Precision, Dimension(:) :: EXX, CS
    Double Precision               :: zs
    ! End of dummy parameters.
    Integer                        :: K

    K = Size(EXX)
    sh%K    = K
    sh%N    = 1
    sh%lmin = 0
    sh%lmax = 0
    Allocate(sh%l(1))
    Allocate(sh%iao(1))
    Allocate(sh%alpha(K))
    Allocate(sh%D(K,1))
    sh%l      = 0
    sh%iao    = 0
    sh%alpha  = EXX * (zs*zs)
    sh%D(:,1) = CS
  End Subroutine SShell



  Subroutine PShell(sh, EXX, CP, zp)
    Implicit None
    Type(Shell)                    :: sh
    Double Precision, Dimension(:) :: EXX, CP
    Double Precision               :: zp
    ! End of dummy parameters.
    Integer                        :: K

    K = Size(EXX)
    sh%K    = K
    sh%N    = 1
    sh%lmin = 1
    sh%lmax = 1
    Allocate(sh%l(1))
    Allocate(sh%iao(1))
    Allocate(sh%alpha(K))
    Allocate(sh%D(K,1))
    sh%l      = 1
    sh%iao    = 1
    sh%alpha  = EXX * (zp*zp)
    sh%D(:,1) = CP
  End Subroutine PShell



  Subroutine DShell(sh, EXXD, CD, zd)
    Implicit None
    Type(Shell)                    :: sh
    Double Precision, Dimension(:) :: EXXD, CD
    Double Precision               :: zd
    ! End of dummy parameters.
    Integer                        :: K

    K = Size(EXXD)
    sh%K    = K
    sh%N    = 1
    sh%lmin = 2
    sh%lmax = 2
    Allocate(sh%l(1))
    Allocate(sh%iao(1))
    Allocate(sh%alpha(K))
    Allocate(sh%D(K,1))
    sh%l      = 2
    sh%iao    = 4
    sh%alpha  = EXXD * (zd*zd)
    sh%D(:,1) = CD
  End Subroutine DShell



  Subroutine SPShell(sh, EXX, CS, CP, zsp)
    Implicit None
    Type(Shell)                    :: sh
    Double Precision, Dimension(:) :: EXX, CS, CP
    Double Precision               :: zsp
    ! End of dummy parameters.
    Integer                        :: K

    K = Size(EXX)
    sh%K    = K
    sh%N    = 2
    sh%lmin = 0
    sh%lmax = 1
    Allocate(sh%l(2))
    Allocate(sh%iao(2))
    Allocate(sh%alpha(K))
    Allocate(sh%D(K,2))
    sh%l(1)   = 0
    sh%l(2)   = 1
    sh%iao(1) = 0
    sh%iao(2) = 1
    sh%alpha  = EXX * (zsp*zsp)
    sh%D(:,1) = CS
    sh%D(:,2) = CP
  End Subroutine SPShell



  Subroutine Renormalize(bs)
    Implicit None
    Type(BasisSet), Dimension(:), Target      :: bs
    ! End of dummy parameters.
    Type(Shell), Pointer                      :: sh
    Double Precision, Dimension(:,:), Pointer :: S
    Double Precision, Parameter               :: PI = 3.1415926535897932385D0
    Double Precision                          :: zeta, x, w
    Integer                                   :: Z, ish, i, j, k, l, n

    Do Z=1, UBound(bs,1)
       Do ish=1, bs(Z)%shells
          sh => bs(Z)%shell(ish)
          If (sh%K == 0)  Exit
          Allocate(S(sh%K, sh%K))

          Do l=sh%lmin, sh%lmax
             Do i=1, sh%K
                Do j=1, i
                   zeta = sh%alpha(i) + sh%alpha(j)
                   x    = PI / zeta
                   S(i,j) = x * Sqrt(x) / ((zeta+zeta) ** l)
                   If (i /= j)  S(j,i) = S(i,j)
                EndDo
             EndDo

             Do n=1, sh%N
                If (sh%l(n) == l) Then
                   Do k=1, sh%K
                      sh%D(k,n) = sh%D(k,n) / Sqrt(S(k,k))
                   EndDo
                   x = 0
                   Do i=1, sh%K
                      Do j=1, sh%K
                         x = x + sh%D(i,n) * sh%D(j,n) * S(i,j)
                      EndDo
                   EndDo

                   w = 1.D0 / Sqrt(x)

                   If (Abs(w-1.D0) >= 1.D-4) Then
                      Print '(A)', 'WARNING: Renormalization factor deviates significantly from 1:'
                      Print '(A,I2,A,I2,A,I2,A,I2,A,E16.8)', &
                           & 'Z=', Z, ', shell #', ish, ', n=', n, ', l=', l, ', W=', w
                   EndIf
                   Do k=1, sh%K
                      sh%D(k,n) = sh%D(k,n) * w
                   EndDo
                EndIf
             EndDo
          EndDo

          Deallocate(S)
       EndDo
    EndDo
  End Subroutine Renormalize



  Subroutine PrintBasisSet(bs, stdout)
    Type(BasisSet), Dimension(:) :: bs
    Integer                      :: stdout
    ! End of dummy parameters.
    Integer                      :: LMZ, Z, i, k, n

    LMZ = Size(bs)

    Write(stdout,'(1X,A)') 'AO basis set for the evaluation of spectroscopic observables:'

    Do Z=1, LMZ
       Do i=1, bs(Z)%shells
          Write(stdout,'(/23X,2I20)') (bs(Z)%shell(i)%l(n), n=1,bs(Z)%shell(i)%N)
          Do k=1, bs(Z)%shell(i)%K
             Write(stdout,'(1X,I2,2I4,3E20.10)') Z, i, k, bs(Z)%shell(i)%alpha(k), (bs(Z)%shell(i)%D(k,n), n=1,bs(Z)%shell(i)%N)
          End Do
       End Do
    End Do

    Write(stdout,'(//)')
  End Subroutine PrintBasisSet


  Subroutine FreeBasisSet(bs)
    Type(BasisSet), Pointer :: bs(:)
    ! End of dummy parameters.
    Integer                 :: Z, i

    Do Z=LBound(bs,1), UBound(bs,1)
       If (bs(Z)%shells > 0) Then
          Do i=1, bs(Z)%shells
             Deallocate(bs(Z)%shell(i)%l)
             Deallocate(bs(Z)%shell(i)%iao)
             Deallocate(bs(Z)%shell(i)%alpha)
             Deallocate(bs(Z)%shell(i)%D)
          End Do
          Deallocate(bs(Z)%shell)
       End If
    End Do

    Deallocate(bs)
  End Subroutine FreeBasisSet


End Module basissets
