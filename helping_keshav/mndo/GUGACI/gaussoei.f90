Module gaussoei

!! Calculate one-electron integrals over Gaussian basis functions,
!! as needed for the evaluation of electric and magnetic dipole
!! transition moments.

!!
!! Written by
!! Axel Koslowski,
!! Max-Planck-Institut fuer Kohlenforschung,
!! April 2000
!!

!! History:
!!
!! 0 Initial version adapted from C++ program
!!   PACI (Program for the Analysis of Chromophoric Interactions)
!!   by Axel Koslowski, Department of Biochemistry, Colorado State
!!   University, 1999

Use basissets
Use sphrharm

Implicit None

Private

Public :: CalcPropOEI


! The following constants must be large enough for the
! Gaussian basis set which will be used:

Integer, Parameter :: MaxL  = 2
Integer, Parameter :: MaxN  = 2

Integer, Parameter :: MaxL1 = MaxL + 1


! Miscellaneous numerical constants:

Double Precision, Parameter :: ovlpconst = 1.5749609947135716530D+01


! Structure that contains the primitive shell pair data:

Type PrimPair
  Double Precision :: oneov2ze, alovze, beovze, kappa
End Type PrimPair


! Structures to store the intermediates:

Type Cube1
  Double Precision, Dimension(0:MaxL1, 0:MaxL1, 0:MaxL1) :: b
End Type Cube1

Type CubeCube1
  Type (Cube1), Dimension(0:MaxL1, 0:MaxL1, 0:MaxL1) :: a
End Type CubeCube1

Type Sqr1
  Double Precision, Dimension(0:MaxL1, 0:MaxL1) :: b
End Type Sqr1

Type SqrSqr1
  Type (Sqr1), Dimension(0:MaxL1, 0:MaxL1) :: a
End Type SqrSqr1

Type Vector
  Double Precision, Dimension(3) :: vec
End Type Vector

Type Sqr3
  Type (Vector), Dimension(0:MaxL, 0:MaxL) :: b
End Type Sqr3

Type SqrSqr3
  Type (Sqr3), Dimension(0:MaxL, 0:MaxL) :: a
End Type SqrSqr3


Contains

  ! The following subroutine is the public driver routine for the calculation of the
  ! overlap and R, DEL, and R x DEL operator matrices in a Gaussian atomic orbital
  ! basis. The overlap integrals are generated first using the procedure by
  ! Head-Gordon and Pople (insert citation here). However, the horizontal
  ! recurrence relation is applied to the uncontracted overlap integrals, because
  ! these are used to compute the R and DEL integrals. The R x DEL integrals are
  ! finally calculated from the contracted overlap and DEL integrals.
  ! All quantities in atomic units.
  !
  Subroutine CalcPropOEI(Smat, R, Del, RxDel, Z, First, Last, Coord, zs, zp, zd, ns, np, nd, atoms, stdout, ecp)
  Implicit None
  Double Precision, Dimension(:,:)      :: Smat
  Double Precision, Dimension(:,:,:)    :: R, Del, RxDel   ! Dimension(:,:,3)
  Integer, Dimension(:)                 :: Z, First, Last
  Double Precision, Dimension(:,:)      :: Coord           ! Dimension(3,:)
  Double Precision, Dimension(:)        :: zs, zp, zd
  Integer,          Dimension(:)        :: ns, np, nd
  Integer                               :: atoms, stdout
  Logical                               :: ecp
  ! End of dummy parameters.
  Type(BasisSet), Dimension(:), Pointer :: bs
  Type(Shell), Pointer                  :: sha, shb
  Type(PrimPair)                        :: pp
  Type(CubeCube1)                       :: PrimS
  Type(SqrSqr1), Dimension(MaxN, MaxN)  :: ContS
  Type(SqrSqr3), Dimension(MaxN, MaxN)  :: ContR
  Type(SqrSqr3), Dimension(MaxN, MaxN)  :: ContDel
  Type(SqrSqr3), Dimension(MaxN, MaxN)  :: ContRxDel
  Double Precision, Dimension(3)        :: ABdiff
  Integer                               :: a, b, s, t, mu, nu
  Integer                               :: ka, kb, na, nb, la, lb
  Double Precision                      :: r2AB, DD

  Call InitYlm

  If (ecp) Then
     Write(stdout,'(1X,"Using basis sets ECP-3G (first-row elements) and ECP-4G (second-row")')
     Write(stdout,'(1X,"elements) for the evaluation spectroscopic observables."//)')
     Call InitECP3G(bs, zs, ns)
  Else
     Write(stdout,'(1X,A//)') 'Using STO-6G basis for the evaluation spectroscopic observables.'
     Call InitSTO6G(bs, zs, zp, zd, ns, np, nd)
  End If

  Call PrintBasisSet(bs, stdout)
  Call Renormalize(bs)

  Do a=1, atoms
    If (bs(Z(a))%shells == 0) Then
      Write(stdout,'(A,I3,A)') ' No Gaussian basis set for element', Z(a), '.'
      Stop 'CalcPropOEI'
    End If
    Do b=1, a
      ABdiff = Coord(:,a) - Coord(:,b)
      r2AB = ABdiff(1) * ABdiff(1) + ABdiff(2) * ABdiff(2) + ABdiff(3) * ABdiff(3)
      Do s=1, bs(Z(a))%shells
        sha => bs(Z(a))%shell(s)
        mu = First(a) + sha%iao(1)
        If (mu > Last(a)) Exit
        Do t=1, bs(Z(b))%shells
          shb => bs(Z(b))%shell(t)
          nu = First(b) + shb%iao(1)
          If (nu > Last(b)) Exit
          Do na=1, sha%N
            Do nb=1, shb%N
              la = sha%l(na)
              lb = shb%l(nb)
              Call ClearSqrSqr1(ContS(na,nb), la, lb)
              Call ClearSqrSqr3(ContR(na,nb), la, lb)
              Call ClearSqrSqr3(ContDel(na,nb), la, lb)
              ! There is no need to clear ContRxDel because it is
              ! computed from the contracted Del and S integrals.
            EndDo
          EndDo
          Do ka=1, sha%K
            Do kb=1, shb%K
              Call InitPrimPair(pp, sha, shb, ka, kb, r2AB)
              Call CalcS(PrimS, pp, ABdiff, sha%lmin, sha%lmax, shb%lmax+1)
              Do na=1, sha%N
                Do nb=1, shb%N
                  la = sha%l(na)
                  lb = shb%l(nb)
                  DD = sha%D(ka,na) * shb%D(kb,nb)
                  Call ContractS(ContS(na,nb), PrimS, la, lb, DD)
                  Call CalcRfromS(ContR(na,nb), PrimS, Coord(:,b), la, lb, DD)
                  Call CalcDelFromS(ContDel(na,nb), PrimS, 2.D0*shb%alpha(kb), la, lb, DD)
                EndDo
              EndDo
            EndDo
          EndDo
          Do na=1, sha%N
            Do nb=1, shb%N
              la = sha%l(na)
              lb = shb%l(nb)
              mu = First(a) + sha%iao(na)
              nu = First(b) + shb%iao(nb)
              Call CalcRxDel(ContRxDel(na,nb), ContDel(na,nb), ContS(na,nb), Coord(:,b), sha%l(na), shb%l(nb))
              Call Normalize(ContS(na,nb), ContR(na,nb), ContDel(na,nb), ContRxDel(na,nb), sha%l(na), shb%l(nb))
              Call Store1(Smat,  ContS(na,nb),     la, lb, mu, nu, +1.D0)
              Call Store3(R,     ContR(na,nb),     la, lb, mu, nu, +1.D0)
              Call Store3(Del,   ContDel(na,nb),   la, lb, mu, nu, -1.D0)
              Call Store3(RxDel, ContRxDel(na,nb), la, lb, mu, nu, -1.D0)
            EndDo
          EndDo
        EndDo
      EndDo
    EndDo
  EndDo

  Call FreeBasisSet(bs)
  End Subroutine CalcPropOEI



  Subroutine ClearSqrSqr1(sqsq1, la, lb)
  Implicit None
  Type(SqrSqr1) :: sqsq1
  Integer       :: la, lb
  ! End of dummy parameters.
  Integer       :: ka, ja, kb, jb

  Do ka=0, la
    Do ja=0, la-ka
      Do kb=0, lb
        Do jb=0, lb-kb
          sqsq1%a(ka,ja)%b(kb,jb) = 0.D0
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine ClearSqrSqr1



  Subroutine ClearSqrSqr3(sqsq3, la, lb)
  Implicit None
  Type(SqrSqr3) :: sqsq3
  Integer       :: la, lb
  ! End of dummy parameters.
  Integer       :: ka, ja, kb, jb

  Do ka=0, la
    Do ja=0, la-ka
      Do kb=0, lb
        Do jb=0, lb-kb
          sqsq3%a(ka,ja)%b(kb,jb)%vec = 0.D0
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine ClearSqrSqr3



  Subroutine InitPrimPair(pp, sha, shb, ka, kb, r2AB)
  Implicit None
  Type(PrimPair)                 :: pp
  Type(Shell)                    :: sha, shb
  Integer                        :: ka, kb
  Double Precision               :: r2ab
  ! End of dummy parameters.
  Double Precision               :: alpha, beta, oneovze

  alpha       = sha%alpha(ka)
  beta        = shb%alpha(kb)
  oneovze     = 1.D0 / (alpha + beta)
  pp%oneov2ze = 0.5D0 * oneovze
  pp%alovze   = alpha * oneovze
  pp%beovze   = beta  * oneovze
  pp%kappa    = Exp(-oneovze * alpha * beta * r2AB)
  End Subroutine InitPrimPair



  Subroutine CalcS(S, pp, ABdiff, Llamin, Llamax, Llblim)
  Implicit None
  Type(CubeCube1), Target        :: S
  Type(PrimPair)                 :: pp
  Double Precision, Dimension(3) :: ABdiff
  Integer                        :: Llamin, Llamax, Llblim
  ! End of dummy parameters.
  Type(Cube1), Pointer           :: So, Sa
  Integer                        :: ia, ja, ka, la
  Integer                        :: ib, jb, kb, lb
  Integer                        :: lblim
  Double Precision               :: help

  So => S%a(0,0,0)

  So%b(0,0,0) = ovlpconst * pp%oneov2ze * sqrt(pp%oneov2ze) * pp%kappa

  lblim = Min(Llblim-Llamin+1, Llblim)

  Do lb=1, lblim
    Do kb=0, lb
      Do jb=0, lb-kb
        ib = lb - kb - jb

        If (ib == 1) Then
          So%b(lb,kb,jb) = pp%alovze * ABdiff(1) * So%b(lb-1, kb, jb);
        Else If (jb == 1) Then
          So%b(lb,kb,jb) = pp%alovze * ABdiff(2) * So%b(lb-1, kb,  0);
        Else If (kb == 1) Then
          So%b(lb,kb,jb) = pp%alovze * ABdiff(3) * So%b(lb-1,  0, jb);
        Else If (ib /= 0) Then
          So%b(lb,kb,jb) = pp%alovze * ABdiff(1) * So%b(lb-1, kb,   jb)   + pp%oneov2ze * (ib-1) * So%b(lb-2, kb,   jb);
        Else If (jb /= 0) Then
          So%b(lb,kb,jb) = pp%alovze * ABdiff(2) * So%b(lb-1, kb,   jb-1) + pp%oneov2ze * (jb-1) * So%b(lb-2, kb,   jb-2);
        Else
          So%b(lb,kb,jb) = pp%alovze * ABdiff(3) * So%b(lb-1, kb-1, jb)   + pp%oneov2ze * (kb-1) * So%b(lb-2, kb-2, jb);
        EndIf
      EndDo
    EndDo
  EndDo

  Do la=1, Llamax
    Do ka=0, la
      Do ja=0, la-ka
        ia = la - ka - ja

        Sa => S%a(la,ka,ja)

        If (ia == 1) Then
      Sa%b(0,0,0) = -pp%beovze * ABdiff(1) * S%a(la-1, ka,   ja  )%b(0,0,0)
        Else If (ja == 1) Then
      Sa%b(0,0,0) = -pp%beovze * ABdiff(2) * S%a(la-1, ka,   0   )%b(0,0,0)
        Else If (ka == 1) Then
      Sa%b(0,0,0) = -pp%beovze * ABdiff(3) * S%a(la-1, 0,    ja  )%b(0,0,0)
        Else If (ia /= 0) Then
      Sa%b(0,0,0) = -pp%beovze * ABdiff(1) * S%a(la-1, ka,   ja  )%b(0,0,0) + pp%oneov2ze * (ia-1) * S%a(la-2, ka, ja)%b(0,0,0)
        Else If (ja /= 0) Then
      Sa%b(0,0,0) = -pp%beovze * ABdiff(2) * S%a(la-1, ka,   ja-1)%b(0,0,0) + pp%oneov2ze * (ja-1) * S%a(la-2, ka, ja-2)%b(0,0,0)
        Else
      Sa%b(0,0,0) = -pp%beovze * ABdiff(3) * S%a(la-1, ka-1, ja  )%b(0,0,0) + pp%oneov2ze * (ka-1) * S%a(la-2, ka-2, ja)%b(0,0,0)
        EndIf

        lblim = Min(Llblim-Llamin+la+1, Llblim)

        Do lb=1, lblim
          Do kb=0, lb
            Do jb=0, lb-kb
              ib = lb - kb - jb

              If (ib == 1) Then
                help = pp%alovze * ABdiff(1) * Sa%b(lb-1, kb, jb)
                If (ia /= 0) help = help + pp%oneov2ze * ia * S%a(la-1, ka,   ja  )%b(lb-1, kb, jb)
                Sa%b(lb,kb,jb) = help
              Else If (jb == 1) Then
                help = pp%alovze * ABdiff(2) * Sa%b(lb-1, kb,  0)
                If (ja /= 0) help = help + pp%oneov2ze * ja * S%a(la-1, ka,   ja-1)%b(lb-1, kb,  0)
                Sa%b(lb,kb,jb) = help
              Else If (kb == 1) Then
                help = pp%alovze * ABdiff(3) * Sa%b(lb-1,  0, jb)
                If (ka /= 0) help = help + pp%oneov2ze * ka * S%a(la-1, ka-1, ja  )%b(lb-1,  0, jb)
                Sa%b(lb,kb,jb) = help
              Else If (ib /= 0) Then
                help = pp%alovze * ABdiff(1) * Sa%b(lb-1, kb, jb  ) + pp%oneov2ze * (ib-1) * Sa%b(lb-2, kb,   jb)
                If (ia /= 0) help = help + pp%oneov2ze * ia * S%a(la-1, ka,   ja  )%b(lb-1, kb,   jb)
                Sa%b(lb,kb,jb) = help
              Else If (jb /= 0) Then
                help = pp%alovze * ABdiff(2) * Sa%b(lb-1, kb, jb-1) + pp%oneov2ze * (jb-1) * Sa%b(lb-2, kb,   jb-2)
                If (ja /= 0) help = help + pp%oneov2ze * ja * S%a(la-1, ka,   ja-1)%b(lb-1, kb,   jb-1)
                Sa%b(lb,kb,jb) = help
              Else
                help = pp%alovze * ABdiff(3) * Sa%b(lb-1, kb-1, jb) + pp%oneov2ze * (kb-1) * Sa%b(lb-2, kb-2, jb)
                If (ka /= 0) help = help + pp%oneov2ze * ka * S%a(la-1, ka-1, ja  )%b(lb-1, kb-1, jb)
                Sa%b(lb,kb,jb) = help
              EndIf
            EndDo
          EndDo
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine CalcS



  Subroutine ContractS(ContS, PrimS, la, lb, DD)
  Implicit None
  Type(SqrSqr1),   Target :: ContS
  Type(CubeCube1), Target :: PrimS
  Integer                 :: la, lb
  Double Precision        :: DD
  ! End of dummy parameters.
  Type(Sqr1),  Pointer    :: ContSa
  Type(Cube1), Pointer    :: PrimSa
  Integer                 :: ja, ka
  Integer                 :: jb, kb

  Do ka=0, la
    Do ja=0, la-ka
      ContSa => ContS%a(ka,ja)
      PrimSa => PrimS%a(la,ka,ja)

      Do kb=0, lb
        Do jb=0, lb-kb
          ContSa%b(kb,jb) = ContSa%b(kb,jb) + DD * PrimSa%b(lb,kb,jb)
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine ContractS



  Subroutine CalcRfromS(R, S, B, la, lb, DD)
  Implicit None
  Type(SqrSqr3),   Target        :: R
  Type(CubeCube1), Target        :: S
  Double Precision, Dimension(3) :: B
  Integer                        :: la, lb
  Double Precision               :: DD
  ! End of dummy parameters.
  Type(Sqr3),   Pointer          :: Ra
  Type(Cube1),  Pointer          :: Sa
  Type(Vector), Pointer          :: Rab
  Double Precision               :: Sab
  Integer                        :: ja, ka
  Integer                        :: jb, kb

  Do ka=0, la
    Do ja=0, la-ka
      Ra => R%a(ka,ja)
      Sa => S%a(la,ka,ja)

      Do kb=0, lb
        Do jb=0, lb-kb
          Rab => Ra%b(kb,jb)
          Sab =  Sa%b(lb,kb,jb)

          Rab%vec(1) = Rab%vec(1) + DD * (Sa%b(lb+1, kb,   jb)   + B(1) * Sab)
          Rab%vec(2) = Rab%vec(2) + DD * (Sa%b(lb+1, kb,   jb+1) + B(2) * Sab)
          Rab%vec(3) = Rab%vec(3) + DD * (Sa%b(lb+1, kb+1, jb)   + B(3) * Sab)
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine CalcRfromS



  Subroutine CalcDelFromS(Del, S, twobeta, la, lb, DD)
  Implicit None
  Type(SqrSqr3),   Target        :: Del
  Type(CubeCube1), Target        :: S
  Double Precision               :: twobeta
  Integer                        :: la, lb
  Double Precision               :: DD
  ! End of dummy parameters.
  Type(Sqr3),   Pointer          :: Dela
  Type(Cube1),  Pointer          :: Sa
  Type(Vector), Pointer          :: Delab
  Double Precision               :: help
  Integer                        :: ja, ka
  Integer                        :: ib, jb, kb

  Do ka=0, la
    Do ja=0, la-ka
      Sa   => S%a(la,ka,ja)
      Dela => Del%a(ka,ja)

      Do kb=0, lb
        Do jb=0, lb-kb
          ib = lb - kb - jb
          Delab => Dela%b(kb,jb)

          help = twobeta * Sa%b(lb+1, kb, jb)
          If (ib /= 0)  help = help - ib * Sa%b(lb-1, kb, jb)
          Delab%vec(1) = Delab%vec(1) - DD * help

          help = twobeta * Sa%b(lb+1, kb, jb+1)
          if (jb /= 0)  help = help - jb * Sa%b(lb-1, kb, jb-1)
          Delab%vec(2) = Delab%vec(2) - DD * help

          help = twobeta * Sa%b(lb+1, kb+1, jb)
          if (kb /= 0)  help = help - kb * Sa%b(lb-1, kb-1, jb)
          Delab%vec(3) = Delab%vec(3) - DD * help
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine CalcDelFromS



  Subroutine CalcRxDel(RxDel, Del, S, B, La, Lb)
  Implicit None
  Type(SqrSqr3), Target          :: RxDel
  Type(SqrSqr3), Target          :: Del
  Type(SqrSqr1), Target          :: S
  Double Precision, Dimension(3) :: B
  Integer                        :: La, Lb
  ! End of dummy parameters.
  Type(Sqr3),   Pointer          :: RxDela
  Type(Sqr3),   Pointer          :: Dela
  Type(Sqr1),   Pointer          :: Sa
  Type(Vector), Pointer          :: RxDelab
  Type(Vector), Pointer          :: Delab
  Integer                        :: ja, ka
  Integer                        :: ib, jb, kb

  Do ka=0, La
    Do ja=0, La-ka
      RxDela => RxDel%a(ka,ja)
      Dela   => Del%a(ka,ja)
      Sa     => S%a(ka,ja)

      Do kb=0, Lb
        Do jb=0, Lb-kb
          ib = Lb - kb - jb

          RxDelab => RxDela%b(kb,jb)
          Delab   => Dela%b(kb,jb)

          RxDelab%vec(1) = B(2) * Delab%vec(3) - B(3) * Delab%vec(2)
          RxDelab%vec(2) = B(3) * Delab%vec(1) - B(1) * Delab%vec(3)
          RxDelab%vec(3) = B(1) * Delab%vec(2) - B(2) * Delab%vec(1)

          If (ib /= 0) Then
            RxDelab%vec(2) = RxDelab%vec(2) + ib * Sa%b(kb+1, jb)
            RxDelab%vec(3) = RxDelab%vec(3) - ib * Sa%b(kb,   jb+1)
          EndIf

          If (jb /= 0) Then
            RxDelab%vec(3) = RxDelab%vec(3) + jb * Sa%b(kb,   jb-1)
            RxDelab%vec(1) = RxDelab%vec(1) - jb * Sa%b(kb+1, jb-1)
          EndIf

          If (kb /= 0) Then
            RxDelab%vec(1) = RxDelab%vec(1) + kb * Sa%b(kb-1, jb+1)
            RxDelab%vec(2) = RxDelab%vec(2) - kb * Sa%b(kb-1, jb)
          EndIf
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine CalcRxDel



  Subroutine Normalize(S, R, Del, RxDel, La, Lb)
  Implicit None
  Type(SqrSqr1), Target            :: S
  Type(SqrSqr3), Target            :: R
  Type(SqrSqr3), Target            :: Del
  Type(SqrSqr3), Target            :: RxDel
  Integer                          :: La, Lb
  ! End of dummy parameters.
  Type(Sqr1), Pointer              :: Sa
  Type(Sqr3), Pointer              :: Ra
  Type(Sqr3), Pointer              :: Dela
  Type(Sqr3), Pointer              :: RxDela
  Double Precision                 :: N, NN
  Integer                          :: ja, ka
  Integer                          :: jb, kb
  Double Precision, Dimension(0:6) :: ofac
  Data ofac /1.D0, 1.D0, 3.D0, 15.D0, 105.D0, 945.D0, 10395.D0/

  Do ka=0, La
    Do ja=0, La-ka
      N = ofac(La-ka-ja) * ofac(ja) * ofac(ka)

      Sa => S%a(ka,ja)
      Ra => R%a(ka,ja)
      Dela => Del%a(ka,ja)
      RxDela => RxDel%a(ka,ja)

      Do kb=0, Lb
        Do jb=0, Lb-kb
          NN = 1.D0 / Sqrt(N * ofac(Lb-kb-jb) * ofac(jb) * ofac(kb))

              Sa%b(kb,jb)     =     Sa%b(kb,jb)     * NN
              Ra%b(kb,jb)%vec =     Ra%b(kb,jb)%vec * NN
            Dela%b(kb,jb)%vec =   Dela%b(kb,jb)%vec * NN
          RxDela%b(kb,jb)%vec = RxDela%b(kb,jb)%vec * NN
        EndDo
      EndDo
    EndDo
  EndDo
  End Subroutine Normalize



  Subroutine Store1(X, ContX, la, lb, mu, nu, factor)
  Implicit None
  Double Precision, Dimension(:,:) :: X
  Type(SqrSqr1)                    :: ContX
  Integer                          :: la, lb, mu, nu
  Double Precision                 :: factor
  ! End of dummy parameters.
  Double Precision                 :: ca, cb, y, z
  Integer                          :: ja, ka, jb, kb, i, j, ii, jj

  Do i=mu, mu+2*la
    Do j=nu, nu+2*lb
      y = 0.D0
      Do ii=1, Ylm(la, i-mu+1)%n
        ka = Ylm(la, i-mu+1)%cart(ii)%k
        ja = Ylm(la, i-mu+1)%cart(ii)%j
        ca = Ylm(la, i-mu+1)%cart(ii)%c
        z = 0.D0
        Do jj=1, Ylm(lb, j-nu+1)%n
          kb = Ylm(lb, j-nu+1)%cart(jj)%k
          jb = Ylm(lb, j-nu+1)%cart(jj)%j
          cb = Ylm(lb, j-nu+1)%cart(jj)%c
          z = z + cb * ContX%a(ka,ja)%b(kb,jb)
        EndDo
        y = y + ca * z
      EndDo
      X(i,j) = y
      if (mu /= nu)  X(j,i) = y * factor
    EndDo
  EndDo
  End Subroutine Store1



  Subroutine Store3(X, ContX, la, lb, mu, nu, factor)
  Implicit None
  Double Precision, Dimension(:,:,:) :: X
  Type(SqrSqr3)                      :: ContX
  Integer                            :: la, lb, mu, nu
  Double Precision                   :: factor
  ! End of dummy parameters.
  Double Precision, Dimension(3)     :: y, z
  Double Precision                   :: ca, cb
  Integer                            :: ja, ka, jb, kb, i, j, ii, jj

  Do i=mu, mu+2*la
    Do j=nu, nu+2*lb
      y = 0.D0
      Do ii=1, Ylm(la, i-mu+1)%n
        ka = Ylm(la, i-mu+1)%cart(ii)%k
        ja = Ylm(la, i-mu+1)%cart(ii)%j
        ca = Ylm(la, i-mu+1)%cart(ii)%c
        z = 0.D0
        Do jj=1, Ylm(lb, j-nu+1)%n
          kb = Ylm(lb, j-nu+1)%cart(jj)%k
          jb = Ylm(lb, j-nu+1)%cart(jj)%j
          cb = Ylm(lb, j-nu+1)%cart(jj)%c
          z = z + cb * ContX%a(ka,ja)%b(kb,jb)%vec
        EndDo
        y = y + ca * z
      EndDo
      X(i,j,:) = y
      if (mu /= nu)  X(j,i,:) = y * factor
    EndDo
  EndDo
  End Subroutine Store3

End Module gaussoei
