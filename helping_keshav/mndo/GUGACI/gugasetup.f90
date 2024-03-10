!
! Initialize ShavittControl structure and construct distinct row table.
!
! First version by Michael E. Beck in 1998 at Zurich University.
!
! Final version by Axel Koslowski in 2000-2005 at MPI Muelheim.
!
! Shape-driven algorithm added by Axel Koslowski in 2012-2013.
!

Module gugasetup

  Use gugaglobal
  Use gugaplots
  Use gugashape
  Use gugautils

  Implicit None

  Private

  Public :: MakeShavittControl  ! Initialize ShavittControl structure and construct DRT.
  Public :: Cleanup             ! Free all dynamic memory.


  ! Special type:

  Type NewPaldusRow
     Integer :: a, b, d, j
  End Type NewPaldusRow


Contains

  !
  ! Public driver routine to initialize ShavittControl structure from
  ! MNDO input information and to construct the distinct row table.
  !
  Subroutine MakeShavittControl(Flags, DRT, integerInputOptions, norbs, nocc, ICIREF, IMOCI,  &
                                IROOTA, NSYM, SymLabel, GrpLabel, NumGPU, NVCapa)
    Implicit None
    Type(ShavittControl),                     Intent(InOut) :: Flags
    Type(PaldusRow),  Dimension(:),   Pointer               :: DRT
    Integer,          Dimension(:),           Intent(In)    :: integerInputOptions  ! Input options (IN2).
    Integer,                                  Intent(In)    :: norbs, nocc          ! Number of SCF MOs, occ. MOs.
    Integer,          Dimension(:,:),         Intent(In)    :: ICIREF               ! Reference configurations.
    Integer,          Dimension(:),   Target, Intent(In)    :: IMOCI                ! List of active MO indices.
    Integer,          Dimension(:),           Intent(In)    :: IROOTA               ! Requested number of CI roots of each irrep.
    Integer,          Dimension(:),           Intent(In)    :: NSYM                 ! Irreps of active MOs.
    Character(Len=4), Dimension(:),           Intent(In)    :: SymLabel             ! List of symmetry labels of all irreps.
    Character(Len=3)                                        :: GrpLabel             ! Point group label.
    Integer,                                  Intent(In)    :: NumGPU               ! Number of CUDA-capable GPUs.
    Integer,          Dimension(NumGPU),      Intent(In)    :: NVCapa               ! Compute capabilities of GPUs.
    ! End of dummy parameters.
    Type(StepVector), Dimension(:),   Pointer               :: RefConf
    Character(Len=4), Dimension(0:3)                        :: OccStr
    Integer                                                 :: numberOfActiveOrbitals, nDRT, i, j, k, l
    Integer                                                 :: firstActiveIndex, lastActiveIndex
    Integer                                                 :: iroot, ncisym, total
    Double Precision                                        :: tuser, tsys, twall

    Data OccStr / '   -', '   a', '   b', '  ab' /

    ! Set "easy" flags:
    Flags%HOMOindex         = nocc
    firstActiveIndex        = nocc - integerInputOptions(131) + 1  ! ici1
    lastActiveIndex         = nocc + integerInputOptions(132)      ! ici2
    Flags%nciref            = integerInputOptions(136)  ! nciref
    Flags%mciref            = integerInputOptions(137)  ! mciref
    Flags%ExcitationLevel   = integerInputOptions(138)  ! levexc
    iroot                   = integerInputOptions(139)  ! iroot
    Flags%LRoot             = integerInputOptions(140)  ! lroot
    Flags%iState            = integerInputOptions(140)  ! lroot
    Flags%jState            = integerInputOptions(140)  ! lroot
    Flags%Charge            = integerInputOptions(141)  ! cichg
    Flags%Multiplicity      = integerInputOptions(142)  ! multci
    ncisym                  = integerInputOptions(143)  ! ncisym
    Flags%Algorithm         = integerInputOptions(144)  ! cidir
    Flags%WhichDiagonalizer = integerInputOptions(145)  ! cidiag
    Flags%MinDav            = integerInputOptions(161)  ! mindav
    Flags%MaxDav            = integerInputOptions(162)  ! maxdav
    Flags%KitDav            = integerInputOptions(163)  ! kitdav
    Flags%QNorm             = 1.D1**(-Dble(integerInputOptions(164)))  ! nrmdav
    Flags%iop               = integerInputOptions(2)    ! iop
    Flags%Leading           = Dble(integerInputOptions(150))*1.D-4  ! cilead
    Flags%Plot              = integerInputOptions(149)  ! ciplot
    Flags%PrintLevel        = integerInputOptions(133)  ! ioutci
    Flags%NumSym            = UBound(SymLabel,1)
    Flags%NRoots            = IROOTA

    ! Save the compute capability of the first NVIDIA GPU.
    ! To use a different device, set the environment variable
    ! CUDA_VISIBLE_DEVICES.
    !
    If (NumGPU > 0) Then
       Flags%NVCapa = NVCapa(1)
    Else
       Flags%NVCapa = 0
    End If

    ! The print level is defined as follows:
    !  -9: no output except error messages
    !   0: standard (i.e. relatively small) output
    !   1: print energies and leading CSFs of all states that have been calculated
    !   2: print Davidson iterations and other information related to the diagonalizers
    !   3: print CSFs as linear combinations of Slater determinants and complete eigensystem
    !   4: print CI Hamiltonian
    !   5: some debug print
    !   6: extensive debug print (HUGE output)

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//A)') ' Setting up...'
       Call GetTime(tuser, tsys, twall)
    End If

    If (ncisym > Flags%NumSym) Then
       Write(Stdout,'(/1X,A,I2,A,A,A)') 'Invalid symmetry number (ncisym =', ncisym, &
                                        ') for point group ', Trim(GrpLabel), '.'
       Stop 'MakeShavittControl'
    End If

    ! Preliminary settings concerning the requested CI roots,
    ! not considering the actual number of CSFs.
    If (iroot > 0) Then
       Flags%TotalRoots         = iroot
       If (ncisym > 0) Then
          Flags%NumberOfCIRoots = 0
          Flags%NRoots(ncisym)  = iroot
       Else
          Flags%NumberOfCIRoots = iroot
       End If
    Else
       Flags%NumberOfCIRoots    = 0
       Flags%TotalRoots         = Sum(IROOTA(1:Flags%NumSym))
    End If

    ! Determine number of active electrons from the first
    ! reference configuration.
    Flags%numberOfActiveElectrons = Sum(ICIREF(:,1))

    ! Initialize list of active orbital indices.
    numberOfActiveOrbitals =  lastActiveIndex - firstActiveIndex + 1
    Flags%WhoIsWho         => IMOCI(1:numberOfActiveOrbitals)

    ! Read symmetries of the active orbitals.
    Allocate(Flags%MOIrreps(numberOfActiveOrbitals))
    Flags%MOIrreps(1:numberOfActiveOrbitals) = NSYM(IMOCI(1:numberOfActiveOrbitals))

    ! Construct DRT.
    Call SetupDistinctRowTable(DRT, iciref, Flags)

    If (ncisym < 0)  Flags%WalkIrrep = 1

    ! Print DRT.
    If(Flags%PrintLevel >= 5)  Call PrintDRT(DRT)

    ! Number of Gelfand states described by DRT.
    nDRT = DRT(1)%yd(4)

    ! Allocate and initialize IndVec.
    Allocate(Flags%IndVec(nDRT))
    Call InitIndexVector(Flags%IndVec, Flags)

    ! Allocate and initialize SymVec.
    ! Element with index 0 is a dummy element that
    ! simplifies the comparisons in module gugashape2.
    Allocate(Flags%SymVec(0:Flags%TotalCSF))
    Flags%SymVec(0) = -1
    Do i=1, nDRT
       j = Flags%IndVec(i)
       If (j /= 0) Then
          Flags%SymVec(j) = Flags%WalkIrrep(i)
       End If
    End Do

    ! Print statistics on CI basis functions:
    If (Flags%PrintLevel >= 1) Then
       Write(Stdout,'(/)')
       Do i=1, Flags%NumSym
          If (Flags%NumCSF(i) > 0) Then
             Write(Stdout,'(A,I9,A,A,A,I1,A)')  " There are ", Flags%NumCSF(i),  &
                " CI basis functions of symmetry ", SymLabel(i), "(", i, ")."
          End If
       End Do
       Write(Stdout,'(A,I9,A)')  &
          " There are ", Flags%TotalCSF, " CI basis functions in the CI problem."
       Write(Stdout,'(A,I9,A)')  &
          " There are ", nDRT, " CI basis functions represented by the DRT."
    End If

    ! Find indices of Gelfand states corresponding to the reference configurations:
    Allocate(Flags%iFirstRef(UBound(iciref,2)))
    Allocate(Flags%iLastRef (UBound(iciref,2)))
    Allocate(Flags%iRefWalk (UBound(iciref,2)))
    Flags%nRefCSF = 0
    Do i=1, Flags%nciref
       Flags%iFirstRef(i) = Flags%nRefCSF + 1
       Call FindGelfandStates(DRT, iciref(:,i), Flags%iRefWalk, Flags%nRefCSF, 1, 1)
       Flags%iLastRef(i)  = Flags%nRefCSF
    End Do
    Allocate(Flags%iRefConf (Flags%nRefCSF))
    Flags%iRefConf = Flags%IndVec(Flags%iRefWalk(1:Flags%nRefCSF))

    ! Find corresponding step vectors:
    Allocate(RefConf(Flags%nRefCSF))
    Do i=1, Flags%nRefCSF
       Allocate(RefConf(i)%d(numberOfActiveOrbitals))
       Call FindStepVector(DRT, RefConf(i)%d, 1, Flags%iRefWalk(i))
    End Do

    ! Print references:
    If (Flags%PrintLevel >= 1) Then
       Write(Stdout,'(//1X,A)') 'Reference configurations and corresponding GUGA-CI many-electron basis functions:'
       Do i=1, Flags%nciref
          Do k=1, (numberOfActiveOrbitals+29)/30
             l = Min(30*k, numberOfActiveOrbitals)
             Write(Stdout,'(/12X,30I4)') Flags%WhoIsWho(30*(k-1)+1:l)
             Write(Stdout,'(1X,A,I7,30I4)') 'Ref.', i, iciref(30*(k-1)+1:l,i)
             If (Flags%iLastRef(i) < Flags%iFirstRef(i)) Then
                Write(Stdout,'(1X,"CSF",8X,A)') 'Not represented by the distinct row table.'
                Cycle
             End If
             Do j=Flags%iFirstRef(i), Flags%iLastRef(i)
                Write(Stdout,'(1X,"CSF",I8,30A)') Flags%iRefConf(j), OccStr(RefConf(j)%d(30*(k-1)+1:l))
             End Do
          End Do
       End Do
    End If

    ! Plot Shavitt graph and loops if requested:
    If (Flags%Plot > 0)  Call PlotShavittGraph(DRT, RefConf, Flags)
    If (Flags%Plot > 1)  Call PlotLoops(DRT, Flags)

    ! Deallocate step vectors:
    Do i=1, Flags%nRefCSF
       Deallocate(RefConf(i)%d)
    End Do
    Deallocate(RefConf)

    ! Final settings concerning the requested CI roots:
    If (Flags%NumberOfCIRoots > 0) Then
       If (Flags%NumberOfCIRoots > Flags%TotalCSF) Then
          Write(Stdout,'(//1X,A)') 'WARNING: Too many CI roots requested. Resetting NumberOfCIRoots.'
          Flags%NumberOfCIRoots = Flags%TotalCSF
          Flags%TotalRoots      = Flags%TotalCSF
          Flags%FirstRoot       = 0
          Flags%LastRoot        = 0
       End If
    Else
       total = 0
       Do i=1, Flags%NumSym
          If (Flags%NRoots(i) > Flags%NumCSF(i)) Then
             Write(Stdout,'(//1X,A,I1,A)') 'WARNING: Too many CI roots requested. Resetting NRoots(',i,').'
             Flags%NRoots(i) = Flags%NumCSF(i)
          End If
          Flags%FirstRoot(i) = total + 1
          total              = total + Flags%NRoots(i)
          Flags%LastRoot(i)  = total
       End Do
       Flags%TotalRoots      = total
    End If

    If (Flags%NumberOfCIRoots > 0) Then
       If (Flags%Algorithm == 2) Then
          If (Flags%PrintLevel > -9) Then
             Write(Stdout,'(//1X,"WARNING: cidir=2 cannot be used if roots of any symmetry are requested. Using cidir=1.")')
          End If
          Flags%Algorithm = 1
       Else If (Flags%Algorithm == -2) Then
          If (Flags%PrintLevel > -9) Then
             Write(Stdout,'(//1X,"WARNING: cidir=-2 cannot be used if roots of any symmetry are requested. Using cidir=-1.")')
          End If
          Flags%Algorithm = -1
       End If
    End If

    Allocate(Flags%RootSym(Flags%TotalRoots))
    Flags%RootSym = 0

    ! Allocate arrays for CI energy and CI coefficients:
    Allocate(Flags%CIE(Flags%TotalRoots))
    Allocate(Flags%CIC(Flags%TotalCSF, Flags%TotalRoots))
    Allocate(Flags%CIESAV(Flags%TotalRoots))
    Allocate(Flags%CICSAV(Flags%TotalCSF, Flags%TotalRoots))

    ! Print timings:
    If (Flags%PrintLevel >= 2) then
       Write(Stdout,'(//A)') ' Setting up done.'
       Call PrintTime(tuser, tsys, twall)
    End If

    ! Convert Gelfand states into linear combinations of Slater determinants:
    If (Flags%PrintLevel >= 3)  &
       Call PrintLinearCombinationsOfSDs(DRT, Flags, SymLabel)

    Return
  End Subroutine MakeShavittControl



  ! Construct the distinct row table.
  Subroutine SetupDistinctRowTable(DRT, iciref, Flags)
    Implicit None
    Type (PaldusRow),    Dimension(:),  Pointer :: DRT
    Integer,             Dimension(:,:)         :: iciref
    Type (ShavittControl)                       :: Flags
    ! End of dummy parameters.
    Type (PaldusRow),    Dimension(:),  Pointer :: oldDRT
    Type (NewPaldusRow), Dimension(:),  Pointer :: newRow
    Integer,             Dimension(:),  Pointer :: null
    Integer                                     :: Mult, Nelec, Norbs
    Integer                                     :: twoa, a, b, k
    Integer                                     :: rows, newrows, irow1, irow2, i, j, d
    Logical                                     :: rowflag, arcflag

    Mult  = Flags%Multiplicity
    Nelec = Flags%numberOfActiveElectrons
    Norbs = UBound(iciref,1)

    twoa  = Nelec - Mult + 1
    a     = twoa / 2
    b     = Mult - 1

    If (2*a /= twoa) Then
       Write(Stdout,'(1X,A)') 'Multiplicity inconsistent with number of electrons.'
       Stop 'SetupDistinctRowTable'
    End If

    ! This is a generous upper bound to the maximum number of rows in the DRT:
    ! At each orbital level there may be at most a*(a+b+1)+b+1 rows, because for
    ! each a there may be at most a+b+1 different values for b, except for the
    ! highest value of a, in which case there may be at most b+1 different
    ! values for b. At orbital levels 0 and Norbs, however, there is exactly
    ! one row. This estimate may be improved considerably.
    Allocate(DRT((Norbs-1)*(a*(a+b+1)+b+1)+2))
    Allocate(newRow(4*(a*(a+b+1)+b+1)))
    Allocate(null(UBound(iciref,2)))

    null  = 0
    rows  = 1
    irow1 = 1
    irow2 = 1
    Call InitRow(DRT(1), Norbs, a, b)

    Do k=Norbs-1, 0, -1
       newrows = 0
       Call FindNewRows(DRT, irow1, irow2, newRow, newrows)
       Call SortNewRows(newRow, newrows)

       i = 1
       Do While (i <= newrows)
          Call InitRow(DRT(rows+1), k, newRow(i)%a, newRow(i)%b)
          rowflag = .False.
          Do While (i <= newrows)
             If (newRow(i)%a /= DRT(rows+1)%a .Or. newRow(i)%b /= DRT(rows+1)%b) Exit
             j = newRow(i)%j
             d = newRow(i)%d
             arcflag = CheckArc(DRT, d, j, iciref, null, null, Flags%ExcitationLevel)
             If (arcflag) Then
                DRT(rows+1)%ju(d) = j
                DRT(j)%jd(d)      = rows+1
                rowflag = .True.
             End If
             i = i + 1
          End Do
          If (rowflag)  rows = rows + 1
       End Do

       irow1 = irow2 + 1
       irow2 = rows
    End Do

    Deallocate(newRow)
    oldDRT => DRT
    Allocate(DRT(rows))
    DRT = oldDRT(1:rows)
    Deallocate(oldDRT)

    Call ArcWeights(DRT, rows)

    Allocate(Flags%ExLev(DRT(1)%yd(4)))
    Allocate(Flags%WalkIrrep(DRT(1)%yd(4)))

    Call ExcLevelAndSymmetry(DRT, 1, iciref, null, 1, 1, Flags)

    Deallocate(null)

    Return
  End Subroutine SetupDistinctRowTable



  ! Initialize new row in the DRT.
  Subroutine InitRow(row, k, a, b)
    Type(PaldusRow) :: row
    Integer         :: k, a, b
    ! End of dummy parameters.

    row%k   = k
    row%a   = a
    row%b   = b
    row%ju  = 0
    row%jd  = 0
    row%yu  = 0
    row%yd  = 0

    Return
  End Subroutine InitRow



  ! Find all allowed rows at the next orbital level.
  Subroutine FindNewRows(DRT, irow1, irow2, newRow, newrows)
    Implicit None
    Type (PaldusRow),    Dimension(:) :: DRT
    Type (NewPaldusRow), Dimension(:) :: newRow
    Integer                           :: irow1, irow2, newrows
    ! End of dummy parameters.
    Integer                           :: irow, a, b, c

    Do irow=irow1, irow2
       a = DRT(irow)%a
       b = DRT(irow)%b
       c = DRT(irow)%k - a - b

       ! d = 0:
       If (c > 0) Then
          newrows = newrows + 1
          newRow(newrows)%a = a
          newRow(newrows)%b = b
          newRow(newrows)%d = 0
          newRow(newrows)%j = irow
       End If

       ! d = 1:
       If (b > 0) Then
          newrows = newrows + 1
          newRow(newrows)%a = a
          newRow(newrows)%b = b-1
          newRow(newrows)%d = 1
          newRow(newrows)%j = irow
       End If

       ! d = 2:
       If (a > 0 .And. c > 0) Then
          newrows = newrows + 1
          newRow(newrows)%a = a-1
          newRow(newrows)%b = b+1
          newRow(newrows)%d = 2
          newRow(newrows)%j = irow
       End If

       ! d = 3:
       If (a > 0) Then
          newrows = newrows + 1
          newRow(newrows)%a = a-1
          newRow(newrows)%b = b
          newRow(newrows)%d = 3
          newRow(newrows)%j = irow
       End If
    End Do

    Return
  End Subroutine FindNewRows



  ! Sort new rows according to their a and b values.
  Subroutine SortNewRows(newRow, newrows)
    Implicit None
    Type (NewPaldusRow), Dimension(:) :: newRow
    Integer                           :: newrows
    ! End of dummy parameters.
    Type (NewPaldusRow)               :: help
    Integer                           :: i, j, k

    Do i=1, newrows-1
       k = i
       Do j=i+1, newrows
          If (newRow(j)%a > newRow(k)%a .Or. (newRow(j)%a == newRow(k)%a) .And. newRow(j)%b > newRow(k)%b)  k = j
       End Do
       If (k /= i) Then
          help      = newRow(i)
          newRow(i) = newRow(k)
          newRow(k) = help
       End If
    End Do

    Return
  End Subroutine SortNewRows



  ! Check if an arc is part of a valid walk (with respect to the excitation level).
  Recursive Logical Function CheckArc(DRT, d, irow, iciref, itmpexc, icumexc, levexc) Result (CheckArcRes)
    Implicit None
    Type (PaldusRow), Dimension(:)                :: DRT
    Integer,          Dimension(:,:)              :: iciref
    Integer,          Dimension(:)                :: itmpexc, icumexc
    Integer                                       :: d, irow, levexc
    ! End of dummy parameters.
    Integer,          Dimension(UBound(iciref,2)) :: jtmpexc, jcumexc
    Integer,          Dimension(1)                :: ivec1, ivec2
    Integer                                       :: k, kocc

    If (irow == 0) Then
       CheckArcRes = .False.
       Return
    End If

    k       = DRT(irow)%k
    kocc    = d
    If (d > 1)  kocc = kocc - 1
    jtmpexc = itmpexc +     iciref(k,:) - kocc
    jcumexc = icumexc + Abs(iciref(k,:) - kocc)
    ivec1   = MinLoc(Abs(jtmpexc))
    ivec2   = MinLoc(jcumexc)
    If (Abs(jtmpexc(ivec1(1))) > levexc .Or. jcumexc(ivec2(1)) > 2*levexc) Then
       CheckArcRes = .False.
       Return
    End If

    If (irow == 1) Then
       CheckArcRes = .True.
       Return
    End If

    CheckArcRes = CheckArc(DRT, 0, DRT(irow)%ju(0), iciref, jtmpexc, jcumexc, levexc) .Or. &
                  CheckArc(DRT, 1, DRT(irow)%ju(1), iciref, jtmpexc, jcumexc, levexc) .Or. &
                  CheckArc(DRT, 2, DRT(irow)%ju(2), iciref, jtmpexc, jcumexc, levexc) .Or. &
                  CheckArc(DRT, 3, DRT(irow)%ju(3), iciref, jtmpexc, jcumexc, levexc)
    Return
  End Function CheckArc



  ! Calculate arc weights.
  Subroutine ArcWeights(DRT, rows)
    Implicit None
    Type(PaldusRow), Dimension(:) :: DRT
    Integer                       :: rows
    ! End of dummy parameters.
    Integer                       :: i, j, d, x

    ! Downward arc weights.
    DRT(rows)%yd(4) = 1
    Do i=rows-1, 1, -1
       DRT(i)%yd(0) = 0
       Do d=0, 3
          j = DRT(i)%jd(d)
          If (j /= 0) Then
             DRT(i)%yd(d+1) = DRT(i)%yd(d) + DRT(j)%yd(4)
             DRT(j)%yu(d)   = DRT(i)%yd(d)
          Else
             DRT(i)%yd(d+1) = DRT(i)%yd(d)
          End If
       End Do

       If (DRT(i)%yd(4) < 0) Then
          Write(Stdout,'(1X,A)') 'Number of Gelfand states exceeds range of type INTEGER.'
          Stop 'ArcWeights'
       End If
    End Do

    ! Upward vertex weights.
    DRT(1)%xu = 1
    Do i=2, rows
       x = 0
       Do d=0, 3
          j = DRT(i)%ju(d)
          If (j /= 0)  x = x + DRT(j)%xu
       End Do
       DRT(i)%xu = x
    End Do

    Return
  End Subroutine ArcWeights



  ! Compute excitation level and symmetry of each walk in the DRT.
  Recursive Subroutine ExcLevelAndSymmetry(DRT, irow, iciref, icumexc, iwalk, isym, flags)
    Implicit None
    Type(PaldusRow), Dimension(:)                :: DRT
    Integer,         Dimension(:,:)              :: iciref
    Integer,         Dimension(:)                :: icumexc
    Integer                                      :: irow, iwalk, isym
    Type(ShavittControl)                         :: flags
    ! End of dummy parameters.
    Integer,         Dimension(8,8)              :: MultiplicationTable
    Integer,         Dimension(UBound(iciref,2)) :: jcumexc
    Integer,         Dimension(1)                :: ivec
    Integer                                      :: d, j, kocc, jsym
    Data MultiplicationTable /  1, 2, 3, 4, 5, 6, 7, 8,  &
                             &  2, 1, 4, 3, 6, 5, 8, 7,  &
                             &  3, 4, 1, 2, 7, 8, 5, 6,  &
                             &  4, 3, 2, 1, 8, 7, 6, 5,  &
                             &  5, 6, 7, 8, 1, 2, 3, 4,  &
                             &  6, 5, 8, 7, 2, 1, 4, 3,  &
                             &  7, 8, 5, 6, 3, 4, 1, 2,  &
                             &  8, 7, 6, 5, 4, 3, 2, 1   /

    If (DRT(irow)%k == 0) Then
       ivec                   = MinLoc(icumexc)
       flags%ExLev(iwalk)     = icumexc(ivec(1))/2
       flags%WalkIrrep(iwalk) = isym
       Return
    End If

    Do d=0, 3
       j = DRT(irow)%jd(d)
       If (j == 0)  Cycle
       kocc = d
       If (d > 1)  kocc = kocc - 1
       jcumexc = icumexc + Abs(iciref(DRT(irow)%k,:) - kocc)
       If (d == 1 .Or. d == 2) Then
          jsym = MultiplicationTable(isym, flags%MOIrreps(DRT(irow)%k))
       Else
          jsym = isym
       End If
       Call ExcLevelAndSymmetry(DRT, j, iciref, jcumexc, iwalk+DRT(irow)%yd(d), jsym, flags)
    End Do

    Return
  End Subroutine ExcLevelAndSymmetry



  ! Print the distinct row table.
  Subroutine PrintDRT(DRT)
    Implicit None
    Type (PaldusRow), Dimension(:), Intent(In) :: DRT  ! Distinct Row Table
    ! End of dummy parameters.
    Integer  :: j

    Write(Stdout,'(//A/)') " Distinct row table:"
    Write(Stdout,fmt=('("                   Paldus row     downward chaining index            arc weights            lower")'))
    Write(Stdout,fmt=('(" orbital   row    -------------    ---------------------     ---------------------------    walks")'))
    Write(Stdout,fmt=('("  level   index   a     b     c    d=0   d=1   d=2   d=3     d=0     d=1     d=2     d=3      y")'))
    Write(Stdout,*)
    Do j=1,Ubound(DRT,1)
       Write(Stdout,fmt=('(I6,I7,7I6,5I8)')) &
         &   DRT(j)%k, j, DRT(j)%a, DRT(j)%b, DRT(j)%k - DRT(j)%a - DRT(j)%b, DRT(j)%jd, DRT(j)%yd
    End Do
  End Subroutine PrintDRT



  ! Check if the Gelfand states are in the reference space and have the desired
  ! symmetry. Write a zero into the index vector otherwise and set the correct
  ! values for the related fields in the ShavittControl structure.
  !
  Subroutine InitIndexVector(IndexVector, Flags)
    Implicit None
    Integer, Dimension(:)     :: IndexVector
    Type(ShavittControl)      :: Flags
    ! End of dummy parameters.
    Integer                   :: ndrt, n, i, j

    IndexVector    = 0
    Flags%FirstCSF = 0
    Flags%LastCSF  = 0
    Flags%NumCSF   = 0
    ndrt           = Ubound(IndexVector,1)
    n              = 0

    Do i=1, Flags%NumSym
       If (Flags%NumberOfCIRoots > 0 .Or. Flags%NRoots(i) > 0) Then
          Flags%FirstCSF(i) = n + 1
          Do j=1, ndrt
             If (Flags%WalkIrrep(j) == i .And. Flags%ExLev(j) <= Flags%ExcitationLevel) Then
                n = n + 1
                IndexVector(j) = n
             End If
          End Do
          Flags%LastCSF(i) = n
          Flags%NumCSF(i)  = n - Flags%FirstCSF(i) + 1
       End If
    End Do

    ! Determine lexical index of each CSF:
    Allocate(Flags%WalkIndex(n))
    Do i=1, ndrt
       j = IndexVector(i)
       If (j > 0)  Flags%WalkIndex(j) = i
    End Do

    Flags%TotalCSF = n

    Return
  End Subroutine InitIndexVector



  ! Public routine to deallocate all dynamic memory.
  Subroutine Cleanup(Flags)
    Implicit None
    Type(ShavittControl), Intent(InOut) :: Flags
    ! End of dummy parameters.

    If (Flags%PrintLevel >= 2) Then
       Write(Stdout,'(//1X,"Cleaning up.")')
    End If

    ! Arrays are deallocated in the order in which
    ! they are declared in Type(ShavittControl).
    Deallocate(Flags%DRT)
    Deallocate(Flags%IndVec)
    Deallocate(Flags%SymVec)

    If (Flags%Algorithm == 1) Then
       Call FreeGenList(Flags%OneList)
       Call FreeGenList(Flags%TwoList)
    Else If (Flags%Algorithm == -1) Then
       Deallocate(Flags%OneGen)
       Deallocate(Flags%TwoGen)
       Deallocate(Flags%OneIndex)
       Deallocate(Flags%TwoIndex)
    End If

    If (Flags%Algorithm > 0) Then
       Deallocate(Flags%UpperPart)
       Deallocate(Flags%LowerPart)
       Deallocate(Flags%OneLoop)
       Deallocate(Flags%TwoLoop)
       Deallocate(Flags%lexUpper)
       Call CleanupCombinations(Flags)
    End If

    Deallocate(Flags%CIE)
    Deallocate(Flags%CIC)
    Deallocate(Flags%CIESAV)
    Deallocate(Flags%CICSAV)
    Deallocate(Flags%WalkIndex)
    Deallocate(Flags%RootSym)
!   Flags%WhoIsWho simply points to IMOCI, so don't deallocate.
!   Deallocate(Flags%WhoIsWho)
    Deallocate(Flags%MOIrreps)
    Deallocate(Flags%WalkIrrep)
    Deallocate(Flags%ExLev)
    Deallocate(Flags%iFirstRef)
    Deallocate(Flags%iLastRef)
    Deallocate(Flags%iRefWalk)
    Deallocate(Flags%iRefConf)

    Return
  End Subroutine Cleanup


End Module gugasetup
