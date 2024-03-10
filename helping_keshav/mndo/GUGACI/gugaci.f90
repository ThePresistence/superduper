!
! This is the interface between the GUGACI module and MNDO99.
!
! First version by Michael E. Beck in 1998 at Zurich University.
! Final version by Axel Koslowski in 2000-2005 at MPI Muelheim.
! Shape-driven algorithm added by Axel Koslowski in 2012-2013.
! Reorganized argument list and removed memory leaks in 2014.
!

Subroutine GugaCI(ICALL, ICALLG, LM2, LM3, LM4, LM6, NUMAT, NUMB, NCIO, NCIGAM,      &
                  IUVCD, IMCD, IPOP, INATUR, ISTATE, JSTATE, NUMLAB, NUMGPU, NB6,    &
                  IN2, IMOCI, ICIREF, IROOTA, EE, ENUCLR, CMO, FAO, WAO, OCCR, NSYM, &
                  COORD, NAT, NFIRST, NLAST, CORE, ZS, ZP, ZD, III, IIID, ICISYM,    &
                  ONEDEN, TWODEN, CNO, PCIAO, CIPROP, XNAC, LABELS, GRPLAB, NVCAPA)
  !
  ! Please refer to ../99/GUGA.f for the definition of ICALLG.
  !
  ! The following comments address the connection between the dummy parameters
  ! of the GugaCI subroutine and the variable names in the standard MNDO99 code.
  ! The GUGA-CI code has been designed to be modular and essentially independent
  ! from the calling program. Data is exchanged through the argument list of the
  ! driver routine GugaCI which is called from subroutines GUGA and CLRGUG.
  !
  ! Survey of the argument list (I=input, O=output).
  !
  ! Argument in GugaCI    /     in SCFCAL   I/O   Remarks
  !
  ! ICALL                        (same)      I    See SUBROUTINE SCF  for definition.
  ! ICALLG                       (same)      I    See SUBROUTINE GUGA for definition.
  ! LM2                          (same)      I    LM2=LM3 or LM2=LM3+1.
  ! LM3                          (same)      I    Number of molecular orbitals.
  ! LM4                          (same)      I    Dimension of a matrix triangle.
  ! LM6                          (same)      I    Number of unique one-center AO pairs.
  ! NUMAT                        (same)      I    Number of atoms.
  ! NUMB                         (same)      I    Index of highest occupied MO.
  ! NCIO                         (same)      I    Number of active orbitals.
  ! NCIGAM                       (same)      I    Size of two-particle density matrix.
  ! IUVCD, IMCD, IPOP, INATUR                I    Control flags for the calculation of properties.
  ! ISTATE, JSTATE               (same)      I    States for the non-adiabatic coupling vectors.
  ! NUMLAB                       IRREPN      I    Number of irreducible representations in this point group.
  ! NUMGPU                       (same)      I    Number of CUDA-capable GPU compute devices.
  ! NB6                          (same)      I    Output file number.
  ! IN2                          (same)     I/O   Standard MNDO99 input options; Array (300).
  ! IMOCI                        (same)      I    Indices of orbitals in the active space; Array (LMACT).
  ! ICIREF                       (same)      I    Reference configurations for the construction of the DRT; Array (LMACT,LMREF);
  ! IROOTA                       (same)      I    Requested number of CI roots of each irrep; Array (8).
  ! EE                           (same)     I/O   SCF+CI energy (in eV) for state LROOT.
  ! ENUCLR                       (same)      I    Nuclear energy (in eV).
  ! CMO                          C           I    MO coefficient matrix; Array (LM2,LM3).
  ! FAO                          F           I    AO basis Fock matrix; Array (LM4).
  ! WAO                          W           I    Square array (LM6,LM6), in eV.
  ! OCCR                         (same)      I    Orbital occupation numbers; Array (LM3).
  ! NSYM                         (same)      I    Symmetry labels for MOs; Array (LM3).
  ! COORD                        (same)      I    Atom coordinates; Array (3,NUMAT).
  ! NAT                          (same)      I    Nuclear charges; Array (NUMAT).
  ! NFIRST, NLAST                (same)      I    Indices of first/last AO at each atomic center; Array (NUMAT).
  ! ZS, ZP, ZD                   (same)      I    Atomic orbital exponents per element; Array (LMZ).
  ! CORE                         (same)      I    Core charge (nucleus plus core electrons) per element; Array (LMZ).
  ! III, IIID                    (same)      I    Principle quantum numbers of AOs per element; Array (LMZ).
  ! ICISYM                       (same)      O    Symmetry of the computed CI states; Array (LMSTAT).
  ! ONEDEN                       A(LCI1)     O    One-particle density matrix; Array (LMACT,NCIO).
  ! TWODEN                       A(LCI2)     O    Two-particle density matrix; Array (NCIGAM).
  ! CNO                          A(LCI3)     O    Natural orbital MO coefficients; Array (LM3,LM3).
  ! PCIAO                        A(LCI4)     O    CI-corrected AO-basis density matrix; Array (LM3,LM3).
  ! CIPROP                       (same)      O    Array of calculated properties; see subroutine Properties
  !                                               (in property.f90) for further explanation; Array (LMPROP,LMSTAT).
  ! XNAC                         (same)      O    Numeric non-adiabatic coupling matrix elements; Array (LMGRD,LMGRD).
  ! LABELS                       IRREP       I    List of symmetry labels; Array (1:NUMLAB).
  ! GRPLAB                       NGROUP      I    Point group label.
  ! NVCAPA                       (same)      I    Compute capabilites (10*major+minor) of NVIDIA GPUs; Array (NUMGPU).
  !
  ! The entries IN2, IMOCI, and ICIREF on the above list are documented in more detail below.
  ! They are also described in the standard MNDO99 input description, section 3.8.
  !
  ! IMOCI(LMACT) is a linear array with dimension LMACT=150 which contains the indices
  ! of the orbitals in the active CI space. The numbering refers to the SCF output.
  ! a) Default case (movo=0): There are NACTIV=ICI1+ICI2 active orbitals ranging
  !    from NUMB-ICI1+1 to NUMB+ICI2 (ascending order with regard to their energies).
  !    In this case, IMOCI(1) = NUMB-ICI1+1, ... , IMOCI(NACTIV) = NUMB+ICI2.
  ! b) Explicit input (movo=1,2): The array (IMOCI(I),I=1,NACTIV) is either read
  !    in completely, or is filled by replacing orbital numbers in the default array.
  !    In either case, IMOCI may then contain the numbers of the active orbitals in
  !    arbitrary order.
  ! c) Explicit input according to symmetry (movo=3). Two lines of input are needed.
  !    The first line contains the MO indices counting only MOs of a given symmetry
  !    which is specified in the second line.
  ! The active MOs (and the corresponding reference configurations) may be given in
  ! any order since they are sorted before the construction of the DRT as described
  ! in ../99/GUGA.f.
  !
  ! ICIREF(LMACT,LMREF) is a rectangular array with dimension LMACT=150 and LMREF=100.
  ! It contains the reference configurations which are used to truncate the distinct row
  ! table which describes the GUGA-CI many-electron basis functions (Gelfand states).
  ! ICIREF is not needed for full CI (nciref=0), but for any meaningful truncated CI
  ! (nciref>0). There are five options for defining ICIREF in the latter case.
  ! a) By default (mciref=0). Implicitly defined reference configurations refer to
  !    an active MO list sorted by increasing MO indices.
  ! b) By explicit input of ICIREF (mciref=1).
  ! c) By implicit input of ICIREF (mciref=2) using excitation indices relative
  !    to the SCF determinant.
  ! d) mciref=3 is similar to a), but after the calculation of all requested CI roots
  !    it is made sure by appending additional reference configurations to ICIREF that
  !    the corresponding GUGA-CI many-electron basis functions (Gelfand states) represent
  !    a sufficiently large percentage of each root (specified by option CISELT).
  ! e) mciref=4 is similar to b) and d).
  ! The order of the active orbitals is identical in IMOCI and ICIREF (first index).
  !
  Use Limit,      Only: LMZ, LMACT, LMREF, LMPROP, LMSTAT, LMGRD, LMCONF
  Use gugaglobal, Only: ShavittControl, Generator, Stdout
  Use gugasetup
  Use gugaincore
  Use gugadirect
  Use gugashape
  Implicit None
  Integer                                    :: ICALL, ICALLG, LM2, LM3, LM4, LM6, NUMAT, NUMB, NCIO, NCIGAM
  Integer                                    :: IUVCD, IMCD, IPOP, INATUR, ISTATE, JSTATE, NUMLAB, NUMGPU, NB6
  Integer,          Dimension(300)           :: IN2
  Integer,          Dimension(LMACT)         :: IMOCI
  Integer,          Dimension(LMACT,LMREF)   :: ICIREF
  Integer,          Dimension(8)             :: IROOTA
  Double Precision                           :: EE, ENUCLR
  Double Precision, Dimension(LM2,LM3)       :: CMO
  Double Precision, Dimension(LM4)           :: FAO
  Double Precision, Dimension(LM6,LM6)       :: WAO
  Double Precision, Dimension(LM3)           :: OCCR
  Integer,          Dimension(LM3)           :: NSYM
  Double Precision, Dimension(3,NUMAT)       :: COORD
  Integer,          Dimension(NUMAT)         :: NAT
  Integer,          Dimension(NUMAT)         :: NFIRST, NLAST
  Double Precision, Dimension(LMZ)           :: CORE, ZS, ZP, ZD
  Integer,          Dimension(LMZ)           :: III, IIID
  Integer,          Dimension(LMSTAT)        :: ICISYM
  Double Precision, Dimension(LMACT,NCIO)    :: ONEDEN
  Double Precision, Dimension(NCIGAM)        :: TWODEN
  Double Precision, Dimension(LM3,LM3)       :: CNO, PCIAO
  Double Precision, Dimension(LMPROP,LMSTAT) :: CIPROP
  Double Precision, Dimension(LMGRD,LMGRD)   :: XNAC
  Character(Len=3)                           :: GRPLAB
  Character(Len=4), Dimension(NUMLAB)        :: LABELS
  Integer,          Dimension(NUMGPU)        :: NVCAPA
  ! End of dummy parameters.
  Integer                                    :: ici1, ici2, ioutci
  Integer                                    :: nciref, nactiv, ncigrd, i, j, k
  Double Precision                           :: DDOT, ciselt
  ! The central control structure is kept in memory via the save attribute:
  Type(ShavittControl), Save                 :: Flags

  ! MD INTERFACE
  INTEGER :: nciconf
  REAL*8  :: cicomp
  INTEGER :: icross
  INTEGER :: nsav15

  ! Note warning about bounds-checking below
  COMMON /DYNVA/ cicomp(LMCONF,6),nciconf
  ! END MD INTERFACE

  ! Standard output is redirected to NB6. Its value is stored as Stdout for convenience.
  Stdout = NB6

  ! Some input flags.
  ici1   = IN2(131)
  ici2   = IN2(132)
  ioutci = IN2(133)
  nciref = IN2(136)
  ciselt = 1D-2 * DBLE(IN2(155))
  nactiv = ici1 + ici2

  ! Debug print: Input data provided to GUGACI module.
  IF(ioutci.GE.5) THEN
     WRITE(NB6,500)
     WRITE(NB6,510) ICALL, ICALLG, LM2, LM3, LM4, LM6, NUMAT, NUMB, NCIO, NCIGAM, &
                    IUVCD, IMCD, IPOP, INATUR, ISTATE, JSTATE, NUMLAB, NUMGPU, NB6
     WRITE(NB6,520) IN2(131:170)
     DO I=1, NACTIV, 20
        K = MIN(I+19, NACTIV)
        WRITE(NB6,530) I, K, IMOCI(I:K)
     END DO
     DO J=1, NCIREF
        DO I=1, NACTIV, 20
           K = MIN(I+19, NACTIV)
           WRITE(NB6,540) I, K, J, ICIREF(I:K, J)
        END DO
     END DO
     WRITE(NB6,550) IROOTA(1:8)
     WRITE(NB6,560) EE, ENUCLR
     DO I=1, LM3, 10
        K = MIN(I+ 9, LM3)
        WRITE(NB6,570) I, K, OCCR(I:K)
     END DO
     DO I=1, LM3, 20
        K = MIN(I+19, LM3)
        WRITE(NB6,580) I, K, NSYM(I:K)
     END DO
     WRITE(NB6,600)
     WRITE(NB6,610)
     DO I=1, NUMAT
        WRITE(NB6,620) I, NAT(I), COORD(:,I)
     END DO
  ENDIF

  Select Case (ICALLG)
  Case (1)
     If (ioutci > 0) Then
        Write(NB6,'(//1X,A,A)') 'Point group used in GUGA-CI calculation:  ', Trim(GRPLAB)
     End If

     ! Initialize ShavittControl structure and construct DRT.
     Call MakeShavittControl(Flags, Flags%DRT, IN2, LM3, NUMB, &
                             ICIREF(1:nactiv,1:nciref), IMOCI(1:nactiv), &
                             IROOTA, NSYM, LABELS, GRPLAB, NUMGPU, NVCAPA)

     ! Calculate upper and lower partial loops
     ! if shape-driven algorithm is requested.
     If (Flags%Algorithm > 0) Then
        Call CalculatePartialLoops(Flags, OCCR)
        Call InitSymCombinations(Flags)
     End If

     ! Calculate generator matrix elements to be stored all the time.
     If (Flags%Algorithm == 1) Then
        Call ShapeDrivenGenerators(Flags, 0)
     Else If (Flags%Algorithm == -1) Then
        Call CalculateGenerators(Flags, OCCR, 0)
     End If
  Case (2)
     ! Solve eigenvalue problem and add references if requested (mciref=3,4).
     If (IAbs(Flags%Algorithm) == 1 .Or. IAbs(Flags%Algorithm) == 2) Then
        Call InCoreCI(Flags, EE, ENUCLR, ciselt, ICIREF(1:nactiv,:), CMO, FAO, WAO, &
                      NFIRST, NLAST, OCCR, LABELS, CIPROP, ICISYM, ICALL)
     Else If (IAbs(Flags%Algorithm) == 3) Then
        Call DirectCI(Flags, ICIREF(1:nactiv,:), LABELS, CMO, OCCR, FAO, WAO, &
                      NFIRST, NLAST, EE, ENUCLR, ciselt, CIPROP, ICISYM, ICALL)
     Else
        Write(NB6,'(1X,A,I2,A)') 'Unexpected value for cidir (', Flags%Algorithm, ').'
        Stop 'GugaCI'
     End If
     ! Store new number of references into IN2.
     IN2(136) = Flags%nciref
  Case (3)
     ! Calculate one- and two-particle (transition) densities for the analytic gradient and NAC vectors.
     If (ISTATE > 0) Then
        Flags%iState = ISTATE
        Flags%jState = ISTATE
     End If
     If (JSTATE > 0) Then
        Flags%jState = JSTATE
     End If
     !
     ! Two-electron generator matrix elements for which bra and ket Gelfand states
     ! have different spatial symmetry are never stored in-core. Therefore,
     ! two-particle transition densities for non-adiabatic coupling vectors must
     ! always be calculated in a direct manner for point groups higher than C1.
     ! For point group C1, subroutines for the calculation of the two-particle
     ! density matrix are equivalent for permanent and transition densities.
     ! For higher point groups, transition density computation is more general
     ! than permanent density calculation, but the latter is more efficient.
     !
     If (IAbs(Flags%Algorithm) == 1  .And.  (Flags%iState == Flags%jState  .Or.  NUMLAB == 1)) Then
        Call InCoreGrad(Flags, ONEDEN(1:nactiv,1:nactiv), TWODEN)
     Else
        Call DirectGrad(Flags, ONEDEN(1:nactiv,1:nactiv), TWODEN, OCCR)
     End If
  Case (4)
     ! Calculate spectroscopic properties.
     If (IAbs(Flags%Algorithm) == 1) Then
        Call InCoreProp(Flags, NAT, COORD, NFIRST, NLAST, CORE, ZS, ZP, ZD, III, IIID, CMO, OCCR, &
                        LABELS, CIPROP, ICISYM, PCIAO, IUVCD, IMCD, IPOP, INATUR, ICALL)
     Else
        Call DirectProp(Flags, NAT, COORD, NFIRST, NLAST, CORE, ZS, ZP, ZD, III, IIID, CMO, OCCR, &
                        LABELS, CIPROP, ICISYM, PCIAO, IUVCD, IMCD, IPOP, INATUR, ICALL)
     End If
  Case (5)

     ! MD INTERFACE
     icross = IN2(160)
     IF (icross .EQ. 6) THEN
        nciconf = Flags%TotalCSF

       ! Presumably DYNVA will be removed when this interface is formalised.
       ! In the meantime these tests need to be kept in sync with the hard-coded 
       ! limits of cicomp.
       If (nciconf > LMCONF) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: nciconf exceeds bounds of cicomp in DYNVA'
          Stop 'GugaCI'
       End If
       If (Flags%TotalRoots > 6) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: TotalRoots exceeds bounds of cicomp in DYNVA'
          Stop 'GugaCI'
       End If

        cicomp=0.D0
        do i = 1, Flags%TotalRoots
           do j = 1, Flags%TotalCSF
              cicomp(j,i) = Flags%CIC(j,i)
           end do
        end do
     END IF
     ! END MD INTERFACE

     ! Calculate numeric non-adiabatic coupling element.
     ncigrd = IN2(159)
     Do i=1, ncigrd
        Do j=1, ncigrd
           XNAC(i,j) = DDOT(Flags%TotalCSF, Flags%CICSAV(1,i), 1, Flags%CIC(1,j), 1) - &
                       DDOT(Flags%TotalCSF, Flags%CICSAV(1,j), 1, Flags%CIC(1,i), 1)
        End Do
     End Do
  Case (6)
     ! MULSAV INTERFACE
     ! TODO: Formalise this interface (no common blocks, no hard-coded 
     ! TODO: limit to cicomp) as it now forms part of the external 
     ! TODO: interface to ChemShell through nsav15=6/7
     icross = IN2(160)
     nsav15 = IN2(18)
     IF ((icross .EQ. 1 .OR. icross .EQ. 2) .AND. (nsav15 .EQ. 6 .OR. nsav15 .EQ. 7)) THEN
        nciconf = Flags%TotalCSF

       ! Presumably DYNVA will be removed when this interface is formalised.
       ! In the meantime these tests need to be kept in sync with the hard-coded 
       ! limits of cicomp.
       If (nciconf > 10000) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: nciconf exceeds bounds of cicomp in DYNVA'
          Stop 'GugaCI'
       End If
       If (Flags%TotalRoots > 6) Then
          Write(Stdout,'(/A)') ' FATAL INTERNAL ERROR: TotalRoots exceeds bounds of cicomp in DYNVA'
          Stop 'GugaCI'
       End If

        cicomp=0.D0
        do i = 1, Flags%TotalRoots
           do j = 1, Flags%TotalCSF
              cicomp(j,i) = Flags%CIC(j,i)
           end do
        end do
     END IF
     ! END MULSAV INTERFACE

     ! Save CI coefficients.
     Flags%CICSAV = Flags%CIC
  Case (7)
     ! Clear all CI data.
     Call Cleanup(Flags)
  Case Default
     Write(StdOut,'(1X,"Invalid call number: ",I2)') ICALLG
  End Select

  Return
  500 FORMAT(// 1X,"DEBUG PRINT: INPUT DATA PROVIDED TO GUGACI.")
  510 FORMAT(/  1X,"ICALL               =",I10, &
             /  1X,"ICALLG              =",I10, &
             /  1X,"LM2                 =",I10, &
             /  1X,"LM3                 =",I10, &
             /  1X,"LM4                 =",I10, &
             /  1X,"LM6                 =",I10, &
             /  1X,"NUMAT               =",I10, &
             /  1X,"NUMB                =",I10, &
             /  1X,"NCIO                =",I10, &
             /  1X,"NCIGAM              =",I10, &
             /  1X,"IUVCD               =",I10, &
             /  1X,"IMCD                =",I10, &
             /  1X,"IPOP                =",I10, &
             /  1X,"INATUR              =",I10, &
             /  1X,"ISTATE              =",I10, &
             /  1X,"JSTATE              =",I10, &
             /  1X,"NUMLAB              =",I10, &
             /  1X,"NUMGPU              =",I10, &
             /  1X,"NB6                 =",I10)
  520 FORMAT(   1X,"IN2   (131:140)     =",10I5, &
             /  1X,"IN2   (141:150)     =",10I5, &
             /  1X,"IN2   (151:160)     =",10I5, &
             /  1X,"IN2   (161:170)     =",10I5)
  530 FORMAT(   1X,"IMOCI (",I3,":",I3,")     =",20I5)
  540 FORMAT(   1X,"ICIREF(",I3,":",I3,",",I3,") =",20I5)
  550 FORMAT(   1X,"IROOTA(  1:  8)     =",8I5)
  560 FORMAT(   1X,"EE                  =",F20.10, &
             /  1X,"ENUCLR              =",F20.10)
  570 FORMAT(   1X,"OCCR  (",I3,":",I3,")     =",10F10.5)
  580 FORMAT(   1X,"NSYM  (",I3,":",I3,")     =",20I5)
  600 FORMAT(// 1X,"QM ATOMS.")
  610 FORMAT(/  4X,"I",6X,"NAT",15X,"X",24X,"Y",24X,"Z")
  620 FORMAT(   1X,I4,I8,3F25.15)
End Subroutine GugaCI
