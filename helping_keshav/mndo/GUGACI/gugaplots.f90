!
! Plot Shavitt graphs as postscript files using the pgplot library and
! print a linear combination of Slater determinants for each Gelfand state.
!
! First version by Michael E. Beck in 1998 at Zurich University.
!
! Final version by Axel Koslowski in 2000-2003 at MPI Muelheim.
!

Module gugaplots

  Use gugaglobal

  Implicit None

  Private

  Public :: PrintLinearCombinationsOfSDs  ! Print a linear combination of Slater determinants for each Gelfand state.
  Public :: PlotShavittGraph              ! Write postscript file with a plot of the Shavitt graph.
  Public :: PlotLoops                     ! Generate a plot for each pair of Gelfand states.


  ! Special types (Michael Beck):

  Type ShavittWalk
     Type(PaldusRow), Pointer :: row     ! some row in the DRT
     Integer                  :: vertex  ! the lexical index of that row
     Integer                  :: d       ! step number
  End Type ShavittWalk


  Type MO
     Integer :: index
     Integer :: spin   ! 1 = alpha, 0 = beta
     Integer :: a,b,d
     Logical :: Open ! true for singly occupied MOs
  End Type MO


  ! Constants:

  Double Precision, Parameter :: Tiny = 1.0d-6

Contains

  !
  ! Plotting routines generating postscript files.
  !

  ! Public routine to plot Shavitt graph (Michael Beck).
  Subroutine PlotShavittGraph(DRT, refConf, ShavittFlags)
    Type(PaldusRow),  Dimension(:)         :: DRT
    Type(StepVector), Dimension(:)         :: refConf
    Type(ShavittControl)                   :: ShavittFlags
    ! End of dummy parameters.

    Integer                                :: Nelectrons,excitationLevel,Multiplicity 
    Integer                                :: tailindex, i, j, norbs, help,d, nconf 

    Integer                                :: oldvertex, nextvertex, k, kmin, kmax, sockel
    Double Precision                       :: fonthoehe, scale
    Double Precision, Dimension(Size(DRT)) :: X, Y
    Double Precision                       :: h1,h2
    Character(Len=3)                       :: zahlstring, zahlstring2, zahlstring3
    Character(Len=6)                       :: bigzahlstring

    Nelectrons      = ShavittFlags%numberOfActiveElectrons
    excitationLevel = ShavittFlags%excitationLevel
    Multiplicity    = ShavittFlags%Multiplicity 
    sockel          = 0

    tailindex = Ubound(DRT,1)
    norbs = DRT(1)%k
    nconf = Ubound(refconf,1)

    help = Step(3)*DRT(1)%a + Step(2)*DRT(1)%b
    Call PGBEG(0,'shavitt.graph.ps/VPS',1,1)                  !ifpgplot!
    h1 = help  + 2 * step(3)
    h2 = norbs + 2
    Call PGENV(0,Real(h1),0,Real(h2),0,-1)                    !ifpgplot!

    ! set roman font
    Call PGSCF(2)                                             !ifpgplot!

    ! set fontheight (PGPLOT default: 1/40 of pageheight)
    If(norbs > 8) Then
       scale = 1.d0/(4.d0*h2)
    Else
       scale = 1.d0/40.d0
    End If
    Call PGSCH(40.*Real(scale))                               !ifpgplot!
    fonthoehe = scale * h2
  
    Call numberToString(Nelectrons, zahlstring)
    Call numberToString(norbs, zahlstring2)
    Call numberToString(Multiplicity, zahlstring3)
    Call PGText(0,Real(h2-1.2*fonthoehe),"     "            & !ifpgplot!
         &      // zahlstring // " electrons in "           & !ifpgplot!
         &      // zahlstring2 //                           & !ifpgplot!
         &      " orbitals, multiplicity: " // zahlstring3)   !ifpgplot!
    Call bigNumberToString(DRT(1)%yd(4),bigzahlstring)
    Call numberToString(nconf,zahlstring2)
    Call numberToString(excitationlevel,zahlstring3)
    Call PGText(0,Real(h2-fonthoehe*2.4),"  "               & !ifpgplot!
         &      // bigzahlstring //                         & !ifpgplot!
         &      " Gelfand states with resp. to "            & !ifpgplot!
         &      // zahlstring2 // " reference CSFs")          !ifpgplot!
    Call PGText(0,Real(h2-fonthoehe*3.6),                   & !ifpgplot!
         &      "                              " //         & !ifpgplot!
         &      "max. excitation level: " // zahlstring3)     !ifpgplot!

    !construct coordinates
    X(1) = step(3)
    Y(1) = h2 - 1.d0
    kmin = 2
    Do i=1,tailindex 
       Do d=0,3
          help = DRT(i)%jd(d)
          If(help > 0) Then
             X(help) = X(i) + step(d) 
             Y(help) = Y(i) - 1.d0
          End If
       End Do
       If(i > 1 .And. DRT(i)%k < DRT(i - 1)%k) Then
          kmax = i - 1
          Do j=kmin,kmax
             Do k=kmin, j-1
                If(Abs(X(j)-X(k)) .Lt. Tiny) X(j) = X(j) + 1.d0
             End Do
          End Do
          Call numberToString(kmin,zahlstring)
          Call PGPTXT(Real(X(kmin)-fonthoehe),              & !ifpgplot!
               &      Real(Y(kmin)-fonthoehe),              & !ifpgplot!
               &      0.0,1.0,zahlstring)                     !ifpgplot!
          Call numberToString(kmax,zahlstring)
          Call PGPTXT(Real(X(kmax)+fonthoehe),              & !ifpgplot!
               &      Real(Y(kmax)),0.0,0.0,zahlstring)       !ifpgplot!
          Call numberToString(DRT(i - 1)%k+sockel,zahlstring)
          Call PGPTXT(Real(h1-fonthoehe),Real(Y(kmax)),     & !ifpgplot!
               &      0.0,1.0,zahlstring)                     !ifpgplot!
          kmin = i
       End If
    End Do
    Call numberToString(tailindex,zahlstring)
    Call PGPTXT(Real(X(tailindex)+fonthoehe),               & !ifpgplot!
         &      Real(Y(tailindex)),0.0,0.0,zahlstring)        !ifpgplot!

    !plot all lines
    Do i=2,tailindex
       Do d=0,3
          help = DRT(i)%ju(d)
          If(help > 0) Then
             Call PGLINE(2,(/Real(X(i)),Real(X(help))/),    & !ifpgplot!
                  & (/Real(Y(i)),Real(Y(help))/))             !ifpgplot!
          End If
       End Do
    End Do

    !plot vertixes
    Call PGPT(tailindex,Real(X),Real(Y),9)                    !ifpgplot!

    !plot reference configurations as walks with fat lines
    Call PGSLW(4)                                             !ifpgplot!
    Do i=1,nconf
       oldvertex = 1
       Do j=norbs,1,-1
          nextvertex = DRT(oldvertex)%jd(refconf(i)%d(j))
          Call PGLINE(2,                                    & !ifpgplot!
               &      Real((/X(oldvertex),X(nextvertex)/)), & !ifpgplot!
               &      Real((/Y(oldvertex),Y(nextvertex)/)))   !ifpgplot!
          oldvertex = nextvertex
       End Do
    End Do

    Call PGEND                                                !ifpgplot!
  End Subroutine PlotShavittGraph



  ! Public driver routine to plot pairs of Gelfand states corresponding
  ! to the upper triangle of the CI Hamiltonian.
  Subroutine PlotLoops(DRT, ShavittFlags)
    Implicit None
    Type(PaldusRow),   Dimension(:), Pointer     :: DRT
    Type(ShavittControl)                         :: ShavittFlags
    ! End of dummy parameters.
    Type(ShavittWalk), Dimension(:), Allocatable :: walk1, walk2

    Allocate(walk1(DRT(1)%k+1))
    Allocate(walk2(DRT(1)%k+1))
    Call ConstructWalks1(DRT, ShavittFlags, walk1, walk2, 1, 0)
    Deallocate(walk2)
    Deallocate(walk1)

    Return
  End Subroutine PlotLoops



  ! Construct upper walk by following downwards all existing arcs.
  ! If we reach the tail of the graph, we have a diagonal element.
  ! If plots of off-diagonal elements are requested (ciplot > 2),
  ! open a loop at each vertex if possible.
  !
  Recursive Subroutine ConstructWalks1(DRT, ShavittFlags, walk1, walk2, irow, y)
    Implicit None
    Type(PaldusRow),   Dimension(:), Pointer :: DRT
    Type(ShavittControl)                     :: ShavittFlags
    Type(ShavittWalk), Dimension(:)          :: walk1, walk2
    Integer                                  :: irow, y
    ! End of dummy parameters.
    Integer                                  :: k, m, d, d2

    If (irow == 0) Return

    k = DRT(irow)%k
    walk1(k+1)%row    => DRT(irow)
    walk1(k+1)%vertex =  irow
    walk2(k+1)%row    => DRT(irow)
    walk2(k+1)%vertex =  irow

    If (k == 0) Then
       m = ShavittFlags%IndVec(y+1)
       If (m/=0) Call PlotLoop(DRT, walk1, walk1, m, m)  !ifpgplot!
    Else
       Do d=0, 3
          walk1(k+1)%d = d
          walk2(k+1)%d = d
          Call ConstructWalks1(DRT, ShavittFlags, walk1, walk2, DRT(irow)%jd(d), y+DRT(irow)%yd(d))

          If (ShavittFlags%plot > 2) Then
             Do d2=d+1, 3
                walk2(k+1)%d = d2
                Call ConstructWalks2(DRT, ShavittFlags, walk1, walk2,       &
                                     DRT(irow)%jd(d),   DRT(irow)%jd(d2),   &
                                     DRT(irow)%yd(d)+y, DRT(irow)%yd(d2)+y)
             End Do
          End If
       End Do
    End If

    Return
  End Subroutine ConstructWalks1



  ! At each orbital level we consider two vertices, one for the bra walk and one for the ket walk.
  ! We follow downwards all existing combinations of arcs emerging from these vertices.
  Recursive Subroutine ConstructWalks2(DRT, ShavittFlags, walk1, walk2, ibra, iket, ybra, yket)
    Implicit None
    Type(PaldusRow),   Dimension(:), Pointer :: DRT
    Type(ShavittControl)                     :: ShavittFlags
    Type(ShavittWalk), Dimension(:)          :: walk1, walk2
    Integer                                  :: ibra, iket, ybra, yket
    ! End of dummy parameters.
    Integer                                  :: k, dbra, dket, m, n

    If (ibra == 0 .Or. iket == 0) Return

    k = DRT(iket)%k
    walk1(k+1)%row    => DRT(ibra)
    walk1(k+1)%vertex =  ibra
    walk2(k+1)%row    => DRT(iket)
    walk2(k+1)%vertex =  iket

    If (k == 0) Then
       m = ShavittFlags%IndVec(ybra+1)
       n = ShavittFlags%IndVec(yket+1)
       If (m/=0 .And. n/=0) Call PlotLoop(DRT, walk2, walk1, m, n)  !ifpgplot!
    Else
       Do dbra=0, 3
          walk1(k+1)%d = dbra
          Do dket=0, 3
             walk2(k+1)%d = dket
             Call ConstructWalks2(DRT, ShavittFlags, walk1, walk2,                  &
                                  DRT(ibra)%jd(dbra),      DRT(iket)%jd(dket),      &
                                  DRT(ibra)%yd(dbra)+ybra, DRT(iket)%yd(dket)+yket)
          End Do
       End Do
    End If

    Return
  End Subroutine ConstructWalks2



  ! Plot pair of walks (Michael Beck).
  Subroutine Plotloop(DRT,walk1,walk2,index,jndex)
    Type(PaldusRow), Dimension(:), Intent(in)    :: DRT
    Type(ShavittWalk), Dimension(:), Intent(in)  :: walk1, walk2
    Integer , Intent(in) :: index,jndex
    ! End of dummy parameters.

    Character(Len =10), Parameter    :: digitstring = '0123456789'
    Integer                          :: tailindex, i, j,k,kmax,kmin,norbs, help,d,dd 
    Integer                          :: old, help2
    Double Precision, Dimension(:), Allocatable  :: X, Y
    Double Precision                             :: h1,h2, fonthoehe, scale
    Character (Len=3)                :: b, bb
    Character (Len=6)                :: indexstring,jndexstring


    tailindex = Ubound(DRT,1)
    norbs = DRT(1)%k
    Allocate(X(tailindex))
    Allocate(Y(tailindex))

    Call bignumberToString(index, indexstring)
    Call bignumberToString(jndex, jndexstring)

    Call PGBEG(0,'loop.' // indexstring // '.' //           & !ifpgplot!
         &     jndexstring // '.ps/VPS',1,1)                  !ifpgplot!
    tailindex = Ubound(DRT,1)

    help = Step(3)*DRT(1)%a + Step(2)*DRT(1)%b
    h1 = help  + 2 * step(3)
    h2 = norbs + 2
    Call PGENV(0,Real(h1),0,Real(h2),0,-1)                    !ifpgplot!

    ! set roman font
    Call PGSCF(2)                                             !ifpgplot!

    ! set fontheight (PGPLOT default: 1/40 of pageheight)
    If(norbs > 8) Then
       scale = 1.d0/(4.d0*h2)
    Else
       scale = 1.d0/40.d0
    End If
    Call PGSCH(40.*Real(scale))                               !ifpgplot!
    fonthoehe = scale * h2
  
    Call PGText(0,Real(h2-1.5*fonthoehe),                   & !ifpgplot!
         &      "  matrixelement < " // indexstring //      & !ifpgplot!
         &      " | H | " // jndexstring // " >")             !ifpgplot!

    !construct coordinates
    X(1) = step(3)
    Y(1) = h2 - 1.d0
    kmin = 2
    Do i=1,tailindex 
       Do d=0,3
          help = DRT(i)%jd(d)
          If(help > 0) Then
             X(help) = X(i) + step(d) 
             Y(help) = Y(i) - 1.d0
          End If
       End Do
       If(i > 1 .And. DRT(i)%k < DRT(i - 1)%k) Then
          kmax = i - 1
          Do j=kmin,kmax
             Do k=kmin, j-1
                If(Abs(X(j)-X(k)) .Lt. Tiny) X(j) = X(j) + 1.d0
             End Do
          End Do
          kmin = i
       End If
    End Do

    !plot vertixes
    Call PGPT(tailindex,Real(X),Real(Y),9)                    !ifpgplot!

    !plot bra configuration with fat lines
    Call PGSLW(4)                                             !ifpgplot!
    old = walk1(1)%vertex
    Do i=1,norbs
       help = walk1(i+1)%vertex
       Call PGLINE(2,Real((/X(old),X(help)/)),              & !ifpgplot!
            &      Real((/Y(old),Y(help)/)))                  !ifpgplot!
       old = help
    End Do

    !plot ket configuration with very fat, broken lines
    old = walk1(1)%vertex
    Call PGSLW(4)                                             !ifpgplot!
    Call PGSLS(2)                                             !ifpgplot!
    Do i=1,norbs
       help = walk2(i+1)%vertex
       Call PGLINE(2,Real((/X(old),X(help)/)),              & !ifpgplot!
            &      Real((/Y(old),Y(help)/)))                  !ifpgplot!
       old = help
    End Do

    !plot ket configuration with very fat, broken lines
    old = walk1(1)%vertex
    Call PGSLW(1)                                             !ifpgplot!
    Call PGSLS(1)                                             !ifpgplot!
    Do i=1,norbs
       help  = walk1(i+1)%vertex
       help2 = walk2(i+1)%vertex
       d  = walk2(i+1)%d + 1
       dd = walk1(i+1)%d + 1
       Call numberToString(walk2(i+1)%row%b, b)
       Call numberToString(walk1(i+1)%row%b, bb)
       If(help == help2) Then
          Call PGPTXT(Real(h1),i+1.,0.0,1.0,"b' = b =" // b)  !ifpgplot!
       Else
          Call PGPTXT(Real(h1),i+1.,0.0,1.0,"b'b=" // b //  & !ifpgplot!
               &      " " // bb)                              !ifpgplot!
       End If
       Call PGPTXT(0.0,i+0.5,0.0,0.0,"d'd="//               & !ifpgplot!
            &      digitstring(d:d) // digitstring(dd:dd))    !ifpgplot!
       old = help2
    End Do
    Call PGEND                                                !ifpgplot!
    !-AK Added deallocation of X and Y:
    Deallocate(Y)
    Deallocate(X)
    Return
  End Subroutine Plotloop



  ! (Michael Beck)
  Subroutine numberToString(number, string)
    Implicit None 
    Integer, Intent(in) :: number 
    Character (Len=3) :: string
    ! end of dummies
    Integer             :: i, kopie
    Character(Len =10), Parameter :: digitstring = '0123456789'
    Integer, Dimension(3)    :: z
    kopie = number

    Do i = 2,0,-1
       z(i+1)= Int(kopie/(10**i))
       kopie  = kopie - z(i+1)*10**i
    End Do
    string = digitstring(z(3)+1:z(3)+1) // digitstring(z(2)+1:z(2)+1) // digitstring(z(1)+1:z(1)+1)
    
    Return
  End Subroutine numberToString



  ! (Michael Beck)
  Subroutine  bignumberToString(number, string)
    Implicit None 
    Integer, Intent(in) :: number 
    Character (Len=6) :: string
    ! end of dummies
    Integer             :: i, kopie
    Character(Len =10), Parameter :: digitstring = '0123456789'
    Integer, Dimension(6)    :: z
    kopie = number

    Do i = 5,0,-1
       z(i+1)= Int(kopie/(10**i))
       kopie  = kopie - z(i+1)*10**i
    End Do
    string = digitstring(z(6)+1:z(6)+1) // digitstring(z(5)+1:z(5)+1) // digitstring(z(4)+1:z(4)+1) //&
         &   digitstring(z(3)+1:z(3)+1) // digitstring(z(2)+1:z(2)+1) // digitstring(z(1)+1:z(1)+1)
    
    Return
  End Subroutine bignumberToString



  ! (Michael Beck)
  Function step(d)
    Integer step
    Integer d
    If(d==0) step = 0
    If(d==1) step = 1
    If(d==2) step = 4
    If(d==3) step = 5
    Return
  End Function step



  !
  ! Subroutines for printing a linear combination of Slater determinants
  ! corresponding to each Gelfand state.
  !

  ! Public driver routine.
  Subroutine PrintLinearCombinationsOfSDs(DRT, ShavittFlags, SymLabel)
    Type(PaldusRow),   Dimension(:), Pointer     :: DRT
    Type(ShavittControl)                         :: ShavittFlags
    Character(Len=4),  Dimension(:)              :: SymLabel
    ! End of dummy parameters.
    Type(ShavittWalk), Dimension(:), Allocatable :: walk
    Integer                                      :: ilex

    ilex = 0
    Allocate(walk(DRT(1)%k+1))
    write(Stdout,*)
    write(Stdout,fmt=('(" =====  List of Gelfand states represented by the DRT  =====")'))
    write(Stdout,*)
    Call ConstructWalks(DRT, ShavittFlags, SymLabel, walk, 1, ilex)
    Deallocate(walk)

    Return
  End Subroutine PrintLinearCombinationsOfSDs



  ! At each vertex of the Shavitt graph follow downwards all existing arcs.
  Recursive Subroutine ConstructWalks(DRT, ShavittFlags, SymLabel, walk, irow, ilex)
    Implicit None
    Type(PaldusRow),   Dimension(:), Pointer :: DRT
    Type(ShavittControl)                     :: ShavittFlags
    Character(Len=4),  Dimension(:)          :: SymLabel
    Type(ShavittWalk), Dimension(:)          :: walk
    Integer                                  :: irow, ilex
    ! End of dummy parameters.
    Integer                                  :: k, i, d

    If (irow == 0) Return

    k = DRT(irow)%k
    walk(k+1)%row    => DRT(irow)
    walk(k+1)%vertex =  irow

    If (k == 0) Then
       ilex  = ilex + 1
       i     = ShavittFlags%IndVec(ilex)
       If (i /= 0) Then
          Write(Stdout,'(2(A,I6))')         ' Gelfand state #: ', i, &
                                            ', lexical index:              ', ilex
          Write(Stdout,'(A,I6,A,A,A,I1,A)') ' ExcitationLevel: ',ShavittFlags%ExLev(ilex),                            &
                                            ', irreducible representation: ', SymLabel(ShavittFlags%WalkIrrep(ilex)), &
                                            '(', ShavittFlags%WalkIrrep(ilex), ')'
          Call ConstructLinearCombinationOfSDs(walk, ShavittFlags%Multiplicity,                             &
               &                               ShavittFlags%numberOfActiveElectrons, ShavittFlags%WhoIsWho)
          Write(Stdout,*) '---------------------------------------------------------------'
       End If
    Else
       Do d=0, 3
          walk(k+1)%d = d
          Call ConstructWalks(DRT, ShavittFlags, SymLabel, walk, DRT(irow)%jd(d), ilex)
       End Do
    End If

    Return
  End Subroutine ConstructWalks



  ! Print a linear combination of Slater determinants for a given walk (Michael Beck).
  ! THIS ROUTINE DOES NOTHING BUT PRINT THESE LINEAR COMBINATIONS OF DETERMINANTS.
  ! THE REMAINING PROGRAM DOES NOT USE SLATER DETERMINANTS AT ALL!
  !
  Subroutine ConstructLinearCombinationOfSDs(walk, Multiplicity, Nelec, whoiswho)
    ! This routine expands a given Gelfand state (i.e. a walk) into a linear
    ! combination of Slater determinants (SD). See eqs. (11), (12) in Shavitt81.
    Implicit None
    Type(ShavittWalk), Dimension(:), Intent(In)  :: walk
    Integer,                         Intent(In)  :: Multiplicity, Nelec
    Integer,           Dimension(:), Intent(In)  :: whoiswho ! labels of active MOs
    ! End of dummy parameters.
    Double Precision                             :: f, norm
    Type(MO),          Dimension(:), Allocatable :: SD
    Integer                                      :: i, j, k
    Integer                                      :: alphas, betas, doublecount, halfoccupied
    Character, Dimension(0:1), Parameter         :: spinstring = (/'b','a'/)

    Allocate(SD(Nelec))
    j = 0
    doublecount = 0
    Do i=2,Ubound(walk,1)
       If(walk(i)%d .Ne. 0) Then
          If(walk(i)%d == 3) Then
             ! this orbital is doubly occupied. We assign an alpha and 
             ! an beta electron once and for all
             doublecount = doublecount + 1
             j=j+1
             SD(j)%index = walk(i)%row%k
             SD(j)%a     = walk(i)%row%a
             SD(j)%b     = walk(i)%row%b
             SD(j)%d     = walk(i)%d
             SD(j)%spin  = 1                 ! alpha
             SD(j)%Open  = .False.
             j=j+1
             SD(j)%index = walk(i)%row%k
             SD(j)%a     = walk(i)%row%a
             SD(j)%b     = walk(i)%row%b
             SD(j)%d     = walk(i)%d
             SD(j)%spin  = 0                 ! beta
             SD(j)%Open  = .False.
          Else
             ! this is an open shell
             j=j+1
             SD(j)%index = walk(i)%row%k
             SD(j)%a     = walk(i)%row%a
             SD(j)%b     = walk(i)%row%b
             SD(j)%d     = walk(i)%d
             SD(j)%Open  = .True.
          End If
       End If
    End Do

    halfoccupied = Nelec-2*doublecount              ! number of open shells in the SD
    alphas = (Nelec+Multiplicity-1)/2 - doublecount ! number of unpaired alpha spins
    betas  = halfoccupied - alphas                  ! number of unpaired beta spins

    Write(Stdout,*)
    Write(Stdout,fmt=('(" decomposition into Slater determinants:")'))
    Write(Stdout,fmt=('("                 # of active MO (a = alpha spin, b = beta spin)")'))
    Write(Stdout,fmt=('("       coeff.  ",40i3)')) (whoiswho(SD(j)%index), j=1,Nelec)
    norm = 0.
    Do i=2**alphas-1, (2**alphas-1)*2**betas
       ! This is a little tricky code for distributing (alphas) alpha-spins 
       ! among (halfoccupied) open shell MOs. The idea is to use all binary 
       ! numbers with (halfoccupied) digits and a trace equal (alphas). 
       ! the lowest such number has the decimal value 2**alphas-1, the highest
       ! equals (2**alphas-1)*2**betas. Advantage of this procedure: the SD
       ! are constructed in something like a lexical order. Disadvantage:
       ! each binary number has to be checked for the correct trace. This 
       ! has been simplified in f90 by the btest intrinsic. I've done a little
       ! search in number theory literatur, but havn't found any analytical
       ! procedure for finding a gapless series of numbers with equal traces.
       j=0
       Do k=1,Nelec
          If(SD(k)%Open) Then
             If(Btest(i,k-1-j)) Then
                SD(k)%spin = 1
             Else
                SD(k)%spin = 0
             End If
          Else
             j=j+1
          End If
       End Do
       f = SDCoefficientInCSF(SD)
       If((Sum(SD%spin) .Eq. alphas+doublecount) .And. (Abs(f) .Ge. Tiny)) Then
          norm = norm + f**2
          Write(Stdout,fmt=('(F15.6,40a3)')) f,(spinstring(SD(j)%spin), j=1,Nelec)
       End If
    End Do
    Write(Stdout,fmt=('(" Norm:",F9.6)')) norm
    
    Deallocate(SD)
    Return
  End Subroutine ConstructLinearCombinationOfSDs



  ! (Michael Beck)
  Function SDCoefficientInCSF(SD)
    ! Auxiliary function for ConstructLinearCombinationOfSDs 
    Implicit None 
    Double Precision         :: SDCoefficientInCSF
    Type(MO), Dimension(:)   :: SD
    ! End of dummy parameters.
    Integer                  :: k
    Integer,  Dimension(2)   :: gamma ! gamma(1/0) gives the number of spins 0/1 preceding a orbital

    SDCoefficientInCSF = 1.D0
    gamma = 0
    Do k=1,Ubound(SD,1)
       If (SD(k)%spin == 1) Then
          gamma(1) = gamma(1) + 1
       Else
          gamma(2) = gamma(2) + 1
       End If
       Select Case(SD(k)%d)
          Case(1)
            SDCoefficientInCSF =  SDCoefficientInCSF &
                 &               * Sqrt(Dble((SD(k)%a+SD(k)%b-gamma(SD(k)%spin+1)))/Dble(SD(k)%b))
          Case(2)
            SDCoefficientInCSF =  SDCoefficientInCSF            &
                 &               * ((-1)**(SD(k)%b+SD(k)%spin)) &
                 &               * Sqrt(Dble((gamma(SD(k)%spin+1)-SD(k)%a+1))/Dble((SD(k)%b+2)))
          Case(3)
            If (SD(k)%spin == 1)  SDCoefficientInCSF = SDCoefficientInCSF * (-1)**SD(k)%b
       End Select
       If (Abs(SDCoefficientInCSF) .Lt. Tiny) Return
    End Do
    Return
  End Function SDCoefficientInCSF


End Module gugaplots
