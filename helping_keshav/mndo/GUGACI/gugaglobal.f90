!
! Definition of global data types and variables.
!
! First version by Michael E. Beck in 1998 at Zurich University.
!
! Final version by Axel Koslowski in 2000-2005 at MPI Muelheim.
!
! Shape-driven algorithm added by Axel Koslowski in 2012-2013.
!

!
! Some words of caution:
!
! GUGA-CI (actually: UGA) many-electron basis functions are not configuration
! state functions. The latter term refers to spin-adapted linear combinations
! of Slater determinants. GUGA-CI basis functions are Gelfand states. This
! term originates from the representation theory of lie groups. These lack
! the specification of the magnetic total spin quantum number which must be
! specified in addition to convert a GUGA-CI basis function into a linear
! combination of Slater determinants. The names of some variables in this
! module include the misleading acronym CSF although they refer to Gelfand
! states. These should be renamed accordingly. Moreover, variables concerning
! the reference Gelfand states (the Gelfand states that correspond to the
! reference configurations) should be named such that it is obvious whether
! the index is a lexical index (counting all Gelfand states represented by
! the Shavitt graph) or whether it refers to the Gelfand states being arranged
! by spatial symmetry.
!

Module gugaglobal

  Implicit None

  Public

  !
  ! Parameters used in the definition of data types.
  !

  ! Maximum number of generators per buffer during loop construction.
  Integer, Parameter :: genBufSize   = 100
  ! Parameters corresponding to partial loops.
  Integer, Parameter :: loopBufSize  = 10000
  Integer, Parameter :: numUpTypes   = 10
  Integer, Parameter :: numLoTypes   = 13
  ! Parameters corresponding to partial loop combinations.
  Integer, Parameter :: combiBufSize = 1000
  ! Upper partial loop type with one index (l).
  Integer, Parameter :: iUpperRt     =  1 ! Cut below an Rt segment.
  ! Upper partial loop types with two indices (k,l).
  Integer, Parameter :: iUpperW      =  2 ! Cut below a W segment.
  Integer, Parameter :: iUpperOneR   =  3 ! Cut below an RtRb one-body loop.
  Integer, Parameter :: iUpperRLdg   =  4 ! Cut below an RtLt segment, diagonal partial loop.
  Integer, Parameter :: iUpperRLco   =  5 ! Cut below an RtLt or RLt segment, walks converged at
                                          ! the border between upper and lower part of the DRT.
  Integer, Parameter :: iUpperRL     =  6 ! Cut below an RtLt or RLt segment, general shape.
  Integer, Parameter :: iUpperRtRt   =  7 ! Cut below an RtRt segment.
  Integer, Parameter :: iUpperRRt    =  8 ! Cut below an RRt segment.
  ! Upper partial loop types with three indices (j,k,l).
  Integer, Parameter :: iUpperR      =  9 ! Cut above an Rb segment.
  Integer, Parameter :: iUpperL      = 10 ! Cut above an Lb segment.
  ! Lower partial loop types with one index (i).
  Integer, Parameter :: iLowerRb     =  1 ! Cut above an Rb segment.
  Integer, Parameter :: iLowerLb     =  2 ! Cut above an Lb segment.
  ! Lower partial loop types with two indices (i,j).
  Integer, Parameter :: iLowerW      =  3 ! Cut above a W segment.
  Integer, Parameter :: iLowerOneR   =  4 ! Cut above an RtRb one-body loop.
  Integer, Parameter :: iLowerOneL   =  5 ! Cut above an LtLb one-body loop.
  Integer, Parameter :: iLowerRLdg   =  6 ! Cut above an RbLb segment, diagonal partial loop.
  Integer, Parameter :: iLowerRLcp   =  7 ! Cut above an RbLb, RLb, or RbL segment, walks converged
                                          ! at the border between upper and lower part of the DRT,
                                          ! lexbra < lexket.
  Integer, Parameter :: iLowerRLcm   =  8 ! Cut above an RbLb, RLb, or RbL segment, walks converged
                                          ! at the border between upper and lower part of the DRT,
                                          ! lexbra > lexket.
  Integer, Parameter :: iLowerRL     =  9 ! Cut above an RbLb, RLb, or RbL segment, general shape.
  Integer, Parameter :: iLowerRbRb   = 10 ! Cut above an RbRb segment.
  Integer, Parameter :: iLowerRbR    = 11 ! Cut above an RbR segment.
  Integer, Parameter :: iLowerRRb    = 12 ! Cut above an RRb segment.
  ! Lower partial loop type with three indices (i,j,k).
  Integer, Parameter :: iLowerR      = 13 ! Cut below an Rt segment.

  ! Combinations of partial loop types associated with the upper
  ! triangle of the CI Hamitonian and corresponding integrals.
  ! Active orbital indices i, j, k, and l associated with the
  ! critical segments are rearranged accordingly during
  ! construction of the partial loops.
  !
  ! Index l is always associated with an upper, index i with
  ! a lower partial loop. Moreover it is guaranteed that index l
  ! corresponds to the highest and index i to the lowest orbital
  ! level (in the DRT) of the four critical orbitals. Indices
  ! j and k may be associated with an upper or lower partial
  ! loop depending on the partial loop type.
  !
  ! iUpperRt   + iLowerRb   => <i|O|l>
  ! iUpperRt   + iLowerR    => (ij|kl)
  ! iUpperW    + iLowerW    => (ii|ll)
  ! iUpperW    + iLowerOneR => (ij|ll)
  ! iUpperOneR + iLowerW    => (ii|kl)
  ! iUpperOneR + iLowerOneR => (ij|kl)
  ! iUpperOneR + iLowerOneL => (ij|kl)
  ! iUpperRLdg + iLowerRLdg => (il|il)
  ! iUpperRLdg + iLowerRLcp => (il|jl)
  ! iUpperRLco + iLowerRLdg => (ik|il)
  ! iUpperRLco + iLowerRLcp => (ik|jl)
  ! iUpperRLco + iLowerRLcm => (ik|jl)
  ! iUpperRL   + iLowerRL   => (ik|jl)
  ! iUpperRtRt + iLowerRbRb => (il|il)*0.5
  ! iUpperRtRt + iLowerRbR  => (il|jl)
  ! iUpperRRt  + iLowerRbRb => (ik|il)
  ! iUpperRRt  + iLowerRbR  => (ik|jl)
  ! iUpperRRt  + iLowerRRb  => (ik|jl)
  ! iUpperR    + iLowerRb   => (ij|kl)
  ! iUpperL    + iLowerLb   => (ij|kl)

  !
  ! Definition of data types.
  !

  Type ShavittControl
     Integer  :: numberOfActiveElectrons
     Integer  :: Charge
     Integer  :: Multiplicity
     Integer  :: ExcitationLevel   ! default: 2 (that is SDCI)
     Integer  :: HOMOindex
     Integer  :: NumberOfCIRoots   ! Requested number of CI roots of unspecified symmetry.
     Integer  :: LRoot             ! IN2(140).
     Integer  :: iState, jState    ! Indices of CI roots of interest.
     Integer  :: WhichDiagonalizer
     Integer  :: MinDav            ! Minimum number of valid columns in B matrix
     Integer  :: MaxDav            ! Maximum number of valid columns in B matrix
     Integer  :: KitDav            ! Maximum number of Davidson iterations
     Double Precision :: QNorm     ! Convergence criterion
     Integer  :: iop               ! see MNDO99.f, IN2(2)
                                   !
     Integer  :: Algorithm         ! Algorithm < 0:  Loop-driven algorithm.
                                   ! Algorithm > 0:  Shape-driven algorithm.
                                   !
                                   ! IABS(Algorithm)
                                   !  = 1  All generator matrix elements between states of
                                   !       the same symmetry are kept in core all the time.
                                   !  = 2  Only the generator matrix elements connecting
                                   !       states of the current symmetry are kept in core.
                                   !  = 3  Direct algorithm.
                                   !
     Integer  :: PrintLevel        !  -9 = no print at all
                                   !   0 = default print
                                   ! 1-6 = gradually growing debug print
                                   !
     Integer  :: Plot              !   0 = no plots
                                   !   1 = plot Shavitt graph (default)
                                   !   2 = like 1 plus Gelfand states
                                   !   3 = like 2 plus loops
                                   !
     Integer  :: NVCapa            ! Compute capability of the first NVIDIA GPU.
                                   ! The value of NVCapa is 10 * Major + Minor or 0.
                                   !
     Type(PaldusRow),    Dimension(:),   Pointer  :: DRT                   ! Distinct row table.
     Integer,            Dimension(:),   Pointer  :: IndVec                ! Projection of the lexical indices on the CSF indices.
                                                                           ! CSF index (or 0) for each lexical index.
     Integer,            Dimension(:),   Pointer  :: SymVec                ! Symmetry of each CSF (see also WalkIrrep below).
     Type(Generator),    Dimension(:),   Pointer  :: OneGen,    TwoGen     ! Buffers for one- and two-electron
                                                                           ! generator matrix elements.
     Integer,            Dimension(:),   Pointer  :: OneIndex,  TwoIndex   ! OneIndex(i) and TwoIndex(i) are the indices
                                                                           ! of the LAST generator matrix element in row i
                                                                           ! of the upper triangle of the CI Hamiltonian
                                                                           ! (unlike in the CSR sparse matrix format).
                                                                           ! OneIndex(0) and TwoIndex(0) are 0 unless the
                                                                           ! corresponding CI Hamiltonian represents a block
                                                                           ! of a larger matrix.
     Type(GenBufList),   Dimension(:),   Pointer  :: OneList,   TwoList    ! Linked lists of buffers with generator matrix
                                                                           ! elements for each diagonal of the CI Hamiltonian,
                                                                           ! for use with the shape-driven algorithm.
     Type(PartialLoop),  Dimension(:),   Pointer  :: UpperPart, LowerPart  ! Upper and lower partial loops.
     Type(CompleteLoop), Dimension(:),   Pointer  :: OneLoop,   TwoLoop    ! Complete diagonal and off-diagonal
     Integer,            Dimension(numUpTypes)    :: iFirstUp,  iLastUp    ! Index of first and last upper partial loop
                                                                           ! with the corresponding type.
     Integer,            Dimension(numLoTypes)    :: iFirstLo,  iLastLo    ! Index of first and last lower partial loop
                                                                           ! with the corresponding type.
     Integer,            Dimension(8)             :: iFirst1,   iFirst2    ! Index of first one-body and two-body loop
                                                                           ! with the corresponding symmetry.
     Integer,            Dimension(8)             :: iDiag1,    iDiag2     ! Index of last diagonal one-body and two-body loop
                                                                           ! with the corresponding symmetry.
     Integer,            Dimension(8)             :: iOff1,     iOff2      ! Index of first off-diagonal one-body and two-body loop
                                                                           ! with the corresponding symmetry.
     Integer,            Dimension(8)             :: iLast1,    iLast2     ! Index of last one-body and two-body loop
                                                                           ! with the corresponding symmetry.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListOne          ! One-body combinations of upper and lower partial loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym13ijkl    ! Two-body combinations of upper and lower partial loops
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym22ijkl    ! considering pairs of equal iProd(sybra, syket) only,
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym22ikjl    ! for CI Hamiltonian and gradient.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym22iill    ! Corresponding to W-W loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym22ilil    ! Corresponding to RtLt-RbLb loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym22half    ! Corresponding to RtRt-RbRb loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListSym31ijkl
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll13ijkl    ! Two-body combinations of upper and lower partial loops
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll22ijkl    ! of any symmetry, for non-adiabatic coupling elements.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll22ikjl    !
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll22iill    ! Corresponding to W-W loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll22ilil    ! Corresponding to RtLt-RbLb loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll22half    ! Corresponding to RtRt-RbRb loops.
     Type(Combination),  Dimension(:),   Pointer  :: CombiListAll31ijkl
     Integer                                      :: NumOne                ! Numbers of such combinations.
     Integer                                      :: NumSym13ijkl
     Integer                                      :: NumSym22ijkl
     Integer                                      :: NumSym22ikjl
     Integer                                      :: NumSym22iill
     Integer                                      :: NumSym22ilil
     Integer                                      :: NumSym22half
     Integer                                      :: NumSym31ijkl
     Integer                                      :: NumAll13ijkl
     Integer                                      :: NumAll22ijkl
     Integer                                      :: NumAll22ikjl
     Integer                                      :: NumAll22iill
     Integer                                      :: NumAll22ilil
     Integer                                      :: NumAll22half
     Integer                                      :: NumAll31ijkl
     Integer,            Dimension(:),   Pointer  :: lexUpper              ! Lexical indices of upper walks from vertices.
     Double Precision,   Dimension(:),   Pointer  :: CIE, CIESAV           ! CI eigenvalues and backup copy.
     Double Precision,   Dimension(:,:), Pointer  :: CIC, CICSAV           ! CI eigenvectors and backup copy.

     Integer                        :: NumSym      ! Number of irreducible representations.
     Integer, Dimension(8)          :: NRoots      ! Requested number of CI roots of each irrep.
     Integer, Dimension(8)          :: FirstRoot   ! Index of first root of each symmetry.
     Integer, Dimension(8)          :: LastRoot    ! Index of last root of each symmetry.
                                                   ! (NRoots, FirstRoot, LastRoot only relevant
                                                   !  if NumberOfCIRoots == 0)
     Integer                        :: TotalRoots  ! Total number of requested CI roots.
     Integer, Dimension(8)          :: NumCSF      ! Number of CSFs of each symmetry in CI Hamiltonian.
     Integer, Dimension(8)          :: FirstCSF    ! Index of first CSF of each symmetry in CI Hamiltonian.
     Integer, Dimension(8)          :: LastCSF     ! Index of last CSF of each symmetry in CI Hamiltonian.
     Integer                        :: TotalCSF    ! Total number of relevant CSFs.
     Integer, Dimension(:), Pointer :: WalkIndex   ! Lexical index corresponding to each CSF.
     Integer, Dimension(:), Pointer :: RootSym     ! Symmetry of each CI root.
     Integer, Dimension(:), Pointer :: WhoIsWho    ! Orbital indices of active orbitals.
     Integer, Dimension(:), Pointer :: MOIrreps    ! Symmetries of active orbitals.
     Integer, Dimension(:), Pointer :: WalkIrrep   ! Holds the irreducible representation of each walk.
     Integer, Dimension(:), Pointer :: ExLev       ! Excitation level of each walk with respect to
                                                   ! the nearest reference configuration.
     Double Precision               :: Leading     ! Threshold for the printing of leading configurations.
     Integer, Dimension(:), Pointer :: iFirstRef   ! Index of first Gelfand state for each reference configuration.
     Integer, Dimension(:), Pointer :: iLastRef    ! Index of last  Gelfand state for each reference configuration.
     Integer, Dimension(:), Pointer :: iRefWalk    ! List of walks corresponding to the reference configurations.
     Integer, Dimension(:), Pointer :: iRefConf    ! List of Gelfand states corresponding to the reference configurations.
     Integer                        :: mciref      ! How the references were specified in the input.
     Integer                        :: nciref      ! Number of reference configurations.
     Integer                        :: nRefCSF     ! Number of corresponding Gelfand states.
     Integer                        :: kBorder     ! Orbital level that separates upper and lower part of the DRT.
  End Type ShavittControl



  Type PaldusRow
     Integer                  ::   k   ! Orbital level.
     Integer                  ::   a   ! Number of double occupancies.
     Integer                  ::   b   ! Number of single occupancies.
     !                             c   ! Number of empty MOs, c=k-a-b.
     Integer, Dimension(0:3)  ::   jd  ! Downward chaining indices.
     Integer, Dimension(0:3)  ::   ju  ! Upward chaining indices.
     Integer, Dimension(0:4)  ::   yd  ! Weights of arcs going down.
     Integer, Dimension(0:3)  ::   yu  ! Weights of arcs coming from above.
     Integer                  ::   xu  ! Upward vertex weight.
     Integer                  ::   iu  ! Index of first upper walk.
  End Type PaldusRow



  Type StepVector
     Integer, Dimension(:), Pointer :: d
  End Type StepVector



  Type Generator
     Sequence
     Double Precision :: value  ! value of the generator matrix element
     Integer(4)       :: icol   ! column index in the CI Hamiltonian
     Integer(2)       :: p, q   ! MO indices (canonical indices for two-electron generators)
  End Type Generator



  Type GenBuffer
     Type(Generator), Dimension(genBufSize) :: gen    ! buffer contents
     Type(GenBuffer), Pointer               :: next   ! pointer to next buffer
     Integer                                :: count  ! number of generators in this buffer
  End Type GenBuffer



  Type GenBufList
     Type(GenBuffer), Pointer               :: first  ! pointer to first buffer
     Integer                                :: count  ! number of buffers in list
  End Type GenBufList



  Type PartialLoop
     Double Precision :: singlet     ! Partial singlet coupling value.
     Double Precision :: triplet     ! Partial triplet coupling value.
     Integer          :: itype       ! Code for the upper partial loop type.
     Integer          :: itip        ! Index of the row at which bra and ket walk coincide.
     Integer          :: ibra        ! Row index of the bra walk at the cut.
     Integer          :: iket        ! Row index of the ket walk at the cut.
     Integer          :: lexbra      ! Partial lexical index of the partial bra walk connecting itop and ibra.
     Integer          :: lexket      ! Partial lexical index of the partial ket walk connecting itop and iket.
     Integer          :: sybra       ! Symmetry of the partial bra walk connecting itop and ibra.
     Integer          :: syket       ! Symmetry of the partial ket walk connecting itop and iket.
     Integer          :: i, j, k, l  ! For upper partial loops, j, k, and l are the relevant indices and
                                     ! i may contain a compound index, depending on the loop type.
                                     ! For lower partial loops, i, j, and k are the relevant indices and
                                     ! l may contain a compound index, depending on the loop type.
  End Type PartialLoop



  Type CompleteLoop
     Double Precision :: value       ! Loop value.
     Integer          :: itop        ! Index of the top row.
     Integer          :: ibot        ! Index of the bottom row.
     Integer          :: lexbra      ! Partial lexical index of the partial bra walk connecting itop and ibot.
     Integer          :: lexket      ! Partial lexical index of the partial ket walk connecting itop and ibot.
     Integer          :: sybra       ! Symmetry of the partial bra walk connecting itop and ibot.
     Integer          :: syket       ! Symmetry of the partial ket walk connecting itop and ibot.
     Integer          :: p, q        ! Indices corresponding to the one-electron integral <i|O|j> or
                                     ! compound indices for the two-electron integral (ij|kl).
  End Type CompleteLoop



  Type Combination
     Integer          :: iUpB, iUpE  ! Index of first and last upper partial loop.
     Integer          :: iLoB, iLoE  ! Index of first and last lower partial loop.
  End Type Combination



  Type PartialLoopBuffer
     Type(PartialLoop),        Dimension(loopBufSize)  :: part    ! Memory for the partial loops in this buffer.
     Type(PartialLoopBuffer),  Pointer                 :: next    ! Pointer to the next buffer in this list.
     Integer                                           :: ncount  ! Number of partial loops in this buffer.
  End Type PartialLoopBuffer



  Type CompleteLoopBuffer
     Type(CompleteLoop),       Dimension(loopBufSize)  :: loop    ! Memory for the loops in this buffer.
     Type(CompleteLoopBuffer), Pointer                 :: next    ! Pointer to the next buffer in this list.
     Integer                                           :: ncount  ! Number of loops in this buffer.
  End Type CompleteLoopBuffer



  Type CombinationBuffer
     Type(Combination),        Dimension(combiBufSize) :: combi   ! Memory for the partial loop combinations in this buffer.
     Type(CombinationBuffer),  Pointer                 :: next    ! Pointer to the next buffer in this list.
     Integer                                           :: ncount  ! Number of elements in this buffer.
  End Type CombinationBuffer



  Type LoopBufLists
     Type(PartialLoopBuffer),  Pointer :: firstU   ! First buffer of linked list with upper partial loops.
     Type(PartialLoopBuffer),  Pointer :: firstL   ! First buffer of linked list with lower partial loops.
     Type(CompleteLoopBuffer), Pointer :: first1   ! First buffer of linked list with complete one-body loops.
     Type(CompleteLoopBuffer), Pointer :: first2   ! First buffer of linked list with complete two-body loops.
     Integer                           :: ncountU  ! Total number of upper partial loops.
     Integer                           :: ncountL  ! Total number of lower partial loops.
     Integer                           :: ncount1  ! Total number of complete one-body loops.
     Integer                           :: ncount2  ! Total number of complete two-body loops.
  End Type LoopBufLists



  !
  ! Global variables (not included in ShavittControl for more efficient access).
  !

  ! Unit number of standard output.
  Integer, Save :: Stdout

  ! Cayley table valid for all point groups in this program.
  ! We access a common block of the core program to get the values.
  Integer       :: iProd(8,8)
  Common /MulTab/ iProd

End Module gugaglobal
