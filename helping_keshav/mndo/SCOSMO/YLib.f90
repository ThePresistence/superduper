MODULE DataTypes
!!!****h* YorkLib/DataTypes
!!!
!!! NAME 
!!!     DataType -- Module containing general purpose derived datatypes and 
!!!     associated KIND PARAMETERs.
!!! AUTHOR
!!!     Darrin M. York
!!!     Department of Chemistry, Harvard University
!!! DATE
!!!     1997
!!! DESCRIPTION
!!!     The library is built such that 3 real data types are defined: MP, 
!!!     SP and DP. "Machine precision" (MP) is the default precision  
!!!     on the architecture being used. MP should only be used when it is
!!!     known that calculations will not need more than 6 to 7 decimal digits
!!!     of precision, or when backward compatability is needed for reading 
!!!     certain binary data files. "double precision" (DP) is used only 
!!!     in programming units where extended precision, if available, is 
!!!     specially desired (rare). All other cases should use the "single 
!!!     precision" KIND parameter. 
!!!
!!!     NOTE: With this precision convention, there can be no overloaded
!!!     programming units that are differentiated by the SP and DP derived
!!!     types, since these can in fact be the same; e.g. you CANNOT have
!!!     a situation like:
!!!
!!!     INTERFACE swap
!!!        MODULE PROCEDURE swap_SP, swap_DP
!!!     END INTERFACE
!!!     ...
!!!     SUBROUTINE swap_SP(asp, bsp)
!!!     REAL(SP), INTENT(INOUT) :: asp, bsp
!!!     ...
!!!     SUBROUTINE swap_DP(adp, bdp)
!!!     REAL(DP), INTENT(INOUT) :: adp, bdp
!!!     ...
!!!     So ONLY use REAL(DP) in subroutine and function dummy arguments
!!!     that are NOT overloaded with REAL(SP) and that ESPECIALLY NEED 
!!!     the increased precision (this is rare).
!!!     
!!!***

  PUBLIC

  INTEGER, PARAMETER :: I8B = SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: MP = KIND(1.0E0)
  INTEGER, PARAMETER :: SP = KIND(1.0D0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0_SP,1.0_SP))
  INTEGER, PARAMETER :: DPC = KIND((1.0_DP,1.0_DP))
  INTEGER, PARAMETER :: LGD = KIND(.true.)


  INTEGER(I4B),PARAMETER  :: LEN_FILENAME = 256
  INTEGER(I4B),PARAMETER  :: LEN_TOKEN = 20
  INTEGER(I4B),PARAMETER  :: LEN_LINE  = 500
  CHARACTER(LEN=6),PARAMETER :: FMT_LINE = '(A500)'

!!!****s* DataTypes/TokValT
!!!
!!! NAME
!!!     TokValT
!!! SOURCE
!!!
  TYPE TokValT
     CHARACTER(LEN=LEN_TOKEN)   :: Tok          =  ""
     LOGICAL(LGD)               :: MultiTok     = .FALSE.
     INTEGER(I4B)               :: LEN_TRUNC    =   4
     INTEGER(I4B)               :: NVal         =   0
     CHARACTER(LEN=1)           :: Type         =  "C" ! C,I,R,L => Character, Integer, Real, Logical

     CHARACTER(LEN=LEN_LINE),POINTER :: Line(:) => NULL()
     INTEGER(I4B),POINTER       :: Col(:)       => NULL()
     INTEGER(I4B),POINTER       :: LineN(:)     => NULL()

     LOGICAL(LGD),POINTER       :: IsInt(:)     => NULL()
     LOGICAL(LGD),POINTER       :: IsReal(:)    => NULL()
     LOGICAL(LGD),POINTER       :: IsLogic(:)   => NULL()

     LOGICAL(LGD),POINTER       :: LogicVal(:)  => NULL()
     REAL(SP),POINTER           :: RealVal(:)   => NULL()
     INTEGER(I4B),POINTER       :: IntVal(:)    => NULL()
     CHARACTER(LEN=200),POINTER :: CharVal(:)   => NULL()
  END TYPE TokValT
!!!***


!!!****s* DataTypes/InpFileT
!!!
!!! NAME
!!!     InpFileT
!!! SOURCE
!!!
  TYPE InpFileT
     CHARACTER(LEN=LEN_LINE),POINTER :: Line(:)  => NULL()
     INTEGER(I4B),POINTER            :: LineN(:) => NULL()
     LOGICAL(LGD)                    :: ERROR = .FALSE.
  END TYPE InpFileT
!!!***





!!!****s* DataTypes/Sparse_Matrix_Data
!!!
!!! NAME
!!!     sprs2_SP and sprs2_DP
!!! SOURCE
!!!
  TYPE sprs2_SP
     INTEGER(I4B) :: n = 0
     INTEGER(I4B) :: len = 0
     REAL(SP), DIMENSION(:), POINTER :: val      => NULL()
     INTEGER(I4B), DIMENSION(:), POINTER :: irow => NULL()
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol => NULL()
  END TYPE sprs2_SP

  TYPE sprs2_DP
     INTEGER(I4B) :: n = 0
     INTEGER(I4B) :: len = 0
     REAL(DP), DIMENSION(:), POINTER :: val      => NULL()
     INTEGER(I4B), DIMENSION(:), POINTER :: irow => NULL()
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol => NULL()
  END TYPE sprs2_DP
!!!***

!!!****s* DataTypes/Vector_Data
!!!
!!! NAME
!!!     Vec3 
!!! SOURCE
!!!
  TYPE Vec3
     REAL(SP) :: v(3) = 0.0_SP
  END TYPE Vec3
!!!***



!!!****s* DataTypes/Unit_Cell_Data
!!!
!!! NAME
!!!     UnitCell
!!! SOURCE
!!!
  TYPE UnitCell
     REAL(SP), DIMENSION(3,3) :: RealSpaceLvec  = 0.0_SP
     REAL(SP), DIMENSION(3,3) :: RecipSpaceLvec = 0.0_SP
     REAL(SP) :: abc(3)                         = 0.0_SP
     REAL(SP) :: AlpBetGam(3)                   = 0.0_SP
     REAL(SP) :: Volume                         = 0.0_SP
     INTEGER(I4B) :: Nop
     REAL(SP), POINTER :: PermOp(:,:,:) =>NULL()
     REAL(SP), POINTER :: TransOp(:,:) =>NULL()
     CHARACTER(LEN=20) :: SpaceGroup = 'UNKNOWN             '
     CHARACTER(LEN=20) :: Type       = 'UNKNOWN             '
  END TYPE UnitCell
!!!***

!!!****s* DataTypes/ZmatrixT
!!!
!!! NAME
!!!     ZmatrixT
!!! SOURCE
!!!
  TYPE ZmatrixT
    REAL(SP), DIMENSION(3) :: IntCrd
    INTEGER(I4B), DIMENSION(3) :: Index
  END TYPE ZmatrixT
!!!***

!!!****s* DataTypes/IO_Data
!!!
!!! NAME
!!!     IO
!!! SOURCE
!!!
  TYPE IO
     INTEGER(I4B) :: unit = 0
     INTEGER(I4B) :: recl = 0
     INTEGER(I4B) :: iostat = 0
     INTEGER(I4B) :: LineNumber = 0
     CHARACTER(LEN=LEN_FILENAME) :: file
     CHARACTER(LEN=20) :: format
     CHARACTER(LEN=20) :: status
     CHARACTER(LEN=20) :: access
     CHARACTER(LEN=20) :: position
     CHARACTER(LEN=20) :: action
     TYPE(IO), POINTER :: nextIO     => NULL()
     TYPE(IO), POINTER :: previousIO => NULL()
  END TYPE IO
!!!***

!!!****s* DataTypes/Pointer_Arrays
!!!
!!! NAME
!!!     I4Bptr, I4B1Dptr, I4B2Dptr, I4B3Dptr, SPptr, SP1Dptr, SP2Dptr, SP3Dptr
!!! SOURCE
!!!
  TYPE I4Bptr
     INTEGER(I4B), POINTER :: I4Bptr => NULL()
  END TYPE I4Bptr

  TYPE I4B1Dptr
     INTEGER(I4B), POINTER :: I4B1Dptr(:) => NULL()
  END TYPE I4B1Dptr

  TYPE I4B2Dptr
     INTEGER(I4B), POINTER :: I4B2Dptr(:,:) => NULL()
  END TYPE I4B2Dptr

  TYPE I4B3Dptr
     INTEGER(I4B), POINTER :: I4B3Dptr(:,:,:) => NULL()
  END TYPE I4B3Dptr

  TYPE SPptr
     REAL(SP), POINTER :: SPptr => NULL()
  END TYPE SPptr

  TYPE SP1Dptr
     REAL(SP), POINTER :: SP1Dptr(:) => NULL()
  END TYPE SP1Dptr

  TYPE SP2Dptr
     REAL(SP), POINTER :: SP2Dptr(:,:) => NULL()
  END TYPE SP2Dptr

  TYPE SP3Dptr
     REAL(SP), POINTER :: SP3Dptr(:,:,:) => NULL()
  END TYPE SP3Dptr
!!!***


!!!***
END MODULE DataTypes

MODULE ConstantsMod
!!!****h* YorkLib/ConstantsMod
!!!
!!! NAME
!!!     ConstantsMod -- Module which includes frequently used 
!!!                     PARAMETER constants.
!!! COPYRIGHT
!!!     Prof. Darrin M. York
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Darrin M. York
!!! CREATION DATE
!!!     1997
!!! DESCRIPTION
!!!     Contains fundamental physical constants and unit conversions to
!!!     very high accuracy with built-in dependencies.  See, for example:
!!!     http://physics.nist.gov/PhysRefData/contents.html
!!!
  USE DataTypes

  PUBLIC

  TYPE, PRIVATE :: ConstantT
     CHARACTER(LEN=25) :: Name
     CHARACTER(LEN=25) :: Unit
     CHARACTER(LEN=25) :: Abbrev
     REAL(SP) :: Value
  END TYPE ConstantT

  ! Numbers
  ! THERE CAN BE NO THING CALLED TINY!  TINY IS A FUNCTION
  REAL(SP), PARAMETER :: TINY_float=1.0E-10_SP
  REAL(SP), PARAMETER :: CCELEC=332.0716_SP !for CHARMM comparisons

  ! Fundamental mathematical constants (no units)
  !
  ! Single precision (for almost all purposes)
  !
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_SP
  REAL(SP), PARAMETER :: PI_D_180=(PI/180.0_SP)
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_SP
  REAL(SP), PARAMETER :: TWO_PI=2*PI
  REAL(SP), PARAMETER :: FOUR_PI=4*PI
  REAL(SP), PARAMETER :: SQRT_PI=1.772453850905516027298167483341_SP
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_SP
  REAL(SP), PARAMETER :: SQRT3=1.73205080756887729352744634150587236694281_SP
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_SP
  REAL(SP), PARAMETER :: RAD2DEG=(180.0_SP/PI) 
  REAL(SP), PARAMETER :: DEG2RAD=(PI/180.0_SP)
  !
  ! Double precision (for special applications requiring extended precision)
  !
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_DP
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_DP
  REAL(DP), PARAMETER :: TWO_PI_D=6.283185307179586476925286766559005768394_DP
  REAL(DP), PARAMETER :: FOUR_PI_D=4*PI_D
  REAL(DP), PARAMETER :: SQRT_PI_D=1.772453850905516027298167483341_DP
  REAL(DP), PARAMETER :: SQRT2_D=1.41421356237309504880168872420969807856967_DP
  REAL(DP), PARAMETER :: SQRT3_D=1.73205080756887729352744634150587236694281_DP
  REAL(DP), PARAMETER :: EULER_D=0.5772156649015328606065120900824024310422_DP
  REAL(DP), PARAMETER :: RAD2DEG_D=(180.0_DP/PI_D) 
  REAL(DP), PARAMETER :: DEG2RAD_D=(PI_D/180.0_DP)

  ! Fundamental physical constants (SI units) CODATA 1998 values
  REAL(SP), PARAMETER :: AVOGADRO_CONSTANT  = 6.02214199E+23_sp
  REAL(SP), PARAMETER :: SPEED_OF_LIGHT     = 299792458.0_sp
  REAL(SP), PARAMETER :: BOLTZMANN_CONSTANT = 1.3806503E-23_SP
  REAL(SP), PARAMETER :: PLANCK_CONSTANT = 6.62606876E-34_SP
  REAL(SP), PARAMETER :: HBAR = PLANCK_CONSTANT/(2*PI)
  REAL(SP), PARAMETER :: ELECTRON_CHARGE = 1.602176462E-19_sp
  REAL(SP), PARAMETER :: ELECTRON_REST_MASS = 9.10938188E-31
  REAL(SP), PARAMETER :: PROTON_REST_MASS = 1.67262158E-27_sp
  REAL(SP), PARAMETER :: NEUTRON_REST_MASS = 1.67492716E-27_sp
  REAL(SP), PARAMETER :: PERMITTIVITY_FREE_SPACE = 8.854187817E-12
  REAL(SP), PARAMETER :: FARADAY = 96485.3415_SP

  !
  ! Atomic Units in terms of SI
  !
  ! Independent quantities
  !
  REAL(SP), PARAMETER :: AU_MASS   = ELECTRON_REST_MASS
  REAL(SP), PARAMETER :: AU_CHARGE = ELECTRON_CHARGE
  REAL(SP), PARAMETER :: AU_ANGULAR_MOMENTUM = HBAR
  REAL(SP), PARAMETER :: AU_PERMITTIVITY = 4*PI*PERMITTIVITY_FREE_SPACE
  !
  ! Dependent quantities
  !
  REAL(SP), PARAMETER :: AU_DISTANCE = &
       & (AU_PERMITTIVITY*AU_ANGULAR_MOMENTUM**2)/(AU_MASS*AU_CHARGE**2)
  REAL(SP), PARAMETER :: AU_ENERGY = &
       & (AU_CHARGE**2)/(AU_PERMITTIVITY*AU_DISTANCE)
  REAL(SP), PARAMETER :: AU_VELOCITY = &
       & (AU_CHARGE**2)/(AU_PERMITTIVITY*AU_ANGULAR_MOMENTUM)
  REAL(SP), PARAMETER :: AU_TIME = &
       & (AU_DISTANCE/AU_VELOCITY)
  REAL(SP), PARAMETER :: AU_ELECTROSTATIC_POTENTIAL = &
       & (AU_CHARGE)/(AU_PERMITTIVITY*AU_DISTANCE)
  REAL(SP), PARAMETER :: AU_FORCE = &
       & (AU_CHARGE**2)/(AU_PERMITTIVITY*AU_DISTANCE**2)
  REAL(SP), PARAMETER :: AU_MAGNETIC_DIPOLE = &
       & (AU_CHARGE*AU_ANGULAR_MOMENTUM)/AU_MASS

  ! Fundamental physical constants (AU)
  REAL(SP), PARAMETER :: SPEED_OF_LIGHT_AU     = &
       SPEED_OF_LIGHT / AU_VELOCITY
  REAL(SP), PARAMETER :: BOLTZMANN_CONSTANT_AU = &
       BOLTZMANN_CONSTANT/AU_ENERGY
  REAL(SP), PARAMETER :: PLANCK_CONSTANT_AU = &
       PLANCK_CONSTANT/(AU_TIME*AU_ENERGY)

  REAL(SP), PARAMETER :: PROTON_REST_MASS_AU = &
       PROTON_REST_MASS / AU_MASS
  REAL(SP), PARAMETER :: NEUTRON_REST_MASS_AU = &
       NEUTRON_REST_MASS / AU_MASS

  ! Units (in terms of AU)
  REAL(SP), PARAMETER :: ELECTRON_VOLT = AU_CHARGE/AU_ENERGY
  REAL(SP), PARAMETER :: JOULE_PER_MOL = &
       & 1.0_SP/(AU_ENERGY*AVOGADRO_CONSTANT)
  REAL(SP), PARAMETER :: KJOULE_PER_MOL = &
       & 1000.0_SP * JOULE_PER_MOL
  REAL(SP), PARAMETER :: KCAL_PER_MOL = 4.184_SP * KJOULE_PER_MOL
  REAL(SP), PARAMETER :: KCAL_PER_ANG_MOL = 4.184E13_SP / (AU_FORCE * AVOGADRO_CONSTANT)
  REAL(SP), PARAMETER :: INVERSE_CM = &
       & (100*SPEED_OF_LIGHT*PLANCK_CONSTANT)/AU_ENERGY
  REAL(SP), PARAMETER :: ANGSTROM = 1.0E-10_SP/AU_DISTANCE
  REAL(SP), PARAMETER :: DEBYE = 3.335641E-30_SP/(AU_CHARGE*AU_DISTANCE)
  REAL(SP), PARAMETER :: GRAMS_PER_MOL = 1.00_SP / &
       & (1000.0_SP * AVOGADRO_CONSTANT  * AU_MASS)
  REAL(SP), PARAMETER :: CM3_PER_MOL =  AVOGADRO_CONSTANT * &
       & (100.0_SP * AU_DISTANCE)**3


  ! 1 au of distance is this many angstroms ~ 0.529177249
  REAL(SP), PARAMETER :: AU_DISTANCE_IN_ANGSTROM = AU_DISTANCE*1.0E10_SP
!!!***

  INTEGER(I4B), PARAMETER, PRIVATE :: &
       & NConvFact = 10
  TYPE(ConstantT), PRIVATE, SAVE :: &
       & ConvFact(NConvFact)

 DATA ConvFact / & 
      & ConstantT('ELECTRON_VOLT',  'au / eV',          &
      &           'EV',             ELECTRON_VOLT),     &
      & ConstantT('JOULE_PER_MOL',  'au / (J/mol)',     &
      &           'J/MOL',          JOULE_PER_MOL),     &
      & ConstantT('KJOULE_PER_MOL', 'au / (kJ/mol)',    &
      &           'KJ/MOL',         KJOULE_PER_MOL),    &
      & ConstantT('KCAL_PER_MOL',   'au / (kcal/mol)',  &
      &           'KCAL/MOL',       KCAL_PER_MOL),      &
      & ConstantT('INVERSE_CM',     'au / cm^{-1}',     &
      &           'CM^{-1}',        INVERSE_CM),        &
      & ConstantT('ANGSTROM',       'au / Ang',         &
      &           'ANG',            ANGSTROM),          &
      & ConstantT('GRAMS_PER_MOL',  'au / (g/mol)',     &
      &           'G/MOL',          GRAMS_PER_MOL),     &
      & ConstantT('DEBYE',          'au / D',           &
      &           'D',              DEBYE),             &
      & ConstantT('ANGSTROM3',      '(au / Ang)^{3}',   &
      &           'ANG^{3}',        ANGSTROM**3),       &
      & ConstantT('CM3_PER_MOL',    'au / (cm**3/mol)', &
      &           'CM^{3}/MOL',     CM3_PER_MOL) &
      & /


CONTAINS

  SUBROUTINE OutputConstants(FILENAME)
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FILENAME

    INTEGER(I4B) :: iunit

    IF ( PRESENT(FILENAME) ) THEN
       iunit=29
       OPEN(FILE=FILENAME,UNIT=iunit)
    ELSE
       iunit=6
    END IF

    WRITE(iunit,110)
    WRITE(iunit,100)
    WRITE(iunit,110)

    WRITE(iunit,*)
    WRITE(iunit,*)'Fundamental mathematical constants (no units)'

    WRITE(iunit,150)'PI...............................', &
         &           PI
    WRITE(iunit,150)'PI_D_180.........................', &
         &           PI_D_180
    WRITE(iunit,150)'PIO2.............................', &
         &           PIO2
    WRITE(iunit,150)'TWO_PI...........................', &
         &           TWO_PI
    WRITE(iunit,150)'FOUR_PI..........................', &
         &           FOUR_PI
    WRITE(iunit,150)'SQRT_PI..........................', &
         &           SQRT_PI
    WRITE(iunit,150)'SQRT2............................', &
         &           SQRT2
    WRITE(iunit,150)'EULER............................', &
         &           EULER
    WRITE(iunit,150)'DEG2RAD..........................', &
         &           DEG2RAD
    WRITE(iunit,150)'RAD2DEG..........................', &
         &           RAD2DEG

    WRITE(iunit,*)
    WRITE(iunit,*)'Fundamental physical constants (SI units)'

    WRITE(iunit,150)'AVOGADRO_CONSTANT................', &
         &           AVOGADRO_CONSTANT,' mol^{-1}'
    WRITE(iunit,150)'SPEED_OF_LIGHT...................', &
         &           SPEED_OF_LIGHT,' m s^{-1}'
    WRITE(iunit,150)'PLANCK_CONSTANT..................', &
         &           PLANCK_CONSTANT,' J s'
    WRITE(iunit,150)'HBAR.............................', &
         &           HBAR,' J s'
    WRITE(iunit,150)'ELECTRON_CHARGE..................', &
         &           ELECTRON_CHARGE,' C'
    WRITE(iunit,150)'ELECTRON_REST_MASS...............', &
         &           ELECTRON_REST_MASS,' kg'
    WRITE(iunit,150)'NEUTRON_REST_MASS................', &
         &           NEUTRON_REST_MASS,' kg'
    WRITE(iunit,150)'PERMITTIVITY_FREE_SPACE..........', &
         &           PERMITTIVITY_FREE_SPACE,' C^2 N^{-1} m^{-2}'
    WRITE(iunit,150)'FARADAY..........................', &
         &           FARADAY,' C mol{-1}'

    WRITE(iunit,*)
    WRITE(iunit,*)'Atomic Units in terms of SI'

    WRITE(iunit,*)
    WRITE(iunit,*)'Dependent quantities'

    WRITE(iunit,150)'AU_MASS..........................', &
         &           AU_MASS, ' kg / au'
    WRITE(iunit,150)'AU_CHARGE........................', &
         &           AU_CHARGE,' C / au'
    WRITE(iunit,150)'AU_ANGULAR_MOMENTUM..............', &
         &           AU_ANGULAR_MOMENTUM,' kg m s^{-1} / au'
    WRITE(iunit,150)'AU_PERMITTIVITY..................', &
         &           AU_PERMITTIVITY,' C^2 N^{-1} m^{-2} / au'


    WRITE(iunit,*)
    WRITE(iunit,*)'Dependent quantities'

    WRITE(iunit,150)'AU_DISTANCE......................', &
         &           AU_DISTANCE,' m / au'
    WRITE(iunit,150)'AU_DISTANCE_IN_ANGSTROM..........', &
         &           AU_DISTANCE_IN_ANGSTROM,' Ang / au'
    WRITE(iunit,150)'AU_ENERGY........................', &
         &           AU_ENERGY,' J / au'
    WRITE(iunit,150)'AU_VELOCITY......................', &
         &           AU_VELOCITY,' m s^{-1} / au'
    WRITE(iunit,150)'AU_TIME..........................', &
         &           AU_TIME,' s / au'
    WRITE(iunit,150)'AU_ELECTROSTATIC_POTENTIAL.......', &
         &           AU_ELECTROSTATIC_POTENTIAL,' J C^{-1} / au'
    WRITE(iunit,150)'AU_FORCE.........................', & 
         &           AU_FORCE,' N / au'
    WRITE(iunit,150)'AU_MAGNETIC_DIPOLE...............', &
         &           AU_MAGNETIC_DIPOLE,' J T^{-1} / au'
    WRITE(iunit,*)
    WRITE(iunit,*)'Fundamental physical constants (AU)'

    WRITE(iunit,150)'SPEED_OF_LIGHT_AU................', &
         &           SPEED_OF_LIGHT_AU,' au'
    WRITE(iunit,150)'BOLTZMANN_CONSTANT_AU............', &
         &           BOLTZMANN_CONSTANT_AU,' au'
    WRITE(iunit,150)'PLANCK_CONSTANT_AU...............', &
         &           PLANCK_CONSTANT_AU,' au'
    WRITE(iunit,150)'PROTON_REST_MASS_AU..............', &
         &           PROTON_REST_MASS_AU,' au'
    WRITE(iunit,150)'NEUTRON_REST_MASS_AU.............', &
         &           NEUTRON_REST_MASS_AU,' au'

    WRITE(iunit,*)
    WRITE(iunit,*)'Conversion Factors'
    WRITE(iunit,*)
    WRITE(iunit,*)"(Example:If you have a length L_Ang in Angstroms"
    WRITE(iunit,*)"and you want it in AU, L_AU = L_Ang * ANGSTROM)"
    WRITE(iunit,*)

    WRITE(iunit,150)'ELECTRON_VOLT....................', &
         &           ELECTRON_VOLT,' au / eV'
    WRITE(iunit,150)'JOULE_PER_MOL....................', &
         &           JOULE_PER_MOL,' au / (J/mol)'
    WRITE(iunit,150)'KJOULE_PER_MOL...................', &
         &           KJOULE_PER_MOL,' au / (kJ/mol)'
    WRITE(iunit,150)'KCAL_PER_MOL.....................', &
         &           KCAL_PER_MOL,' au / (kcal/mol)'
    WRITE(iunit,150)'KCAL_PER_ANG_MOL.................', &
         &           KCAL_PER_ANG_MOL, ' au / (kcal/(Ang*mol))'
    WRITE(iunit,150)'INVERSE_CM.......................', &
         &           INVERSE_CM,' au / cm^{-1}'
    WRITE(iunit,150)'ANGSTROM.........................', &
         &           ANGSTROM,' au / Ang'
    WRITE(iunit,150)'GRAMS_PER_MOL....................', &
         &           GRAMS_PER_MOL,' au / (g/mol)'
    WRITE(iunit,150)'DEBYE............................', &
         &           DEBYE,' au / D'
    WRITE(iunit,150)'CM3_PER_MOL......................', &
         &           CM3_PER_MOL,' au / (cm**3/mol)'




100 FORMAT('************************** CONSTANTS **************************')
110 FORMAT('***************************************************************')
150 FORMAT(a,E18.10,a)

    IF ( PRESENT(FILENAME) ) THEN
       CLOSE(iunit)
    END IF

  END SUBROUTINE OutputConstants

!!$===========================================================================
  FUNCTION GetUnitConversionFactor(string) RESULT(cf)
!!!****f* ConstantsMod/GetUnitConversionFactor
!!!
!!! NAME
!!!     GetUnitConversionFactor -- gets the unit conversion factor in au.
!!! USAGE
!!!     Dipole_au = Dipole_Debye * GetUnitConversionFactor('DEBYE')
!!!     Dipole_Debye = Dipole_au / GetUnitConversionFactor('DEBYE')
!!!     Dipole_Debye = Dipole_au * GetUnitConversionFactor('1/DEBYE')
!!! DESCRIPTION
!!!     Gets the appropriate unit conversion factor from ConstantsMod.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     REAL(SP) :: cf
!!! SEE ALSO
!!!     OutputConstants
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string

    REAL(SP) :: cf
    CHARACTER(LEN=LEN_TRIM(string)) :: tmp_string
    INTEGER(I4B) :: i,l
    LOGICAL(LGD) :: Linvert

    l = LEN_TRIM(string)
    tmp_string = TRIM(string)

    IF ( tmp_string(1:2) == '1/' ) THEN
       Linvert = .TRUE.
       IF ( l > 2 ) THEN
          tmp_string = tmp_string(3:l)//'  '
          l=l-2
       END IF
    ELSE
       Linvert = .FALSE.
    END IF

    cf = 1.0_SP

    DO i=1,NConvFact
       IF (   tmp_string == ConvFact(i)%Name(1:l) .OR. &
            & tmp_string == ConvFact(i)%Unit(1:l) .OR. &
            & tmp_string == ConvFact(i)%Abbrev(1:l) ) THEN
          cf = ConvFact(i)%Value
          EXIT
       END IF
    END DO

    IF ( Linvert ) cf = 1.0_SP / cf

  END FUNCTION GetUnitConversionFactor

END MODULE ConstantsMod










MODULE errormod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

    INTEGER(I4B), SAVE :: BOMBLev = 5
    LOGICAL(LGD), SAVE :: BOMBCore = .FALSE.




CONTAINS
   SUBROUTINE error(string,i1,i2,i3,r1,r2,r3,l1,l2,l3,a1,a2,a3,WARNLev)
      USE DataTypes
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER(I4B), OPTIONAL, INTENT(IN) :: i1,i2,i3
      REAL(SP), OPTIONAL, INTENT(IN) :: r1,r2,r3
      LOGICAL(LGD), OPTIONAL, INTENT(IN) :: l1,l2,l3
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: a1,a2,a3
      INTEGER(I4B), OPTIONAL, INTENT(IN) :: WARNLev
      REAL(SP) :: BIG_NUM
      IF ( PRESENT(WARNLev) ) THEN
         WRITE(*,10)WARNLev
10       FORMAT('WARNING LEVEL ',i2,' DETECTED BY SUBROUTINE error with TAG:')
      ELSE
         WRITE(*,20)
20       FORMAT('FATAL ERROR DETECTED BY SUBROUTINE error with TAG:')
      END IF
      WRITE(*,'(a)')string
      IF ( PRESENT(i1) ) WRITE(*,'(a,i10)')'     i1 = ',i1
      IF ( PRESENT(i2) ) WRITE(*,'(a,i10)')'     i2 = ',i2
      IF ( PRESENT(i3) ) WRITE(*,'(a,i10)')'     i3 = ',i3
      IF ( PRESENT(r1) ) WRITE(*,'(a,e16.8)')'     r1 = ',r1
      IF ( PRESENT(r2) ) WRITE(*,'(a,e16.8)')'     r2 = ',r2
      IF ( PRESENT(r3) ) WRITE(*,'(a,e16.8)')'     r3 = ',r3
      IF ( PRESENT(l1) ) WRITE(*,'(a,1l1)')'     l1 = ',l1
      IF ( PRESENT(l2) ) WRITE(*,'(a,1l1)')'     l2 = ',l2
      IF ( PRESENT(l3) ) WRITE(*,'(a,1l1)')'     l3 = ',l3
      IF ( PRESENT(a1) ) WRITE(*,'(a,a)')'     a1 = ',a1
      IF ( PRESENT(a2) ) WRITE(*,'(a,a)')'     a2 = ',a2
      IF ( PRESENT(a3) ) WRITE(*,'(a,a)')'     a3 = ',a3
      IF (PRESENT(WARNLev)) THEN
         IF (WARNLev < BOMBLev) RETURN
      WRITE(*,'(a,i2)')'BOMBLev is ',BOMBLev
      ENDIF
      IF (BOMBCore) THEN
         WRITE(*,*)'Dumping Core for Back Trace'
         BIG_NUM = 3000.0_SP
         BIG_NUM = exp(BIG_NUM)
      ELSE 
         STOP 'Program terminated by error.'
      ENDIF
   END SUBROUTINE error

END MODULE errormod


!doc=INCOMPLETE/CHECK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!   Written by Darrin M. York, Department of Chemistry, Harvard University,  !
!      Cambridge, MA 02138.                                                  !
!   Copyright by Dr. Darrin M. York (1997), all rights reserved.             !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE UtilitiesMod
!!!****h* YorkLib/UtilitiesMod
!!!
!!! NAME
!!!    UtilitiesMod
!!!
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!!
!!! AUTHOR
!!!   Darrin M. York
!!!
!!! CREATION DATE
!!!   1997
!!!
!!! DESCRIPTION
!!!   Module contains general purpose utility tools.
!!!   Many utilities have been adapted from:
!!!      Numerical Recipes in Fortran 90 by W. H. Press, S. A. Teukolsky,
!!!      W. T. Vetterling, and B. P. Flannery.
!!!
!!! TODO
!!!
!!!
!!! USES
!!!   DataTypes, ConstantsMod, (ErrorMod, ONLY : error)
!!!
!!! SEE ALSO
!!!
!!!
!!! USED BY
!!!
!!!
!!!***

  USE DataTypes
  USE ConstantsMod
  USE ErrorMod, ONLY : error
  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
  
!!!****f* UtilitiesMod/array_copy
!!!
!!! NAME
!!!    array_copy
!!!     
!!! USAGE
!!!    array_copy(src,dest,n_copied,n_not_copied)
!!!    
!!! DESCRIPTION
!!!    Determines the dimension of src and dest, takes the minimum value
!!!    of that set (n_copied) and copies src into dest from 1:n_copied 
!!!
!!! INPUTS
!!!    src(:)
!!!    
!!! OUTPUTS
!!!    dest(:), n_copied, n_not_copied
!!!    
!!!*** 
  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_i
  END INTERFACE


!!!****f* UtilitiesMod/swap
!!!
!!! NAME
!!!    swap
!!!     
!!! USAGE
!!!    swap(a,b,mask)
!!!    
!!! DESCRIPTION
!!!    lets a take the value of b, while b takes the value of a
!!!    under the optional condition mask
!!!
!!! INPUTS
!!!    a,b  (can be of any type)
!!!    OPTIONALLY: LOGICAL(LGD), INTENT(IN) :: mask
!!! OUTPUTS
!!!    a,b
!!!    
!!!*** 
  INTERFACE swap
     MODULE PROCEDURE swap_i, &
          swap_r,swap_rv,swap_rm, &
          swap_c,swap_cv,swap_cm, &
          masked_swap_i4bs,masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE



!!!****f* UtilitiesMod/reallocate
!!!
!!! NAME
!!!    reallocate
!!!     
!!! USAGE
!!!    reallocate(p,n,M)
!!!    
!!! DESCRIPTION
!!!    reallocates p to size n if p is a vector
!!!    optionally reallocates p to size n by m if p is a matrix
!!!
!!! INPUTS
!!!    p can be of any type
!!!    n,M ->  integer(I4B)
!!!    
!!!*** 
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm, &
          & reallocate_iv,reallocate_im,           &
          & reallocate_lv,reallocate_lm,           &
          & reallocate_hv,reallocate_hnv
  END INTERFACE




!!!****f* UtilitiesMod/AllocateIfNeeded
!!!
!!! NAME
!!!    AllocateIfNeeded
!!!     
!!! USAGE
!!!    AllocateIfNeeded(p,n,M)
!!!    
!!! DESCRIPTION
!!!    calls function reallocate if p is not associated
!!!    
!!!
!!! INPUTS
!!!    p can be of any type
!!!    n,M ->  integer(I4B)
!!!    
!!! SEE ALSO
!!!    reallocate
!!!
!!!*** 
  INTERFACE AllocateIfNeeded
     MODULE PROCEDURE AllocateIfNeeded_rv,AllocateIfNeeded_rm,    &
          & AllocateIfNeeded_iv,AllocateIfNeeded_im,              &
          & AllocateIfNeeded_lv,AllocateIfNeeded_lm,              &
          & AllocateIfNeeded_hv,AllocateIfNeeded_hnv
  END INTERFACE




!!!****f* UtilitiesMod/minval
!!!
!!! NAME
!!!    minval
!!!     
!!! USAGE
!!!    minval
!!!    
!!! DESCRIPTION
!!!    returns the minimum value element of a matrix of rank 2 or 3
!!!    
!!!
!!! INPUTS
!!!    a -> real or complex, rank = 2 or 3
!!!    
!!! OUTPUTS
!!!    result is the minimum value
!!!
!!!*** 
  INTERFACE minvalue
     MODULE PROCEDURE minval_rm, minval_rt, minval_im, minval_it
  END INTERFACE





!!!****f* UtilitiesMod/maxval
!!!
!!! NAME
!!!    maxval
!!!     
!!! USAGE
!!!    maxval
!!!    
!!! DESCRIPTION
!!!    returns the maximum value element of a matrix of rank 2 or 3
!!!    
!!!
!!! INPUTS
!!!    a -> real or integer, rank = 2 or 3
!!!    
!!! OUTPUTS
!!!    result is the maximum value
!!!
!!!*** 
  INTERFACE maxvalue
     MODULE PROCEDURE maxval_rm, maxval_rt, maxval_im, maxval_it
  END INTERFACE







!!!****f* UtilitiesMod/imaxloc
!!!
!!! NAME
!!!    imaxloc
!!!     
!!! USAGE
!!!    imaxloc(arr)
!!!    
!!! DESCRIPTION
!!!    returns the integer value of the element corresponding to the maximum 
!!!    value in the vector, arr.
!!!    
!!!
!!! INPUTS
!!!    arr -> real or integer, must be a vector
!!!    
!!! OUTPUTS
!!!   returns an integer
!!!*** 
  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE







!!!****f* UtilitiesMod/iminloc
!!!
!!! NAME
!!!    iminloc
!!!     
!!! USAGE
!!!    iminloc(arr)
!!!    
!!! DESCRIPTION
!!!    returns the integer value of the element containing the minimum 
!!!    value in the vector, arr.
!!!    
!!!
!!! INPUTS
!!!    arr -> real or integer, must be a vector
!!!    
!!! OUTPUTS
!!!   returns an integer
!!!*** 
  INTERFACE iminloc
     MODULE PROCEDURE iminloc_r,iminloc_i
  END INTERFACE





!!!****f* UtilitiesMod/assert
!!!
!!! NAME
!!!    assert
!!!     
!!! USAGE
!!!    assert(n1,N2,N3,N4,string)
!!!      or
!!!    assert(n,string)
!!!    
!!! DESCRIPTION
!!!    Stops the current program if (.NOT. (n1 .AND. n2 .AND. n3 .AND. n4))
!!!    and then prints string to screen, where n2, n3, and n4 are optional
!!!
!!!    if more than four logical descriptors are needed, they can be
!!!    passed into the routine as a vector
!!!
!!! INPUTS
!!!        CHARACTER(LEN=*), INTENT(IN) :: string
!!!        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
!!!        LOGICAL, DIMENSION(:), INTENT(IN) :: n
!!!*** 
  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE







!!!****f* UtilitiesMod/assert
!!!
!!! NAME
!!!    assert
!!!     
!!! USAGE
!!!    assert(n1,N2,N3,N4,string)
!!!      or
!!!    assert(n,string)
!!!    
!!! DESCRIPTION
!!!    Stops the current program if any 2 of the 4 n*'s are not equal
!!!    and then prints string to screen, where n2, n3, and n4 are optional
!!!    if all four n*'s are equal the result of the function is n1
!!!
!!!    if more than four logical descriptors are needed, they can be
!!!    passed into the routine as a vector, in which case the result of the
!!!    function is the first element of the vector
!!!  
!!!
!!! INPUTS
!!!        CHARACTER(LEN=*), INTENT(IN) :: string
!!!        INTEGER, INTENT(IN) :: n1,n2,n3,n4
!!!        INTEGER, DIMENSION(:), INTENT(IN) :: n
!!!*** 
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE




!!!****f* UtilitiesMod/arth
!!!
!!! NAME
!!!    arth
!!!     
!!! USAGE
!!!    arth(first,increment,n)
!!!    
!!! DESCRIPTION
!!!    returns an array of length n containing an arithmetic progression whose first
!!!    value is first and whose increment is increment. If first and increment 
!!!    have rank greater than zero, returns an array of one larger rank, with 
!!!    the last subscript having size n and indexing the progressions.
!!!    
!!!
!!! INPUTS
!!!        first,increment  -> real or integer
!!!        INTEGER(I4B), INTENT(IN) :: n
!!!    
!!!*** 
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_i
  END INTERFACE
  !   INTERFACE arth_DP
  !      MODULE PROCEDURE arth_DP
  !   END INTERFACE







!!!****f* UtilitiesMod/geop
!!!
!!! NAME
!!!    geop
!!!     
!!! USAGE
!!!    geop(first,factor,r)
!!!    
!!! DESCRIPTION
!!!   Returns an array of length n containing a geometric progression whose first
!!!   value is first and whose multiplier is factor.  If first and factor have
!!!   rank greter than zero, returns an array of one larger rank, with the last
!!!   subscript having size n and indexing progression.
!!!    
!!!
!!! INPUTS
!!!        first,factor -> real, integer, or complex
!!!        INTEGER(I4B), INTENT(IN) :: n
!!!    
!!!*** 
  INTERFACE geop
     MODULE PROCEDURE geop_r, geop_i, geop_c
  END INTERFACE




!!!****f* UtilitiesMod/cumsum
!!!
!!! NAME
!!!    cumsum
!!!     
!!! USAGE
!!!    cumsum(arr,seed)
!!!    
!!! DESCRIPTION
!!!   Given the rank 1 array arr of type real or integer, returns an array of
!!!   identical type and size containing the cumulative sums of arr.  If the optional
!!!   argument seed is present, it is added to the first component of the result.
!!!
!!! INPUTS
!!!        arr,seed -> real or integer
!!!    
!!!*** 
  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE




!!!****f* UtilitiesMod/poly
!!!
!!! NAME
!!!    poly
!!!     
!!! USAGE
!!!    poly(x,coeffs,mask)
!!!    
!!! DESCRIPTION
!!!    Returns a scalar value or array with the same type and shape as x,
!!!    containing the result of evaluating the polynomial with one 
!!!    dimensional coefficient vector coeffs on each component of x.  
!!!    The optional argument mask, if present, has the same thape as x, 
!!!    and suppresses evaluations of the polynomial where its components are .false..
!!! INPUTS
!!!        x,coeffs-> real or integer
!!!         LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: mask
!!!    
!!!*** 
  INTERFACE poly
     MODULE PROCEDURE poly_rr,poly_rrv, &
          poly_rc,poly_cc,poly_msk_rrv
  END INTERFACE




!!!****f* UtilitiesMod/poly_term
!!!
!!! NAME
!!!    poly_term
!!!     
!!! USAGE
!!!    poly_term(a,x)
!!!    
!!! DESCRIPTION
!!!    Returns an array of type and size as the one-dimensional array a,
!!!    containing the partial cumulants of the polynomial with coefficients a
!!!    (arranged from highest-order to lowest-order coefficients, n.b.) evaluated
!!!    at x.  This is equivalent to synthetic division, and can be parallelized.
!!!    
!!! INPUTS
!!!        a,x-> real or integer
!!!    
!!!*** 
  INTERFACE poly_term
     MODULE PROCEDURE poly_term_rr,poly_term_cc
  END INTERFACE







!!!****f* UtilitiesMod/EXTENDED_PRODUCT
!!!
!!! NAME
!!!    EXTENDED_PRODUCT
!!!     
!!! USAGE
!!!    EXTENDED_PRODUCT(a)
!!!    
!!! DESCRIPTION
!!!    Returns the extended product of an array a.
!!!    
!!! INPUTS
!!!        a -> real
!!!    
!!!*** 
  INTERFACE EXTENDED_PRODUCT
     MODULE PROCEDURE extendedprod_r
  END INTERFACE

!!!****f* UtilitiesMod/OUTER_PRODUCT
!!!
!!! NAME
!!!    OUTER_PRODUCT
!!!     
!!! USAGE
!!!    OUTER_PRODUCT(a,b)
!!!    
!!! DESCRIPTION
!!!    returns the direct product (tensor product) of arrays a and b.
!!!    
!!! INPUTS
!!!        a,b -> real
!!!    
!!!*** 
  INTERFACE OUTER_PRODUCT
     MODULE PROCEDURE outerprod_r
  END INTERFACE



!!!****f* UtilitiesMod/outerdiff
!!!
!!! NAME
!!!    outerdiff
!!!     
!!! USAGE
!!!    outerdiff(a,b)
!!!    
!!! DESCRIPTION
!!!    returns the outer difference of two vectors
!!!    
!!! INPUTS
!!!        a,b -> real or integer
!!!    
!!!*** 
  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_i
  END INTERFACE




!!!****f* UtilitiesMod/scatter_add
!!!
!!! NAME
!!!    scatter_add
!!!     
!!! USAGE
!!!    scatter_add(dest,source,dest_index)
!!!    
!!! DESCRIPTION
!!!    Adds each component of the array source into a component of dest 
!!!    specified by the index dest_index. (The user will usually have zeroed 
!!!    dest before the call to this routine.)  Note that dest_index has the 
!!!    size of source, but must contain vaules in the range from 1 to 
!!!    size(dest), inclusive.  Out-of-range values are ignored.
!!!    
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: source
!!!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
!!! OUTPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
!!!*** 
  INTERFACE scatter_add
     MODULE PROCEDURE scatter_add_r
  END INTERFACE




!!!****f* UtilitiesMod/scatter_max
!!!
!!! NAME
!!!    scatter_max
!!!     
!!! USAGE
!!!    scatter_max(dest,source,dest_index)
!!!    
!!! DESCRIPTION
!!!    Takes the max operation between each component of the array source 
!!!    and a component of dest specified by the index array dest_index, 
!!!    replacing that component of dest with the value obtained 
!!!    ("maxing into" operation).  (The user will often want to fill the 
!!!    array dest with the value -huge before the call to this subroutine.)  
!!!    Note that dest_index has the size of source, but must contain values 
!!!    in the range from 1 to size(dest), inclusive.  Out-of-range 
!!!    values are ignored.
!!!    
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: source
!!!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
!!! OUTPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
!!!*** 
  INTERFACE scatter_max
     MODULE PROCEDURE scatter_max_r
  END INTERFACE




!!!****f* UtilitiesMod/diagadd
!!!
!!! NAME
!!!    diagadd
!!!     
!!! USAGE
!!!    diagadd(mat,diag)
!!!    
!!! DESCRIPTION
!!!    The argument diag, either a scalar or else a vector whose size 
!!!    must be the smaller of the two dimensions of matrix mat, is 
!!!    added to the diagonal of the matrix mat.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
!!!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
!!! OUTPUTS
!!!    mat
!!!*** 
  INTERFACE diagadd
     MODULE PROCEDURE diagadd_rv,diagadd_r
  END INTERFACE




!!!****f* UtilitiesMod/diagmult
!!!
!!! NAME
!!!    diagmult
!!!     
!!! USAGE
!!!    diagmult(mat,diag)
!!!    
!!! DESCRIPTION
!!!    The argument diag, either a scalar or else a vector whose size 
!!!    must be the smaller of the two dimensions of matrix mat, is 
!!!    multiplied to the diagonal of the matrix mat.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
!!!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
!!! OUTPUTS
!!!    mat
!!!*** 
  INTERFACE diagmult
     MODULE PROCEDURE diagmult_rv,diagmult_r
  END INTERFACE




!!!****f* UtilitiesMod/get_diag
!!!
!!! NAME
!!!    get_diag
!!!     
!!! USAGE
!!!    get_diag(mat)
!!!    
!!! DESCRIPTION
!!!    Returns a vector containing the diagonal values of the matrix mat.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
!!!*** 
  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv
  END INTERFACE




!!!****f* UtilitiesMod/put_diag
!!!
!!! NAME
!!!    put_diag
!!!     
!!! USAGE
!!!    put_diag(mat)
!!!    
!!! DESCRIPTION
!!!    Sets the diagonal of matrix mat equal to the argument diag, 
!!!    either a scalar or else a vector whose size must be the 
!!!    smaller of the two dimensions of matrix mat.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
!!!    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
!!! OUTPUTS
!!!    mat
!!!*** 
  INTERFACE put_diag
     MODULE PROCEDURE put_diag_rv, put_diag_r
  END INTERFACE




!!!****f* UtilitiesMod/SUM
!!!
!!! NAME
!!!    SUM
!!!     
!!! USAGE
!!!    SUM(mat)
!!!    
!!! DESCRIPTION
!!!    Returns the sum of the elements of a real tensor of rank 2 or 3
!!! INPUTS
!!!    mat -> real
!!!*** 
  INTERFACE SUM
     MODULE PROCEDURE SUM_2ind_MATRIX, SUM_3ind_MATRIX
  END INTERFACE




!!!****f* UtilitiesMod/ALMOST_EQUAL
!!!
!!! NAME
!!!    ALMOST_EQUAL
!!!     
!!! USAGE
!!!    ALMOST_EQUAL(a,B,C)
!!!    
!!! DESCRIPTION
!!!    Determines of the difference between a and b is tiny, optionally,
!!!    if the difference between three numbers is tiny (a, b, and c).
!!!    If more than 3 numbers are required, then pass a vector a, and
!!!    the difference between the maximum value and the minimum value 
!!!    are compared.  You can also pass 2 vectors and compare the vectors.
!!! INPUTS
!!!    real numbers... maybe scalar maybe vectors
!!! OUTPUTS
!!!    returns .TRUE. or .FALSE.
!!!*** 
  INTERFACE ALMOST_EQUAL
     MODULE PROCEDURE ALMOST_EQUAL_S2, ALMOST_EQUAL_S3, &
          ALMOST_EQUAL_V1, ALMOST_EQUAL_V2
  END INTERFACE


!!!****f* UtilitiesMod/IndexSearch
!!!
!!! NAME
!!!    IndexSearch
!!!     
!!! USAGE
!!!    IndexSearch(avec,a)
!!!    
!!! DESCRIPTION
!!!    Returns the integer index value that the number/character a apears in avec
!!! INPUTS
!!!    a,avec -> real,integer,character,complex
!!! OUTPUTS
!!!    returns .TRUE. or .FALSE.
!!!*** 
  INTERFACE IndexSearch
     MODULE PROCEDURE IndexSearch_string, IndexSearch_iv, &
                    & IndexSearch_rv, IndexSearch_cv
  END INTERFACE


  INTERFACE push
     MODULE PROCEDURE Push_real_v, Push_real_m, Push_int_v, Push_int_m, Push_char_v, Push_logic_v
  END INTERFACE

  INTERFACE shift
     MODULE PROCEDURE shift_real_v, shift_int_v, shift_logic_v, shift_char_v
  END INTERFACE


!!!****f* UtilitiesMod/estimate_real_size
!!!
!!! NAME
!!!    estimate_real_size
!!!     
!!! USAGE
!!!    estimate_real_size(n, THRESH)
!!!    estimate_real_size(n, m, THRESH)
!!!    
!!! DESCRIPTION
!!!    Prints the size of the square or rectangluar matrix (if it is larger
!!!    than THRESH megabytes) in megabytes to stdout.  The default is to
!!!    print if the matrix is more than 100 megabytes.
!!!
!!! INPUTS
!!!    integer(i4b), intent(in) :: n,m
!!!    real(sp), optional, intent(in) :: THRESH
!!! OUTPUTS
!!!    prints to standard out
!!!*** 
  interface estimate_real_size
     module procedure estimate_real_size_sq, estimate_real_size_rect
  end interface

CONTAINS

!!$===========================================================================
  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(SP), DIMENSION(:), INTENT(IN) :: src
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=MIN(SIZE(src),SIZE(dest))
    n_not_copied=SIZE(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r

  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=MIN(SIZE(src),SIZE(dest))
    n_not_copied=SIZE(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i

  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r

  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv

  SUBROUTINE swap_rm(a,b)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a,1),SIZE(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rm

  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c

  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv

  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a,1),SIZE(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm

  SUBROUTINE masked_swap_i4bs(a,b,mask)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    LOGICAL(LGD), INTENT(IN) :: mask
    INTEGER(I4B) :: swp
    IF (mask) THEN
       swp=a
       a=b
       b=swp
    END IF
  END SUBROUTINE masked_swap_i4bs

  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGD), INTENT(IN) :: mask
    REAL(SP) :: swp
    IF (mask) THEN
       swp=a
       a=b
       b=swp
    END IF
  END SUBROUTINE masked_swap_rs

  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(SIZE(a)) :: swp
    WHERE (mask)
       swp=a
       a=b
       b=swp
    END WHERE
  END SUBROUTINE masked_swap_rv

  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGD), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(SIZE(a,1),SIZE(a,2)) :: swp
    WHERE (mask)
       swp=a
       a=b
       b=swp
    END WHERE
  END SUBROUTINE masked_swap_rm
!!$===========================================================================

!!$===========================================================================
  FUNCTION reallocate_rv(p,n)
    REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    ALLOCATE(reallocate_rv(n),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_rv: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p)
    reallocate_rv(1:MIN(nold,n))=p(1:MIN(nold,n))
    IF ( nold < n ) reallocate_rv(nold+1:n) = 0.0_SP 
    DEALLOCATE(p)
  END FUNCTION reallocate_rv

  FUNCTION reallocate_iv(p,n)
    INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    ALLOCATE(reallocate_iv(n),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_iv: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p)
    reallocate_iv(1:MIN(nold,n))=p(1:MIN(nold,n))
    IF ( nold < n ) reallocate_iv(nold+1:n) = 0_I4B
    DEALLOCATE(p)
  END FUNCTION reallocate_iv

  FUNCTION reallocate_rm(p,n,m)
    REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    ALLOCATE(reallocate_rm(n,m),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_rm: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p,1)
    mold=SIZE(p,2)
    reallocate_rm(1:MIN(nold,n),1:MIN(mold,m))=&
         p(1:MIN(nold,n),1:MIN(mold,m))
    IF ( nold < n ) reallocate_rm(nold+1:n,:) = 0.0_SP
    IF ( mold < m ) reallocate_rm(:,mold+1:m) = 0.0_SP
    DEALLOCATE(p)
  END FUNCTION reallocate_rm

  FUNCTION reallocate_im(p,n,m)
    INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    ALLOCATE(reallocate_im(n,m),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_im: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p,1)
    mold=SIZE(p,2)
    reallocate_im(1:MIN(nold,n),1:MIN(mold,m))=&
         p(1:MIN(nold,n),1:MIN(mold,m))
    IF ( nold < n ) reallocate_im(nold+1:n,:) = 0_I4B
    IF ( mold < m ) reallocate_im(:,mold+1:m) = 0_I4B
    DEALLOCATE(p)
  END FUNCTION reallocate_im

  FUNCTION reallocate_lv(p,n)
    LOGICAL(LGD), DIMENSION(:), POINTER :: p, reallocate_lv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    IF (.NOT. ASSOCIATED(p)) RETURN
    ALLOCATE(reallocate_lv(n),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_lv: problem in attempt to ALLOCATE memory')
    nold=SIZE(p)
    reallocate_lv(1:MIN(nold,n))=p(1:MIN(nold,n))
    IF ( nold < n ) reallocate_lv(nold+1:n) = .FALSE.
    DEALLOCATE(p)
  END FUNCTION reallocate_lv

  FUNCTION reallocate_lm(p,n,m)
    LOGICAL(LGD), DIMENSION(:,:), POINTER :: p, reallocate_lm
    INTEGER(I4B), INTENT(IN) :: n,m
    INTEGER(I4B) :: nold,mold,ierr
    ALLOCATE(reallocate_lm(n,m),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_lm: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p,1)
    mold=SIZE(p,2)
    reallocate_lm(1:MIN(nold,n),1:MIN(mold,m))=&
         p(1:MIN(nold,n),1:MIN(mold,m))
    IF ( nold < n ) reallocate_lm(nold+1:n,:) = .FALSE.
    IF ( mold < m ) reallocate_lm(:,mold+1:m) = .FALSE.
    DEALLOCATE(p)
  END FUNCTION reallocate_lm

  FUNCTION reallocate_hv(p,n)
    CHARACTER(LEN=*), DIMENSION(:), POINTER :: p
    CHARACTER(LEN=1), DIMENSION(:), POINTER :: reallocate_hv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
    ALLOCATE(reallocate_hv(n),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_hv: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p)
    reallocate_hv(1:MIN(nold,n))=p(1:MIN(nold,n))
    IF ( nold < n ) reallocate_hv(nold+1:n) = ' '
    DEALLOCATE(p)
  END FUNCTION reallocate_hv

  FUNCTION reallocate_hnv(p,n,LEN)
    INTEGER(I4B), INTENT(IN) :: LEN
    CHARACTER(LEN=*), DIMENSION(:), POINTER :: p
    CHARACTER(LEN=LEN), DIMENSION(:), POINTER :: reallocate_hnv
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B) :: nold,ierr
!    WRITE(6,*)"IN REALLOCATE_HNV"
    ALLOCATE(reallocate_hnv(n),stat=ierr)
    IF (ierr /= 0) CALL &
         error('reallocate_hnv: problem in attempt to ALLOCATE memory')
    IF (.NOT. ASSOCIATED(p)) RETURN
    nold=SIZE(p)
!    WRITE(6,*)"REALLOCATE_HNV NOLD = ",nold,n
    reallocate_hnv(1:MIN(nold,n))=p(1:MIN(nold,n))
!    WRITE(6,*)"LEAVING"
    IF ( nold < n ) reallocate_hnv(nold+1:n) = ' '
    DEALLOCATE(p)
!    WRITE(6,*)"OUT"
  END FUNCTION reallocate_hnv
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE AllocateIfNeeded_iv(p,n,TRIM)
    INTEGER(I4B), DIMENSION(:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 ) THEN
       p => reallocate(p,n)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p) /= n ) THEN
          p => reallocate(p,n)
       END IF
    ELSE IF ( SIZE(p) < n ) THEN
       p => reallocate(p,n)
    END IF

  END SUBROUTINE AllocateIfNeeded_iv

  SUBROUTINE AllocateIfNeeded_im(p,n,m,TRIM)
    INTEGER(I4B), DIMENSION(:,:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n, m
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 .AND. m > 0 ) THEN
       p => reallocate(p,n,m)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p,1) /= n .OR. SIZE(p,2) /= m ) THEN
          p => reallocate(p,n,m)
       END IF
    ELSE IF ( SIZE(p,1) < n .OR. SIZE(p,2) < m ) THEN
       p => reallocate(p,n,m)
    END IF

  END SUBROUTINE AllocateIfNeeded_im

  SUBROUTINE AllocateIfNeeded_rv(p,n,TRIM)
    REAL(SP), DIMENSION(:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 ) THEN
       p => reallocate(p,n)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p) /= n ) THEN
          p => reallocate(p,n)
       END IF
    ELSE IF ( SIZE(p) < n ) THEN
       p => reallocate(p,n)
    END IF

  END SUBROUTINE AllocateIfNeeded_rv

 SUBROUTINE AllocateIfNeeded_rm(p,n,m,TRIM)
    REAL(SP), DIMENSION(:,:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n, m
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 .AND. m > 0 ) THEN
       p => reallocate(p,n,m)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p,1) /= n .OR. SIZE(p,2) /= m ) THEN
          p => reallocate(p,n,m)
       END IF
    ELSE IF ( SIZE(p,1) < n .OR. SIZE(p,2) < m ) THEN
       p => reallocate(p,n,m)
    END IF

  END SUBROUTINE AllocateIfNeeded_rm

  SUBROUTINE AllocateIfNeeded_lv(p,n,TRIM)
    LOGICAL(LGD), DIMENSION(:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 ) THEN
       p => reallocate(p,n)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p) /= n ) THEN
          p => reallocate(p,n)
       END IF
    ELSE IF ( SIZE(p) < n ) THEN
       p => reallocate(p,n)
    END IF

  END SUBROUTINE AllocateIfNeeded_lv

 SUBROUTINE AllocateIfNeeded_lm(p,n,m,TRIM)
    LOGICAL(LGD), DIMENSION(:,:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n, m
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 .AND. m > 0 ) THEN
       p => reallocate(p,n,m)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p,1) /= n .OR. SIZE(p,2) /= m ) THEN
          p => reallocate(p,n,m)
       END IF
    ELSE IF ( SIZE(p,1) < n .OR. SIZE(p,2) < m ) THEN
       p => reallocate(p,n,m)
    END IF

  END SUBROUTINE AllocateIfNeeded_lm

  SUBROUTINE AllocateIfNeeded_hv(p,n,TRIM)
    CHARACTER(LEN=1), DIMENSION(:), POINTER :: p
    INTEGER(I4B), INTENT(IN) :: n
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 ) THEN
       p => reallocate(p,n)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p) /= n ) THEN
          p => reallocate(p,n)
       END IF
    ELSE IF ( SIZE(p) < n ) THEN
       p => reallocate(p,n)
    END IF

  END SUBROUTINE AllocateIfNeeded_hv

  SUBROUTINE AllocateIfNeeded_hnv(p,n,LEN,TRIM)
    INTEGER(I4B), INTENT(IN) :: n
    INTEGER(I4B), INTENT(IN) :: LEN
    CHARACTER(LEN=LEN), DIMENSION(:), POINTER :: p
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRIM

    IF ( .NOT. ASSOCIATED(p) .AND. n > 0 ) THEN
       p => reallocate(p,n,LEN)
    ELSE IF ( OPTIONAL_FLAG(TRIM) ) THEN
       IF ( SIZE(p) /= n ) THEN
          p => reallocate(p,n,LEN)
       END IF
    ELSE IF ( SIZE(p) < n ) THEN
       p => reallocate(p,n,LEN)
    END IF

  END SUBROUTINE AllocateIfNeeded_hnv
!!$===========================================================================

!!$===========================================================================
  FUNCTION IFirstloc(mask)
    LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: mask
    INTEGER(I4B) :: IFirstloc
    INTEGER(I4B), DIMENSION(1) :: loc
    loc=MAXLOC(MERGE(1,0,mask))
    IFirstloc=loc(1)
    IF (.NOT. mask(ifirstloc)) IFirstloc=SIZE(mask)+1
  END FUNCTION IFirstloc
!!$===========================================================================

!!$===========================================================================
  FUNCTION minval_rm(rm) RESULT(mv)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: rm
    REAL(SP) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(rm,2) > 0 ) THEN
       mv = MINVAL( rm(1:SIZE(rm,1),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(rm,2) > 1 ) THEN
       DO i=2,SIZE(rm,2)
          mv = MIN( mv, MINVAL( rm(1:SIZE(rm,1),i) ) )
       END DO
    END IF
  END FUNCTION minval_rm

  FUNCTION minval_rt(rt) RESULT(mv)
    REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: rt
    REAL(SP) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(rt,3) > 0 ) THEN
       mv = MINVAL( rt(1:SIZE(rt,1),1:SIZE(rt,2),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(rt,3) > 1 ) THEN
       DO i=2,SIZE(rt,3)
          mv = MIN( mv, MINVAL( rt(1:SIZE(rt,1),1:SIZE(rt,2),i) ) )
       END DO
    END IF
  END FUNCTION minval_rt

  FUNCTION minval_im(im) RESULT(mv)
    INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: im
    INTEGER(I4B) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(im,2) > 0 ) THEN
       mv = MINVAL( im(1:SIZE(im,1),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(im,2) > 1 ) THEN
       DO i=2,SIZE(im,2)
          mv = MIN( mv, MAXVAL( im(1:SIZE(im,1),i) ) )
       END DO
    END IF
  END FUNCTION minval_im

  FUNCTION minval_it(it) RESULT(mv)
    INTEGER(I4B), DIMENSION(:,:,:), INTENT(IN) :: it
    INTEGER(I4B) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(it,3) > 0 ) THEN
       mv = MINVAL( it(1:SIZE(it,1),1:SIZE(it,2),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(it,3) > 1 ) THEN
       DO i=2,SIZE(it,3)
          mv = MIN( mv, MAXVAL( it(1:SIZE(it,1),1:SIZE(it,2),i) ) )
       END DO
    END IF
  END FUNCTION minval_it
!!$===========================================================================

!!$===========================================================================
  FUNCTION maxval_rm(rm) RESULT(mv)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: rm
    REAL(SP) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(rm,2) > 0 ) THEN
       mv = MAXVAL( rm(1:SIZE(rm,1),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(rm,2) > 1 ) THEN
       DO i=2,SIZE(rm,2)
          mv = MAX( mv, MAXVAL( rm(1:SIZE(rm,1),i) ) )
       END DO
    END IF
  END FUNCTION maxval_rm

  FUNCTION maxval_rt(rt) RESULT(mv)
    REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: rt
    REAL(SP) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(rt,3) > 0 ) THEN
       mv = MAXVAL( rt(1:SIZE(rt,1),1:SIZE(rt,2),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(rt,3) > 1 ) THEN
       DO i=2,SIZE(rt,3)
          mv = MAX( mv, MAXVAL( rt(1:SIZE(rt,1),1:SIZE(rt,2),i) ) )
       END DO
    END IF
  END FUNCTION maxval_rt

  FUNCTION maxval_im(im) RESULT(mv)
    INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: im
    INTEGER(I4B) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(im,2) > 0 ) THEN
       mv = MAXVAL( im(1:SIZE(im,1),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(im,2) > 1 ) THEN
       DO i=2,SIZE(im,2)
          mv = MAX( mv, MAXVAL( im(1:SIZE(im,1),i) ) )
       END DO
    END IF
  END FUNCTION maxval_im

  FUNCTION maxval_it(it) RESULT(mv)
    INTEGER(I4B), DIMENSION(:,:,:), INTENT(IN) :: it
    INTEGER(I4B) :: mv
    INTEGER(I4B) :: i
    IF ( SIZE(it,3) > 0 ) THEN
       mv = MAXVAL( it(1:SIZE(it,1),1:SIZE(it,2),1) )
    ELSE
       mv = 0.0_SP
    END IF
    IF ( SIZE(it,3) > 1 ) THEN
       DO i=2,SIZE(it,3)
          mv = MAX( mv, MAXVAL( it(1:SIZE(it,1),1:SIZE(it,2),i) ) )
       END DO
    END IF
  END FUNCTION maxval_it
!!$===========================================================================

!!$===========================================================================
  FUNCTION imaxloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_r
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=MAXLOC(arr(:))
    imaxloc_r=iMAX(1)
  END FUNCTION imaxloc_r

  FUNCTION imaxloc_i(iarr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(1) :: imax
    INTEGER(I4B) :: imaxloc_i
    imax=MAXLOC(iarr(:))
    imaxloc_i=iMAX(1)
  END FUNCTION imaxloc_i

  FUNCTION iminloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc_r
    imin=MINLOC(arr(:))
    iminloc_r=iMIN(1)
  END FUNCTION iminloc_r

  FUNCTION iminloc_i(arr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc_i
    imin=MINLOC(arr(:))
    iminloc_i=iMIN(1)
  END FUNCTION iminloc_i
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    IF (.NOT. n1) THEN
       WRITE (*,*) 'error: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    END IF
  END SUBROUTINE assert1

  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    IF (.NOT. (n1 .AND. n2)) THEN
       WRITE (*,*) 'error: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    END IF
  END SUBROUTINE assert2

  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    IF (.NOT. (n1 .AND. n2 .AND. n3)) THEN
       WRITE (*,*) 'error: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    END IF
  END SUBROUTINE assert3

  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    IF (.NOT. (n1 .AND. n2 .AND. n3 .AND. n4)) THEN
       WRITE (*,*) 'error: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    END IF
  END SUBROUTINE assert4

  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    IF (.NOT. ALL(n)) THEN
       WRITE (*,*) 'error: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    END IF
  END SUBROUTINE assert_v

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    IF (n1 == n2) THEN
       assert_eq2=n1
    ELSE
       WRITE (*,*) 'error: an assert_eq failed with this tag:', &
            string
       WRITE(*,*)'n1, n2 = ',n1, n2
       STOP 'program terminated by assert_eq2'
    END IF
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    IF (n1 == n2 .AND. n2 == n3) THEN
       assert_eq3=n1
    ELSE
       WRITE (*,*) 'error: an assert_eq failed with this tag:', &
            string
       WRITE(*,*)'n1, n2, n3 = ',n1, n2, n3
       STOP 'program terminated by assert_eq3'
    END IF
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    IF (n1 == n2 .AND. n2 == n3 .AND. n3 == n4) THEN
       assert_eq4=n1
    ELSE
       WRITE (*,*) 'error: an assert_eq failed with this tag:', &
            string
       WRITE(*,*)'n1, n2, n3, n4 = ',n1, n2, n3, n4
       STOP 'program terminated by assert_eq4'
    END IF
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    IF (ALL(nn(2:) == nn(1))) THEN
       assert_eqn=nn(1)
    ELSE
       WRITE (*,*) 'error: an assert_eq failed with this tag:', &
            string
       WRITE(*,*)'nn ',nn
       STOP 'program terminated by assert_eqn'
    END IF
  END FUNCTION assert_eqn
!!$===========================================================================

!!$===========================================================================
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    IF (n > 0) arth_r(1)=first
    IF (n <= NPAR_ARTH) THEN
       DO k=2,n
          arth_r(k)=arth_r(k-1)+increment
       END DO
    ELSE
       DO k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       END DO
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       DO
          IF (k >= n) EXIT
          k2=k+k
          arth_r(k+1:MIN(k2,n))=temp+arth_r(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
       END DO
    END IF
  END FUNCTION arth_r

  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    IF (n > 0) arth_i(1)=first
    IF (n <= NPAR_ARTH) THEN
       DO k=2,n
          arth_i(k)=arth_i(k-1)+increment
       END DO
    ELSE
       DO k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       END DO
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       DO
          IF (k >= n) EXIT
          k2=k+k
          arth_i(k+1:MIN(k2,n))=temp+arth_i(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
       END DO
    END IF
  END FUNCTION arth_i
!!$===========================================================================

!!$===========================================================================
  FUNCTION arth_DP(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_DP
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    IF (n > 0) arth_DP(1)=first
    IF (n <= NPAR_ARTH) THEN
       DO k=2,n
          arth_DP(k)=arth_DP(k-1)+increment
       END DO
    ELSE
       DO k=2,NPAR2_ARTH
          arth_DP(k)=arth_DP(k-1)+increment
       END DO
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       DO
          IF (k >= n) EXIT
          k2=k+k
          arth_DP(k+1:MIN(k2,n))=temp+arth_DP(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
       END DO
    END IF
  END FUNCTION arth_DP
!!$===========================================================================

!!$===========================================================================
  FUNCTION geop_r(first,factor,n)
    REAL(SP), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: geop_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    IF (n > 0) geop_r(1)=first
    IF (n <= NPAR_GEOP) THEN
       DO k=2,n
          geop_r(k)=geop_r(k-1)*factor
       END DO
    ELSE
       DO k=2,NPAR2_GEOP
          geop_r(k)=geop_r(k-1)*factor
       END DO
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       DO
          IF (k >= n) EXIT
          k2=k+k
          geop_r(k+1:MIN(k2,n))=temp*geop_r(1:MIN(k,n-k))
          temp=temp*temp
          k=k2
       END DO
    END IF
  END FUNCTION geop_r

  FUNCTION geop_i(first,factor,n)
    INTEGER(I4B), INTENT(IN) :: first,factor,n
    INTEGER(I4B), DIMENSION(n) :: geop_i
    INTEGER(I4B) :: k,k2,temp
    IF (n > 0) geop_i(1)=first
    IF (n <= NPAR_GEOP) THEN
       DO k=2,n
          geop_i(k)=geop_i(k-1)*factor
       END DO
    ELSE
       DO k=2,NPAR2_GEOP
          geop_i(k)=geop_i(k-1)*factor
       END DO
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       DO
          IF (k >= n) EXIT
          k2=k+k
          geop_i(k+1:MIN(k2,n))=temp*geop_i(1:MIN(k,n-k))
          temp=temp*temp
          k=k2
       END DO
    END IF
  END FUNCTION geop_i

  FUNCTION geop_c(first,factor,n)
    COMPLEX(SP), INTENT(IN) :: first,factor
    INTEGER(I4B), INTENT(IN) :: n
    COMPLEX(SP), DIMENSION(n) :: geop_c
    INTEGER(I4B) :: k,k2
    COMPLEX(SP) :: temp
    IF (n > 0) geop_c(1)=first
    IF (n <= NPAR_GEOP) THEN
       DO k=2,n
          geop_c(k)=geop_c(k-1)*factor
       END DO
    ELSE
       DO k=2,NPAR2_GEOP
          geop_c(k)=geop_c(k-1)*factor
       END DO
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       DO
          IF (k >= n) EXIT
          k2=k+k
          geop_c(k+1:MIN(k2,n))=temp*geop_c(1:MIN(k,n-k))
          temp=temp*temp
          k=k2
       END DO
    END IF
  END FUNCTION geop_c
!!$===========================================================================

!!$===========================================================================
  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(SIZE(arr)) :: ans
    INTEGER(I4B) :: n,j
    REAL(SP) :: sd
    n=SIZE(arr)
    IF (n == 0_i4b) RETURN
    sd=0.0_SP
    IF (PRESENT(seed)) sd=seed
    ans(1)=arr(1)+sd
    IF (n < NPAR_CUMSUM) THEN
       DO j=2,n
          ans(j)=ans(j-1)+arr(j)
       END DO
    ELSE
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    END IF
  END FUNCTION cumsum_r
  !BL
  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
    INTEGER(I4B), DIMENSION(SIZE(arr)) :: ans
    INTEGER(I4B) :: n,j,sd
    n=SIZE(arr)
    IF (n == 0_i4b) RETURN
    sd=0_i4b
    IF (PRESENT(seed)) sd=seed
    ans(1)=arr(1)+sd
    IF (n < NPAR_CUMSUM) THEN
       DO j=2,n
          ans(j)=ans(j-1)+arr(j)
       END DO
    ELSE
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    END IF
  END FUNCTION cumsum_i
!!$===========================================================================

!!$===========================================================================
!!!****f* UtilitiesMod/cumprod
!!!
!!! NAME
!!!    cumprod
!!!     
!!! USAGE
!!!    cumprod(arr,seed)
!!!    
!!! DESCRIPTION
!!!    Given the rank 1 array arr, returns an array of identical type and 
!!!    size containing the cumulative products of arr.  If the optional 
!!!    argument seed is present, it is multiplied into the first component 
!!!    of the result.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
!!!    REAL(SP), OPTIONAL, INTENT(IN) :: seed
!!!*** 


  RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(SIZE(arr)) :: ans
    INTEGER(I4B) :: n,j
    REAL(SP) :: sd
    n=SIZE(arr)
    IF (n == 0_i4b) RETURN
    sd=1.0_SP
    IF (PRESENT(seed)) sd=seed
    ans(1)=arr(1)*sd
    IF (n < NPAR_CUMPROD) THEN
       DO j=2,n
          ans(j)=ans(j-1)*arr(j)
       END DO
    ELSE
       ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
    END IF
  END FUNCTION cumprod
!!$===========================================================================

!!$===========================================================================
  FUNCTION poly_rr(x,coeffs)
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(SP) :: poly_rr
    REAL(SP) :: pow
    REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=SIZE(coeffs)
    IF (n <= 0) THEN
       poly_rr=0.0_SP
    ELSE IF (n < NPAR_POLY) THEN
       poly_rr=coeffs(n)
       DO i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       END DO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_SP
          nn=ISHFT(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF (nn == 1) EXIT
          pow=pow*pow
          n=nn
       END DO
       poly_rr=vec(1)
       DEALLOCATE(vec)
    END IF
  END FUNCTION poly_rr
!!$===========================================================================

!!$===========================================================================
  FUNCTION poly_rc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_rc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=SIZE(coeffs)
    IF (n <= 0) THEN
       poly_rc=0.0_SP
    ELSE IF (n < NPAR_POLY) THEN
       poly_rc=coeffs(n)
       DO i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       END DO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_SP
          nn=ISHFT(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF (nn == 1) EXIT
          pow=pow*pow
          n=nn
       END DO
       poly_rc=vec(1)
       DEALLOCATE(vec)
    END IF
  END FUNCTION poly_rc
!!$===========================================================================

!!$===========================================================================
  FUNCTION poly_cc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_cc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4B) :: i,n,nn
    n=SIZE(coeffs)
    IF (n <= 0) THEN
       poly_cc=0.0_SP
    ELSE IF (n < NPAR_POLY) THEN
       poly_cc=coeffs(n)
       DO i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       END DO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_SP
          nn=ISHFT(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF (nn == 1) EXIT
          pow=pow*pow
          n=nn
       END DO
       poly_cc=vec(1)
       DEALLOCATE(vec)
    END IF
  END FUNCTION poly_cc
!!$===========================================================================

!!$===========================================================================
  FUNCTION poly_rrv(x,coeffs)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(SP), DIMENSION(SIZE(x)) :: poly_rrv
    INTEGER(I4B) :: i,n,m
    m=SIZE(coeffs)
    n=SIZE(x)
    IF (m <= 0) THEN
       poly_rrv=0.0_SP
    ELSE IF (m < n .OR. m < NPAR_POLY) THEN
       poly_rrv=coeffs(m)
       DO i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       END DO
    ELSE
       DO i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       END DO
    END IF
  END FUNCTION poly_rrv
!!$===========================================================================

!!$===========================================================================
  FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(SIZE(x)) :: poly_msk_rrv
    poly_msk_rrv=UNPACK(poly_rrv(PACK(x,mask),coeffs),mask,0.0_SP)
  END FUNCTION poly_msk_rrv
!!$===========================================================================

!!$===========================================================================
  RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a
    REAL(SP), INTENT(IN) :: b
    REAL(SP), DIMENSION(SIZE(a)) :: u
    INTEGER(I4B) :: n,j
    n=SIZE(a)
    IF (n <= 0) RETURN
    u(1)=a(1)
    IF (n < NPAR_POLYTERM) THEN
       DO j=2,n
          u(j)=a(j)+b*u(j-1)
       END DO
    ELSE
       u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    END IF
  END FUNCTION poly_term_rr
!!$===========================================================================

!!$===========================================================================
  RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(SPC), INTENT(IN) :: b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: u
    INTEGER(I4B) :: n,j
    n=SIZE(a)
    IF (n <= 0) RETURN
    u(1)=a(1)
    IF (n < NPAR_POLYTERM) THEN
       DO j=2,n
          u(j)=a(j)+b*u(j-1)
       END DO
    ELSE
       u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    END IF
  END FUNCTION poly_term_cc
!!$===========================================================================

!!$===========================================================================
!!!****f* UtilitiesMod/zroots_unity
!!!
!!! NAME
!!!    zroots_unity
!!!     
!!! USAGE
!!!    zroots_unity(n,nn)
!!!    
!!! DESCRIPTION
!!!    Returns a complex array containing nn consecutive powers of the
!!!    nth complex root of unity.  
!!! INPUTS
!!!    INTEGER(I4B), INTENT(IN) :: n,nn
!!!*** 

  FUNCTION zroots_unity(n,nn)
    INTEGER(I4B), INTENT(IN) :: n,nn
    COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4B) :: k
    REAL(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWO_PI/n
    k=1
    DO
       IF (k >= nn) EXIT
       zroots_unity(k+1)=CMPLX(COS(k*theta),SIN(k*theta),SPC)
       zroots_unity(k+2:MIN(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:MIN(k,nn-k))
       k=2*k
    END DO
  END FUNCTION zroots_unity
!!$===========================================================================

!!$===========================================================================
  FUNCTION extendedprod_r(a)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a
    REAL(SP) :: extendedprod_r

    INTEGER(I4B) :: i

    extendedprod_r = 1.0_SP
    DO i=1,SIZE(a)
       extendedprod_r = extendedprod_r * a(i)
    END DO
  END FUNCTION extendedprod_r
!!$===========================================================================

!!$===========================================================================
  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_r
    outerprod_r = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outerprod_r
!!!****f* UtilitiesMod/outerdiv
!!!
!!! NAME
!!!    outerdiv
!!!     
!!! USAGE
!!!    outerdiv(a,b)
!!!    
!!! DESCRIPTION
!!!    Returns a matrix that is the outer quotient of two vectors
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
!!!*** 
  FUNCTION outerdiv(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerdiv
    outerdiv = SPREAD(a,dim=2,ncopies=SIZE(b)) / &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outerdiv
!!!****f* UtilitiesMod/outersum
!!!
!!! NAME
!!!    outersum
!!!     
!!! USAGE
!!!    outersum(a,b)
!!!    
!!! DESCRIPTION
!!!    Returns a matrix that is the outer sum of the two vectors  
!!! INPUTS
!!!    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
!!!*** 
  FUNCTION outersum(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outersum
    outersum = SPREAD(a,dim=2,ncopies=SIZE(b)) + &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outersum

  FUNCTION outerdiff_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(SIZE(a),SIZE(b)) :: outerdiff_r
    outerdiff_r = SPREAD(a,dim=2,ncopies=SIZE(b)) - &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outerdiff_r

  FUNCTION outerdiff_i(a,b)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    INTEGER(I4B), DIMENSION(SIZE(a),SIZE(b)) :: outerdiff_i
    outerdiff_i = SPREAD(a,dim=2,ncopies=SIZE(b)) - &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outerdiff_i
!!!****f* UtilitiesMod/outerand
!!!
!!! NAME
!!!    outerand
!!!     
!!! USAGE
!!!    outerand(a,b)
!!!    
!!! DESCRIPTION
!!!    Returns a matrix that is the outer logical of two vectors.
!!! INPUTS
!!!    LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: a,b
!!!*** 
  FUNCTION outerand(a,b)
    LOGICAL(LGD), DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL(LGD), DIMENSION(SIZE(a),SIZE(b)) :: outerand
    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)) .AND. &
         SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION outerand
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE scatter_add_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_add_r')
    m=SIZE(dest)
    DO j=1,n
       i=dest_index(j)
       IF (i > 0 .AND. i <= m) dest(i)=dest(i)+source(j)
    END DO
  END SUBROUTINE scatter_add_r

  SUBROUTINE scatter_max_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4B) :: m,n,j,i
    n=assert_eq2(SIZE(source),SIZE(dest_index),'scatter_max_r')
    m=SIZE(dest)
    DO j=1,n
       i=dest_index(j)
       IF (i > 0 .AND. i <= m) dest(i)=MAX(dest(i),source(j))
    END DO
  END SUBROUTINE scatter_max_r
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE diagadd_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(SIZE(diag),MIN(SIZE(mat,1),SIZE(mat,2)),'diagadd_rv')
    DO j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    END DO
  END SUBROUTINE diagadd_rv

  SUBROUTINE diagadd_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = MIN(SIZE(mat,1),SIZE(mat,2))
    DO j=1,n
       mat(j,j)=mat(j,j)+diag
    END DO
  END SUBROUTINE diagadd_r

  SUBROUTINE diagmult_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(SIZE(diag),MIN(SIZE(mat,1),SIZE(mat,2)),'diagmult_rv')
    DO j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    END DO
  END SUBROUTINE diagmult_rv

  SUBROUTINE diagmult_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = MIN(SIZE(mat,1),SIZE(mat,2))
    DO j=1,n
       mat(j,j)=mat(j,j)*diag
    END DO
  END SUBROUTINE diagmult_r

  FUNCTION get_diag_rv(mat)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(SP), DIMENSION(SIZE(mat,1)) :: get_diag_rv
    INTEGER(I4B) :: j
    j=assert_eq2(SIZE(mat,1),SIZE(mat,2),'get_diag_rv')
    DO j=1,SIZE(mat,1)
       get_diag_rv(j)=mat(j,j)
    END DO
  END FUNCTION get_diag_rv

  SUBROUTINE put_diag_rv(diagv,mat)
    REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4B) :: j,n
    n=assert_eq2(SIZE(diagv),MIN(SIZE(mat,1),SIZE(mat,2)),'put_diag_rv')
    DO j=1,n
       mat(j,j)=diagv(j)
    END DO
  END SUBROUTINE put_diag_rv

  SUBROUTINE put_diag_r(scal,mat)
    REAL(SP), INTENT(IN) :: scal
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4B) :: j,n
    n = MIN(SIZE(mat,1),SIZE(mat,2))
    DO j=1,n
       mat(j,j)=scal
    END DO
  END SUBROUTINE put_diag_r
!!!****f* UtilitiesMod/unit_matrix
!!!
!!! NAME
!!!    unit_matrix
!!!     
!!! USAGE
!!!    unit_matrix(mat)
!!!    
!!! DESCRIPTION
!!!    Sets the diagonal components of mat to unity, all other 
!!!    components are set to zero.
!!! INPUTS
!!!    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
!!!*** 
  SUBROUTINE unit_matrix(mat)
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER(I4B) :: i,n
    n=MIN(SIZE(mat,1),SIZE(mat,2))
    mat(:,:)=0.0_SP
    DO i=1,n
       mat(i,i)=1.0_SP
    END DO
  END SUBROUTINE unit_matrix
!!!****f* UtilitiesMod/upper_triangle
!!!
!!! NAME
!!!    upper_triangle
!!!     
!!! USAGE
!!!    upper_triangle(j,k,extra)
!!!    
!!! DESCRIPTION
!!!    When the optional argument extra is zero or absent, 
!!!    returns a logical mask of shape (j,k) whose values are 
!!!    true above and right of the diagonal, false elsewhere 
!!!    (including the diagonal).  When the extra is present and 
!!!    positive, a corresponding number of additional (sub-)diagonals 
!!!    are returned as true. (extra=1 makes the main diagonal return 
!!!    true.)  When extra is present and negative, it suppresses a 
!!!    corresponding number of superdiagonals.
!!! INPUTS
!!!    INTEGER(I4B), INTENT(IN) :: j,k
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
!!!***
  FUNCTION upper_triangle(j,k,extra)
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGD), DIMENSION(j,k) :: upper_triangle
    INTEGER(I4B) :: n
    n=0
    IF (PRESENT(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle
!!!****f* UtilitiesMod/lower_triangle
!!!
!!! NAME
!!!    lower_triangle
!!!     
!!! USAGE
!!!    lower_triangle(j,k,extra)
!!!    
!!! DESCRIPTION
!!!    When the optional argument extra is zero or absent, 
!!!    returns a logical mask of shape (j,k) whose values are 
!!!    true below and left of the diagonal, false elsewhere 
!!!    (including the diagonal).  When the extra is present and 
!!!    positive, a corresponding number of additional (super-)diagonals 
!!!    are returned as true. (extra=1 makes the main diagonal return 
!!!    true.)  When extra is present and negative, it suppresses a 
!!!    corresponding number of subdiagonals.
!!! INPUTS
!!!    INTEGER(I4B), INTENT(IN) :: j,k
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
!!!***
  FUNCTION lower_triangle(j,k,extra)
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGD), DIMENSION(j,k) :: lower_triangle
    INTEGER(I4B) :: n
    n=0
    IF (PRESENT(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  END FUNCTION lower_triangle
!!$===========================================================================

!!$===========================================================================
!!!****f* UtilitiesMod/vABS
!!!
!!! NAME
!!!    vABS
!!!     
!!! USAGE
!!!    vABS(v)
!!!    
!!! DESCRIPTION
!!!    Returns the length of the vector v.
!!! INPUTS
!!!        REAL(SP), DIMENSION(:), INTENT(IN) :: v
!!!***
  FUNCTION vABS(v)
    REAL(SP), DIMENSION(:), INTENT(IN) :: v
    REAL(SP) :: vabs
    vabs=SQRT(DOT_PRODUCT(v,v))
  END FUNCTION vABS
!!$===========================================================================

!!!****f* UtilitiesMod/optional_flag
!!!
!!! NAME
!!!      optional_flag -- function for processing optional logical flags
!!! USAGE
!!!      optional_flag(FLAG,DEFAULT)
!!! DESCRIPTION
!!!      Often you will want to do flow control based on the value of an
!!!      optional logical flag.  This can be tedious because one must usually
!!!      check the variable with PRESENT, and use a dummy local logical
!!!      to hold the value.  This function does away with those tedious,
!!!      hard to read, multi-line solutions.  Just use it like:
!!!
!!!      if ( OPTIONAL_FLAG(FLAG,DEFAULT) ) then
!!!          blah, blah, blah
!!!      else
!!!          and so on.
!!!      
!!!      If the flag is present, optional_flag takes that value.  If not, it
!!!      takes on the value of DEFAULT.  If DEFAULT is not present, 
!!!      optional_flag is false.
!!! INPUTS
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: FLAG
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DEFAULT 
!!! OUTPUTS
!!!    LOGICAL(LGD) :: OPT_FLAG
!!!***

PURE FUNCTION OPTIONAL_FLAG(FLAG,DEFAULT) RESULT(OPT_FLAG)
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: FLAG
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DEFAULT
    LOGICAL(LGD) :: OPT_FLAG

    IF ( PRESENT(FLAG) ) THEN
       OPT_FLAG = FLAG
    ELSE
       IF ( PRESENT(DEFAULT) ) THEN
          OPT_FLAG = DEFAULT
       ELSE
          OPT_FLAG = .FALSE.
       END IF
    END IF

  END FUNCTION OPTIONAL_FLAG
!!$===========================================================================

!!$===========================================================================
  FUNCTION SUM_3ind_MATRIX(A) RESULT(SUM)
    REAL(SP), INTENT(IN) :: A(:,:,:)
    REAL(SP) :: SUM

    INTEGER(I4B) :: i,j,k

    SUM = 0.0_SP
    DO i=1,SIZE(A,3)
       DO j=1,SIZE(A,2)
          DO k=1,SIZE(A,1)
             SUM = SUM + A(k,j,i)
          END DO
       END DO
    END DO

  END FUNCTION SUM_3ind_MATRIX
!!$===========================================================================

!!$===========================================================================
  FUNCTION SUM_2ind_MATRIX(A) RESULT(SUM)
    REAL(SP), INTENT(IN) :: A(:,:)
    REAL(SP) :: SUM

    INTEGER(I4B) :: i,j

    SUM = 0.0_SP
    DO i=1,SIZE(A,2)
       DO j=1,SIZE(A,1)
          SUM = SUM + A(j,i)
       END DO
    END DO

  END FUNCTION SUM_2ind_MATRIX
!!$===========================================================================

!!$===========================================================================
  FUNCTION ALMOST_EQUAL_S2(A,B) RESULT(EQUAL)
    REAL(SP), INTENT(IN) :: A, B
    LOGICAL(LGD) :: EQUAL

    IF ( ABS(A-B) < TINY_FLOAT ) THEN
       EQUAL = .TRUE.
    ELSE
       EQUAL = .FALSE.
    END IF

  END FUNCTION ALMOST_EQUAL_S2

  FUNCTION ALMOST_EQUAL_S3(A,B,C) RESULT(EQUAL)
    REAL(SP), INTENT(IN) :: A, B, C
    LOGICAL(LGD) :: EQUAL

    IF ( ABS(A-B) < TINY_FLOAT .AND. ABS(B-C) < TINY_FLOAT ) THEN
       EQUAL = .TRUE.
    ELSE
       EQUAL = .FALSE.
    END IF

  END FUNCTION ALMOST_EQUAL_S3

  FUNCTION ALMOST_EQUAL_V1(A) RESULT(EQUAL)
    REAL(SP), INTENT(IN) :: A(:)
    LOGICAL(LGD) :: EQUAL

    IF ( ABS(MAXVAL(A)-MINVAL(A)) < TINY_FLOAT ) THEN
       EQUAL = .TRUE.
    ELSE
       EQUAL = .FALSE.
    END IF

  END FUNCTION ALMOST_EQUAL_V1

  FUNCTION ALMOST_EQUAL_V2(A,B) RESULT(EQUAL)
    REAL(SP), INTENT(IN) :: A(:),B(:)
    LOGICAL(LGD) :: EQUAL

    IF ( ALMOST_EQUAL(A-B) .AND. ABS(MAXVAL(A-B)) < TINY_FLOAT ) THEN
       EQUAL = .TRUE.
    ELSE
       EQUAL = .FALSE.
    END IF

  END FUNCTION ALMOST_EQUAL_V2
!!$===========================================================================

!!$===========================================================================
  FUNCTION IndexSearch_string(avec,a,case_sensitive) RESULT(ind)
    CHARACTER(LEN=*), INTENT(IN) :: avec(:)
    CHARACTER(LEN=*), INTENT(IN) :: a
    LOGICAL(LGD), OPTIONAL :: case_sensitive

    INTEGER(I4B) :: ind

    INTEGER(I4B) :: i

    ind = 0

    DO i = 1, SIZE(avec)
       IF (.NOT. OPTIONAL_FLAG(case_sensitive,DEFAULT=.TRUE.)) THEN
          IF ( LOWERCASE(a) == LOWERCASE(avec(i)) ) THEN
             ind = i
             EXIT
          END IF
       ELSE
          IF ( a == avec(i) ) THEN
             ind = i
             EXIT
          END IF
       END IF
    END DO

  END FUNCTION IndexSearch_string

  FUNCTION IndexSearch_iv(avec,a) RESULT(ind)
    INTEGER(I4B), INTENT(IN) :: avec(:)
    INTEGER(I4B), INTENT(IN) :: a

    INTEGER(I4B) :: ind

    INTEGER(I4B) :: i

    ind = 0

    DO i = 1, SIZE(avec)
       IF ( a == avec(i) ) THEN
          ind = i
          EXIT
       END IF
    END DO

  END FUNCTION IndexSearch_iv

  FUNCTION IndexSearch_rv(avec,a) RESULT(ind)
    REAL(SP), INTENT(IN) :: avec(:)
    REAL(SP), INTENT(IN) :: a

    INTEGER(I4B) :: ind

    INTEGER(I4B) :: i

    ind = 0

    DO i = 1, SIZE(avec)
       IF ( a == avec(i) ) THEN
          ind = i
          EXIT
       END IF
    END DO

  END FUNCTION IndexSearch_rv

  FUNCTION IndexSearch_cv(avec,a) RESULT(ind)
    COMPLEX(SPC), INTENT(IN) :: avec(:)
    COMPLEX(SPC), INTENT(IN) :: a

    INTEGER(I4B) :: ind

    INTEGER(I4B) :: i

    ind = 0

    DO i = 1, SIZE(avec)
       IF ( a == avec(i) ) THEN
          ind = i
          EXIT
       END IF
    END DO

  END FUNCTION IndexSearch_cv
!!$===========================================================================


!!!****f* UtilitiesMod/IsOdd
!!!
!!! NAME
!!!    IsOdd - a function
!!!     
!!! USAGE
!!!    IsOdd(i)
!!!    
!!! DESCRIPTION
!!!    Function returns true if integer is odd.
!!!
!!! INPUTS
!!!    integer(i4b) :: i
!!!    
!!! OUTPUTS
!!!    logical(lgd) :: IsOdd
!!!    
!!!*** 

  logical(lgd) function IsOdd(i)

    integer(i4b), intent(in) :: i

    if ( mod(i,2) /= 0 ) then
       IsOdd = .true.
    else
       IsOdd = .false.
    end if

  end function IsOdd

!!!****f* UtilitiesMod/IsEven
!!!
!!! NAME
!!!    IsEven - a function
!!!     
!!! USAGE
!!!    IsEven(i)
!!!    
!!! DESCRIPTION
!!!    Function returns true if integer is even.
!!!
!!! INPUTS
!!!    integer(i4b) :: i
!!!    
!!! OUTPUTS
!!!    logical(lgd) :: IsEven
!!!    
!!!*** 

  logical(lgd) function IsEven(i)

    integer(i4b), intent(in) :: i

    if ( mod(i,2) == 0 ) then
       IsEven = .true.
    else
       IsEven = .false.
    end if

  end function IsEven

  FUNCTION LOWERCASE(string) RESULT(lstring)
!!!****f* UtilitiesMod/LOWERCASE
!!!
!!! NAME
!!!     LOWERCASE -- returns a lowercased string
!!! USAGE
!!!     LOWERCASE(string)
!!! DESCRIPTION
!!!     Returns a string with all of the letters lowercased
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: lstring
!!! NOTES
!!!     Should be encoding independent, if the platform encoding is Fortran90
!!!     compliant.
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=LEN(string)) :: lstring

    INTEGER(I4B) :: iau, ial, izu, iac, ishift, i

    iau = IACHAR('A'); ial = IACHAR('a'); izu = IACHAR('Z')
    ishift = ial - iau

    lstring = string

    DO i=1,LEN(string)
       iac = IACHAR(string(i:i))
       IF ( iac >= iau .AND. iac <= izu ) THEN
          lstring(i:i) = ACHAR(iac+ishift)
       END IF
    END DO

  END FUNCTION LOWERCASE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION UPPERCASE(string) RESULT(ustring)
!!!****f* UtilitiesMod/UPPERCASE
!!!
!!! NAME
!!!     UPPERCASE -- returns a uppercased string
!!! USAGE
!!!     UPPERCASE(string)
!!! DESCRIPTION
!!!     Returns a string with all of the letters uppercased
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: lstring
!!! NOTES
!!!     Should be encoding independent, if the platform encoding is Fortran90
!!!     compliant.
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=LEN(string)) :: ustring

    INTEGER(I4B) :: iau, ial, izl, iac, ishift, i

    iau = IACHAR('A'); ial = IACHAR('a'); izl = IACHAR('z')
    ishift = ial - iau

    ustring = string

    DO i=1,LEN(string)
       iac = IACHAR(string(i:i))
       IF ( iac >= ial .AND. iac <= izl ) THEN
          ustring(i:i) = ACHAR(iac-ishift)
       END IF
    END DO

  END FUNCTION UPPERCASE


!!!****f* UtilitiesMod/push
!!!
!!! NAME
!!!    push
!!!     
!!! USAGE
!!!    CALL push(ptr(:,:),val(:))
!!!    or
!!!    CALL push(ptr(:),val)
!!!    
!!! DESCRIPTION
!!!    reallocates the last index of ptr by 1
!!!    Stores val in the last index of the array.
!!!
!!! INPUTS
!!!    ptr and val must be the same type
!!!    They can be I4B, SP, LGD, or character
!!!  USES
!!!    YorkLib/DataTypes
!!!    UtilitiesMod/reallocate
!!!    UtilitiesMod/assert_eq
!!!    
!!!*** 

  SUBROUTINE Push_real_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    REAL(SP),POINTER :: ptr(:)
    REAL(SP),INTENT(IN) :: val
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       ALLOCATE( ptr(1) )
       ptr(1) = val
    ELSE
       ptr => reallocate( ptr, SIZE(ptr)+1 )
       ptr(SIZE(ptr)) = val
    END IF
  END SUBROUTINE Push_real_v

  SUBROUTINE Push_real_m(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate,assert_eq
    REAL(SP),POINTER :: ptr(:,:)
    REAL(SP),INTENT(IN) :: val(:)
    INTEGER(I4B) :: n
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       ALLOCATE( ptr(SIZE(val),1) )
       ptr(:,1) = val(:)
    ELSE
       n = assert_eq(SIZE(val),SIZE(ptr,1)," SIZE(val),SIZE(ptr,2) are not equal")
       ptr => reallocate( ptr, n , SIZE(ptr,2)+1)
       ptr(:,SIZE(ptr,2)) = val(:)
    END IF
  END SUBROUTINE Push_real_m

  SUBROUTINE Push_int_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    INTEGER(I4B),POINTER :: ptr(:)
    INTEGER(I4B),INTENT(IN) :: val
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       ALLOCATE( ptr(1) )
       ptr(1) = val
    ELSE
       ptr => reallocate( ptr, SIZE(ptr)+1 )
       ptr(SIZE(ptr)) = val
    END IF
  END SUBROUTINE Push_int_v

  SUBROUTINE Push_logic_v(ptr,val)
    USE DataTypes
    !    USE UtilitiesMod, ONLY : reallocate
    LOGICAL(LGD),POINTER    :: ptr(:)
    LOGICAL(LGD),INTENT(IN) :: val
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       ALLOCATE( ptr(1) )
       ptr(1) = val
    ELSE
       ptr => reallocate( ptr, SIZE(ptr)+1 )
       ptr(SIZE(ptr)) = val
    END IF
  END SUBROUTINE Push_logic_v

  SUBROUTINE Push_int_m(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate,assert_eq
    INTEGER(I4B),POINTER :: ptr(:,:)
    INTEGER(I4B),INTENT(IN) :: val(:)
    INTEGER(I4B) :: n
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       ALLOCATE( ptr(SIZE(val),1) )
       ptr(:,1) = val(:)
    ELSE
       n = assert_eq(SIZE(val),SIZE(ptr,1)," SIZE(val),SIZE(ptr,2) are not equal")
       ptr => reallocate( ptr, n , SIZE(ptr,2)+1)
       ptr(:,SIZE(ptr,2)) = val(:)
    END IF
  END SUBROUTINE Push_int_m
  
  SUBROUTINE Push_char_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
!    INTEGER(I4B),INTENT(IN) :: length
    CHARACTER(LEN=*),POINTER :: ptr(:)
    CHARACTER(LEN=*),INTENT(IN) :: val
 !   WRITE(6,*)"IN Push_char_v"
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
!       WRITE(6,*)"ALLOCATE"
       ALLOCATE( ptr(1) )
 !      WRITE(6,*)"PTR1"
       ptr(1) = val(1:MIN(LEN(ptr(1)),LEN(val)))
!       WRITE(6,*)"RETURN"
    ELSE
!       WRITE(6,*)"REALLOCATE"
!       WRITE(6,*)"SIZE = ",SIZE(ptr)
       ptr => reallocate( ptr, SIZE(ptr)+1,LEN=LEN(ptr(1)) )
!       WRITE(6,*)"PTR(SIZE)"
       ptr(SIZE(ptr)) = val(1:MIN(LEN(ptr),LEN(val)))
!       WRITE(6,*)"RETURN"
    END IF
  END SUBROUTINE Push_char_v



!!!****f* UtilitiesMod/shift
!!!
!!! NAME
!!!    shift
!!!     
!!! USAGE
!!!    CALL shift(ptr(:),val)
!!!    
!!! DESCRIPTION
!!!    reallocates the last index of ptr by -1
!!!    Stores the first index of the array in val.
!!!
!!! INPUTS
!!!    ptr and val must be the same type
!!!    They can be I4B, SP, LGD, or character
!!!  USES
!!!    YorkLib/DataTypes
!!!    UtilitiesMod/reallocate
!!!    UtilitiesMod/assert_eq
!!!    
!!!*** 
  SUBROUTINE Shift_real_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    REAL(SP),POINTER :: ptr(:)
    REAL(SP),INTENT(OUT) :: val
    REAL(SP),POINTER :: tmp(:) => NULL()
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       RETURN
    ELSE
       IF ( SIZE(ptr) < 1 ) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
          RETURN
       END IF
       val = ptr(1)
       IF ( SIZE(ptr) == 1 ) THEN
          DEALLOCATE(ptr)
          RETURN
       END IF
       ALLOCATE(tmp(SIZE(ptr)-1))
       tmp = ptr(2:SIZE(ptr))
!       ptr => reallocate( ptr, SIZE(ptr)-1 )
       DEALLOCATE(ptr)
       ALLOCATE(ptr(SIZE(tmp)))
       ptr = tmp
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE Shift_real_v


  SUBROUTINE Shift_int_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    INTEGER(I4B),POINTER     :: ptr(:)
    INTEGER(I4B),INTENT(OUT) :: val
    INTEGER(I4B),POINTER     :: tmp(:) => NULL()
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       RETURN
    ELSE
       IF ( SIZE(ptr) < 1 ) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
          RETURN
       END IF
       val = ptr(1)
       IF ( SIZE(ptr) == 1 ) THEN
          DEALLOCATE(ptr)
          RETURN
       END IF
       ALLOCATE(tmp(SIZE(ptr)-1))
       tmp = ptr(2:SIZE(ptr))
!       ptr => reallocate( ptr, SIZE(ptr)-1 )
       DEALLOCATE(ptr)
       ALLOCATE(ptr(SIZE(tmp)))
       ptr = tmp
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE Shift_int_v


  SUBROUTINE Shift_logic_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    LOGICAL(LGD),POINTER     :: ptr(:)
    LOGICAL(LGD),INTENT(OUT) :: val
    LOGICAL(LGD),POINTER     :: tmp(:) => NULL()
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       RETURN
    ELSE
       IF ( SIZE(ptr) < 1 ) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
          RETURN
       END IF
       val = ptr(1)
       IF ( SIZE(ptr) == 1 ) THEN
          DEALLOCATE(ptr)
          RETURN
       END IF
       ALLOCATE(tmp(SIZE(ptr)-1))
       tmp = ptr(2:SIZE(ptr))
!       ptr => reallocate( ptr, SIZE(ptr)-1 )
       DEALLOCATE(ptr)
       ALLOCATE(ptr(SIZE(tmp)))
       ptr = tmp
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE Shift_logic_v


  SUBROUTINE Shift_char_v(ptr,val)
    USE DataTypes
!    USE UtilitiesMod, ONLY : reallocate
    CHARACTER(LEN=*),POINTER        :: ptr(:)
    CHARACTER(LEN=*),INTENT(OUT)    :: val
    CHARACTER(LEN=LEN(ptr)),POINTER :: tmp(:)
!    INTEGER(I4B) :: i,n
!    tmp => NULL()
    NULLIFY(tmp)
!    WRITE(6,*)"----ENTERED SHIFT----"
    IF ( .NOT. ASSOCIATED( ptr ) ) THEN
       RETURN
    ELSE
       IF ( SIZE(ptr) < 1 ) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
          RETURN
       END IF
       IF ( LEN_TRIM(ptr(1)) < LEN(val) ) THEN
          val = TRIM(ptr(1))
       ELSE
          val = ptr(1)(1:LEN(val))
       END IF
!       WRITE(6,*)"PTR=",TRIM(ptr(1))
!       WRITE(6,*)"VAL=",TRIM(val)
!       DO i=1,(SIZE(ptr))
!          WRITE(6,'(I3,2A)')i," PTR =",TRIM(ptr(i))
!       END DO

       IF ( SIZE(ptr) == 1 ) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
          RETURN
       END IF
       ALLOCATE(tmp(SIZE(ptr)-1))
       tmp = ptr(2:SIZE(ptr))
!       DO i=1,(SIZE(tmp))
!          WRITE(6,'(I3,2A)')i," TMP =",TRIM(tmp(i))
!       END DO
!       ptr => reallocate( ptr, SIZE(ptr)-1 )
       DEALLOCATE(ptr)
       ALLOCATE(ptr(SIZE(tmp)))
       ptr = tmp
!       DO i=1,(SIZE(ptr))
!          WRITE(6,'(I3,2A)')i," PTR2=",TRIM(ptr(i))
!       END DO
       DEALLOCATE(tmp)
       NULLIFY(tmp)
    END IF
!    WRITE(6,*)"EXITING SHIFT"
  END SUBROUTINE Shift_char_v






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine estimate_real_size_sq(n,THRESH)

    ! Estimates the storage requirements of nxn matrix of 
    ! 8 byte (64 bit) reals and print a message if the storage is estimated to
    ! be greater than THRESH

    integer(i4b), intent(in) :: n
    real(sp), optional, intent(in) :: THRESH

    real(sp) :: t = 100.0
    real(sp) :: mb

    if ( present(THRESH) ) then
       t = THRESH
    end if

    mb = n**2 * 8.0_sp / (2**20)
    if ( mb > t ) then
       write(6,'(A,F0.0,A)')&
               'Estimated size of matrix: ',mb,' M'
    end if
 
  end subroutine estimate_real_size_sq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine estimate_real_size_rect(n,m,THRESH)

    ! Estimates the storage requirements of nxn matrix of 
    ! 8 byte (64 bit) reals and print a message if the storage is estimated to
    ! be greater than THRESH

    integer(i4b), intent(in) :: n,m
    real(sp), optional, intent(in) :: THRESH

    real(sp) :: t = 100.0
    real(sp) :: mb

    if ( present(THRESH) ) then
       t = THRESH
    end if

    mb = n*m * 8.0_sp / (2**20)
    if ( mb > t ) then
       write(6,'(A,F0.0,A)')&
               'Estimated size of matrix: ',mb,' M'
    end if
 
  end subroutine estimate_real_size_rect


  SUBROUTINE Spherical2Cartesian(r,theta,phi,x,y,z,DEG)
!!!****f* UtilitiesMod/Spherical2Cartesian
!!!
!!! NAME
!!!     Spherical2Cartesian
!!! USAGE
!!!     CALL Spherical2Cartesian(r,theta,phi,x,y,z,DEG)
!!! DESCRIPTION
!!!     Given r, theta, and phi
!!!     (theta has a range of pi and phi has a range of 2*pi),
!!!     this routine converts the spherical coordinates
!!!     to the cartesian coordinates x,y,z.
!!!     If DEG is present and true, then degrees are assumed
!!!     (theta has a range of 180 and phi has a range of 360).
!!!
!!! INPUTS
!!!     REAL(SP),INTENT(IN) :: r,theta,phi
!!!     logical(lgd),INTENT(IN),OPTIONAL :: DEG
!!! OUTPUTS
!!!     REAL(SP),INTENT(OUT) :: x,y,z
!!! USES
!!!     MathMod, ONLY : factrl,AssociatedLegendrePolynomial
!!!     USE ConstantsMod, ONLY : FOUR_PI,SQRT2
!!!
!!!***
    USE DataTypes
    USE ConstantsMod, ONLY : RAD2DEG
    ! theta goes from 0 to PI
    ! phi goes from -PI to PI
    REAL(SP),INTENT(OUT) :: x,y,z
    REAL(SP),INTENT(IN) :: r,theta,phi
    LOGICAL(LGD),INTENT(IN),OPTIONAL :: DEG
    REAL(SP) :: th,ph

    IF ( OPTIONAL_FLAG( DEG ) ) THEN
       th = theta/RAD2DEG
       ph = phi/RAD2DEG
    ELSE
       th = theta
       ph = phi
    END IF
    x = r*COS(ph)*SIN(th)
    y = r*SIN(ph)*SIN(th)
    z = r*COS(th)
  END SUBROUTINE Spherical2Cartesian

  SUBROUTINE Cartesian2Spherical(x,y,z,r,theta,phi,DEG)
!!!****f* UtilitiesMod/Cartesian2Spherical
!!!
!!! NAME
!!!     Cartesian2Spherical
!!! USAGE
!!!     CALL Cartesian2Spherical(x,y,z,r,theta,phi,DEG)
!!! DESCRIPTION
!!!     Given cartesian coordinates x,y,z, this routine
!!!     converts the point to spherical coordinates
!!!     r, theta, phi, where the ranges of theta and phi
!!!     are
!!!     theta = (0,PI) and phi = (-PI,+PI)
!!!                (If radians are used)
!!!
!!!     theta = (0,180) and phi = (-180,+180)
!!!                (If degrees are used)
!!!
!!!     Radians are used by default.  If DEG is present
!!!     and true, then angles are output in degrees.
!!!
!!! INPUTS
!!!     REAL(SP),INTENT(IN) :: x,y,z
!!!     logical(lgd),INTENT(IN),OPTIONAL :: DEG
!!! OUTPUTS
!!!     REAL(SP),INTENT(OUT) :: r,theta,phi
!!! USES
!!!     MathMod, ONLY : factrl,AssociatedLegendrePolynomial
!!!     USE ConstantsMod, ONLY : FOUR_PI,SQRT2
!!!
!!!***
    USE DataTypes
    USE ConstantsMod, ONLY : RAD2DEG
    ! theta goes from 0 to PI
    ! phi goes from -PI to PI
    REAL(SP),INTENT(IN) :: x,y,z
    REAL(SP),INTENT(OUT) :: r,theta,phi
    LOGICAL(LGD),INTENT(IN),OPTIONAL :: DEG
    r = SQRT(x*x+y*y+z*z)
    IF ( r .EQ. 0.0_SP ) THEN
       theta = 0.0_SP
       phi = 0.0_SP
       RETURN
    END IF
    IF ( x == 0 .AND. y == 0 ) THEN
       phi = 0.0_SP
    ELSE
       phi = ATAN2(y,x)
    END IF
    theta = ACOS(z/r)
    IF ( OPTIONAL_FLAG(DEG) ) THEN
       phi   = RAD2DEG * phi
       theta = RAD2DEG * theta
    END IF
  END SUBROUTINE Cartesian2Spherical




END MODULE UtilitiesMod




MODULE StringMod
!!!****h* YorkLib/StringMod
!!!
!!! NAME
!!!     StringMod - string functions and parameter constants
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Darrin York
!!! CREATION DATE
!!!     1997
!!! DESCRIPTION
!!!     string functions and parameter constants
!!! USES
!!!     DataTypes, ErrorMod, UtilitiesMod
!!!***

  USE DataTypes
  USE ErrorMod, ONLY : error
  USE UtilitiesMod, ONLY: OPTIONAL_FLAG, assert_eq,UPPERCASE,LOWERCASE
  IMPLICIT NONE

  PUBLIC

!!!****v* StringMod/PublicVariables
!!!
!!! DESCRIPTION
!!!     Character parameters
!!! NOTES
!!!     The value of these paramters is dependent on platform encoding
!!!     These values are valid for ASCII encoding
!!! SOURCE
!!!
  CHARACTER(LEN=1), PARAMETER :: LTICKC        = ACHAR(39)  ! '
  CHARACTER(LEN=1), PARAMETER :: RTICKC        = ACHAR(96)  ! `
  CHARACTER(LEN=1), PARAMETER :: DBLQUOTEC     = ACHAR(34)  ! "
  CHARACTER(LEN=1), PARAMETER :: POUNDC        = ACHAR(35)  ! #
  CHARACTER(LEN=1), PARAMETER :: DOLLARC       = ACHAR(36)  ! $
  CHARACTER(LEN=1), PARAMETER :: PERCENTC      = ACHAR(37)  ! %
  CHARACTER(LEN=1), PARAMETER :: AMPERSANDC    = ACHAR(38)  ! &
  CHARACTER(LEN=1), PARAMETER :: STARC         = ACHAR(42)  ! *
  CHARACTER(LEN=1), PARAMETER :: BACKSLASHC    = ACHAR(47)  ! /
  CHARACTER(LEN=1), PARAMETER :: FORWARDSLASHC = ACHAR(92)  ! \
  CHARACTER(LEN=1), PARAMETER :: ATSYMBLC      = ACHAR(64)  ! @
  CHARACTER(LEN=1), PARAMETER :: BARC          = ACHAR(124) ! |
  CHARACTER(LEN=1), PARAMETER :: TILDC         = ACHAR(126) ! ~

  CHARACTER(LEN=2), PARAMETER :: DBLEPOUND     = POUNDC//POUNDC  ! ##

  ! Private default/setup values
  CHARACTER(LEN=5), PARAMETER, PRIVATE :: SEPERATOR = ' ,;:='

  TYPE, PRIVATE :: ReadLineIO_t
     INTEGER(I4B) :: UNIT    = 5
     LOGICAL(LGD) :: ECHO    = .FALSE.
     INTEGER(I4B) :: OUTUNIT = 6
     LOGICAL(LGD) :: ExitFlag = .FALSE.
     CHARACTER(LEN=120) :: EXITSTRING = ' '
     CHARACTER(LEN=120) :: TAG = ' '
     LOGICAL(LGD) :: PROCESS = .FALSE.
     CHARACTER(LEN=20) :: SEP = SEPERATOR
     LOGICAL(LGD) :: COMMENT = .FALSE.
     CHARACTER(LEN=20) :: COMMENTSYM = '!' 
     LOGICAL(LGD) :: CASESENSITIVE = .TRUE.
     LOGICAL(LGD) :: DIE = .TRUE.
  END TYPE ReadLineIO_t

  TYPE(ReadLineIO_t), PRIVATE, SAVE :: ReadLineIO_Default, ReadLineIO_Setup


!!!***

!!!****f* StringMod/ReadField
!!!
!!! NAME
!!!     ReadField -- read a field from a string
!!! USAGE
!!!     ReadField(string,a,SEP,KEYWORD,SKIP,ADVANCE,IOSTAT, &
!!!       & CASESENSITIVE,TRUNC,WILDCARD,DEFAULT,SUCCESS))
!!! DESCRIPTION
!!!     Reads a field from string into a.
!!!     You can specify a field separator with SEP.
!!!     You can read the field following a KEYWORD using KEYWORD.
!!!     You can skip SKIP fields in string.
!!!     You can advance along the string by setting ADVANCE to the number of
!!!     fields you want to ADVANCE.  Basically this is like SKIPping, but after
!!!     you read the string.
!!!     The OPTIONAL argument KEYWORD allows to read a field
!!!     immediately following the KEYWORD.
!!!     In KEYWORD searches, you can toggle case-sensitivity by 
!!!     setting the optional logical argument CASESENSITIVE (default .TRUE.).
!!!     If you want the KEYWORD comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of KEYWORD that is passed - set TRUNC=.TRUE. (default .TRUE.,
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the refstring.
!!!     You can process your own errors by including IOSTAT.  If you leave out
!!!     IOSTAT, ReadField will abort on errors.
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .TRUE.).
!!!     If you want the comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of "refstring" - set TRUNC=.TRUE. (default .TRUE.,
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the refstring.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(INOUT) :: string
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: KEYWORD
!!!     INTEGER(I4B), OPTIONAL, INTENT(IN) :: ADVANCE
!!!     INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP
!!! OUTPUTS
!!!     CHARACTER(LEN=*), INTENT(INOUT) :: a
!!!                 or
!!!     INTEGER(I4B), INTENT(INOUT) :: a
!!!                 or
!!!     REAL(SP), INTENT(INOUT) :: a
!!!     INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
!!! SOURCE

  INTERFACE ReadField
     MODULE PROCEDURE ReadField_h, &
          & ReadField_i, &
          & ReadField_r
  END INTERFACE

  INTERFACE LookupString
     MODULE PROCEDURE LookupString1, LookupString2, LookupString3
  END INTERFACE
!!!***
!!!****f* StringMod/StringIsaNumber
!!!
!!! NAME
!!!     StringIsaNumber -- is this string a number?
!!! USAGE
!!!     isNum = StringIsaNumber(string)
!!!     isNum = StringIsaNumber(string,real)
!!!     isNum = StringIsaNumber(string,int)
!!! DESCRIPTION
!!!     Uses sleazy trick of using a character string 
!!!     for fortran read. Thus, if fortran can read
!!!     the string into variable [type REAL(SP) or INTEGER(I4B)]
!!!     the string is a valid fortran number.
!!!     On exit, variable will contain numerical value of string
!!!     and result will be true (if string was a number)
!!!     else result is false.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     
!!! OUTPUT
!!!     LOGICAL(LGD) :: IsItaNumber
!!!     INTEGER(I4B), INTENT(OUT), variable
!!!       or
!!!     REAL(SP), INTENT(OUT), variable  
!!!***
  INTERFACE StringIsaNumber
     MODULE PROCEDURE StringIsaNumber_real
     MODULE PROCEDURE StringIsaNumber_int
     MODULE PROCEDURE StringIsaNumber_any
  END INTERFACE


CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION IntegerToString(i,LEN) RESULT(string)
!!!****f* StringMod/IntegerToString
!!!
!!! NAME
!!!     IntegerToString -- converts an integer to a string with padded 0's.
!!! USAGE
!!!     IntegerToString(i,LEN)
!!! INPUTS
!!!     INTEGER(I4B), INTENT(IN) :: i
!!!     INTEGER(I4B), INTENT(IN) :: LEN
!!! OUTPUT
!!!     CHARACTER(LEN=LEN) :: string
!!!***
    INTEGER(I4B), INTENT(IN) :: i
    INTEGER(I4B), INTENT(IN) :: LEN
    CHARACTER(LEN=LEN) :: string

    INTEGER(I4B) :: ILEN
    CHARACTER(LEN=80) :: FmtString

    IF ( i < 0 ) THEN
       ILEN = LEN-1
    ELSE
       ILEN = LEN
    END IF

    FmtString = ' '

    IF ( LEN < 10 ) THEN
       WRITE(FmtString,10)LEN,ILEN
    ELSE IF ( LEN >= 10 .AND. ILEN < 10) THEN
       WRITE(FmtString,15)LEN,ILEN
    ELSE
       WRITE(FmtString,20)LEN,ILEN
    END IF

10  FORMAT('(I',i1,'.',i1,')')
15  FORMAT('(I',i2,'.',i1,')')
20  FORMAT('(I',i2,'.',i2,')')

    WRITE(string,FmtString)i

  END FUNCTION IntegerToString
!!$===========================================================================

!!$===========================================================================
  FUNCTION StringIsaNumber_any(string) RESULT(IsItaNumber)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL(LGD) :: IsItaNumber
    REAL(SP) :: temp_r
    INTEGER(I4B) :: temp_i

    IsItaNumber = StringIsaNumber_real(string,temp_r) .OR. &
         & StringIsaNumber_int(string,temp_i) 

  END FUNCTION StringIsaNumber_any

  FUNCTION StringIsaNumber_real(string,variable) RESULT(IsItaNumber)
    CHARACTER(LEN=*), INTENT(IN) :: string
    REAL(SP), INTENT(OUT) :: variable
    LOGICAL(LGD) :: IsItaNumber
    INTEGER(I4B) :: status

    variable=0.0_SP
    READ(string,*,IOSTAT=status)variable
    IF (status /= 0) THEN
       IsItaNumber=.FALSE.
    ELSE
       IsItaNumber=.TRUE.
    ENDIF

  END FUNCTION StringIsaNumber_real

  FUNCTION StringIsaNumber_int(string,variable) RESULT(IsItaNumber)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER(I4B), INTENT(OUT) :: variable
    LOGICAL(LGD) :: IsItaNumber
    INTEGER(I4B) :: status

    variable=0_I4B
    READ(string,*,IOSTAT=status)variable
    IF (status /= 0) THEN
       IsItaNumber=.FALSE.
    ELSE
       IsItaNumber=.TRUE.
    ENDIF

  END FUNCTION StringIsaNumber_int
!!$===========================================================================

!!$===========================================================================
  FUNCTION StringIsNotaNumber(string) RESULT(IsItaNumber)
!!!****f* StringMod/StringIsNotaNumber
!!!
!!! NAME
!!!     StringIsNotaNumber -- is this string NOT a number?
!!! USAGE
!!!     StringIsNotaNumber(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     LOGICAL(LGD) :: IsItaNumber
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL(LGD) :: IsItaNumber

    IsItaNumber = (.NOT. StringIsaNumber(string))

  END FUNCTION StringIsNotaNumber
!!$===========================================================================

!!$===========================================================================
  FUNCTION StringIsaWord(string) RESULT(IsaWord)
!!!****f* StringMod/StringIsaWord
!!!
!!! NAME
!!!     StringIsaWord -- Does the string contain only Letters?
!!! USAGE
!!!     StringIsaWord(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     LOGICAL(LGD) :: IsaWord
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL(LGD) :: IsaWord
    INTEGER(I4B) :: i, NumCharacters, ial,iau,izl,izu, i_temp

    ial = IACHAR("a")
    iau = IACHAR("A")
    izl = IACHAR("z")
    izu = IACHAR("Z")

    IsaWord = .TRUE.

    NumCharacters = LEN(string)
    DO i=1,numCharacters
       i_temp = IACHAR(string(i:i))
       IF( .NOT.(  (i_temp >= iau .AND. i_temp <= izu) .OR. &
            &      (i_temp >= ial .AND. i_temp <= izl) ) ) THEN
          IsaWord = .FALSE.
          EXIT
       ENDIF
    END DO

  END FUNCTION StringIsaWord
!!$===========================================================================

!!$===========================================================================
  FUNCTION FirstCharacter(string) RESULT(FC)
!!!****f* StringMod/FirstCharacter
!!!
!!! NAME
!!!     FirstCharacter -- returns first character
!!! USAGE
!!!     FirstCharacter(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=1) :: FC
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=1) :: FC

    CHARACTER(LEN=LEN(string)) :: tmp_string

    tmp_string = ADJUSTL(string)
    FC = tmp_string(1:1)

  END FUNCTION FirstCharacter
!!$===========================================================================

!!$===========================================================================
  FUNCTION LastCharacter(string) RESULT(LC)
!!!****f* StringMod/FirstCharacter
!!!
!!! NAME
!!!     LastCharacter -- returns last character
!!! USAGE
!!!     LastCharacter(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=1) :: LC
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=1) :: LC

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: l

    tmp_string = ADJUSTL(string)
    l = LEN_TRIM(tmp_string)
    LC = tmp_string(l:l)

  END FUNCTION LastCharacter
!!$===========================================================================

!!$===========================================================================
  SUBROUTINE SetupReadLineIO(DEFAULT,UNIT,ECHO,OUTUNIT,EXITSTRING, &
       & TAG,PROCESS,SEP,COMMENT,CASESENSITIVE,DIE)
!!!****f* StringMod/SetupReadLineIO
!!!
!!! NAME
!!!     SetupReadLineIO -- sets up defaults for ReadLineIO
!!! USAGE
!!!  SUBROUTINE SetupReadLineIO(DEFAULT,UNIT,ECHO,OUTUNIT,EXITSTRING, &
!!!       & TAG,PROCESS,SEP,COMMENT,CASESENSITIVE,DIE)
!!! DESCRIPTION
!!!     Sets up values for future calls to ReadLineIO.  All arguments are
!!!     OPTIONAL.
!!! INPUTS
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DEFAULT
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: UNIT
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: ECHO
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: OUTUNIT
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: EXITSTRING
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TAG
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: PROCESS
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DIE
!!! SEE ALSO
!!!     PublicVariables, ReadLineIO
!!!***
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DEFAULT
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: UNIT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: ECHO
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: OUTUNIT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: EXITSTRING
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TAG
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: PROCESS
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DIE

    IF ( OPTIONAL_FLAG(DEFAULT) ) THEN
       ReadLineIO_Setup = ReadLineIO_Default
    END IF

    IF ( PRESENT(UNIT) ) ReadLineIO_Setup%UNIT = UNIT
    IF ( PRESENT(OUTUNIT) ) ReadLineIO_Setup%OUTUNIT = OUTUNIT
    IF ( PRESENT(ECHO) ) ReadLineIO_Setup%ECHO = ECHO
    IF ( PRESENT(EXITSTRING) ) THEN
       ReadLineIO_Setup%ExitFlag = .TRUE.
       ReadLineIO_Setup%EXITSTRING = EXITSTRING
    END IF
    IF ( PRESENT(TAG) ) ReadLineIO_Setup%TAG = TAG
    IF ( PRESENT(PROCESS) ) ReadLineIO_Setup%PROCESS = PROCESS
    IF ( PRESENT(SEP) ) ReadLineIO_Setup%SEP = SEP
    IF ( PRESENT(COMMENT) ) THEN
       ReadLineIO_Setup%COMMENTSYM = COMMENT
       ReadLineIO_Setup%COMMENT = .TRUE.
    END IF
    IF ( PRESENT(CASESENSITIVE) ) &
         & ReadLineIO_Setup%CASESENSITIVE = CASESENSITIVE
    IF ( PRESENT(DIE) ) ReadLineIO_Setup%DIE = DIE

  END SUBROUTINE SetupReadLineIO

  FUNCTION ReadLineIO(line,UNIT,IOSTAT,ECHO,OUTUNIT,EXITSTRING, &
       & TAG,PROCESS,SEP,COMMENT,CASESENSITIVE,DIE,LINENUM)     &
       & RESULT(SUCCESS)
!!!****f* StringMod/ReadLineIO
!!!
!!! NAME
!!!     ReadLineIO -- reads a non-trivial line of IO
!!! USAGE
!!!  FUNCTION ReadLineIO(line,UNIT,IOSTAT,ECHO,OUTUNIT,EXITSTRING, &
!!!       & TAG,PROCESS,SEP,COMMENT,CASESENSITIVE,LINENUM)         &
!!!       & RESULT(SUCCESS)
!!! DESCRIPTION
!!!     Reads in a line of input from unit UNIT, or else unit 5.
!!!     OPTIONALLY the IOSTAT is returned, and INPUT can be
!!!     echoed using the LOGICAL flag ECHO and (OPTIONALLY)
!!!     specifying an OUTUNIT for the echoed output.  
!!!     ReadLineIO continues to read in lines until a "non-trivial"
!!!     line is read.  A "non-trivial" line is one that actually contains
!!!     something besides white space and COMMENTs.  If such a line
!!!     is found before the OPTIONAL EXITSTRING is found, the function
!!!     returns .TRUE. - else if either the IOSTAT is nonzero or the
!!!     EXITSTRING is found before a non-trivial input line is reached,
!!!     the function returns .FALSE.
!!!     If the DIE logical flag is set to .TRUE. and the EXITSTRING
!!!     is found, the program will stop (DIE).
!!!     A TAG can be specified such that input lines are stripped
!!!     of the TAG in the first field, and echoed lines comtain the TAG.
!!! INPUTS
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: UNIT
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: ECHO
!!!    INTEGER(I4B), OPTIONAL, INTENT(IN) :: OUTUNIT
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: EXITSTRING
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TAG
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: PROCESS
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!!    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DIE
!!!    INTEGER(I4B),  OPTIONAL,INTENT(INOUT) :: LINENUM
!!! OUTPUT
!!!    CHARACTER(LEN=*), INTENT(OUT) :: line
!!!    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
!!!    LOGICAL(LGD) :: SUCCESS
!!! SEE ALSO
!!!     PublicVariables, ProcessLine
!!!***
    CHARACTER(LEN=*), INTENT(OUT) :: line
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: UNIT
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: ECHO
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: OUTUNIT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: EXITSTRING
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: TAG
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: PROCESS
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: DIE
    INTEGER(I4B), OPTIONAL, INTENT(INOUT) :: LINENUM
    LOGICAL(LGD) :: SUCCESS

    INTEGER(I4B) :: iunit, iout, ierr, ltag, lline
    CHARACTER(LEN=10) :: FmtString
    LOGICAL(LGD) :: EchoFlag,DieFlag,ProcessFlag,CommentFlag,TagFlag, &
         & CaseFlag,ExitFlag
    CHARACTER(LEN=80) :: FC,CommentSym,TagString,SepString,Exit_String

    IF ( PRESENT(UNIT) ) THEN
       iunit=UNIT
    ELSE
       iunit=ReadLineIO_Setup%UNIT
    END IF

    IF ( PRESENT(OUTUNIT) ) THEN
       iout = OUTUNIT
    ELSE
       iout=ReadLineIO_Setup%OUTUNIT
    END IF

    IF ( LEN(line) < 10 ) THEN
       FmtString = TRIM('(a'//IntegerToString(LEN(line),1)//')')
    ELSE IF ( LEN(line) < 100 ) THEN
       FmtString = TRIM('(a'//IntegerToString(LEN(line),2)//')')
    ELSE IF ( LEN(line) < 1000 ) THEN
       FmtString = TRIM('(a'//IntegerToString(LEN(line),3)//')')
    ELSE
       CALL error('ReadLineIO - LEN(line) too large!',WARNLev=5)
    END IF

    IF ( PRESENT(ECHO) ) THEN
       EchoFlag = ECHO
    ELSE
       EchoFlag = ReadLineIO_Setup%ECHO
    END IF

    IF ( PRESENT(DIE) ) THEN
       DieFlag = DIE
    ELSE
       DieFlag = ReadLineIO_Setup%DIE
    END IF

    IF ( PRESENT(TAG) ) THEN
       TagFlag = .TRUE.
       TagString = TRIM(TAG)
    ELSE
       IF ( LEN_TRIM(ReadLineIO_Setup%TAG) == 0 ) THEN
          TagFlag = .FALSE.
          TagString = ' '
       ELSE
          TagFlag = .TRUE.
          TagString = TRIM(ReadLineIO_Setup%TAG)
       END IF
    END IF

    IF ( PRESENT(PROCESS) ) THEN
       ProcessFlag = PROCESS
    ELSE
       ProcessFlag = ReadLineIO_Setup%PROCESS
    END IF

    IF ( PRESENT(COMMENT) ) THEN
       CommentFlag = .TRUE.
       CommentSym = TRIM(COMMENT)
    ELSE
       CommentFlag = ReadLineIO_Setup%COMMENT
       CommentSym = TRIM(ReadLineIO_Setup%COMMENTSYM)
    END IF

    IF ( PRESENT(CASESENSITIVE) ) THEN
       CaseFlag = CASESENSITIVE
    ELSE
       CaseFlag =  ReadLineIO_Setup%CASESENSITIVE
    END IF

    IF ( PRESENT(SEP) ) THEN
       SepString = TRIM(SEP)
    ELSE
       SepString = TRIM(ReadLineIO_Setup%SEP)
    END IF

    IF ( PRESENT(EXITSTRING) ) THEN
       ExitFlag = .TRUE.
       Exit_String = TRIM(EXITSTRING)
    ELSE
       ExitFlag = ReadLineIO_Setup%ExitFlag
       Exit_String = TRIM(ReadLineIO_Setup%EXITSTRING)
    END IF

    SUCCESS = .TRUE.

    DO
       ! Read in line
       READ(iunit,FmtString,IOSTAT=ierr)line
       IF ( PRESENT(LINENUM) ) LINENUM = LINENUM+1
       IF ( PRESENT(IOSTAT) ) THEN
          IOSTAT=ierr
       END IF

       ! If line cannot be read - exit, return .FALSE.
       IF ( ierr /= 0 ) THEN
          IF ( DieFlag .AND. .NOT. PRESENT(IOSTAT) ) THEN
             CALL error('ReadLineIO: IOSTAT error!',i1=ierr,WARNLev=5)
          END IF
          SUCCESS = .FALSE.
          EXIT
       END IF

       ! Strip TAG from line if PRESENT
       IF ( TagFlag ) THEN
          ltag = LEN_TRIM(  ADJUSTL(TagString) )
          lline = LEN_TRIM(  ADJUSTL(line) )
          IF ( CompareString(line,TagString, &
               & TRUNC=.TRUE.,CASESENSITIVE=.FALSE.) ) THEN
             IF ( lline > ltag ) THEN
                line = ADJUSTL( TRIM( line ) )
                line = line(ltag+1:lline)
             ELSE
                line = ' '
             END IF
          END IF
       END IF

       ! If ECHO requested, ECHO the line
       IF ( EchoFlag ) THEN
          IF ( TagFlag ) THEN
             ! Add TAG to echoed output
             WRITE(iout,'(a)')ADJUSTL( TRIM(TagString) )// &
                  & ' '//ADJUSTL( TRIM(line) )
             IF ( iout /= 6 ) THEN
                ! Always echo to STDOUT
                WRITE(6,'(a)')ADJUSTL( TRIM(TagString) )// &
                     & ' '//ADJUSTL( TRIM(line) )
             END IF
          ELSE
             WRITE(iout,'(a)')ADJUSTL( TRIM(line) )
             IF ( iout /= 6 ) THEN
                ! Always echo to STDOUT
                WRITE(6,'(a)')ADJUSTL( TRIM(line) )
             END IF
          END IF
       END IF

       ! If just a comment line - cycle and look for non-trivial input
       IF ( CommentFlag ) THEN
          FC = ADJUSTL( line )
          IF ( FC(1:LEN_TRIM(CommentSym)) == TRIM(CommentSym) ) CYCLE
       END IF

       ! If the exit string is found - exit, return .FALSE.
       IF ( ExitFlag ) THEN
          IF ( CompareString(line,TRIM(Exit_String), &
               & CASESENSITIVE=CASESENSITIVE,TRUNC=.TRUE.) ) THEN
             IF ( DieFlag ) THEN
                CALL error('ReadLineIO: EXITSTRING found - program stopped.', &
                     & WARNLev=5)
             END IF
             SUCCESS = .FALSE.
             EXIT
          END IF
       END IF

       ! If PROCESS requested - ProcessLine
       IF ( ProcessFlag ) THEN
          IF ( CommentFlag ) THEN
             line = ProcessLine(line,         &
                  & SEP=TRIM(SepString),      &
                  & COMMENT=TRIM(CommentSym), &
                  & CASESENSITIVE=CaseFlag)
          ELSE
             line = ProcessLine(line, &
                  & SEP=TRIM(SepString),CASESENSITIVE=CaseFlag)
          END IF
       END IF

       IF ( LEN_TRIM(line) > 0 ) THEN
          ! A non-trivial line of input has been read and processed: EXIT.
          EXIT
       ELSE
          ! Trivial (empty) line results after processing: CYCLE.
          CYCLE
       END IF

    END DO

  END FUNCTION ReadLineIO
!!$===========================================================================

!!$===========================================================================
  FUNCTION ProcessLine(string,SEP,COMMENT,CASESENSITIVE) RESULT(a)
!!!****f* StringMod/ProcessLine
!!!
!!! NAME
!!!     ProcessLine -- processes input line for parsing
!!! USAGE
!!!     ProcessLine(string,SEP,COMMENT,CASESENSITIVE)
!!! DESCRIPTION
!!!     Performs processing such as white out separators, left adjustment,
!!!     trimming, and UPPERCASE.  Optionally, if SEP is specified, a
!!!     list of seperators is replaced by white space.  If COMMENT is 
!!!     specified, string is truncated to the right where the first 
!!!     occurence of COMMENT is found in string.
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: a
!!! SEE ALSO
!!!     PublicVariables
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: COMMENT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE

    CHARACTER(LEN=LEN(string)) :: a
    INTEGER(I4B) :: i

    IF ( PRESENT(COMMENT) ) THEN
       i = INDEX(string,COMMENT)
       IF ( i > 1 ) THEN
          a = string(1:i-1)
       ELSE IF ( i == 1 ) THEN
          a = ' '
       ELSE
          a = string
       END IF
    ELSE
       a = string
    END IF

    IF ( .NOT. OPTIONAL_FLAG(CASESENSITIVE,DEFAULT=.TRUE.) ) THEN
       a = UPPERCASE(a)
    END IF

    a = TRIM( ADJUSTL( WhiteOutSeparators( a, SEP=SEP ) ) )

  END FUNCTION ProcessLine
!!$===========================================================================

!!$===========================================================================
  FUNCTION WhiteOutSeparators(string,SEP) RESULT(a)
!!!****f* StringMod/WhiteOutSeparators
!!!
!!! NAME
!!!     WhiteOutSeparators -- replaces separators by spaces
!!! USAGE
!!!     WhiteOutSeparators(string)
!!! DESCRIPTION
!!!     Replaces all instances of separators with spaces.  If optional
!!!     argument SEP is given it used that list of separators, else
!!!     the default parameters SEPARATOR is used.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: a
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!! SEE ALSO
!!!     PublicVariables
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP

    CHARACTER(LEN=LEN(string)) :: a

    INTEGER(I4B) :: i

    a = string
    IF ( PRESENT(SEP) ) THEN
       IF ( LEN(SEP) > 0 ) THEN
          DO i=1,LEN(string)
             IF ( SCAN(SEP,a(i:i)) /= 0 ) a(i:i) = ' '
          END DO
       END IF
    ELSE
       DO i=1,LEN(string)
          IF ( SCAN(SEPERATOR,a(i:i)) /= 0 ) a(i:i) = ' '
       END DO
    END IF

  END FUNCTION WhiteOutSeparators
!!$===========================================================================

!!$===========================================================================
  RECURSIVE SUBROUTINE SkipField(string,SEP,SKIP)
!!!****f* StringMod/SkipField
!!!
!!! NAME
!!!     SkipField -- skips a field
!!! USAGE
!!!     SkipField(string,SEP,SKIP)
!!! DESCRIPTION
!!!     Returns the string with first SKIP fields removed.
!!!     Fields are determined by OPTIONAL argument SEP if PRESENT,
!!!     else by global string SEPARATOR.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(INOUT) :: string
!!!     INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP
!!! OUTPUT
!!!     CHARACTER(LEN=*), INTENT(INOUT) :: string
!!! SEE ALSO
!!!     PublicVariables
!!!***
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: i,iskip

    IF ( PRESENT(SKIP) ) THEN
       iskip=SKIP
       DO i=1,iskip
          CALL SkipField(string,SEP=SEP)
       END DO
    ELSE
       tmp_string = ADJUSTL(string)
       IF ( PRESENT(SEP) ) THEN
          i = SCAN(tmp_string,SEP)
       ELSE
          i = SCAN(tmp_string,SEPERATOR)
       END IF
       IF ( i < LEN(tmp_string) ) THEN
          tmp_string = ADJUSTL( tmp_string(i+1:) )
       ELSE
          tmp_string = ' '
       END IF
       string = TRIM(tmp_string)
    END IF

  END SUBROUTINE SkipField
!!$===========================================================================

!!$===========================================================================
  FUNCTION InvertString(string) RESULT(a)
!!!****f* StringMod/InvertString
!!!
!!! NAME
!!!     InvertString -- Invert a string.
!!! DESCRIPTION
!!!     Returns a string that is the inverted input string.
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=LEN(string)) :: a

    INTEGER(I4B) :: i,j,n

    n = LEN(string)

    j=n
    DO i=1,n
       a(j:j) = string(i:i)
       j=j-1
    END DO

  END FUNCTION InvertString
!!$===========================================================================

!!$===========================================================================
  FUNCTION ReverseFields(string) RESULT(a)
!!!****f* StringMod/ReverseFields
!!!
!!! NAME
!!!     ReverseFields -- Reverse the order of fields in a string.
!!! DESCRIPTION
!!!     Returns a string that has the order of the fields reversed, and
!!!     separated by white spaces.
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=LEN(string)) :: a

    CHARACTER(LEN=LEN(string)) :: tmp_string, tmp_field

    INTEGER(I4B) :: ierr

    tmp_string = InvertString( ADJUSTL( TRIM(string) ) )

    a = ' '

    DO
       tmp_field = ' '
       CALL ReadField(tmp_string, tmp_field, ADVANCE=1, IOSTAT=ierr)
       IF ( ierr == 0 ) THEN
          tmp_field = ADJUSTL( TRIM( InvertString( tmp_field ) ) )
          a = TRIM( a )//' '//TRIM( tmp_field )
       ELSE
          EXIT
       END IF

    END DO

  END FUNCTION ReverseFields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$===========================================================================
  SUBROUTINE ReadField_h(string,a,SEP,KEYWORD,SKIP,ADVANCE,IOSTAT, &
       & CASESENSITIVE,TRUNC,WILDCARD,WARN,WARNLEV,DEFAULT,SUCCESS)
!!!****f* ReadField/ReadField_h
!!!
!!! NAME
!!!     ReadField_h -- Read a character field
!!! DESCRIPTION
!!!     ReadField when a is a character
!!!***
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    CHARACTER(LEN=*), INTENT(INOUT) :: a
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: KEYWORD
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: ADVANCE
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: WARNLEV
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: WARN
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: DEFAULT
    LOGICAL(LGD), OPTIONAL, INTENT(OUT) :: SUCCESS

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: i,ierr,lev
    LOGICAL(LGD) :: lsuccess

    IF ( PRESENT(DEFAULT) ) THEN
       a = TRIM(DEFAULT)
    END IF

    lsuccess = .FALSE.
    IF ( PRESENT(SUCCESS) ) SUCCESS = lsuccess

    lev = 5
    IF ( PRESENT( WARNLEV ) ) lev   = WARNLEV

    tmp_string = TRIM( ADJUSTL(WhiteOutSeparators(string,SEP=SEP) ) )
    !tmp_string = TRIM( ADJUSTL( string ) )

    IF ( PRESENT(KEYWORD) ) THEN
       i = IndexKeyword(tmp_string,KEYWORD, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
       IF ( i == 0 .OR. i == LEN_TRIM(tmp_string) ) THEN
          IF ( PRESENT(IOSTAT) ) THEN
             IOSTAT = 1
          ELSE
             IF ( OPTIONAL_FLAG(WARN,DEFAULT=.TRUE.) ) THEN
                CALL error('ReadField_h: error reading KEYWORD field', &
                     & WARNLev=lev)
             END IF
          END IF
          RETURN
       END IF
       tmp_string = ADJUSTL( tmp_string(i+LEN_TRIM(KEYWORD):) )
       !CALL SkipField( tmp_string, SKIP=1 )
    END IF

    IF ( PRESENT(SKIP) ) THEN
!       CALL SkipField(tmp_string,SEP=SEP,SKIP=SKIP)
       CALL SkipField(tmp_string,SEP=' ',SKIP=SKIP)
    END IF

! Don't need anymore because SEP has been turned to whitespace above
!    IF ( PRESENT(SEP) ) THEN
!       i = SCAN(tmp_string,SEP)
!    ELSE
!       i = SCAN(tmp_string,SEPERATOR)
!    END IF
    i = index(tmp_string,' ')

    ierr = 0
    IF ( i == 0 .AND. LEN_TRIM(tmp_string) > 0 ) THEN
       ! There are no intervening SEParators - pass back the TRIMmed string
       a = TRIM( ADJUSTL(tmp_string) )
    ELSE IF ( i == 1 .OR. LEN_TRIM(tmp_string) == 0 ) THEN
       ! There is nothing BUT SEParators - pass back ierr = 1
       ierr = 1
       a = ' '
    ELSE
       ! Found a SEParator - pass back string field
       a = TRIM(  ADJUSTL(tmp_string(1:i-1)) )
    END IF

    IF ( PRESENT(IOSTAT) ) THEN
       IOSTAT = ierr
    ELSE IF ( ierr /= 0 ) THEN
       IF ( OPTIONAL_FLAG(WARN,DEFAULT=.TRUE.) ) THEN
          CALL error('ReadField_h: error reading string ',WARNLev=lev,i1=ierr)
       END IF
    END IF

    IF ( LEN_TRIM(a) > 0 ) THEN
       lsuccess = .TRUE.
    END IF

    IF ( PRESENT(SUCCESS) ) SUCCESS = lsuccess

    IF ( PRESENT(ADVANCE) ) THEN
!       CALL SkipField(tmp_string,SKIP=ADVANCE,SEP=SEP)
       CALL SkipField(tmp_string,SEP=' ',SKIP=ADVANCE)
       string = TRIM( ADJUSTL(tmp_string) )
    END IF

  END SUBROUTINE ReadField_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ReadField_i(string,a,SEP,KEYWORD,SKIP,ADVANCE,IOSTAT, &
       & CASESENSITIVE,TRUNC,WILDCARD,DEFAULT,SUCCESS)
!!!****f* ReadField/ReadField_i
!!!
!!! NAME
!!!     ReadField_i -- Read an integer field
!!! DESCRIPTION
!!!     ReadField when a is an integer
!!!***    
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    INTEGER(I4B), INTENT(INOUT) :: a
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: KEYWORD
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: ADVANCE
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: DEFAULT
    LOGICAL(LGD), OPTIONAL, INTENT(OUT) :: SUCCESS

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: ierr
    INTEGER(I4B) :: tmp_val
    LOGICAL(LGD) :: lsuccess

    IF ( PRESENT(DEFAULT) ) THEN
       a = DEFAULT
    END IF

    tmp_string = ' '
    CALL ReadField(string,tmp_string,                                       &
         & SEP=SEP,KEYWORD=KEYWORD,SKIP=SKIP,ADVANCE=ADVANCE,IOSTAT=IOSTAT, &
         & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD,       &
         & SUCCESS=lsuccess)

    IF ( lsuccess ) THEN
       READ(tmp_string,*,iostat=ierr)tmp_val
       IF ( PRESENT(IOSTAT) ) THEN
          IOSTAT = ierr
       END IF
       IF ( ierr /= 0 ) THEN
          lsuccess = .FALSE.
          IF ( .NOT. (PRESENT(IOSTAT) .OR. PRESENT(SUCCESS)) ) THEN
             CALL error('ReadField_i: error reading string ',WARNLev=5)
          END IF
       ELSE
          a=tmp_val
       END IF
    END IF

    IF ( PRESENT(SUCCESS) ) SUCCESS = lsuccess

  END SUBROUTINE ReadField_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ReadField_r(string,a,SEP,KEYWORD,SKIP,ADVANCE,IOSTAT, &
       & CASESENSITIVE,TRUNC,WILDCARD,DEFAULT,SUCCESS)
!!!****f* ReadField/ReadField_r
!!!
!!! NAME
!!!     ReadField_r -- Read a real field
!!! DESCRIPTION
!!!     ReadField when a is a real
!!!***
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    REAL(SP), INTENT(INOUT) :: a
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: KEYWORD
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: ADVANCE
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: SKIP
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
    REAL(SP), OPTIONAL, INTENT(IN) :: DEFAULT
    LOGICAL(LGD), OPTIONAL, INTENT(OUT) :: SUCCESS

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: ierr
    REAL(SP):: tmp_val
    LOGICAL(LGD) :: lsuccess

    IF ( PRESENT(DEFAULT) ) THEN
       a = DEFAULT
    END IF

    tmp_string = ' '
    CALL ReadField(string,tmp_string,                                       &
         & SEP=SEP,KEYWORD=KEYWORD,SKIP=SKIP,ADVANCE=ADVANCE,IOSTAT=IOSTAT, &
         & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD,       &
         & SUCCESS=lsuccess)

    IF ( lsuccess ) THEN
       READ(tmp_string,*,iostat=ierr)tmp_val
       IF ( PRESENT(IOSTAT) ) THEN
          IOSTAT = ierr
       END IF
       IF ( ierr /= 0 ) THEN
          lsuccess = .FALSE.
          IF ( .NOT. (PRESENT(IOSTAT) .OR. PRESENT(SUCCESS)) ) THEN
             CALL error('ReadField_r: error reading string ',WARNLev=5)
          END IF
       ELSE
          a=tmp_val
       END IF
    END IF

    IF ( PRESENT(SUCCESS) ) SUCCESS = lsuccess

  END SUBROUTINE ReadField_r
!!$===========================================================================

!!$===========================================================================
  FUNCTION RemoveSpaces(string) RESULT(FS)
!!!****f* StringMod/RemoveSpaces
!!!
!!! NAME
!!!     RemoveSpaces -- returns a string with no spaces
!!! USAGE
!!!     RemoveSpaces(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: FS
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string

    CHARACTER(LEN=LEN(string)) :: FS

    CHARACTER(LEN=LEN(string)) :: tmp_string
    INTEGER(I4B) :: i,j

    j = 0
    DO i=1, LEN(string)
       IF ( string(i:i) /= ' ' ) THEN
          j=j+1
          tmp_string(j:j) = string(i:i)
       END IF
    END DO

    IF ( j > 0 ) THEN
       FS = tmp_string(1:j)
    ELSE
       FS = ''
    END IF

  END FUNCTION RemoveSpaces
!!$===========================================================================

!!$===========================================================================
  FUNCTION RemoveSeparators(string,SEP) RESULT(FS)
!!!****f* StringMod/RemoveSpaces
!!!
!!! NAME
!!!     RemoveSeparators -- returns a string with no separators
!!! USAGE
!!!     RemoveSeparators(string)
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: FS
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: SEP

    CHARACTER(LEN=LEN(string)) :: FS

    FS = TRIM( ADJUSTL ( RemoveSpaces( WhiteOutSeparators(string,SEP=SEP) ) ) )

  END FUNCTION RemoveSeparators
!!$===========================================================================

!!$===========================================================================
  RECURSIVE FUNCTION CompareString(string, refstring, &
       & CASESENSITIVE, TRUNC, WILDCARD) &
       & RESULT(LCOMP)
!!!****f* StringMod/CompareString
!!!
!!! NAME
!!!     CompareString -- compare strings
!!! USAGE
!!!     CompareString(string, refstring, CASESENSITIVE, TRUNC, WILDCARD)
!!! DESCRIPTION
!!!     Returns .TRUE. if string is ostensibly the same as refstring, as
!!!     far as refstring goes.  
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!!     If you want the comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of "refstring" - set TRUNC=.TRUE. (default .TRUE.,
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the refstring.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string, refstring
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     LOGICAL(LGD) :: LCOMP
!!! SEE ALSO
!!!     StringContains, CompareStringFields, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string, refstring
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
    LOGICAL(LGD) :: LCOMP

    CHARACTER(LEN=LEN(string)) :: tmp_string
    CHARACTER(LEN=LEN(refstring)) :: tmp_refstring

    INTEGER(I4B) :: l,lref

    CHARACTER(LEN=LEN(refstring)) :: refstring1, refstring2
    INTEGER(I4B) :: i1,i2,j1,j2,lrefstring1,lrefstring2
    LOGICAL(LGD) :: FirstIsWild, LastIsWild, ThereIsMoreCheckingToDo

    tmp_refstring = TRIM( ADJUSTL( refstring ) )
    tmp_string    = TRIM( ADJUSTL( string ) )

    lref = LEN_TRIM(tmp_refstring)
    l    = LEN_TRIM(tmp_string)

    IF ( OPTIONAL_FLAG(TRUNC) .OR. &
         & (.NOT. PRESENT(TRUNC) .AND. .NOT. PRESENT(WILDCARD)) ) THEN
       IF ( lref == 0 ) THEN
          ! reference string is empty - LCOMP must is always .TRUE.
          ! if TRUNCation has been specified (by default or otherwise).
          LCOMP = .TRUE.
          RETURN
       END IF
       IF ( lref < l  ) THEN
          tmp_string = TRIM( ADJUSTL( tmp_string(1:lref) ) )
       END IF
    END IF

    IF ( .NOT. OPTIONAL_FLAG(CASESENSITIVE,DEFAULT=.FALSE.) ) THEN
       tmp_refstring = UPPERCASE( tmp_refstring )
       tmp_string    = UPPERCASE( tmp_string )
    END IF

    IF ( PRESENT(WILDCARD) ) THEN

       IF ( StringContains(tmp_refstring,WILDCARD) ) THEN
          ! The tmp_refstring contains a WILDCARD symbol

          ! Process tmp_refstring to remove double WILDCARD symbols
          ! that are immediately adjacent or separated only by white space.
          ! These symbols are by the current definition redundant.
          i1 = INDEX(tmp_refstring(1:lref),WILDCARD)
          DO
             IF ( i1 > 0 .AND. i1 < lref ) THEN
                DO
                   tmp_refstring = TRIM( ADJUSTL( tmp_refstring(1:i1) ) ) // &
                        & TRIM( ADJUSTL( tmp_refstring(i1+1:lref) ) )
                   lref = LEN_TRIM( tmp_refstring )
                   IF ( tmp_refstring(i1+1:i1+1) == WILDCARD ) THEN
                      tmp_refstring(i1+1:i1+1) = ' '
                   ELSE
                      EXIT
                   END IF
                END DO
                i2 = i1 + INDEX(tmp_refstring(i1+1:lref),WILDCARD)
                IF ( i1 == i2 ) THEN
                   EXIT
                ELSE
                   i1 = i2
                END IF
             ELSE
                EXIT
             END IF
          END DO

          ! Find where in tmp_refstring the first WILDCARD appears
          i1 = INDEX(tmp_refstring(1:lref),WILDCARD)
          i1 = MIN(i1,lref)

          IF ( i1 <= 0 ) THEN

             ! How can this be if StringContains(tmp_refstring,WILDCARD)?
             CALL error('CompareString: i1 <= 0!',i1=i1,WARNLev=5)

          ELSE IF ( i1 == lref .AND. lref == 1 ) THEN

             ! The reference string is just a single WILDCARD symbol
             LCOMP = .TRUE.
             RETURN

          END IF

          IF ( i1 == 1 ) THEN
             ! The WILDCARD symbol is the first non-blank
             ! CHARACTER in the reference string.
             FirstIsWild = .TRUE.
          ELSE
             FirstIsWild = .FALSE.
          END IF

          ! Calculate i2 - the position of the next WILDCARD symbol
          i2 = i1 + INDEX(tmp_refstring(i1+1:lref),WILDCARD)
          i2 = MIN(i2,lref)

          IF ( i1 /= i2 .AND. i1+1 > i2-1 ) THEN
             CALL error('CompareString: i1 /= i2 .AND. i1+1 > i2-1', &
                  & WARNLev=5)
          END IF

          IF ( FirstIsWild .AND. i1 == i2 ) THEN
             ! The first character in the reference string is the WILDCARD
             ! symbol, and there is no other WILDCARD symbol.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(i1+1:lref) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             LastIsWild = .FALSE.
             ThereIsMoreCheckingToDo = .FALSE.
          ELSE IF ( FirstIsWild .AND. i2 == lref ) THEN
             ! The first character in the reference string is the WILDCARD
             ! symbol, and there one other WILDCARD symbol that is the
             ! last character in the reference string.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(i1+1:i2-1) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             LastIsWild = .TRUE.
             ThereIsMoreCheckingToDo = .FALSE.
          ELSE IF ( FirstIsWild .AND. i2 < lref ) THEN
             ! The first character in the reference string is the WILDCARD
             ! symbol, and there is another WILDCARD symbol that is NOT the
             ! last character in the reference string.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(i1+1:i2-1) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             refstring2    = TRIM( ADJUSTL( tmp_refstring(i2:lref) ) )
             lrefstring2   = LEN_TRIM( refstring2 )
             LastIsWild = .TRUE.
             ThereIsMoreCheckingToDo = .TRUE.
          ELSE IF ( .NOT. FirstIsWild .AND. i1 == i2 .AND. i2 < lref ) THEN
             ! There is only 1 WILDCARD symbol - but it is not the first
             ! or last character of the reference string.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(1:i1-1) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             refstring2    = TRIM( ADJUSTL( tmp_refstring(i1:lref) ) )
             lrefstring2   = LEN_TRIM( refstring2 )
             LastIsWild = .TRUE.
             ThereIsMoreCheckingToDo = .TRUE.
          ELSE IF ( .NOT. FirstIsWild .AND. i1 == i2 .AND. i2 == lref ) THEN
             ! There is only 1 WILDCARD symbol and it is the last
             ! character of the reference string.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(1:i1-1) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             LastIsWild = .TRUE.
             ThereIsMoreCheckingToDo = .FALSE.
          ELSE IF ( .NOT. FirstIsWild .AND. i1 < i2 ) THEN
             ! There is more than 1 WILDCARD symbol, none of which
             ! are the first and the first character of the reference string.
             refstring1    = TRIM( ADJUSTL( tmp_refstring(1:i1-1) ) )
             lrefstring1   = LEN_TRIM( refstring1 )
             refstring2    = TRIM( ADJUSTL( tmp_refstring(i1:lref) ) )
             lrefstring2   = LEN_TRIM( refstring2 )
             LastIsWild = .TRUE.
             ThereIsMoreCheckingToDo = .TRUE.
          ELSE
             CALL error('CompareString: How did we get here?',WARNLev=5)
          END IF

          j1 = INDEX(tmp_string(1:l),refstring1(1:lrefstring1))
          j2 = j1 + lrefstring1

          IF ( j1 == 0 ) THEN
             ! The first part of the reference string was not found
             ! anywhere in the string.
             LCOMP = .FALSE.
             RETURN                
          ELSE IF ( .NOT. FirstIsWild .AND. j1 > 1 ) THEN
             ! The reference string is not lead by a WILDCARD symbol
             ! (.NOT. FirstIsWild) and therefore the FIRST non-blank part
             ! of the string must match with the FIRST part of the 
             ! reference string (j==1).
             ! It is insufficient for the first part of the reference 
             ! string to merely be contained in the string anywhere.  
             LCOMP = .FALSE.
             RETURN
          ELSE IF ( .NOT. LastIsWild .AND. j2 < l ) THEN
             ! The reference string does not end with a WILDCARD symbol
             ! (.NOT. LastIsWild) and therefore the LAST non-blank part
             ! of the string must match with the LAST part of the 
             ! reference string (j2 < l).
             ! It is insufficient for the last part of the reference 
             ! string to merely be contained in the string anywhere.
             LCOMP = .FALSE.
             RETURN
          ELSE IF ( LastIsWild .AND. ThereIsMoreCheckingToDo ) THEN
             ! The first part of the reference string ends with
             ! second WILDCARD - and there are more non-blank characters
             ! in the reference string to the right of the second
             ! WILDCARD - so we must continue to check.
             IF ( j2 > l ) THEN
                ! Oh-oh.  We have more checking to do, but
                ! we are at the end of the string!
                LCOMP = .FALSE.
                RETURN
             END IF
             tmp_string = TRIM( ADJUSTL( tmp_string(j2:l) ) )
             LCOMP = CompareString(TRIM(tmp_string), &
                  & TRIM(refstring2),TRUNC=.FALSE.,WILDCARD=WILDCARD)
             RETURN
          ELSE
             LCOMP = .TRUE.
             RETURN
          END IF

          CALL error('CompareString - should have returned!',WARNLev=5)

       ELSE IF ( StringContains(tmp_string,WILDCARD) ) THEN
          ! The tmp_string contains the WILDCARD

          LCOMP = CompareString(tmp_refstring, tmp_string, WILDCARD=WILDCARD)
          RETURN

       END IF

    END IF

    IF ( TRIM(tmp_string) == TRIM(tmp_refstring) ) THEN
       LCOMP = .TRUE.
    ELSE
       LCOMP = .FALSE.
    END IF

  END FUNCTION CompareString
!!$===========================================================================

!!$===========================================================================
  FUNCTION IndexKeyword(string,substring,CASESENSITIVE,TRUNC,WILDCARD,BACK) &
       & RESULT(KI)
!!!****f* StringMod/StringContains
!!!
!!! NAME
!!!     IndexKeyword - like INDEX intrinsic, except with OPTIONAL
!!!       CASESENSITIVE,TRUNC,WILDCARD arguments.
!!! USAGE
!!!     IndexKeyword(string,substring,CASESENSITIVE,TRUNC,WILDCARD,BACK)
!!! DESCRIPTION
!!!     Basically a souped-up version of the following statement:
!!!     INDEX(string,substring,BACK)
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!!     If you want the comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of "substring" - set TRUNC=.TRUE. (default .TRUE.,
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the string.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), INTENT(IN) :: substring
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     LOGICAL(LGD) :: TF
!!! SEE ALSO
!!!     CompareString, CompareStringFields, ReplaceString, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), INTENT(IN) :: substring
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: BACK

    CHARACTER(LEN=LEN(string)) :: tmp_string
    CHARACTER(LEN=LEN(substring)) :: tmp_substring

    INTEGER(I4B) :: KI
    INTEGER(I4B) :: i

    KI = 0

    ! If there is no WILDCARD or TRUNC, this will be much faster...
    IF ( .NOT. (PRESENT(WILDCARD) .OR. PRESENT(TRUNC)) ) THEN
       IF ( .NOT. OPTIONAL_FLAG(CASESENSITIVE,DEFAULT=.TRUE.) ) THEN
          tmp_string    = UPPERCASE(string)
          tmp_substring = UPPERCASE(substring)
       ELSE
          tmp_string    = string
          tmp_substring = substring
       END IF
       IF ( PRESENT(BACK) ) THEN
          KI = INDEX(tmp_string,tmp_substring,BACK=BACK)
       ELSE
          KI = INDEX(tmp_string,tmp_substring)
       END IF
       RETURN
    END IF

    ! Ok - here is for the more complicated situation - CompareString will be
    ! the real workhorse.
    tmp_string    = string
    tmp_substring = substring

    IF ( OPTIONAL_FLAG(BACK) ) THEN
       tmp_string    = InvertString(tmp_string)
       tmp_substring = InvertString(tmp_substring)
    END IF

    DO i=1,LEN(tmp_string)
       IF ( CompareString(tmp_string, tmp_substring, &
            & CASESENSITIVE=CASESENSITIVE, WILDCARD=WILDCARD, TRUNC=TRUNC) &
            & ) THEN

          IF ( OPTIONAL_FLAG(BACK) ) THEN
             KI = LEN(tmp_string) - i + 1
          ELSE
             KI = i
          END IF
          EXIT
       END IF
    END DO

  END FUNCTION IndexKeyword
!!$===========================================================================


!!$===========================================================================
  FUNCTION ReplaceString(string, FindString, ReplaceWith, &
       & CASESENSITIVE) &
       & RESULT(NewString)
!!!****f* StringMod/ReplaceString
!!!
!!! NAME
!!!     ReplaceString -- search a string for a substring and replace it
!!!     with another substring.
!!! USAGE
!!!     ReplaceString(string, FindString, ReplaceWith, &
!!!       & CASESENSITIVE)
!!! DESCRIPTION
!!!     Returns a new string the same size as the original with all
!!!     instances of "FindString" replaced by "ReplaceWith".
!!!     If ReplaceWith is longer than FindString, the end of the
!!!     string is truncated, and if FindString is longer than ReplaceWith,
!!!     the end of the string is padded with white space.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string, refstring
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: NewString
!!! SEE ALSO
!!!     StringContains, CompareString, CompareStringFields, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), INTENT(IN) :: FindString
    CHARACTER(LEN=*), INTENT(IN) :: ReplaceWith
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    CHARACTER(LEN=LEN(string)) :: NewString

    CHARACTER(LEN=LEN(string)) :: tmp_string
    CHARACTER(LEN=LEN(FindSTring)) :: tmp_find
    INTEGER(I4B) :: i,i0,j0,jplus,l,lfind,lrep

    IF ( OPTIONAL_FLAG(CASESENSITIVE,DEFAULT=.TRUE.) ) THEN
       tmp_string = string
       tmp_find   = FindString
    ELSE
       tmp_string = UPPERCASE(string)
       tmp_find   = UPPERCASE(FindString)
    END IF

    l     = LEN(string)
    lfind = LEN(FindString)
    lrep  = LEN(ReplaceWith)

    NewString = ' '

    i0 = 1
    j0 = 1
    DO
       IF ( i0 <= l ) THEN
          i = INDEX(tmp_string(i0:l),tmp_find)
       ELSE
          i = 0
       END IF

       IF ( i > 0 ) THEN
          ! String to replace has been found.
          jplus = (i-i0)

          IF ( j0+jplus+lrep > l .OR. i0+jplus-1 > l ) THEN
             CALL error('ReplaceString: End of string reached!', &
                  & WARNLev=5)
          ELSE IF ( jplus > 1 ) THEN
             NewString(j0:j0+jplus-1) = string(i0:i0+jplus-1)
          END IF
          NewString(j0+jplus:j0+jplus+lrep-1) = ReplaceWith
          j0 = j0+jplus+lrep
          i0 = i0+jplus+lfind

       ELSE
          ! String to replace has NOT been found - finish up.

          IF ( i0 <= l .AND. j0 <= l ) THEN
             IF ( i0 < j0 ) THEN
                jplus = l-j0
             ELSE
                jplus = l-i0
             END IF
             NewString(j0:j0+jplus) = string(i0:i0+jplus)

             EXIT

          END IF
       END IF

    END DO

  END FUNCTION ReplaceString
!!$===========================================================================

!!$===========================================================================
  FUNCTION StripPrefix(string, Prefix) &
       & RESULT(NewString)
!!!****f* StringMod/StripPrefix
!!!
!!! NAME
!!!     StripPrefix -- strips off the prefix of a string.
!!! USAGE
!!!     StripPrefix(string, Prefix)
!!! DESCRIPTION
!!!     Returns a new string with the prefix, as defined by all
!!!     the characters to the left of (and including) the last 
!!!     character input string Prefix.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), INTENT(IN) :: Prefix
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: NewString
!!! SEE ALSO
!!!     StringContains, CompareString, CompareStringFields, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), INTENT(IN) :: Prefix

    CHARACTER(LEN=LEN(string)) :: NewString

    INTEGER(I4B) :: i,l,ls

    l = LEN_TRIM(string)
    ls = LEN_TRIM(Prefix)

    i = SCAN(string,TRIM(Prefix))

    IF ( i /= 0 .AND. (i+ls) <= l ) THEN
       NewString = TRIM(ADJUSTL(string(i+ls:l)))
    ELSE
       NewString = TRIM(ADJUSTL(string))
    END IF

  END FUNCTION StripPrefix
!!$===========================================================================

!!$===========================================================================
  FUNCTION StripSuffix(string, Suffix) &
       & RESULT(NewString)
!!!****f* StringMod/StripSuffix
!!!
!!! NAME
!!!     StripSuffix -- strips off the suffix of a string.
!!! USAGE
!!!     StripSuffix(string, Suffix)
!!! DESCRIPTION
!!!     Returns a new string with the suffix, as defined by all
!!!     the characters to the right of (and including) the first 
!!!     character input string "Suffix".
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), INTENT(IN) :: Suffix
!!! OUTPUT
!!!     CHARACTER(LEN=LEN(string)) :: NewString
!!! SEE ALSO
!!!     StringContains, CompareString, CompareStringFields, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), INTENT(IN) :: Suffix

    CHARACTER(LEN=LEN(string)) :: NewString

    CHARACTER(LEN=LEN(string)) :: tmp_string
    CHARACTER(LEN=LEN(Suffix)) :: tmp_Suffix

    tmp_string = InvertString( ADJUSTL( TRIM(string) ) )
    tmp_Suffix = InvertString( ADJUSTL( TRIM(Suffix) ) )

    NewString = InvertString( ADJUSTL( TRIM( &
         & StripPrefix(tmp_string,tmp_Suffix) ) ) )

  END FUNCTION StripSuffix
!!$===========================================================================

!!$===========================================================================
  FUNCTION StringContains(string,substring,CASESENSITIVE,TRUNC,WILDCARD) &
       & RESULT(TF)
!!!****f* StringMod/StringContains
!!!
!!! NAME
!!!     StringContains -- does a string contain a substring?
!!! USAGE
!!!     StringContains(string,substring,CASESENSITIVE)
!!! DESCRIPTION
!!!     Basically a souped-up version of the following statement:
!!!     INDEX(string,substring) == 0
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!!     An optional WILDCARD symbol can be included in the substring.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), INTENT(IN) :: substring
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     LOGICAL(LGD) :: TF
!!! SEE ALSO
!!!     CompareString, CompareStringFields, ReplaceString, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), INTENT(IN) :: substring
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD

!!$    CHARACTER(LEN=LEN(substring)+2) :: newsubstring

    LOGICAL(LGD) :: TF

    IF ( IndexKeyword(string,substring, &
         & CASESENSITIVE=CASESENSITIVE, &
         & TRUNC=TRUNC,                 &
         & WILDCARD=WILDCARD) == 0 ) THEN
       TF = .FALSE.
    ELSE
       TF = .TRUE.
    END IF
!!$
!!$    IF ( PRESENT(WILDCARD) ) THEN
!!$       newsubstring = WILDCARD//TRIM(ADJUSTL(substring))//WILDCARD
!!$       TF = CompareString(string, newsubstring, &
!!$            & TRUNC=.FALSE.,CASESENSITIVE=CASESENSITIVE,WILDCARD=WILDCARD)
!!$    ELSE IF ( OPTIONAL_FLAG(CASESENSITIVE) ) THEN
!!$       IF ( INDEX(string,substring) == 0 ) THEN
!!$          TF = .FALSE.
!!$       ELSE
!!$          TF = .TRUE.
!!$       END IF
!!$    ELSE
!!$       IF ( INDEX(LOWERCASE(string),LOWERCASE(substring)) == 0 ) THEN
!!$          TF = .FALSE.
!!$       ELSE
!!$          TF = .TRUE.
!!$       END IF
!!$    END IF

  END FUNCTION StringContains
!!$===========================================================================

!!$===========================================================================
  FUNCTION CompareStringFields(string, &
       & FIELD1, FIELD2, FIELD3, FIELD4, FIELD5, &
       & CASESENSITIVE, TRUNC, WILDCARD) RESULT(LCOMP)
!!!****f* StringMod/CompareStringFields
!!!
!!! NAME
!!!     CompareStringFields -- compare string fields
!!! USAGE
!!!     CompareStringFields(string,FIELD1='field1',FIELD2='field2',... &
!!!     & CASESENSITIVE)
!!! DESCRIPTION
!!!     Returns .TRUE. if fields in string are ostensibly the same as 
!!!     the reference FIELDS1....  This is used, for example, to compare
!!!     a list of keywords from a line of input.  
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!!     If you want the comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of "refstring" - set TRUNC=.TRUE. (default .TRUE.
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the refstring.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FIELD1, FIELD2, FIELD3, &
!!!          & FIELD4, FIELD5
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     LOGICAL(LGD) :: LCOMP
!!! SEE ALSO
!!!     CompareString, StringContains, LookupString
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FIELD1, FIELD2, FIELD3, &
         & FIELD4, FIELD5
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD

    LOGICAL(LGD) :: LCOMP

    CHARACTER(LEN=LEN(string)) :: tmp_string

    tmp_string = TRIM( ADJUSTL(string) )

    LCOMP = .TRUE.

    IF ( LCOMP .AND. PRESENT(FIELD1) ) THEN
       LCOMP = CompareString(tmp_string,FIELD1, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
    END IF
    CALL SkipField(tmp_string,SKIP=1)

    IF ( LCOMP .AND. PRESENT(FIELD2) ) THEN
       LCOMP = CompareString(tmp_string,FIELD2, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
    END IF
    CALL SkipField(tmp_string,SKIP=1)

    IF ( LCOMP .AND. PRESENT(FIELD3) ) THEN
       LCOMP = CompareString(tmp_string,FIELD3, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
    END IF
    CALL SkipField(tmp_string,SKIP=1)

    IF ( LCOMP .AND. PRESENT(FIELD4) ) THEN
       LCOMP = CompareString(tmp_string,FIELD4, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
    END IF
    CALL SkipField(tmp_string,SKIP=1)

    IF ( LCOMP .AND. PRESENT(FIELD5) ) THEN
       LCOMP = CompareString(tmp_string,FIELD5, &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD)
    END IF
    CALL SkipField(tmp_string,SKIP=1)

  END FUNCTION CompareStringFields
!!$===========================================================================

!!$===========================================================================
  FUNCTION LookupString1(string1, string_list1, &
       & CASESENSITIVE, TRUNC, WILDCARD ) RESULT(il)
!!!****f* StringMod/LookupString
!!!
!!! NAME
!!!     LookupString -- look for a string within a list of strings
!!! USAGE
!!!     LookupString(string1, string_list1, CASESENSITIVE)
!!!     LookupString(string1, ...stringN, & 
!!!                 & string_list1, ...string_listN, &
!!!                 & CASESENSITIVE)
!!! DESCRIPTION
!!!     Looks for the member of string_listX that is identical to
!!!     stringX, and returns the single index for the first match.
!!!     If no matches are found, 0 is returned.  Currently, function
!!!     is overloaded up to N=3.  
!!!     You can toggle case-sensitivity by setting the optional
!!!     logical argument CASESENSITIVE (default .FALSE.).
!!!     If you want the comparison to be "truncated" -
!!!     meaning that comparison only takes place up to the
!!!     length of "refstring" - set TRUNC=.TRUE. (default .TRUE.
!!!     unless the WILDCARD option is PRESENT).  An optional
!!!     WILDCARD symbol can be included in the refstring.
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string1, ...stringN
!!!     CHARACTER(LEN=*), INTENT(IN) :: string_list1(:), ...string_listN
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
!!!     LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
!!!     CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD
!!! OUTPUT
!!!     INTEGER(I4B) :: il
!!! SEE ALSO
!!!     CompareString, StringContains, CompareStringFields
!!!***
    CHARACTER(LEN=*), INTENT(IN) :: string1
    CHARACTER(LEN=*), INTENT(IN) :: string_list1(:)
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD

    INTEGER(I4B) :: il

    INTEGER(I4B) :: i

    il = 0
    DO i=1, SIZE(string_list1)
       IF ( CompareString(string1,string_list1(i), &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC,WILDCARD=WILDCARD) ) THEN
          il = i; EXIT
       END IF
    END DO

  END FUNCTION LookupString1

  FUNCTION LookupString2(string1, string2, &
       & string_list1, string_list2, &
       & CASESENSITIVE, TRUNC, WILDCARD ) &
       & RESULT(il)

    CHARACTER(LEN=*), INTENT(IN) :: string1
    CHARACTER(LEN=*), INTENT(IN) :: string2
    CHARACTER(LEN=*), INTENT(IN) :: string_list1(:)
    CHARACTER(LEN=*), INTENT(IN) :: string_list2(:)
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD

    INTEGER(I4B) :: il

    INTEGER(I4B) :: i,n

    n = assert_eq(SIZE(string_list1),SIZE(string_list2), &
         & 'LookupString2 n')

    il = 0
    DO i=1,n
       IF ( CompareString(string1,string_list1(i),     &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC, &
            & WILDCARD=WILDCARD) .AND. &
            & CompareString(string2,string_list2(i),   &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC, &
            & WILDCARD=WILDCARD) ) THEN
          il = i; EXIT
       END IF
    END DO

  END FUNCTION LookupString2

  FUNCTION LookupString3(string1, string2, string3, &
       & string_list1, string_list2, string_list3, &
       & CASESENSITIVE, TRUNC, WILDCARD ) &
       & RESULT(il)

    CHARACTER(LEN=*), INTENT(IN) :: string1
    CHARACTER(LEN=*), INTENT(IN) :: string2
    CHARACTER(LEN=*), INTENT(IN) :: string3
    CHARACTER(LEN=*), INTENT(IN) :: string_list1(:)
    CHARACTER(LEN=*), INTENT(IN) :: string_list2(:)
    CHARACTER(LEN=*), INTENT(IN) :: string_list3(:)
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: CASESENSITIVE
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: TRUNC
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: WILDCARD

    INTEGER(I4B) :: il

    INTEGER(I4B) :: i,n

    n = assert_eq(SIZE(string_list1),SIZE(string_list2),SIZE(string_list3), &
         & 'LookupString3 n')

    il = 0
    DO i=1,n
       IF (   CompareString(string1,string_list1(i),   &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC, &
            & WILDCARD=WILDCARD) .AND. &
            & CompareString(string2,string_list2(i),   &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC, &
            & WILDCARD=WILDCARD) .AND. &
            & CompareString(string3,string_list3(i),   &
            & CASESENSITIVE=CASESENSITIVE,TRUNC=TRUNC, &
            & WILDCARD=WILDCARD) ) THEN
          il = i; EXIT
       END IF
    END DO

  END FUNCTION LookupString3

!!$===========================================================================
!!!****f* StringMod/Levenshtein
!!!
!!! NAME
!!!     Levenshtein -- computes the Levenshtein (edit) distance between 2 strings.
!!!                    strings with smaller distances are more 'similar'
!!! USAGE
!!!     wt = Levenshtein(string1,string2,delCost,insCost,substCost,mods)
!!!        or simply
!!!     wt = Levenshtein(string1,string2)
!!! DESCRIPTION
!!!     This routine implements a modified version of the Levenshtein (edit) distance
!!!     metric used to determine how similar 2 strings are. 
!!!         V. I. Levenshtein {Soviet Physics Doklady 10(8) p707-710, Feb 1966}
!!!     This algorithm allows non-uniform weighting of the deletion, insertion, and
!!!     substitution operations. Be careful when choosing non-uniform values, as poor
!!!     choices could cause certain operations to always happen, or never happen.
!!!     If delCost==insCost==substCost==1.0, wt is the number of string operations performed.    
!!!       Otherwise, the optional integer mods will provide this information. OPERATIONS
!!!       WILL NOT BE RECORDED FOR ANYTHING WITH ZERO COST!!
!!!
!!!     In reality, this function could be expanded to accept taylored deletion, 
!!!     insertion, and substitution cost routines as a function of character. 
!!!     For example, in our molecular database, if we wanted to find thio subtituted
!!!     molecules we could create a substitution function such that subst('O','S') = 0.0
!!!     and all other subst(char1,char2)=WT. Or to find alternate protonation states
!!!     of a molecule, have del('H')==ins('H')==0.0 and all other del(char)==ins(char)==WT.
!!!
!!!     The algorithm as written is O(n*m) for both space and time.
!!!     Additional algorithms/formulations are available at improved computational cost
!!!      (SEE: D. S. Hirschberg {Comm. A.C.M. 18(6) p341-343, 1975}
!!!            D. R. Powell, L. Allison, T. I. Dix {Inf. Proc. Lett. 70 p127-139, 1999}
!!!            L. Allison, {Inf. Proc. Lett., 43(4), pp207-212, 1992})
!!!
!!! NOTES
!!!     Trailing spaces are NOT SIGNIFICANT in this implementation. ie (assuming
!!!     delCost==insCost==substCost==1.0) 
!!!     The strings
!!!        'Blah1' and 'Blah2 '  will have an edit distance of zero. 
!!!     However
!!!        'Blah1' and 'Blah2 x' will have an edit distance of two.
!!! 
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string1
!!!     CHARACTER(LEN=*), INTENT(IN) :: string2
!!!     REAL(SP),INTENT(IN), optional :: delCost   (default 1.0)
!!!     REAL(SP),INTENT(IN), optional :: insCost   (default 1.0)
!!!     REAL(SP),INTENT(IN), optional :: subCost   (default 1.0)
!!! OUTPUT
!!!     REAL(SP) :: wt
!!!     INTEGER(I4B),optional,INTENT(OUT) :: nmods
!!!***
  function Levenshtein(string1, string2, delCost, insCost, subCost, nmods) result(wt)
    CHARACTER(LEN=*), INTENT(IN) :: string1
    CHARACTER(LEN=*), INTENT(IN) :: string2
    REAL(SP), optional :: delCost, insCost, subCost
    INTEGER(I4B),optional,INTENT(OUT) :: nmods ! total modifications (sub/ins/del) 
    REAL(SP) :: wt

    REAL(SP),allocatable:: d(:,:) ! cost matrix
    INTEGER(I4B) :: n,m ! length of string1,string2
    INTEGER(I4B) :: i,j ! iterates through string1,string2
    REAL(SP) :: subst ! tmp substcost
    REAL(SP) :: dCost, iCost, sCost 

    nmods = 0 
    dCost = 1.0_SP
    iCost = 1.0_SP
    sCost = 1.0_SP
    IF(PRESENT(delCost)) THEN
       dCost=ABS(delCost)
    END IF
    IF(PRESENT(insCost)) THEN
       iCost=ABS(insCost)
    END IF
    IF(PRESENT(subCost)) THEN
       sCost=ABS(subCost)
    END IF

    !Step 1
    n = LEN_TRIM(string1)
    m = LEN_TRIM(string2)
    IF(n == 0) THEN
       nmods=m
       wt=m*insCost
       RETURN
    END IF
    IF (m == 0) THEN
       nmods=n
       wt=n*insCost
       RETURN
    END IF
    allocate(d(n+1,m+1))
    d(1,1)=0.0_SP
    
    ! Step 2    
    DO i = 1, n 
       d(i+1,1) = d(i,1) + dCost
    ENDDO

    DO j = 0, m-1
       d(1,j+1) = d(1,j) + iCost
    END DO

    ! Step 3
    DO i = 1, n
      !Step 4
      DO j = 1, m
         !Step 5
         if (string1(i:i) == string2(j:j)) THEN
            subst = 0;
         else
            subst = sCost;
         END if
        !Step 6
        d(i+1,j+1) = MIN( d(i,j)+subst, d(i,j+1)+dCost, d(i+1,j)+iCost)
        if(d(i+1,j+1)/=0.0_SP) THEN
           nmods=nmods+1
        ENDIF
     END DO
  END DO
    
  ! Step 7
  wt=d(n+1,m+1)
  deallocate(d)
END function Levenshtein




FUNCTION CountFields(line) RESULT(n)
!!$===========================================================================
!!!****f* StringMod/CountFields
!!!
!!! NAME
!!!     CountFields (FUNCTION) Counts the number of fields in a line
!!! USAGE
!!!     n = CountFields(string)
!!! DESCRIPTION
!!!     Counts the number of whitespace separated fields in a line
!!!
!!! INPUTS
!!!     CHARACTER(LEN=*), INTENT(IN) :: string
!!! OUTPUT
!!!     INTEGER(I4B) :: n ! FUNCTION RESULT
!!!***
  CHARACTER(LEN=*),INTENT(IN) :: line
  INTEGER(I4B) :: n
  CHARACTER(LEN=LEN(line)) :: tmp
  INTEGER(I4B) :: istat
  CHARACTER(LEN=LEN_LINE) :: field
  tmp   = line
  n     = 0
  istat = 0
  DO
     CALL ReadField_h(tmp,field,SKIP=n,IOSTAT=istat)
     IF ( istat /= 0 ) EXIT
     n=n+1
  END DO
END FUNCTION CountFields

END MODULE StringMod


MODULE filemod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

   INTEGER(I4B), PARAMETER :: MIN_UNIT=8, MAX_UNIT=800


   INTERFACE OpenFile
      MODULE PROCEDURE OpenFile_list
   END INTERFACE

   INTERFACE CloseFile
      MODULE PROCEDURE CloseFile_I4B!, &
   !        CloseFile1, &
    !       CloseFile2, &
       !    CloseFile3, &
     !      CloseFile4, &
        !   CloseFile5
   END INTERFACE

CONTAINS
  FUNCTION AvailableUnitNumber() RESULT(UNIT)
    USE DataTypes
    IMPLICIT NONE
    INTEGER(I4B) :: UNIT
    INTEGER(I4B) :: iunit
    LOGICAL(LGD) :: LOPE
    UNIT=0
    DO iunit=MIN_UNIT,MAX_UNIT
       INQUIRE(UNIT=iunit,OPENED=LOPE)
       IF ( .NOT. LOPE ) THEN
          UNIT=iunit
          EXIT
       END IF
    END DO
  END FUNCTION AvailableUnitNumber

  SUBROUTINE CloseFile_I4B(unit)
    USE DataTypes
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: unit
    CLOSE(unit)
  END SUBROUTINE CloseFile_I4B

  SUBROUTINE OpenFile_list(FILE,UNIT,STAT,FORM,ACC,RECL,POS,ACT, &
       & IOSTAT,FO,NEWUNIT)
    USE DataTypes
    USE StringMod
    IMPLICIT NONE
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FILE
    INTEGER(I4B), OPTIONAL, INTENT(INOUT) :: UNIT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: STAT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FORM
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ACC
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: RECL
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: POS
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ACT
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: IOSTAT
    TYPE(IO), OPTIONAL, INTENT(OUT) :: FO
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: NEWUNIT
    INTEGER(I4B) :: iunit, i1,i2,i3, ierr,err,testUnit
    LOGICAL(LGD) :: isOpened,GetNewUnit
    INTEGER(I4B) :: rl=80
    CHARACTER(LEN=LEN_FILENAME) :: filename,filename2
    CHARACTER(LEN=20) :: status=   'UNKNOWN             '
    CHARACTER(LEN=20) :: FORMAT=   'FORMATTED           '
    CHARACTER(LEN=20) :: access=   'SEQUENTIAL          '
    CHARACTER(LEN=20) :: position= 'ASIS                '
    CHARACTER(LEN=20) :: action=   'READWRITE           '
    status=   'UNKNOWN             '
    FORMAT=   'FORMATTED           '
    access=   'SEQUENTIAL          '
    position= 'ASIS                '
    action=   'READWRITE           '
    GetNewUnit = .FALSE.
    IF ( PRESENT(NEWUNIT) ) THEN
       GetNewUnit = NEWUNIT
    ENDIF
    ! If we are getting a new unit, then get a new unit
    IF (GetNewUnit) THEN
       iunit = AvailableUnitNumber()
       UNIT=iunit
    end IF
    IF ( PRESENT(UNIT) ) THEN 
       ! If we are getting a new unit, then get a new unit
       IF (GetNewUnit) THEN
          iunit = AvailableUnitNumber()
          UNIT=iunit
       else
          IF ( UNIT > 0 ) THEN
             iunit=UNIT
             ! Check the unit number to see if it is opened
             INQUIRE(UNIT=iunit,OPENED=isOpened,NAME=filename2,IOSTAT=err)
             IF(err /= 0) THEN
                CALL error('OpenFile_list: Unit number Inquire Error: err, unit', &
                     & i1=err,i2=iunit)   
             ELSEIF (isOpened) THEN
                CALL error('OpenFile_list: Unit is already connected to a file', &
                     i1= iunit,a1=filename2)
             ENDIF
          ELSE
             iunit = AvailableUnitNumber()
             UNIT=iunit
          END IF
       end IF
    ELSE
       iunit = AvailableUnitNumber()
    END IF
    IF ( PRESENT(FILE) ) THEN
       filename = FILE
    ELSE
       i1 = iunit/100
       i2 = (iunit-i1*100)/10
       i3 = (iunit-i1*100-i2*10)
       WRITE(filename,10)i1,i2,i3
10     FORMAT('ftn90',i1,i1,i1)
    END IF
    IF ( iunit < MIN_UNIT .OR. iunit > MAX_UNIT ) THEN
       CALL error('OpenFile_list: illegal value for iunit!', &
            i1=iunit,WARNLev=5)
    END IF
    ! Check the file name to see if it is opened
    INQUIRE(FILE=filename,OPENED=isOpened,NUMBER=testUnit,IOSTAT=err)
    IF (err /= 0) THEN
       CALL error('OpenFile_list: Inquire Error on file name', a1=filename)     
    ELSEIF (isOpened) THEN
       CALL error('OpenFile_list: Attempting to open file that is already opened', &
            a1=filename)
    ENDIF
    IF ( PRESENT(STAT) ) THEN
       IF ( UPPERCASE(FirstCharacter(STAT)) == 'N' ) status='NEW'
       IF ( UPPERCASE(FirstCharacter(STAT)) == 'O' ) status='OLD'
       IF ( UPPERCASE(FirstCharacter(STAT)) == 'U' ) status='UNKNOWN'
       IF ( UPPERCASE(FirstCharacter(STAT)) == 'R' ) status='REPLACE'
       IF ( UPPERCASE(FirstCharacter(STAT)) == 'S' ) status='SCRATCH'
    END IF
    IF ( PRESENT(FORM) ) THEN
       IF ( UPPERCASE(FirstCharacter(FORM)) == 'F' ) FORMAT='FORMATTED'
       IF ( UPPERCASE(FirstCharacter(FORM)) == 'U' ) FORMAT='UNFORMATTED'
    END IF
    IF ( PRESENT(ACC) ) THEN
       IF ( UPPERCASE(FirstCharacter(ACC)) == 'S' ) access='SEQUENTIAL'
       IF ( UPPERCASE(FirstCharacter(ACC)) == 'D' ) access='DIRECT'
    END IF
    IF ( PRESENT(RECL) ) rl=RECL
    IF ( PRESENT(POS) ) THEN
       IF ( UPPERCASE(FirstCharacter(POS)) == 'R' ) position='REWIND'
       IF ( UPPERCASE(FirstCharacter(POS)) == 'A' ) position='APPEND'
       IF ( UPPERCASE(POS(1:2)) == 'AS' )           position='ASIS'
    END IF
    IF ( PRESENT(ACT) ) THEN
       IF ( UPPERCASE(FirstCharacter(ACT)) == 'R' ) action='READ'
       IF ( UPPERCASE(FirstCharacter(ACT)) == 'W' ) action='WRITE'
       IF ( StringContains(ACT,'RW',CASESENSITIVE=.FALSE.) .OR. &
            & StringContains(ACT,'READWRITE',CASESENSITIVE=.FALSE.) ) THEN
          action='READWRITE'
       END IF
    END IF
    ! Handle scratch files
    if (status == 'SCRATCH') THEN
       IF ( access(1:1) == 'D' ) THEN
          OPEN(UNIT=iunit, STATUS=status, FORM=FORMAT, &
               ACCESS=access, RECL=rl, POSITION=position, ACTION=action, &
               IOSTAT=ierr)
       ELSE
          OPEN(UNIT=iunit, STATUS=status, FORM=FORMAT, &
               ACCESS=access, POSITION=position, ACTION=action, &
               IOSTAT=ierr)
       END IF
    ELSEIF ( access(1:1) == 'D' ) THEN
       OPEN(UNIT=iunit, FILE=filename, STATUS=status, FORM=FORMAT, &
            ACCESS=access, RECL=rl, POSITION=position, ACTION=action, &
            IOSTAT=ierr)
    ELSE
       OPEN(UNIT=iunit, FILE=filename, STATUS=status, FORM=FORMAT, &
            ACCESS=access, POSITION=position, ACTION=action, &
            IOSTAT=ierr)
    END IF
    IF ( PRESENT(IOSTAT) ) THEN
       IOSTAT=ierr
    ELSE IF ( ierr /= 0 ) THEN
       CALL error('OpenFile_list: error opening file to iunit!', &
            i1=iunit,WARNLev=5)
    END IF
    IF ( PRESENT(FO) ) THEN
       FO%unit=iunit
       FO%iostat=ierr
       FO%file=filename
       FO%format=FORMAT
       FO%status=status
       FO%access=access
       FO%recl=rl
       FO%position=position
       FO%action=action
    END IF
  END SUBROUTINE OpenFile_list

END MODULE filemod





MODULE array_io
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE


  INTERFACE write_array
     MODULE PROCEDURE write_array_real_1Dptr, write_array_real_2Dptr
  END INTERFACE
                                                                                                              
  INTERFACE read_array
     MODULE PROCEDURE read_array_real_1Dptr, read_array_real_2Dptr
  END INTERFACE



CONTAINS
  SUBROUTINE write_array_real_2Dptr(array, unit)
    USE StringMod
    REAL(sp), POINTER :: array(:,:)
    INTEGER(i4b) :: unit
    INTEGER(i4b) :: dim1, dim2, i, j, iostat
    CHARACTER(len=8) :: form
    INQUIRE(UNIT=unit,FORM=form, IOSTAT=iostat)
    IF ( iostat > 0 ) THEN
       WRITE(6,*)'Error inquiring about unit number',unit
       WRITE(6,*)'Refusing to write'
       RETURN
    END IF
    IF ( CompareString(form,'form') ) THEN
       ! Formatted I/O
       IF ( ASSOCIATED(array) ) THEN
          dim1 = SIZE(array,1)
          dim2 = SIZE(array,2)
          WRITE(unit,*)dim1,dim2
          DO i = 1, dim2
             DO j = 1, dim1
                WRITE(unit,*)array(j,i)
             END DO
          END DO
       ELSE
          dim1 = 0
          dim2 = 0
          WRITE(unit,*)dim1,dim2
       END IF
    ELSEIF ( CompareString(form,'unform') ) THEN
       ! Unformatted I/O
       IF ( ASSOCIATED(array) ) THEN
          dim1 = SIZE(array,1)
          dim2 = SIZE(array,2)
          WRITE(unit)dim1,dim2
          DO i = 1, dim2
             DO j = 1, dim1
                WRITE(unit)array(j,i)
             END DO
          END DO
       ELSE
          dim1 = 0
          dim2 = 0
          WRITE(unit)dim1,dim2
       END IF
    ELSE
       CALL error('unable to determine file format in write_array_real2Dptr')
    END IF
  END SUBROUTINE write_array_real_2Dptr

  SUBROUTINE read_array_real_1Dptr(array, unit)
    USE DataTypes
    USE StringMod
    IMPLICIT NONE
    REAL(sp), POINTER :: array(:)
    INTEGER(i4b), INTENT(in) :: unit
    INTEGER(i4b) :: dim1, i, iostat
    CHARACTER(len=8) :: form
    INQUIRE(UNIT=unit,FORM=form, IOSTAT=iostat)
    IF ( iostat > 0 ) THEN
       WRITE(6,*)'Error inquiring about unit number',unit
       WRITE(6,*)'Refusing to read'
       RETURN
    END IF
    IF ( ASSOCIATED(array) ) THEN
       DEALLOCATE(array)
    END IF
    IF ( CompareString(form,'form') ) THEN
       ! Formatted I/O
       READ(unit,*)dim1
       IF ( dim1 > 0 ) THEN
          ALLOCATE(array(1:dim1))
       ELSE
          array => null()
          RETURN
       END IF
       DO i = 1, dim1
          READ(unit,*)array(i)
       END DO
    ELSEIF ( CompareString(form,'unform') ) THEN
       ! Unformatted I/O
       READ(unit)dim1
       IF ( dim1 > 0 ) THEN
          ALLOCATE(array(1:dim1))
       ELSE
          array => null()
          RETURN
       END IF
       DO i = 1, dim1
          READ(unit)array(i)
       END DO
    ELSE
       CALL error('unable to determine file format in read_array_real1Dptr')
    END IF
  END SUBROUTINE read_array_real_1Dptr

  SUBROUTINE write_array_real_1Dptr(array, unit)
    USE StringMod
    REAL(sp), POINTER :: array(:)
    INTEGER(i4b) :: unit
    INTEGER(i4b) :: dim1, i, iostat
    CHARACTER(len=8) :: form
    INQUIRE(UNIT=unit,FORM=form, IOSTAT=iostat)
    IF ( iostat > 0 ) THEN
       WRITE(6,*)'Error inquiring about unit number',unit
       WRITE(6,*)'Refusing to write'
       RETURN
    END IF
    IF ( CompareString(form,'form') ) THEN
       ! Formatted I/O
       IF ( ASSOCIATED(array) ) THEN
          dim1 = SIZE(array)
          WRITE(unit,*)dim1
          DO i = 1, dim1
             WRITE(unit,*)array(i)
          END DO
       ELSE
          dim1 = 0
          WRITE(unit,*)dim1
       END IF
    ELSEIF ( CompareString(form,'unform') ) THEN
       ! Unformatted I.O
       IF ( ASSOCIATED(array) ) THEN
          dim1 = SIZE(array)
          WRITE(unit)dim1
          DO i = 1, dim1
             WRITE(unit)array(i)
          END DO
       ELSE
          dim1 = 0
          WRITE(unit)dim1
       END IF
    ELSE
       CALL error('unable to determine file format in write_array_real1Dptr')
    END IF
  END SUBROUTINE write_array_real_1Dptr

  SUBROUTINE read_array_real_2Dptr(array, unit, form, file)
    USE DataTypes
    USE StringMod
    Use FileMod
    IMPLICIT NONE
    REAL(sp), POINTER :: array(:,:)
    INTEGER(i4b), INTENT(in), optional :: unit
    character(len=*), intent(in), optional :: form, file
    INTEGER(i4b) :: dim1, dim2, i, j, iostat, unit_l
    CHARACTER(len=8) :: form_l
    if (present(unit)) then
       INQUIRE(UNIT=unit,FORM=form_l, IOSTAT=iostat)
       IF ( iostat > 0 ) THEN
          WRITE(6,*)'Error inquiring about unit number',unit
          WRITE(6,*)'Refusing to read'
          RETURN
       END IF
       unit_l = unit
    elseif (present(file)) then
       if (present(form)) then
          form_l = form
       else
          form_l = 'unformatted'
       end if
       unit_l = -1
       call Openfile(UNIT=unit_l,FORM=form_l,FILE=file,&
            STAT='old',ACT='read',IOSTAT=iostat)
    else
       write(6,*)'No unit, and no filename in read_array!'
       write(6,*)'refusing to read'
       return
    end if
    IF ( ASSOCIATED(array) ) THEN
       DEALLOCATE(array)
    END IF
    IF ( CompareString(form_l,'form') ) THEN
       ! Formatted I/O
       READ(unit_l,*)dim1,dim2
       IF ( (dim1 > 0) .AND. (dim2 > 0) ) THEN
          ALLOCATE(array(1:dim1,1:dim2))
       ELSE
          array => null()
          RETURN
       END IF
       DO i = 1, dim2
          DO j = 1, dim1
             READ(unit_l,*)array(j,i)
          END DO
       END DO
    ELSEIF ( CompareString(form_l,'unform') ) THEN
       ! Unformatted I/O
       READ(unit_l)dim1,dim2
       IF ( (dim1 > 0) .AND. (dim2 > 0) ) THEN
          ALLOCATE(array(1:dim1,1:dim2))
       ELSE
          array => null()
          RETURN
       END IF
       DO i = 1, dim2
          DO j = 1, dim1
             READ(unit_l)array(j,i)
          END DO
       END DO
    ELSE
       CALL error('unable to determine file format in read_array_real2Dptr')
    END IF
  END SUBROUTINE read_array_real_2Dptr

END MODULE array_io


!doc=updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!   Written by Darrin M. York, Department of Chemistry, Harvard University,  !
!      Cambridge, MA 02138.                                                  !
!   Copyright by Dr. Darrin M. York (1997), all rights reserved.             !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MathMod
!!!****h* YorkLib/MathMod
!!!
!!! NAME 
!!!     MathMod -- Module containing mathematical utilities.
!!! COPYRIGHT
!!!     Darrin M. York
!!!     Department of Chemistry
!!!     Harvard University
!!! AUTHOR
!!!     Darrin M. York
!!! CREATION DATE
!!!     1997
!!! DESCRIPTION
!!!     Module containing mathematical fuctions and operations.  Many routines
!!!     here were taken with little or no modification from:
!!!       "Numerical Recipes in Fortran 90" by W. H. Press, S. A. Teukolsky,
!!!        W. T. Vetterling, and B. P. Flannery.  
!!!     Where appropriate these routines are referenced with the abbreviation
!!!     NR-f90 p. ----.
!!! USES   
!!!     DataTypes, ErrorMod
!!! USED BY
!!!
!!!***

  USE DataTypes
  USE ErrorMod, ONLY: error
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: &
       RealQuadraticRoots, &
       factln, factrl, &
       delta, &
       gammln, gammp, gammq, gser, gcf, &
       pythag, &
       erf, erfc, &
       IdentityMatrix, CROSS_PRODUCT, VectorNorm, &
       NormalVector, OrthogonalizeVector, OrthoNormalizeVector, &
       Trace, AssociatedLegendrePolynomial, safe_exp, safe_log, Pochhammer
!!!****f* MathMod/factln
!!!
!!! NAME
!!!     factln
!!! USAGE
!!!     factln(n)
!!! DESCRIPTION
!!!     Function returns ln(n!) as a REAL number.
!!!     Contains overloads for INTEGER scalar and vector arguments.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: n  
!!! OUTPUTS
!!!     REAL(SP) :: factln
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE factln
     MODULE PROCEDURE factln_s, factln_v
  END INTERFACE
!!!****f* MathMod/factrl
!!!
!!! NAME
!!!     factrl
!!! USAGE
!!!     factln(n)
!!! DESCRIPTION
!!!     Function returns n! as a REAL number.
!!!     Contains overloads for INTEGER scalar and vector arguments.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: n  
!!! OUTPUTS
!!!     REAL(SP) :: factrl 
!!! USES
!!!     DataTypes, UtilitiesMod
!!!***
  INTERFACE factrl
     MODULE PROCEDURE factrl_s, factrl_v
  END INTERFACE
!!!****f* MathMod/delta
!!!
!!! NAME
!!!     delta
!!! USAGE
!!!     delta(l,lp)
!!! DESCRIPTION
!!!     Function returns the Kronecker delta function.
!!!     Contains overloads for INTEGER and REAL scalar and vector arguments,
!!!     and returns a result of the same size and type of its input arguments.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: l,lp 
!!! OUTPUTS
!!!     REAL(SP) :: delta_res(size(l)) 
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE delta
     MODULE PROCEDURE delta_i, delta_iv, delta_r, delta_rv
  END INTERFACE
!!!****f* MathMod/gammaln
!!!
!!! NAME
!!!     gammln
!!! USAGE
!!!     gammln(x) 
!!! DESCRIPTION
!!!     Function returns ln( GammaFunction(x) ) for x > 0.
!!!     Contains overloads for scalar and vector arguments.
!!!     NR-f90 p. 1085.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: xx 
!!! OUTPUTS
!!!     REAL(SP) :: gammln
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE gammln
     MODULE PROCEDURE gammln_s, gammln_v
  END INTERFACE
!!!****f* MathMod/gammp
!!!
!!! NAME
!!!     gammp
!!! USAGE
!!!     gammp(x)
!!! DESCRIPTION
!!!     Function returns the incomplete gamma function P(a,x) for a,x > 0.
!!!     Contains overloads for scalar and vector arguments.
!!!     NR-f90 p. 1089.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a,x
!!! OUTPUTS
!!!     REAL(SP) :: gammp
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE gammp
     MODULE PROCEDURE gammp_s, gammp_v
  END INTERFACE
!!!****f* MathMod/gammq
!!!
!!! NAME
!!!     gammq
!!! USAGE
!!!     gammq(x)
!!! DESCRIPTION
!!!     Function returns the incomplete gamma function Q(a,x) = 1 - P(a,x) 
!!!     for a,x > 0.    Contains overloads for scalar and vector arguments. 
!!!     NR-f90 p. 1090.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a,x
!!! OUTPUTS
!!!     REAL(SP) :: gammq
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE gammq
     MODULE PROCEDURE gammq_s, gammq_v
  END INTERFACE
!!!****f* MathMod/gser
!!!
!!! NAME
!!!     gser
!!! USAGE
!!!     gser(a,x,gln)
!!! DESCRIPTION
!!!     Function returns the incomplete gamma function P(a,x) for a,x > 0 
!!!     evaluated by its series representation as gser.  Also optionally 
!!!     returns ln( GammaFunction(a) ) as gln.
!!!     Contains overloads for scalar and vector arguments. 
!!!     NR-f90 p. 1090.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a,x
!!! OUTPUTS
!!!     REAL(SP) :: gln
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE gser
     MODULE PROCEDURE gser_s, gser_v
  END INTERFACE
!!!****f* MathMod/gcf
!!!
!!! NAME
!!!     gcf
!!! USAGE
!!!     gcf(a,x,gln)
!!! DESCRIPTION
!!!     Function returns the incomplete gamma function Q(a,x) = 1 - P(a,x) 
!!!     for a,x > 0, evaluated by its continued fraction representation as gcf.
!!!     Also optionally returns ln( GammaFunction(a) ) as gln.
!!!     Contains overloads for scalar and vector arguments.
!!!     ITMAX is the maximum allowed number of iterations; EPS is the relative
!!!     accuracy; FPMIN is a number near the smallest representable floating-
!!!     point number.
!!!     NR-f90 p. 1092.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a,x
!!! OUTPUTS
!!!     REAL(SP) :: gln
!!! USES
!!!     DataTypes, UtilitiesMod
!!!*** 
  INTERFACE gcf
     MODULE PROCEDURE gcf_s, gcf_v
  END INTERFACE
!!!****f* MathMod/erf
!!!
!!! NAME
!!!     erf
!!! USAGE
!!!     erf(x)
!!! DESCRIPTION
!!!     Function returns the standard error function for x > 0.
!!!     Contains overloads for scalar and vector arguments.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: x
!!! OUTPUTS
!!!     REAL(SP) :: erf
!!! USES
!!!     DataTypes
!!!*** 
  INTERFACE erf
     MODULE PROCEDURE erf_s, erf_v
  END INTERFACE
!!!****f* MathMod/erfc
!!!
!!! NAME
!!!     erfc
!!! USAGE
!!!     erfc(x)
!!! DESCRIPTION
!!!     Function returns the complimentary error function for x > 0.
!!!     Contains overloads for scalar and vector arguments.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: x
!!! OUTPUTS
!!!     REAL(SP) :: erfc
!!! USES
!!!     DataTypes
!!!*** 
  INTERFACE erfc
     MODULE PROCEDURE erfc_s, erfc_v
  END INTERFACE
!!!****f* MathMod/Trace
!!!
!!! NAME
!!!     Trace
!!! USAGE
!!!     Trace(a,b)
!!! DESCRIPTION
!!!     Function returns the trace of square matrix a, or matrix product a * b.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a, b
!!! OUTPUTS
!!!     REAL(SP) :: TR
!!! USES
!!!     DataTypes, UtilitiesMod
!!!***  
  INTERFACE Trace
     MODULE PROCEDURE Trace_A, Trace_AB
  END INTERFACE

CONTAINS
  SUBROUTINE RealQuadraticRoots(a,b,c,r1,r2)
!!!****f* MathMod/RealQuadraticRoots
!!!
!!! NAME
!!!     RealQuadraticRoots
!!! USAGE
!!!     RealQuadraticRoots(a,b,c,r1,r2)
!!! DESCRIPTION
!!!     Subroutine to compute real roots of a quadratic equation:
!!!     ax*2 +b*x +c = 0, with real roots r1 and r2 
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a, b, c  
!!! OUTPUTS
!!!     REAL(SP), INTENT(OUT) :: r1, r2 
!!! USES
!!!     DataTypes
!!!***           

    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a, b, c
    REAL(SP), INTENT(OUT) :: r1, r2

    REAL(SP) :: x, y

    x = b*b - 4*a*c
    y = 2*a

    IF ( x < 0.0_SP ) THEN
       CALL error('RealQuadraticRoot: no real roots!',WARNLev=3)
       x = 0.0_SP ! zero the imaginary part
    ELSE
       x = SQRT(x)
    END IF

    r1 = (-b + x)/y
    r2 = (-b - x)/y

    ! return lower root as r1
    IF ( r1 > r2 ) THEN
       x = r1
       r1 = r2
       r2 = x
    END IF

  END SUBROUTINE RealQuadraticRoots

  !!$===========================================================================
  FUNCTION factln_s(n)       
    USE DataTypes; USE UtilitiesMod, ONLY : arth,assert
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP) :: factln_s

    INTEGER(I4B), PARAMETER :: TMAX=100
    REAL(SP), DIMENSION(TMAX), SAVE :: a
    LOGICAL(LGD), SAVE :: init=.TRUE.

    IF (init) THEN
       a(1:TMAX)=gammln(arth(1.0_SP,1.0_SP,TMAX))
       init=.FALSE.
    END IF
    CALL assert(n >= 0, 'factln_s arg')
    IF (n < TMAX) THEN
       factln_s=a(n+1)
    ELSE
       factln_s=gammln(n+1.0_SP)
    END IF
  END FUNCTION factln_s

  FUNCTION factln_v(n)
    USE DataTypes; USE UtilitiesMod, ONLY : arth,assert
    IMPLICIT NONE

    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
    REAL(SP), DIMENSION(SIZE(n)) :: factln_v

    LOGICAL(LGD), DIMENSION(SIZE(n)) :: mask
    INTEGER(I4B), PARAMETER :: TMAX=100
    REAL(SP), DIMENSION(TMAX), SAVE :: a
    LOGICAL(LGD), SAVE :: init=.TRUE.

    IF (init) THEN
       a(1:TMAX)=gammln(arth(1.0_SP,1.0_SP,TMAX))
       init=.FALSE.
    END IF
    CALL assert(ALL(n >= 0), 'factln_v arg')
    mask = (n >= TMAX)
    factln_v=UNPACK(gammln(PACK(n,mask)+1.0_SP),mask,0.0_SP)
    WHERE (.NOT. mask) 
       factln_v=a(n+1)
    END WHERE
  END FUNCTION factln_v
  !!$===========================================================================

  FUNCTION factrl_s(n)
    USE DataTypes; USE UtilitiesMod, ONLY : arth,assert,cumprod
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP) :: factrl_s

    INTEGER(I4B), SAVE :: ntop=0
    INTEGER(I4B), PARAMETER :: NMAX=32
    REAL(SP), DIMENSION(NMAX), SAVE :: a

    CALL assert(n >= 0, 'factrl_s arg')
    IF (n < ntop) THEN
       factrl_s=a(n+1)
    ELSE IF (n < NMAX) THEN
       ntop=NMAX
       a(1)=1.0_SP
       a(2:NMAX)=cumprod(arth(1.0_SP,1.0_SP,NMAX-1))
       factrl_s=a(n+1)
    ELSE
       factrl_s=EXP(gammln(n+1.0_SP))
    END IF
  END FUNCTION factrl_s

  FUNCTION factrl_v(n)
    USE DataTypes; USE UtilitiesMod, ONLY : arth,assert,cumprod
    IMPLICIT NONE

    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
    REAL(SP), DIMENSION(SIZE(n)) :: factrl_v

    LOGICAL(LGD), DIMENSION(SIZE(n)) :: mask
    INTEGER(I4B), SAVE :: ntop=0
    INTEGER(I4B), PARAMETER :: NMAX=32
    REAL(SP), DIMENSION(NMAX), SAVE :: a

    CALL assert(ALL(n >= 0), 'factrl_v arg')
    IF (ntop == 0) THEN
       ntop=NMAX
       a(1)=1.0_SP
       a(2:NMAX)=cumprod(arth(1.0_SP,1.0_SP,NMAX-1))
    END IF
    mask = (n >= NMAX)
    factrl_v=UNPACK(EXP(gammln(PACK(n,mask)+1.0_SP)),mask,0.0_SP)
    WHERE (.NOT. mask) 
       factrl_v=a(n+1)
    END WHERE
  END FUNCTION factrl_v
  !!$===========================================================================

  FUNCTION delta_i(l,lp) RESULT(delta_res)
    USE DataTypes
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: l,lp
    INTEGER(I4B) :: delta_res

    IF ( l == lp ) THEN
       delta_res = 1
    ELSE
       delta_res = 0
    END IF

  END FUNCTION delta_i

  FUNCTION delta_iv(l,lp) RESULT(delta_res)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: l(:),lp(:)
    INTEGER(I4B) :: delta_res(SIZE(l))

    INTEGER(I4B) :: n

    n=assert_eq(SIZE(l),SIZE(lp),'delta_iv')
    WHERE ( l == lp )
       delta_res(1:n) = 1
    ELSEWHERE
       delta_res(1:n) = 0
    END WHERE

  END FUNCTION delta_iv

  FUNCTION delta_r(l,lp) RESULT(delta_res)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: l,lp
    REAL(SP) :: delta_res

    IF ( l == lp ) THEN
       delta_res = 1.0_SP
    ELSE
       delta_res = 0.0_SP
    END IF

  END FUNCTION delta_r

  FUNCTION delta_rv(l,lp) RESULT(delta_res)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: l(:),lp(:)
    REAL(SP) :: delta_res(SIZE(l))

    INTEGER(I4B) :: n

    n=assert_eq(SIZE(l),SIZE(lp),'delta_rv')
    WHERE ( l == lp )
       delta_res(1:n) = 1.0_SP
    ELSEWHERE
       delta_res(1:n) = 0.0_SP
    END WHERE

  END FUNCTION delta_rv
  !!$===========================================================================

  FUNCTION gammln_s(xx)
    USE DataTypes; USE UtilitiesMod, ONLY : arth_DP,assert
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: xx
    REAL(SP) :: gammln_s

    REAL(DP) :: tmp,x
    REAL(DP) :: stp = 2.5066282746310005_DP
    REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_DP,&
         -86.50532032941677_DP,24.01409824083091_DP,&
         -1.231739572450155_DP,0.1208650973866179e-2_DP,&
         -0.5395239384953e-5_DP/)

    CALL assert(xx > 0.0, 'gammln_s arg')
    x=xx
    tmp=x+5.5_DP
    tmp=(x+0.5_DP)*LOG(tmp)-tmp
    gammln_s=tmp+LOG(stp*(1.000000000190015_DP+&
         SUM(coef(:)/arth_DP(x+1.0_DP,1.0_DP,SIZE(coef))))/x)
  END FUNCTION gammln_s

  FUNCTION gammln_v(xx)
    USE DataTypes; USE UtilitiesMod, ONLY: assert
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: xx
    REAL(SP), DIMENSION(SIZE(xx)) :: gammln_v

    REAL(DP), DIMENSION(SIZE(xx)) :: ser,tmp,x,y
    REAL(DP) :: stp = 2.5066282746310005_DP
    REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_DP,&
         -86.50532032941677_DP,24.01409824083091_DP,&
         -1.231739572450155_DP,0.1208650973866179e-2_DP,&
         -0.5395239384953e-5_DP/)
    INTEGER(I4B) :: i

    IF (SIZE(xx) == 0) RETURN
    CALL assert(ALL(xx > 0.0), 'gammln_v arg')
    x=xx
    tmp=x+5.5_DP
    tmp=(x+0.5_DP)*LOG(tmp)-tmp
    ser=1.000000000190015_DP
    y=x
    DO i=1,SIZE(coef)
       y=y+1.0_DP
       ser=ser+coef(i)/y
    END DO
    gammln_v=tmp+LOG(stp*ser/x)
  END FUNCTION gammln_v
  !!$===========================================================================

  FUNCTION gammp_s(a,x)
    USE DataTypes; USE UtilitiesMod, ONLY : assert
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,x
    REAL(SP) :: gammp_s

    CALL assert( x >= 0.0,  a > 0.0, 'gammp_s args')
    IF (x<a+1.0_SP) THEN
       gammp_s=gser(a,x)
    ELSE
       gammp_s=1.0_SP-gcf(a,x)
    END IF
  END FUNCTION gammp_s

  FUNCTION gammp_v(a,x)
    USE DataTypes; USE UtilitiesMod, ONLY : assert,assert_eq
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(SIZE(x)) :: gammp_v

    LOGICAL(LGD), DIMENSION(SIZE(x)) :: mask
    INTEGER(I4B) :: ndum

    ndum=assert_eq(SIZE(a),SIZE(x),'gammp_v')
    CALL assert( ALL(x >= 0.0),  ALL(a > 0.0), 'gammp_v args')
    mask = (x<a+1.0_SP)
    gammp_v=MERGE(gser(a,MERGE(x,0.0_SP,mask)), &
         1.0_SP-gcf(a,MERGE(x,0.0_SP,.NOT. mask)),mask)
  END FUNCTION gammp_v
  !!$===========================================================================

  FUNCTION gammq_s(a,x)
    USE DataTypes; USE UtilitiesMod, ONLY : assert
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,x
    REAL(SP) :: gammq_s

    CALL assert( x >= 0.0,  a > 0.0, 'gammq_s args')
    IF (x<a+1.0_SP) THEN
       gammq_s=1.0_SP-gser(a,x)
    ELSE
       gammq_s=gcf(a,x)
    END IF
  END FUNCTION gammq_s

  FUNCTION gammq_v(a,x)
    USE DataTypes; USE UtilitiesMod, ONLY : assert,assert_eq
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(SIZE(a)) :: gammq_v

    LOGICAL(LGD), DIMENSION(SIZE(x)) :: mask
    INTEGER(I4B) :: ndum

    ndum=assert_eq(SIZE(a),SIZE(x),'gammq_v')
    CALL assert( ALL(x >= 0.0),  ALL(a > 0.0), 'gammq_v args')
    mask = (x<a+1.0_SP)
    gammq_v=MERGE(1.0_SP-gser(a,MERGE(x,0.0_SP,mask)), &
         gcf(a,MERGE(x,0.0_SP,.NOT. mask)),mask)
  END FUNCTION gammq_v
  !!$===========================================================================

  FUNCTION gser_s(a,x,gln)
    USE DataTypes; USE UtilitiesMod, ONLY : error
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,x
    REAL(SP), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP) :: gser_s

    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=EPSILON(x)
    INTEGER(I4B) :: n
    REAL(SP) :: ap,del,summ

    IF (x == 0.0) THEN
       gser_s=0.0
       RETURN
    END IF
    ap=a
    summ=1.0_SP/a
    del=summ
    DO n=1,ITMAX
       ap=ap+1.0_SP
       del=del*x/ap
       summ=summ+del
       IF (ABS(del) < ABS(summ)*EPS) EXIT
    END DO
    IF (n > ITMAX) CALL error('a too large, ITMAX too small in gser_s')
    IF (PRESENT(gln)) THEN
       gln=gammln(a)
       gser_s=summ*EXP(-x+a*LOG(x)-gln)
    ELSE
       gser_s=summ*EXP(-x+a*LOG(x)-gammln(a))
    END IF
  END FUNCTION gser_s

  FUNCTION gser_v(a,x,gln)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq,error
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP), DIMENSION(SIZE(a)) :: gser_v

    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=EPSILON(x)
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(SIZE(a)) :: ap,del,summ
    LOGICAL(LGD), DIMENSION(SIZE(a)) :: converged,zero

    n=assert_eq(SIZE(a),SIZE(x),'gser_v')
    zero=(x == 0.0)
    WHERE (zero) gser_v=0.0
    ap=a
    summ=1.0_SP/a
    del=summ
    converged=zero
    DO n=1,ITMAX
       WHERE (.NOT. converged)
          ap=ap+1.0_SP
          del=del*x/ap
          summ=summ+del
          converged = (ABS(del) < ABS(summ)*EPS)
       END WHERE
       IF (ALL(converged)) EXIT
    END DO
    IF (n > ITMAX) CALL error('a too large, ITMAX too small in gser_v')
    IF (PRESENT(gln)) THEN
       IF (SIZE(gln) < SIZE(a)) CALL &
            error('gser: Not enough space for gln')
       gln=gammln(a)
       WHERE (.NOT. zero) 
          gser_v=summ*EXP(-x+a*LOG(x)-gln)
       END WHERE
    ELSE
       WHERE (.NOT. zero) 
          gser_v=summ*EXP(-x+a*LOG(x)-gammln(a))
       END WHERE
    END IF
  END FUNCTION gser_v
  !!$===========================================================================

  FUNCTION gcf_s(a,x,gln)
    USE DataTypes; USE UtilitiesMod, ONLY : error
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,x
    REAL(SP), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP) :: gcf_s

    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=EPSILON(x),FPMIN=TINY(x)/EPS
    INTEGER(I4B) :: i
    REAL(SP) :: an,b,c,d,del,h

    IF (x == 0.0) THEN
       gcf_s=1.0
       RETURN
    END IF
    b=x+1.0_SP-a
    c=1.0_SP/FPMIN
    d=1.0_SP/b
    h=d
    DO i=1,ITMAX
       an=-i*(i-a)
       b=b+2.0_SP
       d=an*d+b
       IF (ABS(d) < FPMIN) d=FPMIN
       c=b+an/c
       IF (ABS(c) < FPMIN) c=FPMIN
       d=1.0_SP/d
       del=d*c
       h=h*del
       IF (ABS(del-1.0_SP) <= EPS) EXIT
    END DO
    IF (i > ITMAX) CALL error('a too large, ITMAX too small in gcf_s')
    IF (PRESENT(gln)) THEN
       gln=gammln(a)
       gcf_s=EXP(-x+a*LOG(x)-gln)*h
    ELSE
       gcf_s=EXP(-x+a*LOG(x)-gammln(a))*h
    END IF
  END FUNCTION gcf_s

  FUNCTION gcf_v(a,x,gln)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq,error
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
    REAL(SP), DIMENSION(SIZE(a)) :: gcf_v

    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(SP), PARAMETER :: EPS=EPSILON(x),FPMIN=TINY(x)/EPS
    INTEGER(I4B) :: i
    REAL(SP), DIMENSION(SIZE(a)) :: an,b,c,d,del,h
    LOGICAL(LGD), DIMENSION(SIZE(a)) :: converged,zero

    i=assert_eq(SIZE(a),SIZE(x),'gcf_v')
    zero=(x == 0.0)
    WHERE (zero)
       gcf_v=1.0
    ELSEWHERE
       b=x+1.0_SP-a
       c=1.0_SP/FPMIN
       d=1.0_SP/b
       h=d
    END WHERE
    converged=zero
    DO i=1,ITMAX
       WHERE (.NOT. converged)
          an=-i*(i-a)
          b=b+2.0_SP
          d=an*d+b
          d=MERGE(FPMIN,d, ABS(d)<FPMIN )
          c=b+an/c
          c=MERGE(FPMIN,c, ABS(c)<FPMIN )
          d=1.0_SP/d
          del=d*c
          h=h*del
          converged = (ABS(del-1.0_SP)<=EPS)
       END WHERE
       IF (ALL(converged)) EXIT
    END DO
    IF (i > ITMAX) CALL error('a too large, ITMAX too small in gcf_v')
    IF (PRESENT(gln)) THEN
       IF (SIZE(gln) < SIZE(a)) CALL &
            error('gser: Not enough space for gln')
       gln=gammln(a)
       WHERE (.NOT. zero) 
          gcf_v=EXP(-x+a*LOG(x)-gln)*h
       END WHERE
    ELSE
       WHERE (.NOT. zero) 
          gcf_v=EXP(-x+a*LOG(x)-gammln(a))*h
       END WHERE
    END IF
  END FUNCTION gcf_v
  !!$===========================================================================

  !!$===========================================================================
  FUNCTION pythag(a,b) RESULT(pythag_res)
!!!****f* MathMod/pythag
!!!
!!! NAME
!!!     pythag
!!! USAGE
!!!     pythag(a,b)
!!! DESCRIPTION
!!!     Computes SQRT(a*a + b*b) without destructive underflow or overflow.
!!!     This subroutine has been adapted from:
!!!        "Numerical Recipes in Fortran 90" by W. H. Press, S. A. Teukolsky,
!!!         W. T. Vetterling, and B. P. Flannery, p. 1029,
!!!     and is explained in the 2nd Ed. of "Numerical Recipes in Fortran" by
!!!     the same authors pp. 62.          
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a,b
!!! OUTPUTS
!!!     REAL(SP) :: pythag_res
!!! USES
!!!     DataTypes
!!!*** 

    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,b
    REAL(SP) :: pythag_res
    REAL(SP) :: absa,absb

    absa=ABS(a)
    absb=ABS(b)
    IF (absa > absb) THEN
       pythag_res=absa*SQRT(1.0_SP+(absb/absa)**2)
    ELSE
       IF (absb == 0.0) THEN
          pythag_res=0.0_SP
       ELSE
          pythag_res=absb*SQRT(1.0_SP+(absa/absb)**2)
       END IF
    END IF

  END FUNCTION pythag
  !!$===========================================================================

  FUNCTION erf_s(x)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: erf_s

    erf_s=gammp(0.5_SP,x**2)
    IF (x < 0.0) erf_s=-erf_s
  END FUNCTION erf_s

  FUNCTION erf_v(x)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(SIZE(x)) :: erf_v

    erf_v=gammp(SPREAD(0.5_SP,1,SIZE(x)),x**2)
    WHERE (x < 0.0) 
       erf_v=-erf_v
    END WHERE
  END FUNCTION erf_v
  !!$===========================================================================
  FUNCTION erfc_s(x)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: erfc_s

    erfc_s=MERGE(1.0_SP+gammp(0.5_SP,x**2),gammq(0.5_SP,x**2), x < 0.0)
  END FUNCTION erfc_s

  FUNCTION erfc_v(x)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(SIZE(x)) :: erfc_v

    LOGICAL(LGD), DIMENSION(SIZE(x)) :: mask

    mask = (x < 0.0)
    erfc_v=MERGE(1.0_SP+gammp(SPREAD(0.5_SP,1,SIZE(x)), &
         MERGE(x,0.0_SP,mask)**2),gammq(SPREAD(0.5_SP,1,SIZE(x)), &
         MERGE(x,0.0_SP,.NOT. mask)**2),mask)
  END FUNCTION erfc_v
  !!$===========================================================================
!!!****f* MathMod/CROSS_PRODUCT
!!!
!!! NAME
!!!     CROSS_PRODUCT
!!! USAGE
!!!     CROSS_PRODUCT(a,b)
!!! DESCRIPTION
!!!     Function returns the vector cross-product of two 3-vectors.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: a(:),b(:)
!!! OUTPUTS
!!!     REAL(SP) :: vp(3)
!!! USES
!!!     DataTypes, UtilitiesMod
!!!***           
  FUNCTION CROSS_PRODUCT(a,b) RESULT(vp)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a(:),b(:)

    REAL(SP) :: vp(3)

    INTEGER(I4B) :: n

    n=assert_eq(SIZE(a),SIZE(b),3,'CROSS_PRODUCT')

    vp = (/ (a(2)*b(3) - b(2)*a(3)), &
         (a(3)*b(1) - b(3)*a(1)), &
         (a(1)*b(2) - b(1)*a(2)) /)

  END FUNCTION CROSS_PRODUCT
  !!$===========================================================================
!!!****f* MathMod/IdentityMatrix
!!!
!!! NAME
!!!     IdentityMatrix
!!! USAGE
!!!     IdentityMatrix(n)
!!! DESCRIPTION
!!!     Function returns the nxn identity matrix.
!!!
!!! INPUTS
!!!     INTEGER(I4B), INTENT(IN) :: n
!!! OUTPUTS
!!!     REAL(SP) :: M(n,n)
!!! USES
!!!     DataTypes, UtilitiesMod
!!!***           
  FUNCTION IdentityMatrix(n) RESULT(M)
    USE DataTypes
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n

    REAL(SP) :: M(n,n)

    INTEGER(I4B) :: i

    M = 0.0_SP
    DO i=1,n
       M(i,i) = 1.0_SP
    END DO

  END FUNCTION IdentityMatrix
  !!$========================================================================
!!!****f* MathMod/VectorNorm
!!!
!!! NAME
!!!    VectorNorm 
!!! USAGE
!!!     VectorNorm(v)
!!! DESCRIPTION
!!!     Function returns the L2 norm of a vector: SQRT( DOT_PRODUCT(v,v) ).
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: v(:)
!!! OUTPUTS
!!!     REAL(SP) :: vn
!!! USES
!!!     DataTypes
!!!***           

  FUNCTION VectorNorm(v) RESULT(vn)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: v(:)
    REAL(SP) :: vn

    vn = SQRT( DOT_PRODUCT(v,v) )

  END FUNCTION VectorNorm
  !!$=========================================================================
!!!****f* MathMod/NormalVector
!!!
!!! NAME
!!!     NormalVector
!!! USAGE
!!!     NormalVector(v)
!!! DESCRIPTION
!!!     Function returns the vector v normalized.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: v(:) 
!!! OUTPUTS
!!!     REAL(SP) :: vp(SIZE(v))
!!! USES
!!!     DataTypes
!!!***           
  FUNCTION NormalVector(v) RESULT(vp)
    USE DataTypes
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: v(:)
    REAL(SP) :: vp(SIZE(v))

    REAL(SP) :: dnorm

    dnorm = SQRT( DOT_PRODUCT(v,v) )
    IF ( dnorm /= 0.0 ) THEN
       dnorm = 1.0_SP / dnorm
    END IF

    vp = v * dnorm

  END FUNCTION NormalVector
  !!$=====================================================================
!!!****f* MathMod/OrthogonalizeVector
!!!
!!! NAME
!!!     OrthogonalizeVector
!!! USAGE
!!!     OrthogonalizeVector(v1,v2)
!!! DESCRIPTION
!!!     Function returns a vector that is v1 orthogonalized to v2.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: v1(:),v2(:) 
!!! OUTPUTS
!!!     REAL(SP) :: vp(SIZE(v1))
!!! USES
!!!     DataTypes, UtilitiesMod 
!!!***           
  FUNCTION OrthogonalizeVector(v1,v2) RESULT(vp)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: v1(:),v2(:)
    REAL(SP) :: vp(SIZE(v1))

    INTEGER(I4B) :: n

    n=assert_eq(SIZE(v1),SIZE(v2),'OrthogonalizeVector')

    vp = NormalVector(v2)
    vp = v1 - DOT_PRODUCT( v1, vp ) * vp

  END FUNCTION OrthogonalizeVector
  !!$=====================================================================
!!!****f* MathMod/OrthoNormalizeVector
!!!
!!! NAME
!!!     OrthoNormalizeVector
!!! USAGE
!!!     OrthoNormalizeVector(v1,v2)
!!! DESCRIPTION
!!!     Function returns a vector that is v1 orthonormalized to v2.
!!!
!!! INPUTS
!!!     REAL(SP), INTENT(IN) :: v1(:),v2(:)
!!! OUTPUTS
!!!     REAL(SP) :: vp(SIZE(v1))
!!! USES
!!!     DataTypes, UtilitiesMod 
!!!***           
  FUNCTION OrthoNormalizeVector(v1,v2) RESULT(vp)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: v1(:),v2(:)
    REAL(SP) :: vp(SIZE(v1))

    INTEGER(I4B) :: n

    n=assert_eq(SIZE(v1),SIZE(v2),'OrthoNormalizeVector')

    vp = NormalVector( OrthogonalizeVector(v1,v2) )

  END FUNCTION OrthoNormalizeVector
  !!$============================================================================

  FUNCTION Trace_A(a) RESULT(TR)                 
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a(:,:)
    REAL(SP) :: TR

    INTEGER(I4B) :: i, n

    n=assert_eq(SIZE(a,1),SIZE(a,2),'Trace_A')

    TR = 0.0_SP
    DO i=1, n
       TR = TR + a(i,i)
    END DO

  END FUNCTION Trace_A

  FUNCTION Trace_AB(a,b) RESULT(TR)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a(:,:), b(:,:)
    REAL(SP) :: TR

    INTEGER(I4B) :: i, j, n

    n=assert_eq(SIZE(a,1),SIZE(a,2),SIZE(b,1),SIZE(b,2),'Trace_AB')

    TR = 0.0_SP
    DO i=1, n
       DO j=1, n
          TR = TR + a(j,i)*b(i,j)
       END DO
    END DO

  END FUNCTION Trace_AB

  !!$===========================================================================

  FUNCTION AssociatedLegendrePolynomial(l,m,x) RESULT(plgndr)

!!!****f* MathMod/AssociatedLegendrePolynomial
!!!
!!! NAME
!!!     AssociatedLegendrePolynomial  (function)
!!! USAGE
!!!     blah = AssociatedLegendrePolynomial(l,m,x)
!!! DESCRIPTION
!!!     Returns the Associated Legendre Polynomial of some l and m
!!!     at some point x.
!!!     This is just 'plgndr' on page 247 of Numerical Recipes
!!!
!!! INPUTS
!!!     integer(i4b) :: l,m
!!!     real(sp)     :: x
!!! OUTPUTS
!!!     REAL(SP) :: AssociatedLegendrePolynomial
!!! USES
!!!     ConstantsMod
!!!*** 
    REAL(sp),INTENT(in) :: x
    INTEGER(i4b),INTENT(in) :: l,m
    REAL(sp) :: plgndr
    INTEGER(i4b) :: i,ll
    REAL(sp) :: fact,pll,pmm,pmmp1,somx2

    IF(m.GT.l) WRITE(6,*)"AssociatedLegendrePolynomial m =",m," > l =",l
    IF(ABS(x).GT.1.0_SP) WRITE(6,*)"AssociatedLegendrePolynomial x =",x," out of bounds"
    IF(m.LT.0) WRITE(6,*)"AssociatedLegendrePolynomial m =",m," < 0"

    pmm = 1.0_SP
    ! Use analytic formula for Pm^m(x)
    IF(m.GT.0) THEN
       somx2 = SQRT( (1.0_SP-x) * (1.0_SP+x) )
       fact  = 1.0_SP
       DO i = 1 , m
          pmm   = -pmm*fact*somx2
          fact  = fact + 2.0_SP
       END DO
    END IF
    IF(l.EQ.m) THEN
       plgndr = pmm
    ELSE
       pmmp1 = x * ( 2.0_SP*m + 1.0_SP ) * pmm
       IF(l.EQ.m+1) THEN
          plgndr = pmmp1
       ELSE
          DO ll = m+2 , l
             pll   = (x*(2.0_SP*ll-1.0_SP)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm   = pmmp1
             pmmp1 = pll
          END DO
          plgndr=pll
       END IF
    END IF

  END FUNCTION AssociatedLegendrePolynomial

!!!****f* MathMod/safe_exp
!!!
!!! NAME
!!!     safe_exp
!!! USAGE
!!!     blah = safe_exp(x)
!!! DESCRIPTION
!!!     An exp that doesn't blow up, no matter what you give it, 
!!!     but comes as close as the computer will let it
!!!
!!! INPUTS
!!!     real(sp)     :: x
!!! OUTPUTS
!!!     REAL(SP) :: y
!!!*** 
  elemental function safe_exp(x) result(y)
    
    real(sp), intent(in) :: x
    
    real(sp) :: y
    
    real(sp) :: psafe, nsafe
    
    ! this is the biggest number we can feed to exp
    psafe = log(huge(1.0_sp))
    
    ! this is the smallest number we can feed to exp
    nsafe = log(tiny(1.0_sp))

    if ( x > psafe ) then
       y = huge(1.0_sp)
    elseif ( x < nsafe ) then
       y = 0.0_sp
    else
       y = exp(x)
    end if

  end function safe_exp

!!!****f* MathMod/safe_log
!!!
!!! NAME
!!!     safe_log
!!! USAGE
!!!     blah = safe_log(x)
!!! DESCRIPTION
!!!     An log that doesn't blow up, no matter what you give it, 
!!!     but comes as close as the computer will let it
!!!     
!!!     for (x < 0) this function will return huge(x)
!!!     for (x == 0) this function will return -huge(x)
!!!
!!! INPUTS
!!!     real(sp)     :: x
!!! OUTPUTS
!!!     REAL(SP) :: y
!!!*** 
  elemental function safe_log(x) result(y)

    real(sp), intent(in) :: x
    
    real(sp) :: y

    if ( x == 0.0_sp ) then
       y = -huge(x)
    elseif ( x < 0 ) then
       y = huge(x)
    else
       y = log(x)
    end if

  end function safe_log

  FUNCTION Pochhammer(a,n) RESULT(b)
!!!****f* MathMod/Pochhammer
!!!
!!! NAME
!!!     Pochhammer
!!! USAGE
!!!     blah = Pochhammer(a,n)
!!! DESCRIPTION
!!!     Computes the Pochhammer symbol (a)_n
!!!     (a)_n = a*(a+1)*...*(a+n-1) = Gamma(a+n)/Gamma(a)
!!!
!!! INPUTS
!!!     INTEGER(I4B),INTENT(IN) :: a,n
!!! OUTPUTS
!!!     REAL(SP) :: b ! (function result)
!!!*** 
    !  b = a*(a+1)*...*(a+n-1)
    INTEGER(I4B),INTENT(IN) :: a,n
    REAL(SP) :: b
    INTEGER(I4B) :: anm1,i

    IF ( a == 0 ) THEN
       IF ( n == 0 ) THEN
          b = 1.0_SP
       ELSE
          b = 0.0_SP
       END IF
    ELSE IF ( n == 0 ) THEN
       b = 1.0_SP
    ELSE IF ( n == 1 ) THEN
       b = a
    ELSE IF ( a < 0 ) THEN
       b = 0.0_SP
    ELSE IF ( n > 0 .AND. n < 15 ) THEN
       b = 1.0_SP
       anm1 = a+n-1
       DO i=a,anm1
          b = b*i
       END DO
    ELSE
       b = safe_exp( gammln(1.0_SP*(a+n)) ) / safe_exp( gammln(1.0_SP*a) )
    END IF

  END FUNCTION Pochhammer


END MODULE MathMod

MODULE nrMod
contains
  SUBROUTINE gauleg(x1,x2,x,w)
    USE DataTypes; USE ConstantsMod
    USE UtilitiesMod, ONLY : arth,assert_eq
    USE ErrorMod, ONLY : error
    IMPLICIT NONE
    
    REAL(SP), INTENT(IN) :: x1,x2
    REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
    
    REAL(DP), PARAMETER :: EPS=3.0e-14_DP
    INTEGER(I4B) :: its,j,m,n
    INTEGER(I4B), PARAMETER :: MAXIT=10
    REAL(DP) :: xl,xm
    REAL(DP), DIMENSION((SIZE(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL(LGD), DIMENSION((SIZE(x)+1)/2) :: unfinished
    
    n=assert_eq(SIZE(x),SIZE(w),'gauleg')
    m=(n+1)/2
    xm=0.5_DP*(x2+x1)
    xl=0.5_DP*(x2-x1)
    z=COS(PI_D*(arth(1,1,m)-0.25_DP)/(n+0.5_DP))
    unfinished=.TRUE.
    DO its=1,MAXIT
       WHERE (unfinished)
          p1=1.0_DP
          p2=0.0_DP
       END WHERE
       DO j=1,n
          WHERE (unfinished)
             p3=p2
             p2=p1
             p1=((2.0_DP*j-1.0_DP)*z*p2-(j-1.0_DP)*p3)/j
          END WHERE
       END DO
       WHERE (unfinished)
          pp=n*(z*p1-p2)/(z*z-1.0_DP)
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > EPS)
       END WHERE
       IF (.NOT. ANY(unfinished)) EXIT
    END DO
    IF (its == MAXIT+1) CALL error('too many iterations in gauleg')
    x(1:m)=xm-xl*z
    x(n:n-m+1:-1)=xm+xl*z
    w(1:m)=2.0_DP*xl/((1.0_DP-z**2)*pp**2)
    w(n:n-m+1:-1)=w(1:m)
  END SUBROUTINE gauleg

end MODULE nrMod


MODULE eigenmod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE





CONTAINS
  SUBROUTINE cholsl(a,p,b,x)
    USE DataTypes; USE UtilitiesMod, ONLY : assert_eq
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
    REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
    REAL(SP), DIMENSION(:), INTENT(OUT) :: x
!$D---------------------------------------------------------------------------
!$DL SUBROUTINE cholsl(a,p)
!$D
!$D  PURPOSE:
!$DL    Solves A x = b for A a REAL, SYMMETRIC, POSITIVE DEFINITE matrix,
!$DL    with Cholesky decomposition computed with  choldc.
!$DL
!$D
!$D  INPUT:
!$D     A      - Cholesky decomposition stored as lower triangle of a
!$D              as returned from choldc.
!$D     p      - Diagonal elements of Cholesky decomposition as returned 
!$D              from choldc.
!$D     b      - b vector (A x = b).
!$D
!$D  OUTPUT:
!$D     x      - Solution of A x = b
!$D---------------------------------------------------------------------------
    INTEGER(I4B) :: i,n
    n=assert_eq((/SIZE(a,1),SIZE(a,2),SIZE(p),SIZE(b),SIZE(x)/),'cholsl')
    DO i=1,n
       x(i)=(b(i)-DOT_PRODUCT(a(i,1:i-1),x(1:i-1)))/p(i)
    END DO
    DO i=n,1,-1
       x(i)=(x(i)-DOT_PRODUCT(a(i+1:n,i),x(i+1:n)))/p(i)
    END DO
  END SUBROUTINE cholsl

END MODULE eigenmod





MODULE angularquadraturemod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

  PRIVATE


  INTEGER(I4B), PARAMETER :: NLebedev=32
  INTEGER(I4B), PARAMETER :: LebedevNpt(1:NLebedev) = (/ &
       6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, &
       974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)
  INTEGER(I4B), PARAMETER :: LebedevLmax(1:NLebedev) = (/ &
       3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, &
       71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131 /)
  INTEGER(I4B), PARAMETER :: LebedevL2Max(1:NLebedev) = (/ &
       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 20, 23, 26, 29, 32, 35, &
       38, 41, 44, 47, 50, 53, 56, 59, 62, 65 /)
  LOGICAL(LGD), PARAMETER :: T=.TRUE.,F=.FALSE.
  LOGICAL(LGD), PARAMETER :: LebedevValidCOSMO(1:NLebedev) = (/ &
       T, T, T, T, T, F, T, T, T, T, T, F, F, T, T, T, T, T, T, T, T, &
       T, T, T, T, T, T, T, T, T, T, T /)
!!! GaussianProduct minimums
  INTEGER(I4B), PARAMETER :: GaussProdNptmin= 8
  INTEGER(I4B), PARAMETER :: GaussProdLmin  = 3
  INTEGER(I4B), PARAMETER :: GaussProdL2min = 1


PUBLIC :: GetNumAngularQuadPts, GetAngularQuadPtRange, AngularQuad

CONTAINS

  Function ValidGaussProdNpts(npt,rup) Result(vaild_npt)
    INTEGER(I4B), INTENT(IN) :: npt
    logical(lgd),optional :: rup
    INTEGER(I4B) :: vaild_npt,m
    logical(lgd) :: round_up
    IF(.NOT. PRESENT(rup))THEN
       round_up =.TRUE.
    ELSE
       round_up = rup
    ENDIF
    IF(round_up) THEN
       m = CEILING(SQRT(REAL(npt,SP)/2.0_SP) )
    ELSE
       m = FLOOR(SQRT(REAL(npt,SP)/2.0_SP) )
    ENDIF
    vaild_npt = 2*m*m
  END Function ValidGaussProdNpts

  SUBROUTINE GetNumAngularQuadPts(type, npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2, COSMO)
    USE DataTypes
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: type
    INTEGER(I4B), INTENT(INOUT) :: npt
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L1
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L2
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMIN
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMAX
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L1
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L2
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: COSMO
    SELECT CASE(type)
    CASE(0)
       ! Lebedev
       CALL GetNumAngularQuadPts_Lebedev(npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2, COSMO)
    CASE(1)
       ! Gaussian Product  
       CALL GetNumAngularQuadPts_GaussProd(npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2)
    CASE DEFAULT
       ! Lebedev
       CALL GetNumAngularQuadPts_Lebedev(npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2, COSMO)
    END SELECT
  END SUBROUTINE GetNumAngularQuadPts

SUBROUTINE GetNumAngularQuadPts_Lebedev(npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2, COSMO )
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE
    INTEGER(I4B), INTENT(INOUT) :: npt
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L1
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L2
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMIN
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMAX
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L1
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L2
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: COSMO
    INTEGER(I4B) :: i, rule, rule_max, rule_min, rule_l,l
    LOGICAL(LGD) :: INPUT_ARGS
    INPUT_ARGS = ( PRESENT(NMIN) .OR. PRESENT(NMAX) .OR. &
         PRESENT(L1) .OR. PRESENT(L2) )
    IF ( INPUT_ARGS .AND. npt <= 0) THEN
       ! This is the smallest number of 'l' to resonably support
       ! It works out to 6pps for Lebedev and 8pps for GaussProd
       ! Thus 6 and 8 are the smallest discritizations possible
       l = 3
       IF ( PRESENT(L1) ) l = MAX(L1,l)
       ! Convert L2 to L1 value (L1= 2*L2 + 1)
       IF ( PRESENT(L2) ) l = MAX((2*L2+1),l)
          ! defaults for min/max pts
          rule_min = 1
          rule_max = NLebedev
          ! scan best l
          rule_l = NLebedev
          DO i=1,NLebedev
             IF ( l <= LebedevLmax(i) ) THEN
                rule_l = i; EXIT
             END IF
          END DO
          ! scan best nmin
          IF ( PRESENT(NMIN) ) THEN
             IF ( NMIN > LebedevNpt(NLebedev) ) THEN
                CALL error('GetNumAngularQuadPts_lebedev: '//&
                     'NMIN exceeds largest possible Lebedev mesh (nmin,nmesh)',&
                     i1=NMIN,i2=LebedevNpt(NLebedev))
             ENDIF
             DO i=1,NLebedev
                IF ( LebedevNpt(i) >= NMIN ) THEN
                   rule_min = i; EXIT
                END IF
             END DO
          ENDIF
          ! scan best nmax
          IF ( PRESENT(NMAX) ) THEN
             DO i=NLebedev,1,-1
                IF ( LebedevNpt(i) <= NMAX ) THEN
                   rule_max = i; EXIT
                END IF
             END DO
          END IF
          ! select minimal npts
          rule = MIN( rule_max, MAX( rule_min, rule_l ) )
     ELSEIF(npt > 0)THEN
          IF ( NMIN > LebedevNpt(NLebedev) ) THEN
             CALL error('GetNumAngularQuadPts_Lebedev: '//&
                  'NMIN exceeds largest possible Lebedev mesh (nmin,nmesh)',&
                  i1=NMIN,i2=LebedevNpt(NLebedev))
          ENDIF
          DO i=1,NLebedev
            IF ( LebedevNpt(i) >= npt ) THEN
               rule = i; EXIT
            END IF
         END DO
     ELSE
        CALL error('GetNumAngularQuadPts_Lebedev: '//&
                     'Npts < 0 and no optional arguments')
     ENDIF
     IF(Present(COSMO))THEN
       IF(COSMO)THEN
         DO
           IF(LebedevValidCosmo(rule)) EXIT
           rule=rule+1
         ENDDO
       ENDIF
     ENDIF
     npt = LebedevNpt(rule)
     IF(PRESENT(npt_L1)) npt_L1 = LebedevLmax(rule)
     IF(PRESENT(npt_L2)) npt_L2 = LebedevL2max(rule)
END SUBROUTINE GetNumAngularQuadPts_Lebedev

SUBROUTINE GetNumAngularQuadPts_GaussProd(npt, npt_L1, npt_L2, NMIN, NMAX, L1, L2 )
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE
    INTEGER(I4B), INTENT(INOUT) :: npt
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L1
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: npt_L2
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMIN
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: NMAX
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L1
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: L2

    INTEGER(I4B) :: tmp_n,npt_max, npt_min, npt_l,l
    LOGICAL(LGD) :: INPUT_ARGS

    INPUT_ARGS = ( PRESENT(NMIN) .OR. PRESENT(NMAX) .OR. &
         PRESENT(L1) .OR. PRESENT(L2) )
    IF ( INPUT_ARGS .AND. npt <= 0) THEN
       ! This is the smallest number of 'l' to resonably support
       ! It works out to 6pps for Lebedev and 8pps for GaussProd
       ! Thus 6 and 8 are the smallest discritizations possible
       l = 3
       IF ( PRESENT(L1) ) l = MAX(L1,l)
       ! Convert L2 to L1 value (L1= 2*L2 + 1)
       IF ( PRESENT(L2) ) l = MAX((2*L2+1),l)
       ! estimate number of points based on L
       tmp_n = CEILING((l+1)**2/2.0_SP)
       ! get next closest valid value for gauss product
       npt_l = ValidGaussProdNpts(tmp_n)
 
       npt_min=GaussProdNptmin
       IF ( PRESENT(NMIN) ) THEN
          npt_min = ValidGaussProdNpts(MAX(ABS(NMIN),npt_min))
       ENDIF
       npt_max = HUGE(I4B)/10
       IF ( PRESENT(NMAX) ) THEN
          npt_max = ValidGaussProdNpts(MIN(ABS(NMAX),npt_max),.FALSE.)
       ENDIF
       npt = MIN( npt_max, MAX( npt_min, npt_l ) )
     ELSEIF(npt > 0) THEN
        npt = ValidGaussProdNpts(MAX(GaussProdNptmin,npt))
     ELSE
        CALL error('GetNumAngularQuadPts_Lebedev: '//&
                     'Npts < 0 and no optional arguments')
     ENDIF
     l =  FLOOR(SQRT(2.0*npt)-1.0)
     IF(PRESENT(npt_L1)) npt_L1 = l
     IF(PRESENT(npt_L2)) npt_L2 = (l-1)/2
     
END SUBROUTINE GetNumAngularQuadPts_GaussProd


SUBROUTINE GetAngularQuadPtRange(type, pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2, COSMO)
    USE DataTypes
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: type
    INTEGER(I4B), pointer :: pts(:)
    INTEGER(I4B), optional, pointer :: pts_L1(:)
    INTEGER(I4B), optional, pointer :: pts_L2(:)
    INTEGER(I4B), INTENT(IN) :: MinVal, MaxVal
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: NPTS, L1, L2
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: COSMO
    SELECT CASE(type)
    CASE(0)
       ! Lebedev
       CALL GetAngularQuadPtRange_Lebedev(pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2, COSMO)
    CASE(1)
       ! Gaussian Product
       CALL GetAngularQuadPtRange_GaussProd(pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2)
    CASE DEFAULT
       ! Lebedev
       CALL GetAngularQuadPtRange_Lebedev(pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2, COSMO)
    END SELECT
  END SUBROUTINE GetAngularQuadPtRange

SUBROUTINE GetAngularQuadPtRange_GaussProd(pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE
    INTEGER(I4B), pointer :: pts(:)
    INTEGER(I4B), optional, pointer :: pts_L1(:)
    INTEGER(I4B), optional, pointer :: pts_L2(:)
    INTEGER(I4B), INTENT(IN) :: MinVal, MaxVal
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: NPTS, L1, L2
    LOGICAL(LGD) :: lnpts, lL1, lL2
    INTEGER(I4B) :: ntrue,i,m1,m2,lmin,lmax,npt_min,npt_max
    lnpts=.FALSE.
    lL1=.FALSE.
    lL2=.FALSE.
    IF(PRESENT(NPTS)) lnpts=NPTS
    IF(PRESENT(L1)) lL1 = L1
    IF(PRESENT(L2)) lL2 = L2
    IF((MinVal < 0) .OR. (MinVal > MaxVal))THEN
       CALL error('GetAngularQuadPtRange_GaussProd: '//&
            'Minval is negative or bigger than MaxVal')
    ENDIF
    ntrue = COUNT((/lnpts,lL1,lL2/))
    IF(ntrue > 1) THEN
       CALL error('GetAngularQuadPtRange_GaussProd: '//&
            'More than 2 of NPTS,L1,L2 specified and true',l1=lnpts,l2=lL1,l3=lL2)
    ELSEIF(ntrue <= 0) THEN
       lnpts = .TRUE.
       ntrue=1
    ENDIF
    ! convert to range of valid npts
    IF(lL1 .OR. lL2)THEN
       IF(lL1)THEN
          lmin =  MAX(GaussProdLmin,MinVal)
          lmax =  MaxVal
       ELSEIF(lL2)THEN
          lmin= 2*MAX(MinVal,GaussProdL2min) + 1
          lmax= 2*MaxVal+1
       ENDIF
       ! get next closest valid value for gauss product
       npt_min = ValidGaussProdNpts(CEILING((lmin+1)**2/2.0_SP))
       npt_max = ValidGaussProdNpts(CEILING((lmax+1)**2/2.0_SP))
    ELSE
       npt_min = ValidGaussProdNpts(MAX(ABS(MinVal),GaussProdNptMin))
       npt_max = ValidGaussProdNpts(MIN(ABS(MaxVal),HUGE(I4B)/10))
    ENDIF
    m1 = CEILING(SQRT(npt_min/2.0_SP) )
    m2 = CEILING(SQRT(npt_max/2.0_SP) )
    IF(associated(pts)) deallocate(pts)
    allocate(pts(m2-m1+1))
    do i=m1,m2
       pts(i-m1+1) = 2*i*i
    enddo
    IF(present(pts_L1))THEN
       IF(associated(pts_L1)) deallocate(pts_L1)
       allocate(pts_L1(m2-m1+1))
       do i=m1,m2
          pts_L1(i-m1+1) = FLOOR(2.0_SP*i-1.0_SP)
       enddo
    ENDIF
    IF(present(pts_L2))THEN
       IF(associated(pts_L2)) deallocate(pts_L2)
       allocate(pts_L2(m2-m1+1))
       do i=m1,m2
          pts_L2(i-m1+1) = (FLOOR(2.0_SP*i-1.0_SP) - 1) / 2
       enddo
    ENDIF
  END SUBROUTINE GetAngularQuadPtRange_GaussProd

SUBROUTINE GetAngularQuadPtRange_Lebedev(pts, MinVal, MaxVal, pts_L1, pts_L2, NPTS, L1, L2, COSMO)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE
    INTEGER(I4B), pointer :: pts(:)
    INTEGER(I4B), optional, pointer :: pts_L1(:)
    INTEGER(I4B), optional, pointer :: pts_L2(:)
    INTEGER(I4B), INTENT(IN) :: MinVal, MaxVal
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: NPTS, L1, L2
    LOGICAL(LGD), OPTIONAL, INTENT(IN) :: COSMO
    LOGICAL(LGD) :: lnpts, lL1, lL2,cosmo_pts
    INTEGER(I4B) :: ntrue,i,m1,m2,lmin,lmax,npt_min,npt_max,nret,idx

    lnpts=.FALSE.
    lL1=.FALSE.
    lL2=.FALSE.

    IF(PRESENT(NPTS)) lnpts=NPTS
    IF(PRESENT(L1)) lL1 = L1
    IF(PRESENT(L2)) lL2 = L2

    IF((MinVal < 0) .OR. (MinVal > MaxVal))THEN
       CALL error('GetAngularQuadPtRange_GaussProd: '//&
            'Minval is negative or bigger than MaxVal')
    ENDIF

    ntrue = COUNT((/lnpts,lL1,lL2/))
    IF(ntrue > 1) THEN
       CALL error('GetAngularQuadPtRange_GaussProd: '//&
            'More than 2 of NPTS,L1,L2 specified and true',l1=lnpts,l2=lL1,l3=lL2)
    ELSEIF(ntrue <= 0) THEN
       lnpts = .TRUE.
       ntrue=1
    ENDIF

    IF(lL1 .OR. lL2)THEN
       IF(lL1)THEN
          lmin =  MAX(LebedevLmax(1),MinVal)
          lmax =  MIN(MaxVal,LebedevLmax(NLebedev))
       ELSEIF(lL2)THEN
          lmin= 2*MAX(MinVal,LebedevL2max(1)) + 1
          lmax= 2*MIN(MaxVal,LebedevL2max(NLebedev))+1
       ENDIF
       DO i=1,NLebedev
          IF ( lmin <= LebedevLmax(i) ) THEN
             m1=i
             EXIT
          END IF
       ENDDO
       DO i=m1,NLebedev
          IF ( lmax <= LebedevLmax(i) ) THEN
             m2=i
             EXIT
          END IF
       ENDDO
    ELSE
       npt_min = MAX(MinVal,LebedevNpt(1))
       npt_max = MIN(MaxVal,LebedevNpt(NLebedev))
       DO i=1,NLebedev
          IF ( npt_min <= LebedevNpt(i) ) THEN
             m1 = i
             EXIT
          END IF
       END DO
       DO i=m1,NLebedev
          IF ( npt_max <= LebedevNpt(i) ) THEN
             m2 = i
             EXIT
          END IF
       END DO 
    ENDIF
    cosmo_pts=.FALSE.
    nret=m2-m1+1 
    IF(PRESENT(COSMO))THEN
      IF(COSMO)THEN
        cosmo_pts=.TRUE.
         ! make sure m2 point is valid
         DO
           IF(LebedevValidCosmo(m2)) EXIT
           m2=m2+1
         ENDDO
        nret = COUNT(LebedevValidCosmo(m1:m2))
      ENDIF
    ENDIF
    IF(associated(pts)) deallocate(pts)
    allocate(pts(nret))
    IF(present(pts_L1))THEN
       IF(associated(pts_L1)) deallocate(pts_L1)
       allocate(pts_L1(nret))
    ENDIF
    IF(present(pts_L2))THEN
       IF(associated(pts_L2)) deallocate(pts_L2)
       allocate(pts_L2(nret))
    ENDIF

    idx=1
    do i=m1,m2
       IF(cosmo_pts .AND. .NOT. LebedevValidCosmo(i)) CYCLE
       pts(idx) = LebedevNpt(i)
       IF(present(pts_L1)) pts_L1(idx) = LebedevLmax(i)
       IF(present(pts_L2)) pts_L2(idx) = LebedevL2max(i)
       idx=idx+1
    enddo

  END SUBROUTINE GetAngularQuadPtRange_Lebedev


  SUBROUTINE AngularQuad(type, npt, pt, wt, alloc)
    USE DataTypes
    USE constantsMod, only : FOUR_PI
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: type, npt
    REAL(SP), POINTER :: pt(:,:), wt(:)
    LOGICAL(LGD), INTENT(IN) :: alloc

    IF ( alloc ) THEN
       IF(associated(pt)) deallocate(pt)
       IF(associated(wt)) deallocate(wt)
       allocate(pt(3,npt))
       allocate(wt(npt))
    END IF

    IF ( SIZE(pt,2) < npt .OR. SIZE(wt) < npt ) THEN
       CALL error('AngularQuad: '//&
            & 'Size of input arrays insufficient for npt', &
            WARNLev=5,i1=npt)
    END IF

    SELECT CASE (type)
    CASE(1)
       ! Gaussian Product
       CALL GaussianProduct(npt, pt, wt)
    CASE DEFAULT
       ! Lebedev quadrature
       SELECT CASE (npt)
       CASE (6)
          CALL Lebedev_6_pt(npt, pt, wt)
       CASE (14)
          CALL Lebedev_14_pt(npt, pt, wt)
       CASE (26)
          CALL Lebedev_26_pt(npt, pt, wt)
       CASE (38)
          CALL Lebedev_38_pt(npt, pt, wt)
       CASE (50)
          CALL Lebedev_50_pt(npt, pt, wt)
       CASE (74)
          CALL Lebedev_74_pt(npt, pt, wt)
       CASE (86)
          CALL Lebedev_86_pt(npt, pt, wt)
       CASE (110)
          CALL Lebedev_110_pt(npt, pt, wt)
       CASE (146)
          CALL Lebedev_146_pt(npt, pt, wt)
       CASE (170)
          CALL Lebedev_170_pt(npt, pt, wt)
       CASE (194)
          CALL Lebedev_194_pt(npt, pt, wt)
       CASE (230)
          CALL Lebedev_230_pt(npt, pt, wt)
       CASE (266)
          CALL Lebedev_266_pt(npt, pt, wt)
       CASE (302)
          CALL Lebedev_302_pt(npt, pt, wt)
       CASE (350)
          CALL Lebedev_350_pt(npt, pt, wt)
       CASE (434)
          CALL Lebedev_434_pt(npt, pt, wt)
       CASE (590)
          CALL Lebedev_590_pt(npt, pt, wt)
       CASE (770)
          CALL Lebedev_770_pt(npt, pt, wt)
       CASE (974)
          CALL Lebedev_974_pt(npt, pt, wt)
       CASE (1202)
          CALL Lebedev_1202_pt(npt, pt, wt)
       CASE (1454)
          CALL Lebedev_1454_pt(npt, pt, wt)
       CASE (1730)
          CALL Lebedev_1730_pt(npt, pt, wt)
       CASE (2030)
          CALL Lebedev_2030_pt(npt, pt, wt)
       CASE (2354)
          CALL Lebedev_2354_pt(npt, pt, wt)
       CASE (2702)
          CALL Lebedev_2702_pt(npt, pt, wt)
       CASE (3074)
          CALL Lebedev_3074_pt(npt, pt, wt)
       CASE (3470)
          CALL Lebedev_3470_pt(npt, pt, wt)
       CASE (3890)
          CALL Lebedev_3890_pt(npt, pt, wt)
       CASE (4334)
          CALL Lebedev_4334_pt(npt, pt, wt)
       CASE (4802)
          CALL Lebedev_4802_pt(npt, pt, wt)
       CASE (5294)
          CALL Lebedev_5294_pt(npt, pt, wt)
       CASE (5810)
          CALL Lebedev_5810_pt(npt, pt, wt)
       CASE DEFAULT
          CALL error('AngularQuad: '//&
               'Lebedev does not have quadrature for these points.',i1=npt)
       END SELECT
    END SELECT
    wt(1:npt) = FOUR_PI * wt(1:npt) / SUM( wt(1:npt) )
  end SUBROUTINE AngularQuad

!===========================================================================
!             Beginning of  GaussianProduct rules
!===========================================================================
  SUBROUTINE GaussianProduct(n, pt, wt)
    USE DataTypes
    USE ConstantsMod, only : TWO_PI
    USE ErrorMod, only : error
    USE nrMod, only : gauleg
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)

    REAL(SP) :: theta(n),wtheta(n),x1,x2,phi,dphi
    INTEGER(I4B) :: m, ind, ntheta, nphi, itheta, iphi

!    m = INT(SQRT( REAL(n/2,SP) ) )
    m = FLOOR(SQRT(REAL(n,SP)/2.0_SP) )
    IF ( 2*m*m /= n ) THEN
       CALL error('GaussianProduct: n must equal 2*m*m  (n,m)', &
            WARNLev=5,i1=n,i2=m)
    END IF

    x1 = -1.0_SP; x2 = 1.0_SP
    ntheta = m; nphi = 2*m
    CALL gauleg(x1,x2,theta(1:ntheta),wtheta(1:ntheta))

    dphi = TWO_PI / nphi
    ind = 0
    DO iphi = 1,nphi
       DO itheta=1,ntheta
          ind = ind + 1
          phi = (iphi-1)*dphi
          pt(1,ind) = SQRT( 1.0_SP - theta(itheta)**2) * COS(phi)
          pt(2,ind) = SQRT( 1.0_SP - theta(itheta)**2) * SIN(phi)
          pt(3,ind) = theta(itheta)
          wt(ind) = wtheta(itheta)*dphi
       END DO
    END DO

  END SUBROUTINE GaussianProduct
!===========================================================================
!             END of GaussianProduct rules
!===========================================================================


!===========================================================================
!             Beginning of Lebedev rules
!===========================================================================
  Subroutine Lebedev_6_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 3, 1, 0, 0, 0, 0, 0, 6 /)
    REAL(SP),    PARAMETER :: pa(1) = (/ 0.1666666666666667E+00_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_6_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_6_pt

  Subroutine Lebedev_14_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 5, 1, 1, 0, 0, 0, 0, 14 /)
    REAL(SP),    PARAMETER :: pa(2) = (/ &
         0.6666666666666667E-01_SP, 0.7500000000000000E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_14_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_14_pt

  Subroutine Lebedev_26_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 7, 1, 1, 1, 0, 0, 0, 26 /)
    REAL(SP),    PARAMETER :: pa(3) = (/ &
         0.4761904761904762E-01_SP, 0.3214285714285714E-01_SP, 0.3809523809523810E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_26_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_26_pt

  Subroutine Lebedev_38_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 9, 1, 1, 0, 0, 1, 0, 38 /)
    REAL(SP),    PARAMETER :: pa(4) = (/ &
         0.9523809523809524E-02_SP, 0.3214285714285714E-01_SP, 0.4597008433809831E+00_SP, &
         0.2857142857142857E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_38_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_38_pt

  Subroutine Lebedev_50_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 11, 1, 1, 1, 1, 0, 0, 50 /)
    REAL(SP),    PARAMETER :: pa(5) = (/ &
         0.1269841269841270E-01_SP, 0.2109375000000000E-01_SP, 0.2257495590828924E-01_SP, &
         0.3015113445777636E+00_SP, 0.2017333553791887E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_50_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_50_pt

  Subroutine Lebedev_74_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 13, 1, 1, 1, 1, 1, 0, 74 /)
    REAL(SP),    PARAMETER :: pa(7) = (/ &
         0.5130671797338464E-03_SP, -0.2958603896103896E-01_SP, 0.1660406956574204E-01_SP, &
         0.4803844614152614E+00_SP, 0.2657620708215946E-01_SP, 0.3207726489807764E+00_SP, &
         0.1652217099371571E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_74_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_74_pt

  Subroutine Lebedev_86_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 15, 1, 1, 0, 2, 1, 0, 86 /)
    REAL(SP),    PARAMETER :: pa(8) = (/ &
         0.1154401154401154E-01_SP, 0.1194390908585628E-01_SP, 0.3696028464541502E+00_SP, &
         0.1111055571060340E-01_SP, 0.6943540066026664E+00_SP, 0.1187650129453714E-01_SP, &
         0.3742430390903412E+00_SP, 0.1181230374690448E-01_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_86_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_86_pt

  Subroutine Lebedev_110_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 17, 1, 1, 0, 3, 1, 0, 110 /)
    REAL(SP),    PARAMETER :: pa(10) = (/ &
         0.3828270494937162E-02_SP, 0.9793737512487512E-02_SP, 0.1851156353447362E+00_SP, &
         0.8211737283191111E-02_SP, 0.6904210483822922E+00_SP, 0.9942814891178103E-02_SP, &
         0.3956894730559419E+00_SP, 0.9595471336070963E-02_SP, 0.4783690288121502E+00_SP, &
         0.9694996361663028E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_110_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_110_pt

  Subroutine Lebedev_146_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 19, 1, 1, 1, 3, 0, 1, 146 /)
    REAL(SP),    PARAMETER :: pa(12) = (/ &
         0.5996313688621381E-03_SP, 0.7210515360144488E-02_SP, 0.7372999718620756E-02_SP, &
         0.6764410400114264E+00_SP, 0.7116355493117555E-02_SP, 0.4174961227965453E+00_SP, &
         0.6753829486314477E-02_SP, 0.1574676672039082E+00_SP, 0.7574394159054034E-02_SP, &
         0.1403553811713183E+00_SP, 0.4493328323269557E+00_SP, 0.6991087353303262E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_146_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_146_pt

  Subroutine Lebedev_170_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 21, 1, 1, 1, 3, 1, 1, 170 /)
    REAL(SP),    PARAMETER :: pa(14) = (/ &
         0.5544842902037365E-02_SP, 0.6383674773515093E-02_SP, 0.6071332770670752E-02_SP, &
         0.2551252621114134E+00_SP, 0.5183387587747790E-02_SP, 0.6743601460362766E+00_SP, &
         0.6317929009813725E-02_SP, 0.4318910696719410E+00_SP, 0.6201670006589077E-02_SP, &
         0.2613931360335988E+00_SP, 0.5477143385137348E-02_SP, 0.4990453161796037E+00_SP, &
         0.1446630744325115E+00_SP, 0.5968383987681156E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_170_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_170_pt

  Subroutine Lebedev_194_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 23, 1, 1, 1, 4, 1, 1, 194 /)
    REAL(SP),    PARAMETER :: pa(16) = (/ &
         0.1782340447244611E-02_SP, 0.5573383178848738E-02_SP, 0.5716905949977102E-02_SP, &
         0.6712973442695226E+00_SP, 0.5608704082587997E-02_SP, 0.2892465627575439E+00_SP, &
         0.5158237711805383E-02_SP, 0.4446933178717437E+00_SP, 0.5518771467273614E-02_SP, &
         0.1299335447650067E+00_SP, 0.4106777028169394E-02_SP, 0.3457702197611283E+00_SP, &
         0.5051846064614808E-02_SP, 0.1590417105383530E+00_SP, 0.8360360154824589E+00_SP, &
         0.5530248916233094E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_194_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_194_pt

  Subroutine Lebedev_230_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 25, 1, 1, 0, 5, 2, 1, 230 /)
    REAL(SP),    PARAMETER :: pa(19) = (/ &
         -0.5522639919727325E-01_SP, 0.4450274607445226E-02_SP, 0.4492044687397611E+00_SP, &
         0.4496841067921404E-02_SP, 0.2520419490210201E+00_SP, 0.5049153450478750E-02_SP, &
         0.6981906658447242E+00_SP, 0.3976408018051883E-02_SP, 0.6587405243460960E+00_SP, &
         0.4401400650381014E-02_SP, 0.4038544050097660E-01_SP, 0.1724544350544401E-01_SP, &
         0.5823842309715585E+00_SP, 0.4231083095357343E-02_SP, 0.3545877390518688E+00_SP, &
         0.5198069864064399E-02_SP, 0.2272181808998187E+00_SP, 0.4864661535886647E+00_SP, &
         0.4695720972568883E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_230_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_230_pt

  Subroutine Lebedev_266_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 27, 1, 1, 1, 5, 1, 2, 266 /)
    REAL(SP),    PARAMETER :: pa(21) = (/ &
         -0.1313769127326952E-02_SP, 0.4186853881700583E-02_SP, -0.2522728704859336E-02_SP, &
         0.7039373391585475E+00_SP, 0.5315167977810885E-02_SP, 0.1012526248572414E+00_SP, &
         0.4047142377086219E-02_SP, 0.4647448726420539E+00_SP, 0.4112482394406990E-02_SP, &
         0.3277420654971629E+00_SP, 0.3595584899758782E-02_SP, 0.6620338663699974E+00_SP, &
         0.4256131351428158E-02_SP, 0.8506508083520399E+00_SP, 0.4229582700647240E-02_SP, &
         0.3233484542692899E+00_SP, 0.1153112011009701E+00_SP, 0.4080914225780505E-02_SP, &
         0.2314790158712601E+00_SP, 0.5244939240922365E+00_SP, 0.4071467593830964E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_266_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_266_pt

  Subroutine Lebedev_302_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 29, 1, 1, 0, 6, 2, 2, 302 /)
    REAL(SP),    PARAMETER :: pa(24) = (/ &
         0.8545911725128148E-03_SP, 0.3599119285025571E-02_SP, 0.3515640345570105E+00_SP, &
         0.3449788424305883E-02_SP, 0.6566329410219612E+00_SP, 0.3604822601419882E-02_SP, &
         0.4729054132581005E+00_SP, 0.3576729661743367E-02_SP, 0.9618308522614784E-01_SP, &
         0.2352101413689164E-02_SP, 0.2219645236294178E+00_SP, 0.3108953122413675E-02_SP, &
         0.7011766416089545E+00_SP, 0.3650045807677255E-02_SP, 0.2644152887060663E+00_SP, &
         0.2982344963171804E-02_SP, 0.5718955891878961E+00_SP, 0.3600820932216460E-02_SP, &
         0.2510034751770465E+00_SP, 0.8000727494073952E+00_SP, 0.3571540554273387E-02_SP, &
         0.1233548532583327E+00_SP, 0.4127724083168531E+00_SP, 0.3392312205006170E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_302_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_302_pt

  Subroutine Lebedev_350_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 31, 1, 1, 0, 6, 2, 3, 350 /)
    REAL(SP),    PARAMETER :: pa(27) = (/ &
         0.3006796749453936E-02_SP, 0.3050627745650771E-02_SP, 0.7068965463912316E+00_SP, &
         0.1621104600288991E-02_SP, 0.4794682625712025E+00_SP, 0.3005701484901752E-02_SP, &
         0.1927533154878019E+00_SP, 0.2990992529653774E-02_SP, 0.6930357961327123E+00_SP, &
         0.2982170644107595E-02_SP, 0.3608302115520091E+00_SP, 0.2721564237310992E-02_SP, &
         0.6498486161496169E+00_SP, 0.3033513795811141E-02_SP, 0.1932945013230339E+00_SP, &
         0.3007949555218533E-02_SP, 0.3800494919899303E+00_SP, 0.2881964603055307E-02_SP, &
         0.2899558825499574E+00_SP, 0.7934537856582316E+00_SP, 0.2958357626535696E-02_SP, &
         0.9684121455103957E-01_SP, 0.8280801506686862E+00_SP, 0.3036020026407088E-02_SP, &
         0.1833434647041659E+00_SP, 0.9074658265305127E+00_SP, 0.2832187403926303E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_350_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_350_pt

  Subroutine Lebedev_434_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 35, 1, 1, 1, 7, 2, 4, 434 /)
    REAL(SP),    PARAMETER :: pa(33) = (/ &
         0.5265897968224436E-03_SP, 0.2512317418927307E-02_SP, 0.2548219972002607E-02_SP, &
         0.6909346307509111E+00_SP, 0.2530403801186355E-02_SP, 0.1774836054609158E+00_SP, &
         0.2014279020918528E-02_SP, 0.4914342637784746E+00_SP, 0.2501725168402936E-02_SP, &
         0.6456664707424256E+00_SP, 0.2513267174597564E-02_SP, 0.2861289010307638E+00_SP, &
         0.2302694782227416E-02_SP, 0.7568084367178018E-01_SP, 0.1462495621594614E-02_SP, &
         0.3927259763368002E+00_SP, 0.2445373437312980E-02_SP, 0.8818132877794288E+00_SP, &
         0.2417442375638981E-02_SP, 0.9776428111182649E+00_SP, 0.1910951282179532E-02_SP, &
         0.2054823696403044E+00_SP, 0.8689460322872412E+00_SP, 0.2416930044324775E-02_SP, &
         0.5905157048925271E+00_SP, 0.7999278543857286E+00_SP, 0.2512236854563495E-02_SP, &
         0.5550152361076807E+00_SP, 0.7717462626915901E+00_SP, 0.2496644054553086E-02_SP, &
         0.9371809858553722E+00_SP, 0.3344363145343455E+00_SP, 0.2236607760437849E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_434_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_434_pt

  Subroutine Lebedev_590_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 41, 1, 1, 0, 9, 3, 6, 590 /)
    REAL(SP),    PARAMETER :: pa(44) = (/ &
         0.3095121295306187E-03_SP, 0.1852379698597489E-02_SP, 0.7040954938227469E+00_SP, &
         0.1871790639277744E-02_SP, 0.6807744066455243E+00_SP, 0.1858812585438317E-02_SP, &
         0.6372546939258752E+00_SP, 0.1852028828296213E-02_SP, 0.5044419707800358E+00_SP, &
         0.1846715956151242E-02_SP, 0.4215761784010967E+00_SP, 0.1818471778162769E-02_SP, &
         0.3317920736472123E+00_SP, 0.1749564657281154E-02_SP, 0.2384736701421887E+00_SP, &
         0.1617210647254411E-02_SP, 0.1459036449157763E+00_SP, 0.1384737234851692E-02_SP, &
         0.6095034115507196E-01_SP, 0.9764331165051050E-03_SP, 0.6116843442009876E+00_SP, &
         0.1857161196774078E-02_SP, 0.3964755348199858E+00_SP, 0.1705153996395864E-02_SP, &
         0.1724782009907724E+00_SP, 0.1300321685886048E-02_SP, 0.5610263808622060E+00_SP, &
         0.3518280927733519E+00_SP, 0.1842866472905286E-02_SP, 0.4742392842551980E+00_SP, &
         0.2634716655937950E+00_SP, 0.1802658934377451E-02_SP, 0.5984126497885380E+00_SP, &
         0.1816640840360209E+00_SP, 0.1849830560443660E-02_SP, 0.3791035407695563E+00_SP, &
         0.1720795225656878E+00_SP, 0.1713904507106709E-02_SP, 0.2778673190586244E+00_SP, &
         0.8213021581932511E-01_SP, 0.1555213603396808E-02_SP, 0.5033564271075117E+00_SP, &
         0.8999205842074875E-01_SP, 0.1802239128008525E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_590_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_590_pt

  Subroutine Lebedev_770_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 47, 1, 1, 1, 10, 3, 9, 770 /)
    REAL(SP),    PARAMETER :: pa(56) = (/ &
         0.2192942088181184E-03_SP, 0.1421940344335877E-02_SP, 0.1436433617319080E-02_SP, &
         0.5087204410502360E-01_SP, 0.6798123511050502E-03_SP, 0.1228198790178831E+00_SP, &
         0.9913184235294912E-03_SP, 0.2026890814408786E+00_SP, 0.1180207833238949E-02_SP, &
         0.2847745156464294E+00_SP, 0.1296599602080921E-02_SP, 0.3656719078978026E+00_SP, &
         0.1365871427428316E-02_SP, 0.4428264886713469E+00_SP, 0.1402988604775325E-02_SP, &
         0.5140619627249735E+00_SP, 0.1418645563595609E-02_SP, 0.6306401219166803E+00_SP, &
         0.1421376741851662E-02_SP, 0.6716883332022612E+00_SP, 0.1423996475490962E-02_SP, &
         0.6979792685336881E+00_SP, 0.1431554042178567E-02_SP, 0.1446865674195309E+00_SP, &
         0.9254401499865368E-03_SP, 0.3390263475411216E+00_SP, 0.1250239995053509E-02_SP, &
         0.5335804651263506E+00_SP, 0.1394365843329230E-02_SP, 0.6944024393349413E-01_SP, &
         0.2355187894242326E+00_SP, 0.1127089094671749E-02_SP, 0.2269004109529460E+00_SP, &
         0.4102182474045730E+00_SP, 0.1345753760910670E-02_SP, 0.8025574607775339E-01_SP, &
         0.6214302417481605E+00_SP, 0.1424957283316783E-02_SP, 0.1467999527896572E+00_SP, &
         0.3245284345717394E+00_SP, 0.1261523341237750E-02_SP, 0.1571507769824727E+00_SP, &
         0.5224482189696630E+00_SP, 0.1392547106052696E-02_SP, 0.2365702993157246E+00_SP, &
         0.6017546634089558E+00_SP, 0.1418761677877656E-02_SP, 0.7714815866765732E-01_SP, &
         0.4346575516141163E+00_SP, 0.1338366684479554E-02_SP, 0.3062936666210730E+00_SP, &
         0.4908826589037616E+00_SP, 0.1393700862676131E-02_SP, 0.3822477379524787E+00_SP, &
         0.5648768149099500E+00_SP, 0.1415914757466932E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_770_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_770_pt

  Subroutine Lebedev_974_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 53, 1, 1, 0, 12, 4, 12, 974 /)
    REAL(SP),    PARAMETER :: pa(70) = (/ &
         0.1438294190527431E-03_SP, 0.1125772288287004E-02_SP, 0.4292963545341347E-01_SP, &
         0.4948029341949241E-03_SP, 0.1051426854086404E+00_SP, 0.7357990109125470E-03_SP, &
         0.1750024867623087E+00_SP, 0.8889132771304384E-03_SP, 0.2477653379650257E+00_SP, &
         0.9888347838921435E-03_SP, 0.3206567123955957E+00_SP, 0.1053299681709471E-02_SP, &
         0.3916520749849983E+00_SP, 0.1092778807014578E-02_SP, 0.4590825874187624E+00_SP, &
         0.1114389394063227E-02_SP, 0.5214563888415861E+00_SP, 0.1123724788051555E-02_SP, &
         0.6253170244654199E+00_SP, 0.1125239325243814E-02_SP, 0.6637926744523170E+00_SP, &
         0.1126153271815905E-02_SP, 0.6910410398498301E+00_SP, 0.1130286931123841E-02_SP, &
         0.7052907007457760E+00_SP, 0.1134986534363955E-02_SP, 0.1236686762657990E+00_SP, &
         0.6823367927109931E-03_SP, 0.2940777114468387E+00_SP, 0.9454158160447096E-03_SP, &
         0.4697753849207649E+00_SP, 0.1074429975385679E-02_SP, 0.6334563241139567E+00_SP, &
         0.1129300086569132E-02_SP, 0.5974048614181342E-01_SP, 0.2029128752777523E+00_SP, &
         0.8436884500901954E-03_SP, 0.1375760408473636E+00_SP, 0.4602621942484054E+00_SP, &
         0.1075255720448885E-02_SP, 0.3391016526336286E+00_SP, 0.5030673999662036E+00_SP, &
         0.1108577236864462E-02_SP, 0.1271675191439820E+00_SP, 0.2817606422442134E+00_SP, &
         0.9566475323783357E-03_SP, 0.2693120740413512E+00_SP, 0.4331561291720157E+00_SP, &
         0.1080663250717391E-02_SP, 0.1419786452601918E+00_SP, 0.6256167358580814E+00_SP, &
         0.1126797131196295E-02_SP, 0.6709284600738255E-01_SP, 0.3798395216859157E+00_SP, &
         0.1022568715358061E-02_SP, 0.7057738183256172E-01_SP, 0.5517505421423520E+00_SP, &
         0.1108960267713108E-02_SP, 0.2783888477882155E+00_SP, 0.6029619156159187E+00_SP, &
         0.1122790653435766E-02_SP, 0.1979578938917407E+00_SP, 0.3589606329589096E+00_SP, &
         0.1032401847117460E-02_SP, 0.2087307061103274E+00_SP, 0.5348666438135476E+00_SP, &
         0.1107249382283854E-02_SP, 0.4055122137872836E+00_SP, 0.5674997546074373E+00_SP, &
         0.1121780048519972E-02_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_974_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_974_pt

  Subroutine Lebedev_1202_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 59, 1, 1, 1, 13, 4, 16, 1202 /)
    REAL(SP),    PARAMETER :: pa(85) = (/ &
         0.1105189233267572E-03_SP, 0.9133159786443561E-03_SP, 0.9205232738090741E-03_SP, &
         0.3712636449657089E-01_SP, 0.3690421898017899E-03_SP, 0.9140060412262223E-01_SP, &
         0.5603990928680660E-03_SP, 0.1531077852469906E+00_SP, 0.6865297629282609E-03_SP, &
         0.2180928891660612E+00_SP, 0.7720338551145630E-03_SP, 0.2839874532200175E+00_SP, &
         0.8301545958894795E-03_SP, 0.3491177600963764E+00_SP, 0.8686692550179628E-03_SP, &
         0.4121431461444309E+00_SP, 0.8927076285846890E-03_SP, 0.4718993627149127E+00_SP, &
         0.9060820238568219E-03_SP, 0.5273145452842337E+00_SP, 0.9119777254940867E-03_SP, &
         0.6209475332444019E+00_SP, 0.9128720138604181E-03_SP, 0.6569722711857291E+00_SP, &
         0.9130714935691735E-03_SP, 0.6841788309070143E+00_SP, 0.9152873784554116E-03_SP, &
         0.7012604330123631E+00_SP, 0.9187436274321654E-03_SP, 0.1072382215478166E+00_SP, &
         0.5176977312965694E-03_SP, 0.2582068959496968E+00_SP, 0.7331143682101417E-03_SP, &
         0.4172752955306717E+00_SP, 0.8463232836379928E-03_SP, 0.5700366911792503E+00_SP, &
         0.9031122694253992E-03_SP, 0.9827986018263947E+00_SP, 0.1771774022615325E+00_SP, &
         0.6485778453163257E-03_SP, 0.9624249230326228E+00_SP, 0.2475716463426288E+00_SP, &
         0.7435030910982369E-03_SP, 0.9402007994128811E+00_SP, 0.3354616289066489E+00_SP, &
         0.7998527891839054E-03_SP, 0.9320822040143202E+00_SP, 0.3173615246611977E+00_SP, &
         0.8101731497468018E-03_SP, 0.9043674199393299E+00_SP, 0.4090268427085357E+00_SP, &
         0.8483389574594331E-03_SP, 0.8912407560074747E+00_SP, 0.3854291150669224E+00_SP, &
         0.8556299257311812E-03_SP, 0.8676435628462708E+00_SP, 0.4932221184851285E+00_SP, &
         0.8803208679738260E-03_SP, 0.8581979986041619E+00_SP, 0.4785320675922435E+00_SP, &
         0.8811048182425720E-03_SP, 0.8396753624049856E+00_SP, 0.4507422593157064E+00_SP, &
         0.8850282341265444E-03_SP, 0.8165288564022188E+00_SP, 0.5632123020762100E+00_SP, &
         0.9021342299040653E-03_SP, 0.8015469370783529E+00_SP, 0.5434303569693900E+00_SP, &
         0.9010091677105086E-03_SP, 0.7773563069070351E+00_SP, 0.5123518486419871E+00_SP, &
         0.9022692938426915E-03_SP, 0.7661621213900394E+00_SP, 0.6394279634749102E+00_SP, &
         0.9158016174693465E-03_SP, 0.7553584143533510E+00_SP, 0.6269805509024392E+00_SP, &
         0.9131578003189435E-03_SP, 0.7344305757559503E+00_SP, 0.6031161693096310E+00_SP, &
         0.9107813579482705E-03_SP, 0.7043837184021765E+00_SP, 0.5693702498468441E+00_SP, &
         0.9105760258970126E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_1202_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_1202_pt

  Subroutine Lebedev_1454_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 65, 1, 1, 0, 15, 5, 20, 1454 /)
    REAL(SP),    PARAMETER :: pa(102) = (/ &
         0.7777160743261247E-04_SP, 0.7557646413004701E-03_SP, 0.3229290663413854E-01_SP, &
         0.2841633806090617E-03_SP, 0.8036733271462222E-01_SP, 0.4374419127053555E-03_SP, &
         0.1354289960531653E+00_SP, 0.5417174740872172E-03_SP, 0.1938963861114426E+00_SP, &
         0.6148000891358593E-03_SP, 0.2537343715011275E+00_SP, 0.6664394485800705E-03_SP, &
         0.3135251434752570E+00_SP, 0.7025039356923220E-03_SP, 0.3721558339375338E+00_SP, &
         0.7268511789249627E-03_SP, 0.4286809575195696E+00_SP, 0.7422637534208629E-03_SP, &
         0.4822510128282994E+00_SP, 0.7509545035841214E-03_SP, 0.5320679333566263E+00_SP, &
         0.7548535057718401E-03_SP, 0.6172998195394274E+00_SP, 0.7554088969774001E-03_SP, &
         0.6510679849127481E+00_SP, 0.7553147174442808E-03_SP, 0.6777315251687360E+00_SP, &
         0.7564767653292297E-03_SP, 0.6963109410648741E+00_SP, 0.7587991808518730E-03_SP, &
         0.7058935009831749E+00_SP, 0.7608261832033027E-03_SP, 0.9955546194091857E+00_SP, &
         0.4021680447874916E-03_SP, 0.9734115901794209E+00_SP, 0.5804871793945964E-03_SP, &
         0.9275693732388626E+00_SP, 0.6792151955945159E-03_SP, 0.8568022422795103E+00_SP, &
         0.7336741211286294E-03_SP, 0.7623495553719372E+00_SP, 0.7581866300989608E-03_SP, &
         0.5707522908892223E+00_SP, 0.4387028039889501E+00_SP, 0.7538257859800743E-03_SP, &
         0.5196463388403083E+00_SP, 0.3858908414762617E+00_SP, 0.7483517247053123E-03_SP, &
         0.4646337531215351E+00_SP, 0.3301937372343854E+00_SP, 0.7371763661112059E-03_SP, &
         0.4063901697557691E+00_SP, 0.2725423573563777E+00_SP, 0.7183448895756934E-03_SP, &
         0.3456329466643087E+00_SP, 0.2139510237495250E+00_SP, 0.6895815529822191E-03_SP, &
         0.2831395121050332E+00_SP, 0.1555922309786647E+00_SP, 0.6480105801792886E-03_SP, &
         0.2197682022925330E+00_SP, 0.9892878979686097E-01_SP, 0.5897558896594636E-03_SP, &
         0.1564696098650355E+00_SP, 0.4598642910675510E-01_SP, 0.5095708849247346E-03_SP, &
         0.6027356673721295E+00_SP, 0.3376625140173426E+00_SP, 0.7536906428909755E-03_SP, &
         0.5496032320255096E+00_SP, 0.2822301309727988E+00_SP, 0.7472505965575118E-03_SP, &
         0.4921707755234567E+00_SP, 0.2248632342592540E+00_SP, 0.7343017132279698E-03_SP, &
         0.4309422998598483E+00_SP, 0.1666224723456479E+00_SP, 0.7130871582177445E-03_SP, &
         0.3664108182313672E+00_SP, 0.1086964901822169E+00_SP, 0.6817022032112776E-03_SP, &
         0.2990189057758436E+00_SP, 0.5251989784120085E-01_SP, 0.6380941145604121E-03_SP, &
         0.6268724013144998E+00_SP, 0.2297523657550023E+00_SP, 0.7550381377920310E-03_SP, &
         0.5707324144834607E+00_SP, 0.1723080607093800E+00_SP, 0.7478646640144802E-03_SP, &
         0.5096360901960365E+00_SP, 0.1140238465390513E+00_SP, 0.7335918720601220E-03_SP, &
         0.4438729938312456E+00_SP, 0.5611522095882537E-01_SP, 0.7110120527658118E-03_SP, &
         0.6419978471082389E+00_SP, 0.1164174423140873E+00_SP, 0.7571363978689501E-03_SP, &
         0.5817218061802611E+00_SP, 0.5797589531445219E-01_SP, 0.7489908329079234E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_1454_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_1454_pt

  Subroutine Lebedev_1730_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 71, 1, 1, 1, 16, 5, 25, 1730 /)
    REAL(SP),    PARAMETER :: pa(120) = (/ &
         0.6309049437420976E-04_SP, 0.6357185073530720E-03_SP, 0.6398287705571748E-03_SP, &
         0.2860923126194662E-01_SP, 0.2221207162188168E-03_SP, 0.7142556767711522E-01_SP, &
         0.3475784022286848E-03_SP, 0.1209199540995559E+00_SP, 0.4350742443589804E-03_SP, &
         0.1738673106594379E+00_SP, 0.4978569136522127E-03_SP, 0.2284645438467734E+00_SP, &
         0.5435036221998053E-03_SP, 0.2834807671701512E+00_SP, 0.5765913388219542E-03_SP, &
         0.3379680145467339E+00_SP, 0.6001200359226003E-03_SP, 0.3911355454819537E+00_SP, &
         0.6162178172717512E-03_SP, 0.4422860353001403E+00_SP, 0.6265218152438485E-03_SP, &
         0.4907781568726057E+00_SP, 0.6323987160974212E-03_SP, 0.5360006153211468E+00_SP, &
         0.6350767851540569E-03_SP, 0.6142105973596603E+00_SP, 0.6354362775297107E-03_SP, &
         0.6459300387977504E+00_SP, 0.6352302462706235E-03_SP, 0.6718056125089225E+00_SP, &
         0.6358117881417972E-03_SP, 0.6910888533186254E+00_SP, 0.6373101590310117E-03_SP, &
         0.7030467416823252E+00_SP, 0.6390428961368665E-03_SP, 0.8354951166354646E-01_SP, &
         0.3186913449946576E-03_SP, 0.2050143009099486E+00_SP, 0.4678028558591711E-03_SP, &
         0.3370208290706637E+00_SP, 0.5538829697598626E-03_SP, 0.4689051484233963E+00_SP, &
         0.6044475907190476E-03_SP, 0.5939400424557334E+00_SP, 0.6313575103509012E-03_SP, &
         0.1394983311832261E+00_SP, 0.4097581162050343E-01_SP, 0.4078626431855630E-03_SP, &
         0.1967999180485014E+00_SP, 0.8851987391293348E-01_SP, 0.4759933057812725E-03_SP, &
         0.2546183732548967E+00_SP, 0.1397680182969819E+00_SP, 0.5268151186413440E-03_SP, &
         0.3121281074713875E+00_SP, 0.1929452542226526E+00_SP, 0.5643048560507316E-03_SP, &
         0.3685981078502492E+00_SP, 0.2467898337061562E+00_SP, 0.5914501076613073E-03_SP, &
         0.4233760321547856E+00_SP, 0.3003104124785409E+00_SP, 0.6104561257874195E-03_SP, &
         0.4758671236059246E+00_SP, 0.3526684328175033E+00_SP, 0.6230252860707806E-03_SP, &
         0.5255178579796463E+00_SP, 0.4031134861145713E+00_SP, 0.6305618761760796E-03_SP, &
         0.5718025633734589E+00_SP, 0.4509426448342351E+00_SP, 0.6343092767597889E-03_SP, &
         0.2686927772723415E+00_SP, 0.4711322502423248E-01_SP, 0.5176268945737826E-03_SP, &
         0.3306006819904809E+00_SP, 0.9784487303942695E-01_SP, 0.5564840313313692E-03_SP, &
         0.3904906850594983E+00_SP, 0.1505395810025273E+00_SP, 0.5856426671038980E-03_SP, &
         0.4479957951904390E+00_SP, 0.2039728156296050E+00_SP, 0.6066386925777091E-03_SP, &
         0.5027076848919780E+00_SP, 0.2571529941121107E+00_SP, 0.6208824962234458E-03_SP, &
         0.5542087392260217E+00_SP, 0.3092191375815670E+00_SP, 0.6296314297822907E-03_SP, &
         0.6020850887375187E+00_SP, 0.3593807506130276E+00_SP, 0.6340423756791859E-03_SP, &
         0.4019851409179594E+00_SP, 0.5063389934378671E-01_SP, 0.5829627677107342E-03_SP, &
         0.4635614567449800E+00_SP, 0.1032422269160612E+00_SP, 0.6048693376081110E-03_SP, &
         0.5215860931591575E+00_SP, 0.1566322094006254E+00_SP, 0.6202362317732461E-03_SP, &
         0.5758202499099271E+00_SP, 0.2098082827491099E+00_SP, 0.6299005328403779E-03_SP, &
         0.6259893683876795E+00_SP, 0.2618824114553391E+00_SP, 0.6347722390609353E-03_SP, &
         0.5313795124811891E+00_SP, 0.5263245019338556E-01_SP, 0.6203778981238834E-03_SP, &
         0.5893317955931995E+00_SP, 0.1061059730982005E+00_SP, 0.6308414671239979E-03_SP, &
         0.6426246321215801E+00_SP, 0.1594171564034221E+00_SP, 0.6362706466959498E-03_SP, &
         0.6511904367376113E+00_SP, 0.5354789536565540E-01_SP, 0.6375414170333233E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_1730_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_1730_pt

  Subroutine Lebedev_2030_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 77, 1, 1, 0, 18, 6, 30, 2030 /)
    REAL(SP),    PARAMETER :: pa(140) = (/ &
         0.4656031899197431E-04_SP, 0.5421549195295507E-03_SP, 0.2540835336814348E-01_SP, &
         0.1778522133346553E-03_SP, 0.6399322800504915E-01_SP, 0.2811325405682796E-03_SP, &
         0.1088269469804125E+00_SP, 0.3548896312631459E-03_SP, 0.1570670798818287E+00_SP, &
         0.4090310897173364E-03_SP, 0.2071163932282514E+00_SP, 0.4493286134169965E-03_SP, &
         0.2578914044450844E+00_SP, 0.4793728447962723E-03_SP, 0.3085687558169623E+00_SP, &
         0.5015415319164265E-03_SP, 0.3584719706267024E+00_SP, 0.5175127372677937E-03_SP, &
         0.4070135594428709E+00_SP, 0.5285522262081019E-03_SP, 0.4536618626222638E+00_SP, &
         0.5356832703713962E-03_SP, 0.4979195686463577E+00_SP, 0.5397914736175170E-03_SP, &
         0.5393075111126999E+00_SP, 0.5416899441599930E-03_SP, 0.6115617676843916E+00_SP, &
         0.5419308476889938E-03_SP, 0.6414308435160159E+00_SP, 0.5416936902030596E-03_SP, &
         0.6664099412721607E+00_SP, 0.5419544338703164E-03_SP, 0.6859161771214913E+00_SP, &
         0.5428983656630975E-03_SP, 0.6993625593503890E+00_SP, 0.5442286500098193E-03_SP, &
         0.7062393387719380E+00_SP, 0.5452250345057301E-03_SP, 0.7479028168349763E-01_SP, &
         0.2568002497728530E-03_SP, 0.1848951153969366E+00_SP, 0.3827211700292145E-03_SP, &
         0.3059529066581305E+00_SP, 0.4579491561917824E-03_SP, 0.4285556101021362E+00_SP, &
         0.5042003969083574E-03_SP, 0.5468758653496526E+00_SP, 0.5312708889976025E-03_SP, &
         0.6565821978343439E+00_SP, 0.5438401790747117E-03_SP, 0.1253901572367117E+00_SP, &
         0.3681917226439641E-01_SP, 0.3316041873197344E-03_SP, 0.1775721510383941E+00_SP, &
         0.7982487607213301E-01_SP, 0.3899113567153771E-03_SP, 0.2305693358216114E+00_SP, &
         0.1264640966592335E+00_SP, 0.4343343327201309E-03_SP, 0.2836502845992063E+00_SP, &
         0.1751585683418957E+00_SP, 0.4679415262318919E-03_SP, 0.3361794746232590E+00_SP, &
         0.2247995907632670E+00_SP, 0.4930847981631031E-03_SP, 0.3875979172264824E+00_SP, &
         0.2745299257422246E+00_SP, 0.5115031867540091E-03_SP, 0.4374019316999074E+00_SP, &
         0.3236373482441118E+00_SP, 0.5245217148457367E-03_SP, 0.4851275843340022E+00_SP, &
         0.3714967859436741E+00_SP, 0.5332041499895321E-03_SP, 0.5303391803806868E+00_SP, &
         0.4175353646321745E+00_SP, 0.5384583126021542E-03_SP, 0.5726197380596287E+00_SP, &
         0.4612084406355461E+00_SP, 0.5411067210798852E-03_SP, 0.2431520732564863E+00_SP, &
         0.4258040133043952E-01_SP, 0.4259797391468714E-03_SP, 0.3002096800895869E+00_SP, &
         0.8869424306722721E-01_SP, 0.4604931368460021E-03_SP, 0.3558554457457432E+00_SP, &
         0.1368811706510655E+00_SP, 0.4871814878255202E-03_SP, 0.4097782537048887E+00_SP, &
         0.1860739985015033E+00_SP, 0.5072242910074885E-03_SP, 0.4616337666067458E+00_SP, &
         0.2354235077395853E+00_SP, 0.5217069845235350E-03_SP, 0.5110707008417874E+00_SP, &
         0.2842074921347011E+00_SP, 0.5315785966280310E-03_SP, 0.5577415286163795E+00_SP, &
         0.3317784414984102E+00_SP, 0.5376833708758905E-03_SP, 0.6013060431366950E+00_SP, &
         0.3775299002040700E+00_SP, 0.5408032092069521E-03_SP, 0.3661596767261781E+00_SP, &
         0.4599367887164592E-01_SP, 0.4842744917904866E-03_SP, 0.4237633153506581E+00_SP, &
         0.9404893773654421E-01_SP, 0.5048926076188130E-03_SP, 0.4786328454658452E+00_SP, &
         0.1431377109091971E+00_SP, 0.5202607980478373E-03_SP, 0.5305702076789774E+00_SP, &
         0.1924186388843570E+00_SP, 0.5309932388325743E-03_SP, 0.5793436224231788E+00_SP, &
         0.2411590944775190E+00_SP, 0.5377419770895208E-03_SP, 0.6247069017094747E+00_SP, &
         0.2886871491583605E+00_SP, 0.5411696331677717E-03_SP, 0.4874315552535204E+00_SP, &
         0.4804978774953206E-01_SP, 0.5197996293282420E-03_SP, 0.5427337322059053E+00_SP, &
         0.9716857199366665E-01_SP, 0.5311120836622945E-03_SP, 0.5943493747246700E+00_SP, &
         0.1465205839795055E+00_SP, 0.5384309319956951E-03_SP, 0.6421314033564943E+00_SP, &
         0.1953579449803574E+00_SP, 0.5421859504051886E-03_SP, 0.6020628374713980E+00_SP, &
         0.4916375015738108E-01_SP, 0.5390948355046314E-03_SP, 0.6529222529856881E+00_SP, &
         0.9861621540127005E-01_SP, 0.5433312705027845E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_2030_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_2030_pt

  Subroutine Lebedev_2354_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 83, 1, 1, 1, 19, 6, 36, 2354 /)
    REAL(SP),    PARAMETER :: pa(161) = (/ &
         0.3922616270665292E-04_SP, 0.4678202801282136E-03_SP, 0.4703831750854424E-03_SP, &
         0.2290024646530589E-01_SP, 0.1437832228979900E-03_SP, 0.5779086652271284E-01_SP, &
         0.2303572493577644E-03_SP, 0.9863103576375984E-01_SP, 0.2933110752447454E-03_SP, &
         0.1428155792982185E+00_SP, 0.3402905998359838E-03_SP, 0.1888978116601463E+00_SP, &
         0.3759138466870372E-03_SP, 0.2359091682970210E+00_SP, 0.4030638447899798E-03_SP, &
         0.2831228833706171E+00_SP, 0.4236591432242211E-03_SP, 0.3299495857966693E+00_SP, &
         0.4390522656946746E-03_SP, 0.3758840802660796E+00_SP, 0.4502523466626247E-03_SP, &
         0.4204751831009480E+00_SP, 0.4580577727783541E-03_SP, 0.4633068518751051E+00_SP, &
         0.4631391616615899E-03_SP, 0.5039849474507313E+00_SP, 0.4660928953698676E-03_SP, &
         0.5421265793440747E+00_SP, 0.4674751807936953E-03_SP, 0.6092660230557310E+00_SP, &
         0.4676414903932920E-03_SP, 0.6374654204984869E+00_SP, 0.4674086492347870E-03_SP, &
         0.6615136472609892E+00_SP, 0.4674928539483207E-03_SP, 0.6809487285958127E+00_SP, &
         0.4680748979686447E-03_SP, 0.6952980021665196E+00_SP, 0.4690449806389040E-03_SP, &
         0.7041245497695400E+00_SP, 0.4699877075860818E-03_SP, 0.6744033088306065E-01_SP, &
         0.2099942281069176E-03_SP, 0.1678684485334166E+00_SP, 0.3172269150712804E-03_SP, &
         0.2793559049539613E+00_SP, 0.3832051358546523E-03_SP, 0.3935264218057639E+00_SP, &
         0.4252193818146985E-03_SP, 0.5052629268232558E+00_SP, 0.4513807963755000E-03_SP, &
         0.6107905315437531E+00_SP, 0.4657797469114178E-03_SP, 0.1135081039843524E+00_SP, &
         0.3331954884662588E-01_SP, 0.2733362800522836E-03_SP, 0.1612866626099378E+00_SP, &
         0.7247167465436538E-01_SP, 0.3235485368463559E-03_SP, 0.2100786550168205E+00_SP, &
         0.1151539110849745E+00_SP, 0.3624908726013453E-03_SP, 0.2592282009459942E+00_SP, &
         0.1599491097143677E+00_SP, 0.3925540070712828E-03_SP, 0.3081740561320203E+00_SP, &
         0.2058699956028027E+00_SP, 0.4156129781116235E-03_SP, 0.3564289781578164E+00_SP, &
         0.2521624953502911E+00_SP, 0.4330644984623263E-03_SP, 0.4035587288240703E+00_SP, &
         0.2982090785797674E+00_SP, 0.4459677725921312E-03_SP, 0.4491671196373903E+00_SP, &
         0.3434762087235733E+00_SP, 0.4551593004456795E-03_SP, 0.4928854782917489E+00_SP, &
         0.3874831357203437E+00_SP, 0.4613341462749918E-03_SP, 0.5343646791958988E+00_SP, &
         0.4297814821746926E+00_SP, 0.4651019618269806E-03_SP, 0.5732683216530990E+00_SP, &
         0.4699402260943537E+00_SP, 0.4670249536100625E-03_SP, 0.2214131583218986E+00_SP, &
         0.3873602040643895E-01_SP, 0.3549555576441708E-03_SP, 0.2741796504750071E+00_SP, &
         0.8089496256902013E-01_SP, 0.3856108245249010E-03_SP, 0.3259797439149485E+00_SP, &
         0.1251732177620872E+00_SP, 0.4098622845756882E-03_SP, 0.3765441148826891E+00_SP, &
         0.1706260286403185E+00_SP, 0.4286328604268950E-03_SP, 0.4255773574530558E+00_SP, &
         0.2165115147300408E+00_SP, 0.4427802198993945E-03_SP, 0.4727795117058430E+00_SP, &
         0.2622089812225259E+00_SP, 0.4530473511488561E-03_SP, 0.5178546895819012E+00_SP, &
         0.3071721431296201E+00_SP, 0.4600805475703138E-03_SP, 0.5605141192097460E+00_SP, &
         0.3508998998801138E+00_SP, 0.4644599059958017E-03_SP, 0.6004763319352512E+00_SP, &
         0.3929160876166931E+00_SP, 0.4667274455712508E-03_SP, 0.3352842634946949E+00_SP, &
         0.4202563457288019E-01_SP, 0.4069360518020356E-03_SP, 0.3891971629814670E+00_SP, &
         0.8614309758870850E-01_SP, 0.4260442819919195E-03_SP, 0.4409875565542281E+00_SP, &
         0.1314500879380001E+00_SP, 0.4408678508029063E-03_SP, 0.4904893058592484E+00_SP, &
         0.1772189657383859E+00_SP, 0.4518748115548597E-03_SP, 0.5375056138769549E+00_SP, &
         0.2228277110050294E+00_SP, 0.4595564875375116E-03_SP, 0.5818255708669969E+00_SP, &
         0.2677179935014386E+00_SP, 0.4643988774315846E-03_SP, 0.6232334858144959E+00_SP, &
         0.3113675035544165E+00_SP, 0.4668827491646946E-03_SP, 0.4489485354492058E+00_SP, &
         0.4409162378368174E-01_SP, 0.4400541823741973E-03_SP, 0.5015136875933150E+00_SP, &
         0.8939009917748489E-01_SP, 0.4514512890193797E-03_SP, 0.5511300550512623E+00_SP, &
         0.1351806029383365E+00_SP, 0.4596198627347549E-03_SP, 0.5976720409858000E+00_SP, &
         0.1808370355053196E+00_SP, 0.4648659016801781E-03_SP, 0.6409956378989354E+00_SP, &
         0.2257852192301602E+00_SP, 0.4675502017157673E-03_SP, 0.5581222330827514E+00_SP, &
         0.4532173421637160E-01_SP, 0.4598494476455523E-03_SP, 0.6074705984161695E+00_SP, &
         0.9117488031840314E-01_SP, 0.4654916955152048E-03_SP, 0.6532272537379033E+00_SP, &
         0.1369294213140155E+00_SP, 0.4684709779505137E-03_SP, 0.6594761494500487E+00_SP, &
         0.4589901487275583E-01_SP, 0.4691445539106986E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_2354_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_2354_pt

  Subroutine Lebedev_2702_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 89, 1, 1, 0, 21, 7, 42, 2702 /)
    REAL(SP),    PARAMETER :: pa(184) = (/ &
         0.2998675149888161E-04_SP, 0.4077860529495355E-03_SP, 0.2065562538818703E-01_SP, &
         0.1185349192520667E-03_SP, 0.5250918173022379E-01_SP, 0.1913408643425751E-03_SP, &
         0.8993480082038376E-01_SP, 0.2452886577209897E-03_SP, 0.1306023924436019E+00_SP, &
         0.2862408183288702E-03_SP, 0.1732060388531418E+00_SP, 0.3178032258257357E-03_SP, &
         0.2168727084820249E+00_SP, 0.3422945667633690E-03_SP, 0.2609528309173586E+00_SP, &
         0.3612790520235922E-03_SP, 0.3049252927938952E+00_SP, 0.3758638229818521E-03_SP, &
         0.3483484138084404E+00_SP, 0.3868711798859953E-03_SP, 0.3908321549106406E+00_SP, &
         0.3949429933189938E-03_SP, 0.4320210071894814E+00_SP, 0.4006068107541156E-03_SP, &
         0.4715824795890053E+00_SP, 0.4043192149672723E-03_SP, 0.5091984794078453E+00_SP, &
         0.4064947495808078E-03_SP, 0.5445580145650803E+00_SP, 0.4075245619813152E-03_SP, &
         0.6072575796841768E+00_SP, 0.4076423540893566E-03_SP, 0.6339484505755803E+00_SP, &
         0.4074280862251555E-03_SP, 0.6570718257486958E+00_SP, 0.4074163756012244E-03_SP, &
         0.6762557330090709E+00_SP, 0.4077647795071246E-03_SP, 0.6911161696923790E+00_SP, &
         0.4084517552782530E-03_SP, 0.7012841911659961E+00_SP, 0.4092468459224052E-03_SP, &
         0.7064559272410020E+00_SP, 0.4097872687240906E-03_SP, 0.6123554989894765E-01_SP, &
         0.1738986811745028E-03_SP, 0.1533070348312393E+00_SP, 0.2659616045280191E-03_SP, &
         0.2563902605244206E+00_SP, 0.3240596008171533E-03_SP, 0.3629346991663361E+00_SP, &
         0.3621195964432943E-03_SP, 0.4683949968987538E+00_SP, 0.3868838330760539E-03_SP, &
         0.5694479240657952E+00_SP, 0.4018911532693111E-03_SP, 0.6634465430993955E+00_SP, &
         0.4089929432983252E-03_SP, 0.1033958573552305E+00_SP, 0.3034544009063584E-01_SP, &
         0.2279907527706409E-03_SP, 0.1473521412414395E+00_SP, 0.6618803044247135E-01_SP, &
         0.2715205490578897E-03_SP, 0.1924552158705967E+00_SP, 0.1054431128987715E+00_SP, &
         0.3057917896703976E-03_SP, 0.2381094362890328E+00_SP, 0.1468263551238858E+00_SP, &
         0.3326913052452555E-03_SP, 0.2838121707936760E+00_SP, 0.1894486108187886E+00_SP, &
         0.3537334711890037E-03_SP, 0.3291323133373415E+00_SP, 0.2326374238761579E+00_SP, &
         0.3700567500783129E-03_SP, 0.3736896978741460E+00_SP, 0.2758485808485768E+00_SP, &
         0.3825245372589122E-03_SP, 0.4171406040760013E+00_SP, 0.3186179331996921E+00_SP, &
         0.3918125171518296E-03_SP, 0.4591677985256915E+00_SP, 0.3605329796303794E+00_SP, &
         0.3984720419937579E-03_SP, 0.4994733831718418E+00_SP, 0.4012147253586509E+00_SP, &
         0.4029746003338211E-03_SP, 0.5377731830445096E+00_SP, 0.4403050025570692E+00_SP, &
         0.4057428632156627E-03_SP, 0.5737917830001331E+00_SP, 0.4774565904277483E+00_SP, &
         0.4071719274114857E-03_SP, 0.2027323586271389E+00_SP, 0.3544122504976147E-01_SP, &
         0.2990236950664119E-03_SP, 0.2516942375187273E+00_SP, 0.7418304388646328E-01_SP, &
         0.3262951734212878E-03_SP, 0.3000227995257181E+00_SP, 0.1150502745727186E+00_SP, &
         0.3482634608242413E-03_SP, 0.3474806691046342E+00_SP, 0.1571963371209364E+00_SP, &
         0.3656596681700892E-03_SP, 0.3938103180359209E+00_SP, 0.1999631877247100E+00_SP, &
         0.3791740467794218E-03_SP, 0.4387519590455703E+00_SP, 0.2428073457846535E+00_SP, &
         0.3894034450156905E-03_SP, 0.4820503960077787E+00_SP, 0.2852575132906155E+00_SP, &
         0.3968600245508371E-03_SP, 0.5234573778475101E+00_SP, 0.3268884208674639E+00_SP, &
         0.4019931351420050E-03_SP, 0.5627318647235282E+00_SP, 0.3673033321675939E+00_SP, &
         0.4052108801278599E-03_SP, 0.5996390607156954E+00_SP, 0.4061211551830290E+00_SP, &
         0.4068978613940934E-03_SP, 0.3084780753791947E+00_SP, 0.3860125523100059E-01_SP, &
         0.3454275351319704E-03_SP, 0.3589988275920223E+00_SP, 0.7928938987104867E-01_SP, &
         0.3629963537007920E-03_SP, 0.4078628415881973E+00_SP, 0.1212614643030087E+00_SP, &
         0.3770187233889873E-03_SP, 0.4549287258889735E+00_SP, 0.1638770827382693E+00_SP, &
         0.3878608613694378E-03_SP, 0.5000278512957279E+00_SP, 0.2065965798260176E+00_SP, &
         0.3959065270221274E-03_SP, 0.5429785044928199E+00_SP, 0.2489436378852235E+00_SP, &
         0.4015286975463570E-03_SP, 0.5835939850491711E+00_SP, 0.2904811368946891E+00_SP, &
         0.4050866785614717E-03_SP, 0.6216870353444856E+00_SP, 0.3307941957666609E+00_SP, &
         0.4069320185051913E-03_SP, 0.4151104662709091E+00_SP, 0.4064829146052554E-01_SP, &
         0.3760120964062763E-03_SP, 0.4649804275009218E+00_SP, 0.8258424547294755E-01_SP, &
         0.3870969564418064E-03_SP, 0.5124695757009662E+00_SP, 0.1251841962027289E+00_SP, &
         0.3955287790534055E-03_SP, 0.5574711100606224E+00_SP, 0.1679107505976331E+00_SP, &
         0.4015361911302668E-03_SP, 0.5998597333287227E+00_SP, 0.2102805057358715E+00_SP, &
         0.4053836986719548E-03_SP, 0.6395007148516600E+00_SP, 0.2518418087774107E+00_SP, &
         0.4073578673299117E-03_SP, 0.5188456224746252E+00_SP, 0.4194321676077518E-01_SP, &
         0.3954628379231406E-03_SP, 0.5664190707942778E+00_SP, 0.8457661551921499E-01_SP, &
         0.4017645508847530E-03_SP, 0.6110464353283153E+00_SP, 0.1273652932519396E+00_SP, &
         0.4059030348651293E-03_SP, 0.6526430302051563E+00_SP, 0.1698173239076354E+00_SP, &
         0.4080565809484880E-03_SP, 0.6167551880377548E+00_SP, 0.4266398851548864E-01_SP, &
         0.4063018753664651E-03_SP, 0.6607195418355383E+00_SP, 0.8551925814238349E-01_SP, &
         0.4087191292799671E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_2702_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_2702_pt

  Subroutine Lebedev_3074_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 95, 1, 1, 1, 22, 7, 49, 3074 /)
    REAL(SP),    PARAMETER :: pa(208) = (/ &
         0.2599095953754734E-04_SP, 0.3586067974412447E-03_SP, 0.3603134089687541E-03_SP, &
         0.1886108518723392E-01_SP, 0.9831528474385880E-04_SP, 0.4800217244625303E-01_SP, &
         0.1605023107954450E-03_SP, 0.8244922058397242E-01_SP, 0.2072200131464099E-03_SP, &
         0.1200408362484023E+00_SP, 0.2431297618814187E-03_SP, 0.1595773530809965E+00_SP, &
         0.2711819064496707E-03_SP, 0.2002635973434064E+00_SP, 0.2932762038321116E-03_SP, &
         0.2415127590139982E+00_SP, 0.3107032514197368E-03_SP, 0.2828584158458477E+00_SP, &
         0.3243808058921213E-03_SP, 0.3239091015338138E+00_SP, 0.3349899091374030E-03_SP, &
         0.3643225097962194E+00_SP, 0.3430580688505218E-03_SP, 0.4037897083691802E+00_SP, &
         0.3490124109290343E-03_SP, 0.4420247515194127E+00_SP, 0.3532148948561955E-03_SP, &
         0.4787572538464938E+00_SP, 0.3559862669062833E-03_SP, 0.5137265251275234E+00_SP, &
         0.3576224317551411E-03_SP, 0.5466764056654611E+00_SP, 0.3584050533086076E-03_SP, &
         0.6054859420813535E+00_SP, 0.3584903581373224E-03_SP, 0.6308106701764562E+00_SP, &
         0.3582991879040586E-03_SP, 0.6530369230179584E+00_SP, 0.3582371187963125E-03_SP, &
         0.6718609524611158E+00_SP, 0.3584353631122350E-03_SP, 0.6869676499894013E+00_SP, &
         0.3589120166517785E-03_SP, 0.6980467077240748E+00_SP, 0.3595445704531601E-03_SP, &
         0.7048241721250522E+00_SP, 0.3600943557111074E-03_SP, 0.5591105222058232E-01_SP, &
         0.1456447096742039E-03_SP, 0.1407384078513916E+00_SP, 0.2252370188283782E-03_SP, &
         0.2364035438976309E+00_SP, 0.2766135443474897E-03_SP, 0.3360602737818170E+00_SP, &
         0.3110729491500851E-03_SP, 0.4356292630054665E+00_SP, 0.3342506712303391E-03_SP, &
         0.5321569415256174E+00_SP, 0.3491981834026860E-03_SP, 0.6232956305040554E+00_SP, &
         0.3576003604348932E-03_SP, 0.9469870086838469E-01_SP, 0.2778748387309470E-01_SP, &
         0.1921921305788564E-03_SP, 0.1353170300568141E+00_SP, 0.6076569878628364E-01_SP, &
         0.2301458216495632E-03_SP, 0.1771679481726077E+00_SP, 0.9703072762711040E-01_SP, &
         0.2604248549522893E-03_SP, 0.2197066664231751E+00_SP, 0.1354112458524762E+00_SP, &
         0.2845275425870697E-03_SP, 0.2624783557374927E+00_SP, 0.1750996479744100E+00_SP, &
         0.3036870897974840E-03_SP, 0.3050969521214442E+00_SP, 0.2154896907449802E+00_SP, &
         0.3188414832298066E-03_SP, 0.3472252637196021E+00_SP, 0.2560954625740152E+00_SP, &
         0.3307046414722089E-03_SP, 0.3885610219026360E+00_SP, 0.2965070050624096E+00_SP, &
         0.3398330969031360E-03_SP, 0.4288273776062765E+00_SP, 0.3363641488734497E+00_SP, &
         0.3466757899705373E-03_SP, 0.4677662471302948E+00_SP, 0.3753400029836788E+00_SP, &
         0.3516095923230054E-03_SP, 0.5051333589553359E+00_SP, 0.4131297522144286E+00_SP, &
         0.3549645184048486E-03_SP, 0.5406942145810492E+00_SP, 0.4494423776081795E+00_SP, &
         0.3570415969441392E-03_SP, 0.5742204122576457E+00_SP, 0.4839938958841502E+00_SP, &
         0.3581251798496118E-03_SP, 0.1865407027225188E+00_SP, 0.3259144851070796E-01_SP, &
         0.2543491329913348E-03_SP, 0.2321186453689432E+00_SP, 0.6835679505297343E-01_SP, &
         0.2786711051330776E-03_SP, 0.2773159142523882E+00_SP, 0.1062284864451989E+00_SP, &
         0.2985552361083679E-03_SP, 0.3219200192237254E+00_SP, 0.1454404409323047E+00_SP, &
         0.3145867929154039E-03_SP, 0.3657032593944029E+00_SP, 0.1854018282582510E+00_SP, &
         0.3273290662067609E-03_SP, 0.4084376778363622E+00_SP, 0.2256297412014750E+00_SP, &
         0.3372705511943501E-03_SP, 0.4499004945751427E+00_SP, 0.2657104425000896E+00_SP, &
         0.3448274437851510E-03_SP, 0.4898758141326335E+00_SP, 0.3052755487631557E+00_SP, &
         0.3503592783048583E-03_SP, 0.5281547442266309E+00_SP, 0.3439863920645423E+00_SP, &
         0.3541854792663162E-03_SP, 0.5645346989813992E+00_SP, 0.3815229456121914E+00_SP, &
         0.3565995517909428E-03_SP, 0.5988181252159848E+00_SP, 0.4175752420966734E+00_SP, &
         0.3578802078302898E-03_SP, 0.2850425424471603E+00_SP, 0.3562149509862536E-01_SP, &
         0.2958644592860982E-03_SP, 0.3324619433027876E+00_SP, 0.7330318886871096E-01_SP, &
         0.3119548129116835E-03_SP, 0.3785848333076282E+00_SP, 0.1123226296008472E+00_SP, &
         0.3250745225005984E-03_SP, 0.4232891028562115E+00_SP, 0.1521084193337708E+00_SP, &
         0.3355153415935208E-03_SP, 0.4664287050829722E+00_SP, 0.1921844459223610E+00_SP, &
         0.3435847568549328E-03_SP, 0.5078458493735726E+00_SP, 0.2321360989678303E+00_SP, &
         0.3495786831622488E-03_SP, 0.5473779816204180E+00_SP, 0.2715886486360520E+00_SP, &
         0.3537767805534621E-03_SP, 0.5848617133811376E+00_SP, 0.3101924707571355E+00_SP, &
         0.3564459815421428E-03_SP, 0.6201348281584888E+00_SP, 0.3476121052890973E+00_SP, &
         0.3578464061225468E-03_SP, 0.3852191185387871E+00_SP, 0.3763224880035108E-01_SP, &
         0.3239748762836212E-03_SP, 0.4325025061073423E+00_SP, 0.7659581935637135E-01_SP, &
         0.3345491784174287E-03_SP, 0.4778486229734490E+00_SP, 0.1163381306083900E+00_SP, &
         0.3429126177301782E-03_SP, 0.5211663693009000E+00_SP, 0.1563890598752899E+00_SP, &
         0.3492420343097421E-03_SP, 0.5623469504853703E+00_SP, 0.1963320810149200E+00_SP, &
         0.3537399050235257E-03_SP, 0.6012718188659246E+00_SP, 0.2357847407258738E+00_SP, &
         0.3566209152659172E-03_SP, 0.6378179206390117E+00_SP, 0.2743846121244060E+00_SP, &
         0.3581084321919782E-03_SP, 0.4836936460214534E+00_SP, 0.3895902610739024E-01_SP, &
         0.3426522117591512E-03_SP, 0.5293792562683797E+00_SP, 0.7871246819312640E-01_SP, &
         0.3491848770121379E-03_SP, 0.5726281253100033E+00_SP, 0.1187963808202981E+00_SP, &
         0.3539318235231476E-03_SP, 0.6133658776169068E+00_SP, 0.1587914708061787E+00_SP, &
         0.3570231438458694E-03_SP, 0.6515085491865307E+00_SP, 0.1983058575227646E+00_SP, &
         0.3586207335051714E-03_SP, 0.5778692716064976E+00_SP, 0.3977209689791542E-01_SP, &
         0.3541196205164025E-03_SP, 0.6207904288086192E+00_SP, 0.7990157592981152E-01_SP, &
         0.3574296911573953E-03_SP, 0.6608688171046802E+00_SP, 0.1199671308754309E+00_SP, &
         0.3591993279818963E-03_SP, 0.6656263089489130E+00_SP, 0.4015955957805969E-01_SP, &
         0.3595855034661997E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_3074_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_3074_pt

  Subroutine Lebedev_3470_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 101, 1, 1, 0, 24, 8, 56, 3470 /)
    REAL(SP),    PARAMETER :: pa(234) = (/ &
         0.2040382730826330E-04_SP, 0.3178149703889544E-03_SP, 0.1721420832906233E-01_SP, &
         0.8288115128076110E-04_SP, 0.4408875374981770E-01_SP, 0.1360883192522954E-03_SP, &
         0.7594680813878681E-01_SP, 0.1766854454542662E-03_SP, 0.1108335359204799E+00_SP, &
         0.2083153161230153E-03_SP, 0.1476517054388567E+00_SP, 0.2333279544657158E-03_SP, &
         0.1856731870860615E+00_SP, 0.2532809539930247E-03_SP, 0.2243634099428821E+00_SP, &
         0.2692472184211158E-03_SP, 0.2633006881662727E+00_SP, 0.2819949946811885E-03_SP, &
         0.3021340904916283E+00_SP, 0.2920953593973030E-03_SP, 0.3405594048030089E+00_SP, &
         0.2999889782948352E-03_SP, 0.3783044434007372E+00_SP, 0.3060292120496902E-03_SP, &
         0.4151194767407910E+00_SP, 0.3105109167522192E-03_SP, 0.4507705766443257E+00_SP, &
         0.3136902387550312E-03_SP, 0.4850346056573187E+00_SP, 0.3157984652454632E-03_SP, &
         0.5176950817792470E+00_SP, 0.3170516518425422E-03_SP, 0.5485384240820989E+00_SP, &
         0.3176568425633755E-03_SP, 0.6039117238943308E+00_SP, 0.3177198411207062E-03_SP, &
         0.6279956655573113E+00_SP, 0.3175519492394733E-03_SP, 0.6493636169568952E+00_SP, &
         0.3174654952634756E-03_SP, 0.6677644117704504E+00_SP, 0.3175676415467654E-03_SP, &
         0.6829368572115624E+00_SP, 0.3178923417835410E-03_SP, 0.6946195818184121E+00_SP, &
         0.3183788287531909E-03_SP, 0.7025711542057026E+00_SP, 0.3188755151918807E-03_SP, &
         0.7066004767140119E+00_SP, 0.3191916889313849E-03_SP, 0.5132537689946062E-01_SP, &
         0.1231779611744508E-03_SP, 0.1297994661331225E+00_SP, 0.1924661373839880E-03_SP, &
         0.2188852049401307E+00_SP, 0.2380881867403424E-03_SP, 0.3123174824903457E+00_SP, &
         0.2693100663037885E-03_SP, 0.4064037620738195E+00_SP, 0.2908673382834366E-03_SP, &
         0.4984958396944782E+00_SP, 0.3053914619381535E-03_SP, 0.5864975046021365E+00_SP, &
         0.3143916684147777E-03_SP, 0.6686711634580175E+00_SP, 0.3187042244055363E-03_SP, &
         0.8715738780835950E-01_SP, 0.2557175233367578E-01_SP, 0.1635219535869790E-03_SP, &
         0.1248383123134007E+00_SP, 0.5604823383376681E-01_SP, 0.1968109917696070E-03_SP, &
         0.1638062693383378E+00_SP, 0.8968568601900765E-01_SP, 0.2236754342249974E-03_SP, &
         0.2035586203373176E+00_SP, 0.1254086651976279E+00_SP, 0.2453186687017181E-03_SP, &
         0.2436798975293774E+00_SP, 0.1624780150162012E+00_SP, 0.2627551791580541E-03_SP, &
         0.2838207507773806E+00_SP, 0.2003422342683208E+00_SP, 0.2767654860152220E-03_SP, &
         0.3236787502217692E+00_SP, 0.2385628026255263E+00_SP, 0.2879467027765895E-03_SP, &
         0.3629849554840691E+00_SP, 0.2767731148783578E+00_SP, 0.2967639918918702E-03_SP, &
         0.4014948081992087E+00_SP, 0.3146542308245309E+00_SP, 0.3035900684660351E-03_SP, &
         0.4389818379260225E+00_SP, 0.3519196415895088E+00_SP, 0.3087338237298308E-03_SP, &
         0.4752331143674377E+00_SP, 0.3883050984023654E+00_SP, 0.3124608838860167E-03_SP, &
         0.5100457318374018E+00_SP, 0.4235613423908649E+00_SP, 0.3150084294226743E-03_SP, &
         0.5432238388954868E+00_SP, 0.4574484717196220E+00_SP, 0.3165958398598402E-03_SP, &
         0.5745758685072442E+00_SP, 0.4897311639255524E+00_SP, 0.3174320440957372E-03_SP, &
         0.1723981437592809E+00_SP, 0.3010630597881105E-01_SP, 0.2182188909812599E-03_SP, &
         0.2149553257844597E+00_SP, 0.6326031554204694E-01_SP, 0.2399727933921445E-03_SP, &
         0.2573256081247422E+00_SP, 0.9848566980258631E-01_SP, 0.2579796133514652E-03_SP, &
         0.2993163751238106E+00_SP, 0.1350835952384266E+00_SP, 0.2727114052623535E-03_SP, &
         0.3407238005148000E+00_SP, 0.1725184055442181E+00_SP, 0.2846327656281355E-03_SP, &
         0.3813454978483264E+00_SP, 0.2103559279730725E+00_SP, 0.2941491102051334E-03_SP, &
         0.4209848104423343E+00_SP, 0.2482278774554860E+00_SP, 0.3016049492136107E-03_SP, &
         0.4594519699996300E+00_SP, 0.2858099509982883E+00_SP, 0.3072949726175648E-03_SP, &
         0.4965640166185930E+00_SP, 0.3228075659915428E+00_SP, 0.3114768142886460E-03_SP, &
         0.5321441655571562E+00_SP, 0.3589459907204151E+00_SP, 0.3143823673666223E-03_SP, &
         0.5660208438582166E+00_SP, 0.3939630088864310E+00_SP, 0.3162269764661535E-03_SP, &
         0.5980264315964364E+00_SP, 0.4276029922949089E+00_SP, 0.3172164663759821E-03_SP, &
         0.2644215852350733E+00_SP, 0.3300939429072552E-01_SP, 0.2554575398967435E-03_SP, &
         0.3090113743443063E+00_SP, 0.6803887650078501E-01_SP, 0.2701704069135677E-03_SP, &
         0.3525871079197808E+00_SP, 0.1044326136206709E+00_SP, 0.2823693413468940E-03_SP, &
         0.3950418005354029E+00_SP, 0.1416751597517679E+00_SP, 0.2922898463214289E-03_SP, &
         0.4362475663430163E+00_SP, 0.1793408610504821E+00_SP, 0.3001829062162428E-03_SP, &
         0.4760661812145854E+00_SP, 0.2170630750175722E+00_SP, 0.3062890864542953E-03_SP, &
         0.5143551042512103E+00_SP, 0.2545145157815807E+00_SP, 0.3108328279264746E-03_SP, &
         0.5509709026935597E+00_SP, 0.2913940101706601E+00_SP, 0.3140243146201245E-03_SP, &
         0.5857711030329428E+00_SP, 0.3274169910910705E+00_SP, 0.3160638030977130E-03_SP, &
         0.6186149917404392E+00_SP, 0.3623081329317265E+00_SP, 0.3171462882206275E-03_SP, &
         0.3586894569557064E+00_SP, 0.3497354386450040E-01_SP, 0.2812388416031796E-03_SP, &
         0.4035266610019441E+00_SP, 0.7129736739757095E-01_SP, 0.2912137500288045E-03_SP, &
         0.4467775312332510E+00_SP, 0.1084758620193165E+00_SP, 0.2993241256502206E-03_SP, &
         0.4883638346608543E+00_SP, 0.1460915689241772E+00_SP, 0.3057101738983822E-03_SP, &
         0.5281908348434601E+00_SP, 0.1837790832369980E+00_SP, 0.3105319326251432E-03_SP, &
         0.5661542687149311E+00_SP, 0.2212075390874021E+00_SP, 0.3139565514428167E-03_SP, &
         0.6021450102031452E+00_SP, 0.2580682841160985E+00_SP, 0.3161543006806366E-03_SP, &
         0.6360520783610050E+00_SP, 0.2940656362094121E+00_SP, 0.3172985960613294E-03_SP, &
         0.4521611065087196E+00_SP, 0.3631055365867002E-01_SP, 0.2989400336901431E-03_SP, &
         0.4959365651560963E+00_SP, 0.7348318468484350E-01_SP, 0.3054555883947677E-03_SP, &
         0.5376815804038283E+00_SP, 0.1111087643812648E+00_SP, 0.3104764960807702E-03_SP, &
         0.5773314480243768E+00_SP, 0.1488226085145408E+00_SP, 0.3141015825977616E-03_SP, &
         0.6148113245575056E+00_SP, 0.1862892274135151E+00_SP, 0.3164520621159896E-03_SP, &
         0.6500407462842380E+00_SP, 0.2231909701714456E+00_SP, 0.3176652305912204E-03_SP, &
         0.5425151448707213E+00_SP, 0.3718201306118944E-01_SP, 0.3105097161023939E-03_SP, &
         0.5841860556907931E+00_SP, 0.7483616335067346E-01_SP, 0.3143014117890550E-03_SP, &
         0.6234632186851500E+00_SP, 0.1125990834266120E+00_SP, 0.3168172866287200E-03_SP, &
         0.6602934551848843E+00_SP, 0.1501303813157619E+00_SP, 0.3181401865570968E-03_SP, &
         0.6278573968375105E+00_SP, 0.3767559930245720E-01_SP, 0.3170663659156037E-03_SP, &
         0.6665611711264577E+00_SP, 0.7548443301360158E-01_SP, 0.3185447944625510E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_3470_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_3470_pt

  Subroutine Lebedev_3890_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 107, 1, 1, 1, 25, 8, 64, 3890 /)
    REAL(SP),    PARAMETER :: pa(261) = (/ &
         0.1807395252196920E-04_SP, 0.2836065837530581E-03_SP, 0.2848008782238827E-03_SP, &
         0.1587876419858352E-01_SP, 0.7013149266673816E-04_SP, 0.4069193593751206E-01_SP, &
         0.1162798021956766E-03_SP, 0.7025888115257997E-01_SP, 0.1518728583972105E-03_SP, &
         0.1027495450028704E+00_SP, 0.1798796108216934E-03_SP, 0.1371457730893426E+00_SP, &
         0.2022593385972785E-03_SP, 0.1727758532671953E+00_SP, 0.2203093105575464E-03_SP, &
         0.2091492038929037E+00_SP, 0.2349294234299855E-03_SP, 0.2458813281751915E+00_SP, &
         0.2467682058747003E-03_SP, 0.2826545859450066E+00_SP, 0.2563092683572224E-03_SP, &
         0.3191957291799622E+00_SP, 0.2639253896763318E-03_SP, 0.3552621469299578E+00_SP, &
         0.2699137479265108E-03_SP, 0.3906329503406230E+00_SP, 0.2745196420166739E-03_SP, &
         0.4251028614093031E+00_SP, 0.2779529197397593E-03_SP, 0.4584777520111870E+00_SP, &
         0.2803996086684265E-03_SP, 0.4905711358710193E+00_SP, 0.2820302356715842E-03_SP, &
         0.5212011669847385E+00_SP, 0.2830056747491068E-03_SP, 0.5501878488737995E+00_SP, &
         0.2834808950776839E-03_SP, 0.6025037877479342E+00_SP, 0.2835282339078929E-03_SP, &
         0.6254572689549016E+00_SP, 0.2833819267065800E-03_SP, 0.6460107179528248E+00_SP, &
         0.2832858336906784E-03_SP, 0.6639541138154251E+00_SP, 0.2833268235451244E-03_SP, &
         0.6790688515667495E+00_SP, 0.2835432677029253E-03_SP, 0.6911338580371512E+00_SP, &
         0.2839091722743049E-03_SP, 0.6999385956126490E+00_SP, 0.2843308178875841E-03_SP, &
         0.7053037748656896E+00_SP, 0.2846703550533846E-03_SP, 0.4732224387180115E-01_SP, &
         0.1051193406971900E-03_SP, 0.1202100529326803E+00_SP, 0.1657871838796974E-03_SP, &
         0.2034304820664855E+00_SP, 0.2064648113714232E-03_SP, 0.2912285643573002E+00_SP, &
         0.2347942745819741E-03_SP, 0.3802361792726768E+00_SP, 0.2547775326597726E-03_SP, &
         0.4680598511056146E+00_SP, 0.2686876684847025E-03_SP, 0.5528151052155599E+00_SP, &
         0.2778665755515867E-03_SP, 0.6329386307803041E+00_SP, 0.2830996616782929E-03_SP, &
         0.8056516651369069E-01_SP, 0.2363454684003124E-01_SP, 0.1403063340168372E-03_SP, &
         0.1156476077139389E+00_SP, 0.5191291632545936E-01_SP, 0.1696504125939477E-03_SP, &
         0.1520473382760421E+00_SP, 0.8322715736994519E-01_SP, 0.1935787242745390E-03_SP, &
         0.1892986699745931E+00_SP, 0.1165855667993712E+00_SP, 0.2130614510521968E-03_SP, &
         0.2270194446777792E+00_SP, 0.1513077167409504E+00_SP, 0.2289381265931048E-03_SP, &
         0.2648908185093273E+00_SP, 0.1868882025807859E+00_SP, 0.2418630292816186E-03_SP, &
         0.3026389259574136E+00_SP, 0.2229277629776224E+00_SP, 0.2523400495631193E-03_SP, &
         0.3400220296151384E+00_SP, 0.2590951840746235E+00_SP, 0.2607623973449605E-03_SP, &
         0.3768217953335510E+00_SP, 0.2951047291750847E+00_SP, 0.2674441032689209E-03_SP, &
         0.4128372900921884E+00_SP, 0.3307019714169930E+00_SP, 0.2726432360343356E-03_SP, &
         0.4478807131815630E+00_SP, 0.3656544101087634E+00_SP, 0.2765787685924545E-03_SP, &
         0.4817742034089257E+00_SP, 0.3997448951939695E+00_SP, 0.2794428690642224E-03_SP, &
         0.5143472814653344E+00_SP, 0.4327667110812024E+00_SP, 0.2814099002062895E-03_SP, &
         0.5454346213905650E+00_SP, 0.4645196123532293E+00_SP, 0.2826429531578994E-03_SP, &
         0.5748739313170252E+00_SP, 0.4948063555703345E+00_SP, 0.2832983542550884E-03_SP, &
         0.1599598738286342E+00_SP, 0.2792357590048985E-01_SP, 0.1886695565284976E-03_SP, &
         0.1998097412500951E+00_SP, 0.5877141038139065E-01_SP, 0.2081867882748234E-03_SP, &
         0.2396228952566202E+00_SP, 0.9164573914691377E-01_SP, 0.2245148680600796E-03_SP, &
         0.2792228341097746E+00_SP, 0.1259049641962687E+00_SP, 0.2380370491511872E-03_SP, &
         0.3184251107546741E+00_SP, 0.1610594823400863E+00_SP, 0.2491398041852455E-03_SP, &
         0.3570481164426244E+00_SP, 0.1967151653460898E+00_SP, 0.2581632405881230E-03_SP, &
         0.3949164710492144E+00_SP, 0.2325404606175168E+00_SP, 0.2653965506227417E-03_SP, &
         0.4318617293970503E+00_SP, 0.2682461141151439E+00_SP, 0.2710857216747087E-03_SP, &
         0.4677221009931678E+00_SP, 0.3035720116011973E+00_SP, 0.2754434093903659E-03_SP, &
         0.5023417939270955E+00_SP, 0.3382781859197439E+00_SP, 0.2786579932519380E-03_SP, &
         0.5355701836636128E+00_SP, 0.3721383065625942E+00_SP, 0.2809011080679474E-03_SP, &
         0.5672608451328771E+00_SP, 0.4049346360466055E+00_SP, 0.2823336184560987E-03_SP, &
         0.5972704202540162E+00_SP, 0.4364538098633802E+00_SP, 0.2831101175806309E-03_SP, &
         0.2461687022333596E+00_SP, 0.3070423166833368E-01_SP, 0.2221679970354546E-03_SP, &
         0.2881774566286831E+00_SP, 0.6338034669281885E-01_SP, 0.2356185734270703E-03_SP, &
         0.3293963604116978E+00_SP, 0.9742862487067941E-01_SP, 0.2469228344805590E-03_SP, &
         0.3697303822241377E+00_SP, 0.1323799532282290E+00_SP, 0.2562726348642046E-03_SP, &
         0.4090663023135127E+00_SP, 0.1678497018129336E+00_SP, 0.2638756726753028E-03_SP, &
         0.4472819355411712E+00_SP, 0.2035095105326114E+00_SP, 0.2699311157390862E-03_SP, &
         0.4842513377231437E+00_SP, 0.2390692566672091E+00_SP, 0.2746233268403837E-03_SP, &
         0.5198477629962928E+00_SP, 0.2742649818076149E+00_SP, 0.2781225674454771E-03_SP, &
         0.5539453011883145E+00_SP, 0.3088503806580094E+00_SP, 0.2805881254045684E-03_SP, &
         0.5864196762401251E+00_SP, 0.3425904245906614E+00_SP, 0.2821719877004913E-03_SP, &
         0.6171484466668390E+00_SP, 0.3752562294789468E+00_SP, 0.2830222502333124E-03_SP, &
         0.3350337830565727E+00_SP, 0.3261589934634747E-01_SP, 0.2457995956744870E-03_SP, &
         0.3775773224758284E+00_SP, 0.6658438928081572E-01_SP, 0.2551474407503706E-03_SP, &
         0.4188155229848973E+00_SP, 0.1014565797157954E+00_SP, 0.2629065335195311E-03_SP, &
         0.4586805892009344E+00_SP, 0.1368573320843822E+00_SP, 0.2691900449925075E-03_SP, &
         0.4970895714224235E+00_SP, 0.1724614851951608E+00_SP, 0.2741275485754276E-03_SP, &
         0.5339505133960747E+00_SP, 0.2079779381416412E+00_SP, 0.2778530970122595E-03_SP, &
         0.5691665792531440E+00_SP, 0.2431385788322288E+00_SP, 0.2805010567646741E-03_SP, &
         0.6026387682680377E+00_SP, 0.2776901883049853E+00_SP, 0.2822055834031040E-03_SP, &
         0.6342676150163307E+00_SP, 0.3113881356386632E+00_SP, 0.2831016901243473E-03_SP, &
         0.4237951119537067E+00_SP, 0.3394877848664351E-01_SP, 0.2624474901131803E-03_SP, &
         0.4656918683234929E+00_SP, 0.6880219556291447E-01_SP, 0.2688034163039377E-03_SP, &
         0.5058857069185980E+00_SP, 0.1041946859721635E+00_SP, 0.2738932751287636E-03_SP, &
         0.5443204666713996E+00_SP, 0.1398039738736393E+00_SP, 0.2777944791242523E-03_SP, &
         0.5809298813759742E+00_SP, 0.1753373381196155E+00_SP, 0.2806011661660987E-03_SP, &
         0.6156416039447128E+00_SP, 0.2105215793514010E+00_SP, 0.2824181456597460E-03_SP, &
         0.6483801351066604E+00_SP, 0.2450953312157051E+00_SP, 0.2833585216577828E-03_SP, &
         0.5103616577251688E+00_SP, 0.3485560643800719E-01_SP, 0.2738165236962878E-03_SP, &
         0.5506738792580681E+00_SP, 0.7026308631512033E-01_SP, 0.2778365208203180E-03_SP, &
         0.5889573040995292E+00_SP, 0.1059035061296403E+00_SP, 0.2807852940418966E-03_SP, &
         0.6251641589516930E+00_SP, 0.1414823925236026E+00_SP, 0.2827245949674705E-03_SP, &
         0.6592414921570178E+00_SP, 0.1767207908214530E+00_SP, 0.2837342344829828E-03_SP, &
         0.5930314017533384E+00_SP, 0.3542189339561672E-01_SP, 0.2809233907610981E-03_SP, &
         0.6309812253390175E+00_SP, 0.7109574040369549E-01_SP, 0.2829930809742694E-03_SP, &
         0.6666296011353230E+00_SP, 0.1067259792282730E+00_SP, 0.2841097874111479E-03_SP, &
         0.6703715271049922E+00_SP, 0.3569455268820809E-01_SP, 0.2843455206008783E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_3890_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_3890_pt

  Subroutine Lebedev_4334_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 113, 1, 1, 0, 27, 9, 72, 4334 /)
    REAL(SP),    PARAMETER :: pa(290) = (/ &
         0.1449063022537883E-04_SP, 0.2546377329828424E-03_SP, 0.1462896151831013E-01_SP, &
         0.6018432961087496E-04_SP, 0.3769840812493139E-01_SP, 0.1002286583263673E-03_SP, &
         0.6524701904096891E-01_SP, 0.1315222931028093E-03_SP, 0.9560543416134648E-01_SP, &
         0.1564213746876724E-03_SP, 0.1278335898929198E+00_SP, 0.1765118841507736E-03_SP, &
         0.1613096104466031E+00_SP, 0.1928737099311080E-03_SP, 0.1955806225745371E+00_SP, &
         0.2062658534263270E-03_SP, 0.2302935218498028E+00_SP, 0.2172395445953787E-03_SP, &
         0.2651584344113027E+00_SP, 0.2262076188876047E-03_SP, 0.2999276825183209E+00_SP, &
         0.2334885699462397E-03_SP, 0.3343828669718798E+00_SP, 0.2393355273179203E-03_SP, &
         0.3683265013750518E+00_SP, 0.2439559200468863E-03_SP, 0.4015763206518108E+00_SP, &
         0.2475251866060002E-03_SP, 0.4339612026399770E+00_SP, 0.2501965558158773E-03_SP, &
         0.4653180651114582E+00_SP, 0.2521081407925925E-03_SP, 0.4954893331080803E+00_SP, &
         0.2533881002388081E-03_SP, 0.5243207068924930E+00_SP, 0.2541582900848261E-03_SP, &
         0.5516590479041704E+00_SP, 0.2545365737525860E-03_SP, 0.6012371927804176E+00_SP, &
         0.2545726993066799E-03_SP, 0.6231574466449819E+00_SP, 0.2544456197465555E-03_SP, &
         0.6429416514181271E+00_SP, 0.2543481596881064E-03_SP, 0.6604124272943595E+00_SP, &
         0.2543506451429194E-03_SP, 0.6753851470408250E+00_SP, 0.2544905675493763E-03_SP, &
         0.6876717970626160E+00_SP, 0.2547611407344429E-03_SP, 0.6970895061319234E+00_SP, &
         0.2551060375448869E-03_SP, 0.7034746912553310E+00_SP, 0.2554291933816039E-03_SP, &
         0.7067017217542295E+00_SP, 0.2556255710686343E-03_SP, 0.4382223501131123E-01_SP, &
         0.9041339695118195E-04_SP, 0.1117474077400006E+00_SP, 0.1438426330079022E-03_SP, &
         0.1897153252911440E+00_SP, 0.1802523089820518E-03_SP, 0.2724023009910331E+00_SP, &
         0.2060052290565496E-03_SP, 0.3567163308709902E+00_SP, 0.2245002248967466E-03_SP, &
         0.4404784483028087E+00_SP, 0.2377059847731150E-03_SP, 0.5219833154161411E+00_SP, &
         0.2468118955882525E-03_SP, 0.5998179868977553E+00_SP, 0.2525410872966528E-03_SP, &
         0.6727803154548222E+00_SP, 0.2553101409933397E-03_SP, 0.7476563943166086E-01_SP, &
         0.2193168509461185E-01_SP, 0.1212879733668632E-03_SP, 0.1075341482001416E+00_SP, &
         0.4826419281533887E-01_SP, 0.1472872881270931E-03_SP, 0.1416344885203259E+00_SP, &
         0.7751191883575742E-01_SP, 0.1686846601010828E-03_SP, 0.1766325315388586E+00_SP, &
         0.1087558139247680E+00_SP, 0.1862698414660208E-03_SP, 0.2121744174481514E+00_SP, &
         0.1413661374253096E+00_SP, 0.2007430956991861E-03_SP, 0.2479669443408145E+00_SP, &
         0.1748768214258880E+00_SP, 0.2126568125394796E-03_SP, 0.2837600452294113E+00_SP, &
         0.2089216406612073E+00_SP, 0.2224394603372113E-03_SP, 0.3193344933193984E+00_SP, &
         0.2431987685545972E+00_SP, 0.2304264522673135E-03_SP, 0.3544935442438745E+00_SP, &
         0.2774497054377770E+00_SP, 0.2368854288424087E-03_SP, 0.3890571932288154E+00_SP, &
         0.3114460356156915E+00_SP, 0.2420352089461772E-03_SP, 0.4228581214259090E+00_SP, &
         0.3449806851913012E+00_SP, 0.2460597113081295E-03_SP, 0.4557387211304052E+00_SP, &
         0.3778618641248256E+00_SP, 0.2491181912257687E-03_SP, 0.4875487950541643E+00_SP, &
         0.4099086391698978E+00_SP, 0.2513528194205857E-03_SP, 0.5181436529962997E+00_SP, &
         0.4409474925853973E+00_SP, 0.2528943096693220E-03_SP, 0.5473824095600661E+00_SP, &
         0.4708094517711291E+00_SP, 0.2538660368488136E-03_SP, 0.5751263398976174E+00_SP, &
         0.4993275140354637E+00_SP, 0.2543868648299022E-03_SP, 0.1489515746840028E+00_SP, &
         0.2599381993267017E-01_SP, 0.1642595537825183E-03_SP, 0.1863656444351767E+00_SP, &
         0.5479286532462190E-01_SP, 0.1818246659849308E-03_SP, 0.2238602880356348E+00_SP, &
         0.8556763251425254E-01_SP, 0.1966565649492420E-03_SP, 0.2612723375728160E+00_SP, &
         0.1177257802267011E+00_SP, 0.2090677905657991E-03_SP, 0.2984332990206190E+00_SP, &
         0.1508168456192700E+00_SP, 0.2193820409510504E-03_SP, 0.3351786584663333E+00_SP, &
         0.1844801892177727E+00_SP, 0.2278870827661928E-03_SP, 0.3713505522209120E+00_SP, &
         0.2184145236087598E+00_SP, 0.2348283192282090E-03_SP, 0.4067981098954663E+00_SP, &
         0.2523590641486229E+00_SP, 0.2404139755581477E-03_SP, 0.4413769993687534E+00_SP, &
         0.2860812976901373E+00_SP, 0.2448227407760734E-03_SP, 0.4749487182516394E+00_SP, &
         0.3193686757808996E+00_SP, 0.2482110455592573E-03_SP, 0.5073798105075426E+00_SP, &
         0.3520226949547602E+00_SP, 0.2507192397774103E-03_SP, 0.5385410448878654E+00_SP, &
         0.3838544395667890E+00_SP, 0.2524765968534880E-03_SP, 0.5683065353670530E+00_SP, &
         0.4146810037640963E+00_SP, 0.2536052388539425E-03_SP, 0.5965527620663510E+00_SP, &
         0.4443224094681121E+00_SP, 0.2542230588033068E-03_SP, 0.2299227700856157E+00_SP, &
         0.2865757664057584E-01_SP, 0.1944817013047896E-03_SP, 0.2695752998553267E+00_SP, &
         0.5923421684485993E-01_SP, 0.2067862362746635E-03_SP, 0.3086178716611389E+00_SP, &
         0.9117817776057715E-01_SP, 0.2172440734649114E-03_SP, 0.3469649871659077E+00_SP, &
         0.1240593814082605E+00_SP, 0.2260125991723423E-03_SP, 0.3845153566319655E+00_SP, &
         0.1575272058259175E+00_SP, 0.2332655008689523E-03_SP, 0.4211600033403215E+00_SP, &
         0.1912845163525413E+00_SP, 0.2391699681532458E-03_SP, 0.4567867834329882E+00_SP, &
         0.2250710177858171E+00_SP, 0.2438801528273928E-03_SP, 0.4912829319232061E+00_SP, &
         0.2586521303440910E+00_SP, 0.2475370504260665E-03_SP, 0.5245364793303812E+00_SP, &
         0.2918112242865407E+00_SP, 0.2502707235640574E-03_SP, 0.5564369788915756E+00_SP, &
         0.3243439239067890E+00_SP, 0.2522031701054241E-03_SP, 0.5868757697775287E+00_SP, &
         0.3560536787835351E+00_SP, 0.2534511269978784E-03_SP, 0.6157458853519617E+00_SP, &
         0.3867480821242581E+00_SP, 0.2541284914955151E-03_SP, 0.3138461110672113E+00_SP, &
         0.3051374637507278E-01_SP, 0.2161509250688394E-03_SP, 0.3542495872050569E+00_SP, &
         0.6237111233730755E-01_SP, 0.2248778513437852E-03_SP, 0.3935751553120181E+00_SP, &
         0.9516223952401907E-01_SP, 0.2322388803404617E-03_SP, 0.4317634668111147E+00_SP, &
         0.1285467341508517E+00_SP, 0.2383265471001355E-03_SP, 0.4687413842250821E+00_SP, &
         0.1622318931656033E+00_SP, 0.2432476675019525E-03_SP, 0.5044274237060283E+00_SP, &
         0.1959581153836453E+00_SP, 0.2471122223750674E-03_SP, 0.5387354077925727E+00_SP, &
         0.2294888081183837E+00_SP, 0.2500291752486870E-03_SP, 0.5715768898356105E+00_SP, &
         0.2626031152713945E+00_SP, 0.2521055942764682E-03_SP, 0.6028627200136111E+00_SP, &
         0.2950904075286713E+00_SP, 0.2534472785575503E-03_SP, 0.6325039812653463E+00_SP, &
         0.3267458451113286E+00_SP, 0.2541599713080121E-03_SP, 0.3981986708423407E+00_SP, &
         0.3183291458749821E-01_SP, 0.2317380975862936E-03_SP, 0.4382791182133300E+00_SP, &
         0.6459548193880908E-01_SP, 0.2378550733719775E-03_SP, 0.4769233057218166E+00_SP, &
         0.9795757037087952E-01_SP, 0.2428884456739118E-03_SP, 0.5140823911194238E+00_SP, &
         0.1316307235126655E+00_SP, 0.2469002655757292E-03_SP, 0.5496977833862983E+00_SP, &
         0.1653556486358704E+00_SP, 0.2499657574265851E-03_SP, 0.5837047306512727E+00_SP, &
         0.1988931724126510E+00_SP, 0.2521676168486082E-03_SP, 0.6160349566926879E+00_SP, &
         0.2320174581438950E+00_SP, 0.2535935662645334E-03_SP, 0.6466185353209440E+00_SP, &
         0.2645106562168662E+00_SP, 0.2543356743363214E-03_SP, 0.4810835158795404E+00_SP, &
         0.3275917807743992E-01_SP, 0.2427353285201535E-03_SP, 0.5199925041324341E+00_SP, &
         0.6612546183967181E-01_SP, 0.2468258039744386E-03_SP, 0.5571717692207494E+00_SP, &
         0.9981498331474143E-01_SP, 0.2500060956440310E-03_SP, 0.5925789250836378E+00_SP, &
         0.1335687001410374E+00_SP, 0.2523238365420979E-03_SP, 0.6261658523859670E+00_SP, &
         0.1671444402896463E+00_SP, 0.2538399260252846E-03_SP, 0.6578811126669331E+00_SP, &
         0.2003106382156076E+00_SP, 0.2546255927268069E-03_SP, 0.5609624612998100E+00_SP, &
         0.3337500940231335E-01_SP, 0.2500583360048449E-03_SP, 0.5979959659984670E+00_SP, &
         0.6708750335901803E-01_SP, 0.2524777638260203E-03_SP, 0.6330523711054002E+00_SP, &
         0.1008792126424850E+00_SP, 0.2540951193860656E-03_SP, 0.6660960998103972E+00_SP, &
         0.1345050343171794E+00_SP, 0.2549524085027472E-03_SP, 0.6365384364585819E+00_SP, &
         0.3372799460737052E-01_SP, 0.2542569507009158E-03_SP, 0.6710994302899275E+00_SP, &
         0.6755249309678028E-01_SP, 0.2552114127580376E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_4334_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_4334_pt

  Subroutine Lebedev_4802_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 119, 1, 1, 1, 28, 9, 81, 4802 /)
    REAL(SP),    PARAMETER :: pa(320) = (/ &
         0.9687521879420705E-04_SP, 0.2297310852498558E-03_SP, 0.2307897895367918E-03_SP, &
         0.2335728608887064E-01_SP, 0.7386265944001919E-04_SP, 0.4352987836550653E-01_SP, &
         0.8257977698542210E-04_SP, 0.6439200521088801E-01_SP, 0.9706044762057630E-04_SP, &
         0.9003943631993181E-01_SP, 0.1302393847117003E-03_SP, 0.1196706615548473E+00_SP, &
         0.1541957004600968E-03_SP, 0.1511715412838134E+00_SP, 0.1704459770092199E-03_SP, &
         0.1835982828503801E+00_SP, 0.1827374890942906E-03_SP, 0.2165081259155405E+00_SP, &
         0.1926360817436107E-03_SP, 0.2496208720417563E+00_SP, 0.2008010239494833E-03_SP, &
         0.2827200673567900E+00_SP, 0.2075635983209175E-03_SP, 0.3156190823994346E+00_SP, &
         0.2131306638690909E-03_SP, 0.3481476793749115E+00_SP, 0.2176562329937335E-03_SP, &
         0.3801466086947226E+00_SP, 0.2212682262991018E-03_SP, 0.4114652119634011E+00_SP, &
         0.2240799515668565E-03_SP, 0.4419598786519751E+00_SP, 0.2261959816187525E-03_SP, &
         0.4714925949329543E+00_SP, 0.2277156368808855E-03_SP, 0.4999293972879466E+00_SP, &
         0.2287351772128336E-03_SP, 0.5271387221431248E+00_SP, 0.2293490814084085E-03_SP, &
         0.5529896780837761E+00_SP, 0.2296505312376273E-03_SP, 0.6000856099481712E+00_SP, &
         0.2296793832318756E-03_SP, 0.6210562192785175E+00_SP, 0.2295785443842974E-03_SP, &
         0.6401165879934240E+00_SP, 0.2295017931529102E-03_SP, 0.6571144029244334E+00_SP, &
         0.2295059638184868E-03_SP, 0.6718910821718863E+00_SP, 0.2296232343237362E-03_SP, &
         0.6842845591099010E+00_SP, 0.2298530178740771E-03_SP, 0.6941353476269816E+00_SP, &
         0.2301579790280501E-03_SP, 0.7012965242212991E+00_SP, 0.2304690404996513E-03_SP, &
         0.7056471428242644E+00_SP, 0.2307027995907102E-03_SP, 0.4595557643585895E-01_SP, &
         0.9312274696671092E-04_SP, 0.1049316742435023E+00_SP, 0.1199919385876926E-03_SP, &
         0.1773548879549274E+00_SP, 0.1598039138877690E-03_SP, 0.2559071411236127E+00_SP, &
         0.1822253763574900E-03_SP, 0.3358156837985898E+00_SP, 0.1988579593655040E-03_SP, &
         0.4155835743763893E+00_SP, 0.2112620102533307E-03_SP, 0.4937894296167472E+00_SP, &
         0.2201594887699007E-03_SP, 0.5691569694793316E+00_SP, 0.2261622590895036E-03_SP, &
         0.6405840854894251E+00_SP, 0.2296458453435705E-03_SP, 0.7345133894143348E-01_SP, &
         0.2177844081486067E-01_SP, 0.1006006990267000E-03_SP, 0.1009859834044931E+00_SP, &
         0.4590362185775188E-01_SP, 0.1227676689635876E-03_SP, 0.1324289619748758E+00_SP, &
         0.7255063095690877E-01_SP, 0.1467864280270117E-03_SP, 0.1654272109607127E+00_SP, &
         0.1017825451960684E+00_SP, 0.1644178912101232E-03_SP, 0.1990767186776461E+00_SP, &
         0.1325652320980364E+00_SP, 0.1777664890718961E-03_SP, 0.2330125945523278E+00_SP, &
         0.1642765374496765E+00_SP, 0.1884825664516690E-03_SP, 0.2670080611108287E+00_SP, &
         0.1965360374337889E+00_SP, 0.1973269246453848E-03_SP, 0.3008753376294316E+00_SP, &
         0.2290726770542238E+00_SP, 0.2046767775855328E-03_SP, 0.3344475596167860E+00_SP, &
         0.2616645495370823E+00_SP, 0.2107600125918040E-03_SP, 0.3675709724070786E+00_SP, &
         0.2941150728843141E+00_SP, 0.2157416362266829E-03_SP, 0.4001000887587812E+00_SP, &
         0.3262440400919066E+00_SP, 0.2197557816920721E-03_SP, 0.4318956350436028E+00_SP, &
         0.3578835350611916E+00_SP, 0.2229192611835437E-03_SP, 0.4628239056795531E+00_SP, &
         0.3888751854043678E+00_SP, 0.2253385110212775E-03_SP, 0.4927563229773636E+00_SP, &
         0.4190678003222840E+00_SP, 0.2271137107548774E-03_SP, 0.5215687136707969E+00_SP, &
         0.4483151836883852E+00_SP, 0.2283414092917525E-03_SP, 0.5491402346984905E+00_SP, &
         0.4764740676087880E+00_SP, 0.2291161673130077E-03_SP, 0.5753520160126075E+00_SP, &
         0.5034021310998277E+00_SP, 0.2295313908576598E-03_SP, 0.1388326356417754E+00_SP, &
         0.2435436510372806E-01_SP, 0.1438204721359031E-03_SP, 0.1743686900537244E+00_SP, &
         0.5118897057342652E-01_SP, 0.1607738025495257E-03_SP, 0.2099737037950268E+00_SP, &
         0.8014695048539634E-01_SP, 0.1741483853528379E-03_SP, 0.2454492590908548E+00_SP, &
         0.1105117874155699E+00_SP, 0.1851918467519151E-03_SP, 0.2807219257864278E+00_SP, &
         0.1417950531570966E+00_SP, 0.1944628638070613E-03_SP, 0.3156842271975842E+00_SP, &
         0.1736604945719597E+00_SP, 0.2022495446275152E-03_SP, 0.3502090945177752E+00_SP, &
         0.2058466324693981E+00_SP, 0.2087462382438514E-03_SP, 0.3841684849519686E+00_SP, &
         0.2381284261195919E+00_SP, 0.2141074754818308E-03_SP, 0.4174372367906016E+00_SP, &
         0.2703031270422569E+00_SP, 0.2184640913748162E-03_SP, 0.4498926465011892E+00_SP, &
         0.3021845683091309E+00_SP, 0.2219309165220329E-03_SP, 0.4814146229807701E+00_SP, &
         0.3335993355165720E+00_SP, 0.2246123118340624E-03_SP, 0.5118863625734701E+00_SP, &
         0.3643833735518232E+00_SP, 0.2266062766915125E-03_SP, 0.5411947455119144E+00_SP, &
         0.3943789541958179E+00_SP, 0.2280072952230796E-03_SP, 0.5692301500357246E+00_SP, &
         0.4234320144403542E+00_SP, 0.2289082025202583E-03_SP, 0.5958857204139576E+00_SP, &
         0.4513897947419260E+00_SP, 0.2294012695120025E-03_SP, 0.2156270284785766E+00_SP, &
         0.2681225755444491E-01_SP, 0.1722434488736947E-03_SP, 0.2532385054909710E+00_SP, &
         0.5557495747805614E-01_SP, 0.1830237421455091E-03_SP, 0.2902564617771537E+00_SP, &
         0.8569368062950249E-01_SP, 0.1923855349997633E-03_SP, 0.3266979823143256E+00_SP, &
         0.1167367450324135E+00_SP, 0.2004067861936271E-03_SP, 0.3625039627493614E+00_SP, &
         0.1483861994003304E+00_SP, 0.2071817297354263E-03_SP, 0.3975838937548699E+00_SP, &
         0.1803821503011405E+00_SP, 0.2128250834102103E-03_SP, 0.4318396099009774E+00_SP, &
         0.2124962965666424E+00_SP, 0.2174513719440102E-03_SP, 0.4651706555732742E+00_SP, &
         0.2445221837805913E+00_SP, 0.2211661839150214E-03_SP, 0.4974752649620969E+00_SP, &
         0.2762701224322987E+00_SP, 0.2240665257813102E-03_SP, 0.5286517579627517E+00_SP, &
         0.3075627775211328E+00_SP, 0.2262439516632620E-03_SP, 0.5586001195731895E+00_SP, &
         0.3382311089826877E+00_SP, 0.2277874557231869E-03_SP, 0.5872229902021319E+00_SP, &
         0.3681108834741399E+00_SP, 0.2287854314454994E-03_SP, 0.6144258616235123E+00_SP, &
         0.3970397446872839E+00_SP, 0.2293268499615575E-03_SP, 0.2951676508064861E+00_SP, &
         0.2867499538750441E-01_SP, 0.1912628201529828E-03_SP, 0.3335085485472725E+00_SP, &
         0.5867879341903510E-01_SP, 0.1992499672238701E-03_SP, 0.3709561760636381E+00_SP, &
         0.8961099205022284E-01_SP, 0.2061275533454027E-03_SP, 0.4074722861667498E+00_SP, &
         0.1211627927626297E+00_SP, 0.2119318215968572E-03_SP, 0.4429923648839117E+00_SP, &
         0.1530748903554898E+00_SP, 0.2167416581882652E-03_SP, 0.4774428052721736E+00_SP, &
         0.1851176436721877E+00_SP, 0.2206430730516600E-03_SP, 0.5107446539535904E+00_SP, &
         0.2170829107658179E+00_SP, 0.2237186938699523E-03_SP, 0.5428151370542935E+00_SP, &
         0.2487786689026271E+00_SP, 0.2260480075032884E-03_SP, 0.5735699292556964E+00_SP, &
         0.2800239952795016E+00_SP, 0.2277098884558542E-03_SP, 0.6029253794562866E+00_SP, &
         0.3106445702878119E+00_SP, 0.2287845715109671E-03_SP, 0.6307998987073145E+00_SP, &
         0.3404689500841194E+00_SP, 0.2293547268236294E-03_SP, 0.3752652273692719E+00_SP, &
         0.2997145098184479E-01_SP, 0.2056073839852528E-03_SP, 0.4135383879344028E+00_SP, &
         0.6086725898678011E-01_SP, 0.2114235865831876E-03_SP, 0.4506113885153907E+00_SP, &
         0.9238849548435643E-01_SP, 0.2163175629770551E-03_SP, 0.4864401554606072E+00_SP, &
         0.1242786603851851E+00_SP, 0.2203392158111650E-03_SP, 0.5209708076611709E+00_SP, &
         0.1563086731483386E+00_SP, 0.2235473176847839E-03_SP, 0.5541422135830122E+00_SP, &
         0.1882696509388506E+00_SP, 0.2260024141501235E-03_SP, 0.5858880915113817E+00_SP, &
         0.2199672979126059E+00_SP, 0.2277675929329182E-03_SP, 0.6161399390603444E+00_SP, &
         0.2512165482924867E+00_SP, 0.2289102112284834E-03_SP, 0.6448296482255090E+00_SP, &
         0.2818368701871888E+00_SP, 0.2295027954625118E-03_SP, 0.4544796274917948E+00_SP, &
         0.3088970405060312E-01_SP, 0.2161281589879992E-03_SP, 0.4919389072146628E+00_SP, &
         0.6240947677636835E-01_SP, 0.2201980477395102E-03_SP, 0.5279313026985183E+00_SP, &
         0.9430706144280313E-01_SP, 0.2234952066593166E-03_SP, 0.5624169925571135E+00_SP, &
         0.1263547818770374E+00_SP, 0.2260540098520838E-03_SP, 0.5953484627093287E+00_SP, &
         0.1583430788822594E+00_SP, 0.2279157981899988E-03_SP, 0.6266730715339185E+00_SP, &
         0.1900748462555988E+00_SP, 0.2291296918565571E-03_SP, 0.6563363204278871E+00_SP, &
         0.2213599519592567E+00_SP, 0.2297533752536649E-03_SP, 0.5314574716585696E+00_SP, &
         0.3152508811515374E-01_SP, 0.2234927356465995E-03_SP, 0.5674614932298185E+00_SP, &
         0.6343865291465561E-01_SP, 0.2261288012985219E-03_SP, 0.6017706004970264E+00_SP, &
         0.9551503504223951E-01_SP, 0.2280818160923688E-03_SP, 0.6343471270264178E+00_SP, &
         0.1275440099801196E+00_SP, 0.2293773295180159E-03_SP, 0.6651494599127802E+00_SP, &
         0.1593252037671960E+00_SP, 0.2300528767338634E-03_SP, 0.6050184986005704E+00_SP, &
         0.3192538338496105E-01_SP, 0.2281893855065666E-03_SP, 0.6390163550880400E+00_SP, &
         0.6402824353962306E-01_SP, 0.2295720444840727E-03_SP, 0.6711199107088448E+00_SP, &
         0.9609805077002909E-01_SP, 0.2303227649026753E-03_SP, 0.6741354429572275E+00_SP, &
         0.3211853196273233E-01_SP, 0.2304831913227114E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_4802_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_4802_pt

  Subroutine Lebedev_5294_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 125, 1, 1, 0, 30, 10, 90, 5294 /)
    REAL(SP),    PARAMETER :: pa(352) = (/ &
         0.9080510764308163E-04_SP, 0.2084824361987793E-03_SP, 0.2303261686261450E-01_SP, &
         0.5011105657239616E-04_SP, 0.3757208620162394E-01_SP, 0.5942520409683854E-04_SP, &
         0.5821912033821852E-01_SP, 0.9564394826109721E-04_SP, 0.8403127529194872E-01_SP, &
         0.1185530657126338E-03_SP, 0.1122927798060578E+00_SP, 0.1364510114230331E-03_SP, &
         0.1420125319192987E+00_SP, 0.1505828825605415E-03_SP, 0.1726396437341978E+00_SP, &
         0.1619298749867023E-03_SP, 0.2038170058115696E+00_SP, 0.1712450504267789E-03_SP, &
         0.2352849892876508E+00_SP, 0.1789891098164999E-03_SP, 0.2668363354312461E+00_SP, &
         0.1854474955629795E-03_SP, 0.2982941279900452E+00_SP, 0.1908148636673661E-03_SP, &
         0.3295002922087076E+00_SP, 0.1952377405281833E-03_SP, 0.3603094918363593E+00_SP, &
         0.1988349254282232E-03_SP, 0.3905857895173920E+00_SP, 0.2017079807160050E-03_SP, &
         0.4202005758160837E+00_SP, 0.2039473082709094E-03_SP, 0.4490310061597227E+00_SP, &
         0.2056360279288953E-03_SP, 0.4769586160311491E+00_SP, 0.2068525823066865E-03_SP, &
         0.5038679887049750E+00_SP, 0.2076724877534488E-03_SP, 0.5296454286519961E+00_SP, &
         0.2081694278237885E-03_SP, 0.5541776207164850E+00_SP, 0.2084157631219326E-03_SP, &
         0.5990467321921213E+00_SP, 0.2084381531128593E-03_SP, 0.6191467096294587E+00_SP, &
         0.2083476277129307E-03_SP, 0.6375251212901849E+00_SP, 0.2082686194459732E-03_SP, &
         0.6540514381131168E+00_SP, 0.2082475686112415E-03_SP, 0.6685899064391510E+00_SP, &
         0.2083139860289915E-03_SP, 0.6810013009681648E+00_SP, 0.2084745561831237E-03_SP, &
         0.6911469578730340E+00_SP, 0.2087091313375890E-03_SP, 0.6988956915141736E+00_SP, &
         0.2089718413297697E-03_SP, 0.7041335794868720E+00_SP, 0.2092003303479793E-03_SP, &
         0.7067754398018567E+00_SP, 0.2093336148263241E-03_SP, 0.3840368707853623E-01_SP, &
         0.7591708117365267E-04_SP, 0.9835485954117399E-01_SP, 0.1083383968169186E-03_SP, &
         0.1665774947612998E+00_SP, 0.1403019395292510E-03_SP, 0.2405702335362910E+00_SP, &
         0.1615970179286436E-03_SP, 0.3165270770189046E+00_SP, 0.1771144187504911E-03_SP, &
         0.3927386145645443E+00_SP, 0.1887760022988168E-03_SP, 0.4678825918374656E+00_SP, &
         0.1973474670768214E-03_SP, 0.5408022024266935E+00_SP, 0.2033787661234659E-03_SP, &
         0.6104967445752438E+00_SP, 0.2072343626517331E-03_SP, 0.6760910702685738E+00_SP, &
         0.2091177834226918E-03_SP, 0.6655644120217392E-01_SP, 0.1936508874588424E-01_SP, &
         0.9316684484675566E-04_SP, 0.9446246161270182E-01_SP, 0.4252442002115869E-01_SP, &
         0.1116193688682976E-03_SP, 0.1242651925452509E+00_SP, 0.6806529315354374E-01_SP, &
         0.1298623551559414E-03_SP, 0.1553438064846751E+00_SP, 0.9560957491205369E-01_SP, &
         0.1450236832456426E-03_SP, 0.1871137110542670E+00_SP, 0.1245931657452888E+00_SP, &
         0.1572719958149914E-03_SP, 0.2192612628836257E+00_SP, 0.1545385828778978E+00_SP, &
         0.1673234785867195E-03_SP, 0.2515682807206955E+00_SP, 0.1851004249723368E+00_SP, &
         0.1756860118725188E-03_SP, 0.2838535866287290E+00_SP, 0.2160182608272384E+00_SP, &
         0.1826776290439367E-03_SP, 0.3159578817528521E+00_SP, 0.2470799012277111E+00_SP, &
         0.1885116347992865E-03_SP, 0.3477370882791392E+00_SP, 0.2781014208986402E+00_SP, &
         0.1933457860170574E-03_SP, 0.3790576960890540E+00_SP, 0.3089172523515731E+00_SP, &
         0.1973060671902064E-03_SP, 0.4097938317810200E+00_SP, 0.3393750055472244E+00_SP, &
         0.2004987099616311E-03_SP, 0.4398256572859637E+00_SP, 0.3693322470987730E+00_SP, &
         0.2030170909281499E-03_SP, 0.4690384114718480E+00_SP, 0.3986541005609877E+00_SP, &
         0.2049461460119080E-03_SP, 0.4973216048301053E+00_SP, 0.4272112491408562E+00_SP, &
         0.2063653565200186E-03_SP, 0.5245681526132446E+00_SP, 0.4548781735309936E+00_SP, &
         0.2073507927381027E-03_SP, 0.5506733911803888E+00_SP, 0.4815315355023251E+00_SP, &
         0.2079764593256122E-03_SP, 0.5755339829522475E+00_SP, 0.5070486445801855E+00_SP, &
         0.2083150534968778E-03_SP, 0.1305472386056362E+00_SP, 0.2284970375722366E-01_SP, &
         0.1262715121590664E-03_SP, 0.1637327908216477E+00_SP, 0.4812254338288384E-01_SP, &
         0.1414386128545972E-03_SP, 0.1972734634149637E+00_SP, 0.7531734457511935E-01_SP, &
         0.1538740401313898E-03_SP, 0.2308694653110130E+00_SP, 0.1039043639882017E+00_SP, &
         0.1642434942331432E-03_SP, 0.2643899218338160E+00_SP, 0.1334526587117626E+00_SP, &
         0.1729790609237496E-03_SP, 0.2977171599622171E+00_SP, 0.1636414868936382E+00_SP, &
         0.1803505190260828E-03_SP, 0.3307293903032310E+00_SP, 0.1942195406166568E+00_SP, &
         0.1865475350079657E-03_SP, 0.3633069198219073E+00_SP, 0.2249752879943753E+00_SP, &
         0.1917182669679069E-03_SP, 0.3953346955922727E+00_SP, 0.2557218821820032E+00_SP, &
         0.1959851709034382E-03_SP, 0.4267018394184914E+00_SP, 0.2862897925213193E+00_SP, &
         0.1994529548117882E-03_SP, 0.4573009622571704E+00_SP, 0.3165224536636518E+00_SP, &
         0.2022138911146548E-03_SP, 0.4870279559856109E+00_SP, 0.3462730221636496E+00_SP, &
         0.2043518024208592E-03_SP, 0.5157819581450322E+00_SP, 0.3754016870282835E+00_SP, &
         0.2059450313018110E-03_SP, 0.5434651666465393E+00_SP, 0.4037733784993613E+00_SP, &
         0.2070685715318472E-03_SP, 0.5699823887764627E+00_SP, 0.4312557784139123E+00_SP, &
         0.2077955310694373E-03_SP, 0.5952403350947741E+00_SP, 0.4577175367122110E+00_SP, &
         0.2081980387824712E-03_SP, 0.2025152599210369E+00_SP, 0.2520253617719557E-01_SP, &
         0.1521318610377956E-03_SP, 0.2381066653274425E+00_SP, 0.5223254506119000E-01_SP, &
         0.1622772720185755E-03_SP, 0.2732823383651612E+00_SP, 0.8060669688588620E-01_SP, &
         0.1710498139420709E-03_SP, 0.3080137692611118E+00_SP, 0.1099335754081255E+00_SP, &
         0.1785911149448736E-03_SP, 0.3422405614587601E+00_SP, 0.1399120955959857E+00_SP, &
         0.1850125313687736E-03_SP, 0.3758808773890420E+00_SP, 0.1702977801651705E+00_SP, &
         0.1904229703933298E-03_SP, 0.4088458383438932E+00_SP, 0.2008799256601680E+00_SP, &
         0.1949259956121987E-03_SP, 0.4410450550841152E+00_SP, 0.2314703052180836E+00_SP, &
         0.1986161545363960E-03_SP, 0.4723879420561312E+00_SP, 0.2618972111375892E+00_SP, &
         0.2015790585641370E-03_SP, 0.5027843561874343E+00_SP, 0.2920013195600270E+00_SP, &
         0.2038934198707418E-03_SP, 0.5321453674452458E+00_SP, 0.3216322555190551E+00_SP, &
         0.2056334060538251E-03_SP, 0.5603839113834030E+00_SP, 0.3506456615934198E+00_SP, &
         0.2068705959462289E-03_SP, 0.5874150706875146E+00_SP, 0.3789007181306267E+00_SP, &
         0.2076753906106002E-03_SP, 0.6131559381660038E+00_SP, 0.4062580170572782E+00_SP, &
         0.2081179391734803E-03_SP, 0.2778497016394506E+00_SP, 0.2696271276876226E-01_SP, &
         0.1700345216228943E-03_SP, 0.3143733562261912E+00_SP, 0.5523469316960465E-01_SP, &
         0.1774906779990410E-03_SP, 0.3501485810261827E+00_SP, 0.8445193201626464E-01_SP, &
         0.1839659377002642E-03_SP, 0.3851430322303653E+00_SP, 0.1143263119336083E+00_SP, &
         0.1894987462975169E-03_SP, 0.4193013979470415E+00_SP, 0.1446177898344475E+00_SP, &
         0.1941548809452595E-03_SP, 0.4525585960458567E+00_SP, 0.1751165438438091E+00_SP, &
         0.1980078427252384E-03_SP, 0.4848447779622947E+00_SP, 0.2056338306745660E+00_SP, &
         0.2011296284744488E-03_SP, 0.5160871208276894E+00_SP, 0.2359965487229226E+00_SP, &
         0.2035888456966776E-03_SP, 0.5462112185696926E+00_SP, 0.2660430223139146E+00_SP, &
         0.2054516325352142E-03_SP, 0.5751425068101757E+00_SP, 0.2956193664498032E+00_SP, &
         0.2067831033092635E-03_SP, 0.6028073872853596E+00_SP, 0.3245763905312779E+00_SP, &
         0.2076485320284876E-03_SP, 0.6291338275278409E+00_SP, 0.3527670026206972E+00_SP, &
         0.2081141439525255E-03_SP, 0.3541797528439391E+00_SP, 0.2823853479435550E-01_SP, &
         0.1834383015469222E-03_SP, 0.3908234972074657E+00_SP, 0.5741296374713106E-01_SP, &
         0.1889540591777677E-03_SP, 0.4264408450107590E+00_SP, 0.8724646633650199E-01_SP, &
         0.1936677023597375E-03_SP, 0.4609949666553286E+00_SP, 0.1175034422915616E+00_SP, &
         0.1976176495066504E-03_SP, 0.4944389496536006E+00_SP, 0.1479755652628428E+00_SP, &
         0.2008536004560983E-03_SP, 0.5267194884346086E+00_SP, 0.1784740659484352E+00_SP, &
         0.2034280351712291E-03_SP, 0.5577787810220990E+00_SP, 0.2088245700431244E+00_SP, &
         0.2053944466027758E-03_SP, 0.5875563763536670E+00_SP, 0.2388628136570763E+00_SP, &
         0.2068077642882360E-03_SP, 0.6159910016391269E+00_SP, 0.2684308928769185E+00_SP, &
         0.2077250949661599E-03_SP, 0.6430219602956268E+00_SP, 0.2973740761960252E+00_SP, &
         0.2082062440705320E-03_SP, 0.4300647036213646E+00_SP, 0.2916399920493977E-01_SP, &
         0.1934374486546626E-03_SP, 0.4661486308935531E+00_SP, 0.5898803024755659E-01_SP, &
         0.1974107010484300E-03_SP, 0.5009658555287261E+00_SP, 0.8924162698525409E-01_SP, &
         0.2007129290388658E-03_SP, 0.5344824270447704E+00_SP, 0.1197185199637321E+00_SP, &
         0.2033736947471293E-03_SP, 0.5666575997416371E+00_SP, 0.1502300756161382E+00_SP, &
         0.2054287125902493E-03_SP, 0.5974457471404752E+00_SP, 0.1806004191913564E+00_SP, &
         0.2069184936818894E-03_SP, 0.6267984444116886E+00_SP, 0.2106621764786252E+00_SP, &
         0.2078883689808782E-03_SP, 0.6546664713575417E+00_SP, 0.2402526932671914E+00_SP, &
         0.2083886366116359E-03_SP, 0.5042711004437253E+00_SP, 0.2982529203607657E-01_SP, &
         0.2006593275470817E-03_SP, 0.5392127456774380E+00_SP, 0.6008728062339922E-01_SP, &
         0.2033728426135397E-03_SP, 0.5726819437668618E+00_SP, 0.9058227674571398E-01_SP, &
         0.2055008781377608E-03_SP, 0.6046469254207278E+00_SP, 0.1211219235803400E+00_SP, &
         0.2070651783518502E-03_SP, 0.6350716157434952E+00_SP, 0.1515286404791580E+00_SP, &
         0.2080953335094320E-03_SP, 0.6639177679185454E+00_SP, 0.1816314681255552E+00_SP, &
         0.2086284998988521E-03_SP, 0.5757276040972253E+00_SP, 0.3026991752575440E-01_SP, &
         0.2055549387644668E-03_SP, 0.6090265823139755E+00_SP, 0.6078402297870770E-01_SP, &
         0.2071871850267654E-03_SP, 0.6406735344387661E+00_SP, 0.9135459984176636E-01_SP, &
         0.2082856600431965E-03_SP, 0.6706397927793709E+00_SP, 0.1218024155966590E+00_SP, &
         0.2088705858819358E-03_SP, 0.6435019674426665E+00_SP, 0.3052608357660639E-01_SP, &
         0.2083995867536322E-03_SP, 0.6747218676375681E+00_SP, 0.6112185773983089E-01_SP, &
         0.2090509712889637E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_5294_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_5294_pt

  Subroutine Lebedev_5810_pt(n,pt,wt)
    USE DataTypes
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    INTEGER(I4B),PARAMETER :: nk(8) = (/ 131, 1, 1, 1, 31, 10, 100, 5810 /)
    REAL(SP),    PARAMETER :: pa(385) = (/ &
         0.9735347946175486E-05_SP, 0.1901059546737578E-03_SP, 0.1907581241803167E-03_SP, &
         0.1182361662400277E-01_SP, 0.3926424538919212E-04_SP, 0.3062145009138958E-01_SP, &
         0.6667905467294382E-04_SP, 0.5329794036834243E-01_SP, 0.8868891315019135E-04_SP, &
         0.7848165532862220E-01_SP, 0.1066306000958872E-03_SP, 0.1054038157636201E+00_SP, &
         0.1214506743336128E-03_SP, 0.1335577797766211E+00_SP, 0.1338054681640871E-03_SP, &
         0.1625769955502252E+00_SP, 0.1441677023628504E-03_SP, 0.1921787193412792E+00_SP, &
         0.1528880200826557E-03_SP, 0.2221340534690548E+00_SP, 0.1602330623773609E-03_SP, &
         0.2522504912791132E+00_SP, 0.1664102653445244E-03_SP, 0.2823610860679697E+00_SP, &
         0.1715845854011323E-03_SP, 0.3123173966267560E+00_SP, 0.1758901000133069E-03_SP, &
         0.3419847036953789E+00_SP, 0.1794382485256736E-03_SP, 0.3712386456999758E+00_SP, &
         0.1823238106757407E-03_SP, 0.3999627649876828E+00_SP, 0.1846293252959976E-03_SP, &
         0.4280466458648093E+00_SP, 0.1864284079323098E-03_SP, 0.4553844360185711E+00_SP, &
         0.1877882694626914E-03_SP, 0.4818736094437834E+00_SP, 0.1887716321852025E-03_SP, &
         0.5074138709260629E+00_SP, 0.1894381638175673E-03_SP, 0.5319061304570707E+00_SP, &
         0.1898454899533629E-03_SP, 0.5552514978677286E+00_SP, 0.1900497929577815E-03_SP, &
         0.5981009025246183E+00_SP, 0.1900671501924092E-03_SP, 0.6173990192228116E+00_SP, &
         0.1899837555533510E-03_SP, 0.6351365239411131E+00_SP, 0.1899014113156229E-03_SP, &
         0.6512010228227200E+00_SP, 0.1898581257705106E-03_SP, 0.6654758363948120E+00_SP, &
         0.1898804756095753E-03_SP, 0.6778410414853370E+00_SP, 0.1899793610426402E-03_SP, &
         0.6881760887484110E+00_SP, 0.1901464554844117E-03_SP, 0.6963645267094598E+00_SP, &
         0.1903533246259542E-03_SP, 0.7023010617153579E+00_SP, 0.1905556158463228E-03_SP, &
         0.7059004636628753E+00_SP, 0.1907037155663528E-03_SP, 0.3552470312472575E-01_SP, &
         0.5992997844249967E-04_SP, 0.9151176620841283E-01_SP, 0.9749059382456978E-04_SP, &
         0.1566197930068980E+00_SP, 0.1241680804599158E-03_SP, 0.2265467599271907E+00_SP, &
         0.1437626154299360E-03_SP, 0.2988242318581361E+00_SP, 0.1584200054793902E-03_SP, &
         0.3717482419703886E+00_SP, 0.1694436550982744E-03_SP, 0.4440094491758889E+00_SP, &
         0.1776617014018108E-03_SP, 0.5145337096756642E+00_SP, 0.1836132434440077E-03_SP, &
         0.5824053672860230E+00_SP, 0.1876494727075983E-03_SP, 0.6468283961043370E+00_SP, &
         0.1899906535336482E-03_SP, 0.6095964259104373E-01_SP, 0.1787828275342931E-01_SP, &
         0.8143252820767350E-04_SP, 0.8811962270959388E-01_SP, 0.3953888740792096E-01_SP, &
         0.9998859890887728E-04_SP, 0.1165936722428831E+00_SP, 0.6378121797722990E-01_SP, &
         0.1156199403068359E-03_SP, 0.1460232857031785E+00_SP, 0.8985890813745037E-01_SP, &
         0.1287632092635513E-03_SP, 0.1761197110181755E+00_SP, 0.1172606510576162E+00_SP, &
         0.1398378643365139E-03_SP, 0.2066471190463718E+00_SP, 0.1456102876970995E+00_SP, &
         0.1491876468417391E-03_SP, 0.2374076026328152E+00_SP, 0.1746153823011775E+00_SP, &
         0.1570855679175456E-03_SP, 0.2682305474337051E+00_SP, 0.2040383070295584E+00_SP, &
         0.1637483948103775E-03_SP, 0.2989653312142369E+00_SP, 0.2336788634003698E+00_SP, &
         0.1693500566632843E-03_SP, 0.3294762752772209E+00_SP, 0.2633632752654219E+00_SP, &
         0.1740322769393633E-03_SP, 0.3596390887276086E+00_SP, 0.2929369098051601E+00_SP, &
         0.1779126637278296E-03_SP, 0.3893383046398812E+00_SP, 0.3222592785275512E+00_SP, &
         0.1810908108835412E-03_SP, 0.4184653789358347E+00_SP, 0.3512004791195743E+00_SP, &
         0.1836529132600190E-03_SP, 0.4469172319076166E+00_SP, 0.3796385677684537E+00_SP, &
         0.1856752841777379E-03_SP, 0.4745950813276976E+00_SP, 0.4074575378263879E+00_SP, &
         0.1872270566606832E-03_SP, 0.5014034601410262E+00_SP, 0.4345456906027828E+00_SP, &
         0.1883722645591307E-03_SP, 0.5272493404551239E+00_SP, 0.4607942515205134E+00_SP, &
         0.1891714324525297E-03_SP, 0.5520413051846366E+00_SP, 0.4860961284181720E+00_SP, &
         0.1896827480450146E-03_SP, 0.5756887237503077E+00_SP, 0.5103447395342790E+00_SP, &
         0.1899628417059528E-03_SP, 0.1225039430588352E+00_SP, 0.2136455922655793E-01_SP, &
         0.1123301829001669E-03_SP, 0.1539113217321372E+00_SP, 0.4520926166137188E-01_SP, &
         0.1253698826711277E-03_SP, 0.1856213098637712E+00_SP, 0.7086468177864818E-01_SP, &
         0.1366266117678531E-03_SP, 0.2174998728035131E+00_SP, 0.9785239488772918E-01_SP, &
         0.1462736856106918E-03_SP, 0.2494128336938330E+00_SP, 0.1258106396267210E+00_SP, &
         0.1545076466685412E-03_SP, 0.2812321562143480E+00_SP, 0.1544529125047001E+00_SP, &
         0.1615096280814007E-03_SP, 0.3128372276456111E+00_SP, 0.1835433512202753E+00_SP, &
         0.1674366639741759E-03_SP, 0.3441145160177973E+00_SP, 0.2128813258619585E+00_SP, &
         0.1724225002437900E-03_SP, 0.3749567714853510E+00_SP, 0.2422913734880829E+00_SP, &
         0.1765810822987288E-03_SP, 0.4052621732015610E+00_SP, 0.2716163748391453E+00_SP, &
         0.1800104126010751E-03_SP, 0.4349335453522385E+00_SP, 0.3007127671240280E+00_SP, &
         0.1827960437331284E-03_SP, 0.4638776641524965E+00_SP, 0.3294470677216479E+00_SP, &
         0.1850140300716308E-03_SP, 0.4920046410462687E+00_SP, 0.3576932543699155E+00_SP, &
         0.1867333507394938E-03_SP, 0.5192273554861704E+00_SP, 0.3853307059757764E+00_SP, &
         0.1880178688638289E-03_SP, 0.5454609081136522E+00_SP, 0.4122425044452694E+00_SP, &
         0.1889278925654758E-03_SP, 0.5706220661424140E+00_SP, 0.4383139587781027E+00_SP, &
         0.1895213832507346E-03_SP, 0.5946286755181518E+00_SP, 0.4634312536300553E+00_SP, &
         0.1898548277397420E-03_SP, 0.1905370790924295E+00_SP, 0.2371311537781979E-01_SP, &
         0.1349105935937341E-03_SP, 0.2242518717748009E+00_SP, 0.4917878059254806E-01_SP, &
         0.1444060068369326E-03_SP, 0.2577190808025936E+00_SP, 0.7595498960495142E-01_SP, &
         0.1526797390930008E-03_SP, 0.2908724534927187E+00_SP, 0.1036991083191100E+00_SP, &
         0.1598208771406474E-03_SP, 0.3236354020056219E+00_SP, 0.1321348584450234E+00_SP, &
         0.1659354368615331E-03_SP, 0.3559267359304543E+00_SP, 0.1610316571314789E+00_SP, &
         0.1711279910946440E-03_SP, 0.3876637123676956E+00_SP, 0.1901912080395707E+00_SP, &
         0.1754952725601440E-03_SP, 0.4187636705218842E+00_SP, 0.2194384950137950E+00_SP, &
         0.1791247850802529E-03_SP, 0.4491449019883107E+00_SP, 0.2486155334763858E+00_SP, &
         0.1820954300877716E-03_SP, 0.4787270932425445E+00_SP, 0.2775768931812335E+00_SP, &
         0.1844788524548449E-03_SP, 0.5074315153055574E+00_SP, 0.3061863786591120E+00_SP, &
         0.1863409481706220E-03_SP, 0.5351810507738336E+00_SP, 0.3343144718152556E+00_SP, &
         0.1877433008795068E-03_SP, 0.5619001025975381E+00_SP, 0.3618362729028427E+00_SP, &
         0.1887444543705232E-03_SP, 0.5875144035268046E+00_SP, 0.3886297583620408E+00_SP, &
         0.1894009829375006E-03_SP, 0.6119507308734495E+00_SP, 0.4145742277792031E+00_SP, &
         0.1897683345035198E-03_SP, 0.2619733870119463E+00_SP, 0.2540047186389353E-01_SP, &
         0.1517327037467653E-03_SP, 0.2968149743237949E+00_SP, 0.5208107018543989E-01_SP, &
         0.1587740557483543E-03_SP, 0.3310451504860488E+00_SP, 0.7971828470885599E-01_SP, &
         0.1649093382274097E-03_SP, 0.3646215567376676E+00_SP, 0.1080465999177927E+00_SP, &
         0.1701915216193265E-03_SP, 0.3974916785279360E+00_SP, 0.1368413849366629E+00_SP, &
         0.1746847753144065E-03_SP, 0.4295967403772029E+00_SP, 0.1659073184763559E+00_SP, &
         0.1784555512007570E-03_SP, 0.4608742854473447E+00_SP, 0.1950703730454614E+00_SP, &
         0.1815687562112174E-03_SP, 0.4912598858949903E+00_SP, 0.2241721144376724E+00_SP, &
         0.1840864370663302E-03_SP, 0.5206882758945558E+00_SP, 0.2530655255406489E+00_SP, &
         0.1860676785390006E-03_SP, 0.5490940914019819E+00_SP, 0.2816118409731066E+00_SP, &
         0.1875690583743703E-03_SP, 0.5764123302025542E+00_SP, 0.3096780504593238E+00_SP, &
         0.1886453236347225E-03_SP, 0.6025786004213506E+00_SP, 0.3371348366394987E+00_SP, &
         0.1893501123329645E-03_SP, 0.6275291964794956E+00_SP, 0.3638547827694396E+00_SP, &
         0.1897366184519868E-03_SP, 0.3348189479861771E+00_SP, 0.2664841935537443E-01_SP, &
         0.1643908815152736E-03_SP, 0.3699515545855295E+00_SP, 0.5424000066843495E-01_SP, &
         0.1696300350907768E-03_SP, 0.4042003071474669E+00_SP, 0.8251992715430854E-01_SP, &
         0.1741553103844483E-03_SP, 0.4375320100182624E+00_SP, 0.1112695182483710E+00_SP, &
         0.1780015282386092E-03_SP, 0.4699054490335947E+00_SP, 0.1402964116467816E+00_SP, &
         0.1812116787077125E-03_SP, 0.5012739879431952E+00_SP, 0.1694275117584291E+00_SP, &
         0.1838323158085421E-03_SP, 0.5315874883754966E+00_SP, 0.1985038235312689E+00_SP, &
         0.1859113119837737E-03_SP, 0.5607937109622117E+00_SP, 0.2273765660020893E+00_SP, &
         0.1874969220221698E-03_SP, 0.5888393223495521E+00_SP, 0.2559041492849764E+00_SP, &
         0.1886375612681076E-03_SP, 0.6156705979160163E+00_SP, 0.2839497251976899E+00_SP, &
         0.1893819575809276E-03_SP, 0.6412338809078123E+00_SP, 0.3113791060500690E+00_SP, &
         0.1897794748256767E-03_SP, 0.4076051259257167E+00_SP, 0.2757792290858463E-01_SP, &
         0.1738963926584846E-03_SP, 0.4423788125791520E+00_SP, 0.5584136834984293E-01_SP, &
         0.1777442359873466E-03_SP, 0.4760480917328258E+00_SP, 0.8457772087727143E-01_SP, &
         0.1810010815068719E-03_SP, 0.5085838725946297E+00_SP, 0.1135975846359248E+00_SP, &
         0.1836920318248129E-03_SP, 0.5399513637391218E+00_SP, 0.1427286904765053E+00_SP, &
         0.1858489473214328E-03_SP, 0.5701118433636380E+00_SP, 0.1718112740057635E+00_SP, &
         0.1875079342496592E-03_SP, 0.5990240530606021E+00_SP, 0.2006944855985351E+00_SP, &
         0.1887080239102310E-03_SP, 0.6266452685139695E+00_SP, 0.2292335090598907E+00_SP, &
         0.1894905752176822E-03_SP, 0.6529320971415942E+00_SP, 0.2572871512353714E+00_SP, &
         0.1898991061200695E-03_SP, 0.4791583834610126E+00_SP, 0.2826094197735932E-01_SP, &
         0.1809065016458791E-03_SP, 0.5130373952796940E+00_SP, 0.5699871359683649E-01_SP, &
         0.1836297121596799E-03_SP, 0.5456252429628476E+00_SP, 0.8602712528554394E-01_SP, &
         0.1858426916241869E-03_SP, 0.5768956329682385E+00_SP, 0.1151748137221281E+00_SP, &
         0.1875654101134641E-03_SP, 0.6068186944699046E+00_SP, 0.1442811654136362E+00_SP, &
         0.1888240751833503E-03_SP, 0.6353622248024907E+00_SP, 0.1731930321657680E+00_SP, &
         0.1896497383866979E-03_SP, 0.6624927035731797E+00_SP, 0.2017619958756061E+00_SP, &
         0.1900775530219121E-03_SP, 0.5484933508028488E+00_SP, 0.2874219755907391E-01_SP, &
         0.1858525041478814E-03_SP, 0.5810207682142106E+00_SP, 0.5778312123713695E-01_SP, &
         0.1876248690077947E-03_SP, 0.6120955197181352E+00_SP, 0.8695262371439526E-01_SP, &
         0.1889404439064607E-03_SP, 0.6416944284294319E+00_SP, 0.1160893767057166E+00_SP, &
         0.1898168539265290E-03_SP, 0.6697926391731260E+00_SP, 0.1450378826743251E+00_SP, &
         0.1902779940661772E-03_SP, 0.6147594390585488E+00_SP, 0.2904957622341456E-01_SP, &
         0.1890125641731815E-03_SP, 0.6455390026356783E+00_SP, 0.5823809152617197E-01_SP, &
         0.1899434637795751E-03_SP, 0.6747258588365477E+00_SP, 0.8740384899884715E-01_SP, &
         0.1904520856831751E-03_SP, 0.6772135750395347E+00_SP, 0.2919946135808105E-01_SP, &
         0.1905534498734563E-03_SP /)

    IF(n /= nk(8))THEN
       call error('Lebedev_5810_pt: Bad Mesh Request. '//&
            'Expected npts does not equal actual npts',i1=n,i2=nk(8))
    ENDIF
    CALL Lebedev_Oh_Mesh(nk, pa, pt, wt)
  end Subroutine Lebedev_5810_pt


!!$===========================================================================
  !  SUBROUTINE Lebedev_Oh_Mesh(  nk,    pa,   pt,    wt)
  !
  !  Generates integration mesh for integrations over the unit sphere
  !  with octahedral symmetry.
  !
  !  Input:                                       
  !   nk(8)  vector [nk(2:7)] defining mesh rules to use (6,8,12,24a,24b,48)
  !          nk(1)=Lmax, nk(8)=npts
  !   pa(*)  vector defining parameters for meshing rules
  !
  !  Output:
  !   pt(3,*)   cartesian points on the unit sphere
  !   wt(*)     integration weight factors
!$----------------------------------------------------------------------------

!!!****f* AngularQuadratureMod/Lebedev_Oh_Mesh
!!!
!!! NAME
!!!     Lebedev_Oh_Mesh -- Obtain possible quadrature values between [MinVal, MaxVal] 
!!! USAGE
!!!     Lebedev_Oh_Mesh(  nk,    pa,   pt,    wt)
!!! DESCRIPTION
!!!     Generates integration meshs for integrations over the unit sphere
!!!     with octahedral symmetry.
!!! INPUTS
!!!     INTEGER(I4B), INTENT(IN) :: nk(8)  
!!!               nk(1) = Lmax, 
!!!               nk(2:7) mesh rules to use (6,8,12,24a,24b,48), and how many times
!!!               nk(8) = npts
!!!     REAL(SP), INTENT(IN) :: pa(:) - vector defining parameters for meshing rules
!!! OUTPUTS
!!!     pt(3,nk(8)) -  cartesian points on the unit sphere
!!!     wt(nk(8))   -  integration weight factors
!!!***
  SUBROUTINE Lebedev_Oh_Mesh(  nk,    pa,   pt,    wt)
    USE DataTypes
    USE ConstantsMod, only : SQRT2,SQRT3
    USE ErrorMod, only : error
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: nk(:)
    REAL(SP), INTENT(IN) :: pa(:)
    REAL(SP), INTENT(OUT) :: pt(:,:), wt(:)
    REAL(SP), PARAMETER :: c08pt = 1.0_SP/SQRT3
    REAL(SP), PARAMETER :: c12pt = 1.0_SP/SQRT2

    INTEGER(I4B) :: idx, ip, n1,jj
    REAL(SP) :: rr, ss, tt

    idx = 0
    ip = 0
    
    ! There can be only 1 of the 6pt, 8pt and 12 pt
    ! rules since their locations are fixed
    IF (nk(2)  >  0) THEN   ! 6 points
       ip = ip + 1
       pt(1:3,idx+1) = (/ 1.0_SP, 0.0_SP, 0.0_SP /)
       pt(1:3,idx+2) = (/-1.0_SP, 0.0_SP, 0.0_SP /)
       pt(1:3,idx+3) = (/ 0.0_SP, 1.0_SP, 0.0_SP /)
       pt(1:3,idx+4) = (/ 0.0_SP,-1.0_SP, 0.0_SP /)
       pt(1:3,idx+5) = (/ 0.0_SP, 0.0_SP, 1.0_SP /)
       pt(1:3,idx+6) = (/ 0.0_SP, 0.0_SP,-1.0_SP /)
       wt(idx+1:idx+6)   = pa(ip)
       idx=idx+6
    ENDIF

    IF (nk(3)  >  0) THEN   ! 8 points
       ip = ip + 1
       pt(1:3,idx+1) = (/  c08pt,  c08pt,  c08pt /)
       pt(1:3,idx+2) = (/ -c08pt,  c08pt,  c08pt /)
       pt(1:3,idx+3) = (/  c08pt, -c08pt,  c08pt /)
       pt(1:3,idx+4) = (/ -c08pt, -c08pt,  c08pt /)
       pt(1:3,idx+5) = (/  c08pt,  c08pt, -c08pt /)
       pt(1:3,idx+6) = (/ -c08pt,  c08pt, -c08pt /)
       pt(1:3,idx+7) = (/  c08pt, -c08pt, -c08pt /)
       pt(1:3,idx+8) = (/ -c08pt, -c08pt, -c08pt /)
       wt(idx+1:idx+8) =  pa(ip)
       idx=idx+8
    END IF

    IF (nk(4)  >  0) THEN   ! 12 points
       ip = ip + 1
       pt(1:3,idx+1)  = (/ 0.0_SP,  c12pt,  c12pt /)
       pt(1:3,idx+2)  = (/ 0.0_SP, -c12pt,  c12pt /)
       pt(1:3,idx+3)  = (/ 0.0_SP,  c12pt, -c12pt /)
       pt(1:3,idx+4)  = (/ 0.0_SP, -c12pt, -c12pt /)
       pt(1:3,idx+5)  = (/  c12pt, 0.0_SP,  c12pt /)
       pt(1:3,idx+6)  = (/ -c12pt, 0.0_SP,  c12pt /)
       pt(1:3,idx+7)  = (/  c12pt, 0.0_SP, -c12pt /)
       pt(1:3,idx+8)  = (/ -c12pt, 0.0_SP, -c12pt /)
       pt(1:3,idx+9)  = (/  c12pt,  c12pt, 0.0_SP /)
       pt(1:3,idx+10) = (/ -c12pt,  c12pt, 0.0_SP /)
       pt(1:3,idx+11) = (/  c12pt, -c12pt, 0.0_SP /)
       pt(1:3,idx+12) = (/ -c12pt, -c12pt, 0.0_SP /)
       wt(idx+1:idx+12) =  pa(ip)
       idx=idx+12
    END IF

    n1 = nk(5)              ! 24a points
    DO jj=1,n1
       ip = ip + 1
       rr = pa(ip)
       ss = SQRT(1.0_SP - 2.0_SP*rr*rr) 
       ip = ip + 1
       pt(1:3,idx+1) =  (/  rr,  rr,  ss /)
       pt(1:3,idx+2) =  (/ -rr,  rr,  ss /)
       pt(1:3,idx+3) =  (/  rr, -rr,  ss /)
       pt(1:3,idx+4) =  (/ -rr, -rr,  ss /)
       pt(1:3,idx+5) =  (/  rr,  rr, -ss /)
       pt(1:3,idx+6) =  (/ -rr,  rr, -ss /)
       pt(1:3,idx+7) =  (/  rr, -rr, -ss /)
       pt(1:3,idx+8) =  (/ -rr, -rr, -ss /)
       pt(1:3,idx+9) =  (/  rr,  ss,  rr /)
       pt(1:3,idx+10) = (/ -rr,  ss,  rr /)
       pt(1:3,idx+11) = (/  rr, -ss,  rr /)
       pt(1:3,idx+12) = (/ -rr, -ss,  rr /)
       pt(1:3,idx+13) = (/  rr,  ss, -rr /)
       pt(1:3,idx+14) = (/ -rr,  ss, -rr /)
       pt(1:3,idx+15) = (/  rr, -ss, -rr /)
       pt(1:3,idx+16) = (/ -rr, -ss, -rr /)
       pt(1:3,idx+17) = (/  ss,  rr,  rr /)
       pt(1:3,idx+18) = (/ -ss,  rr,  rr /)
       pt(1:3,idx+19) = (/  ss, -rr,  rr /)
       pt(1:3,idx+20) = (/ -ss, -rr,  rr /)
       pt(1:3,idx+21) = (/  ss,  rr, -rr /)
       pt(1:3,idx+22) = (/ -ss,  rr, -rr /)
       pt(1:3,idx+23) = (/  ss, -rr, -rr /)
       pt(1:3,idx+24) = (/ -ss, -rr, -rr /)
       wt(idx+1:idx+24) =  pa(ip)
       idx=idx+24
    ENDDO

    n1 = nk(6)              ! 24b points
    DO jj=1,n1
       ip = ip + 1
       rr = pa(ip)
       ss = SQRT(1.0_SP - rr*rr)  
       ip = ip + 1
       pt(1:3,idx+1) =  (/  rr,  ss,  0.0_SP /)
       pt(1:3,idx+2) =  (/ -rr,  ss,  0.0_SP /)
       pt(1:3,idx+3) =  (/  rr, -ss,  0.0_SP /)
       pt(1:3,idx+4) =  (/ -rr, -ss,  0.0_SP /)
       pt(1:3,idx+5) =  (/  ss,  rr,  0.0_SP /)
       pt(1:3,idx+6) =  (/ -ss,  rr,  0.0_SP /)
       pt(1:3,idx+7) =  (/  ss, -rr,  0.0_SP /)
       pt(1:3,idx+8) =  (/ -ss, -rr,  0.0_SP /)
       pt(1:3,idx+9) =  (/  rr,  0.0_SP,  ss /)
       pt(1:3,idx+10) = (/ -rr,  0.0_SP,  ss /)
       pt(1:3,idx+11) = (/  rr,  0.0_SP, -ss /)
       pt(1:3,idx+12) = (/ -rr,  0.0_SP, -ss /)
       pt(1:3,idx+13) = (/  ss,  0.0_SP,  rr /)
       pt(1:3,idx+14) = (/ -ss,  0.0_SP,  rr /)
       pt(1:3,idx+15) = (/  ss,  0.0_SP, -rr /)
       pt(1:3,idx+16) = (/ -ss,  0.0_SP, -rr /)
       pt(1:3,idx+17) = (/  0.0_SP,  rr,  ss /)
       pt(1:3,idx+18) = (/  0.0_SP, -rr,  ss /)
       pt(1:3,idx+19) = (/  0.0_SP,  rr, -ss /)
       pt(1:3,idx+20) = (/  0.0_SP, -rr, -ss /)
       pt(1:3,idx+21) = (/  0.0_SP,  ss,  rr /)
       pt(1:3,idx+22) = (/  0.0_SP, -ss,  rr /)
       pt(1:3,idx+23) = (/  0.0_SP,  ss, -rr /)
       pt(1:3,idx+24) = (/  0.0_SP, -ss, -rr /)
       wt(idx+1:idx+24) =  pa(ip)
       idx=idx+24
    END DO

    n1 = nk(7)              ! 48 points
    DO jj=1,n1
       ip = ip + 1
       rr = pa(ip)
       ip = ip + 1
       ss = pa(ip)
       tt = SQRT(1.0_SP - rr*rr - ss*ss)
       ip = ip + 1
       pt(1:3,idx+1) =  (/  rr,  ss,  tt /)
       pt(1:3,idx+2) =  (/ -rr,  ss,  tt /)
       pt(1:3,idx+3) =  (/  rr, -ss,  tt /)
       pt(1:3,idx+4) =  (/ -rr, -ss,  tt /)
       pt(1:3,idx+5) =  (/  rr,  ss, -tt /)
       pt(1:3,idx+6) =  (/ -rr,  ss, -tt /)
       pt(1:3,idx+7) =  (/  rr, -ss, -tt /)
       pt(1:3,idx+8) =  (/ -rr, -ss, -tt /)
       pt(1:3,idx+9) =  (/  rr,  tt,  ss /)
       pt(1:3,idx+10) = (/ -rr,  tt,  ss /)
       pt(1:3,idx+11) = (/  rr, -tt,  ss /)
       pt(1:3,idx+12) = (/ -rr, -tt,  ss /)
       pt(1:3,idx+13) = (/  rr,  tt, -ss /)
       pt(1:3,idx+14) = (/ -rr,  tt, -ss /)
       pt(1:3,idx+15) = (/  rr, -tt, -ss /)
       pt(1:3,idx+16) = (/ -rr, -tt, -ss /)
       pt(1:3,idx+17) = (/  ss,  rr,  tt /)
       pt(1:3,idx+18) = (/ -ss,  rr,  tt /)
       pt(1:3,idx+19) = (/  ss, -rr,  tt /)
       pt(1:3,idx+20) = (/ -ss, -rr,  tt /)
       pt(1:3,idx+21) = (/  ss,  rr, -tt /)
       pt(1:3,idx+22) = (/ -ss,  rr, -tt /)
       pt(1:3,idx+23) = (/  ss, -rr, -tt /)
       pt(1:3,idx+24) = (/ -ss, -rr, -tt /)
       pt(1:3,idx+25) = (/  ss,  tt,  rr /)
       pt(1:3,idx+26) = (/ -ss,  tt,  rr /)
       pt(1:3,idx+27) = (/  ss, -tt,  rr /)
       pt(1:3,idx+28) = (/ -ss, -tt,  rr /)
       pt(1:3,idx+29) = (/  ss,  tt, -rr /)
       pt(1:3,idx+30) = (/ -ss,  tt, -rr /)
       pt(1:3,idx+31) = (/  ss, -tt, -rr /)
       pt(1:3,idx+32) = (/ -ss, -tt, -rr /)
       pt(1:3,idx+33) = (/  tt,  rr,  ss /)
       pt(1:3,idx+34) = (/ -tt,  rr,  ss /)
       pt(1:3,idx+35) = (/  tt, -rr,  ss /)
       pt(1:3,idx+36) = (/ -tt, -rr,  ss /)
       pt(1:3,idx+37) = (/  tt,  rr, -ss /)
       pt(1:3,idx+38) = (/ -tt,  rr, -ss /)
       pt(1:3,idx+39) = (/  tt, -rr, -ss /)
       pt(1:3,idx+40) = (/ -tt, -rr, -ss /)
       pt(1:3,idx+41) = (/  tt,  ss,  rr /)
       pt(1:3,idx+42) = (/ -tt,  ss,  rr /)
       pt(1:3,idx+43) = (/  tt, -ss,  rr /)
       pt(1:3,idx+44) = (/ -tt, -ss,  rr /)
       pt(1:3,idx+45) = (/  tt,  ss, -rr /)
       pt(1:3,idx+46) = (/ -tt,  ss, -rr /)
       pt(1:3,idx+47) = (/  tt, -ss, -rr /)
       pt(1:3,idx+48) = (/ -tt, -ss, -rr /)
       wt(idx+1:idx+48) =  pa(ip)
       idx=idx+48
    END DO
    IF(idx /= nk(8))THEN
       call error('Lebedev_Oh_Mesh: Quatrature rule broken. '//&
            'Expected nmesh does not equal actual nmesh',i1=nk(8),i2=idx)
    ENDIF

  END SUBROUTINE Lebedev_Oh_Mesh


END MODULE angularquadraturemod





MODULE periodictablemod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

   INTEGER(I4B), PARAMETER,PUBLIC :: MAXATOM=54
   CHARACTER(LEN=2),PUBLIC,SAVE :: ElementNameData(1:MAXATOM)
   REAL(SP)         :: AtomicWeightData(1:MAXATOM)
   REAL(SP)         :: StandardIsotopeMassData(1:MAXATOM)
   REAL(SP)         :: CovalentRadiusData(1:MAXATOM)
   REAL(SP)         :: vdWRadiusData(1:MAXATOM)
   REAL(SP)         :: ElectronegativityData(1:MAXATOM)
   REAL(SP)         :: HardnessData(1:MAXATOM)
   REAL(SP)         :: IonizationPotentialData(1:MAXATOM)
   REAL(SP)         :: ElectronAffinityData(1:MAXATOM)
   CHARACTER(LEN=15),PUBLIC,SAVE :: ElementPrettyNameData(1:MAXATOM)


  INTERFACE ElementName
     MODULE PROCEDURE ElementName_iZ !, ElementName_name
  END INTERFACE


CONTAINS
  SUBROUTINE RangeCheck(iZ)
    USE DataTypes
    USE ErrorMod
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: iZ
    IF ( iZ < 1 .OR. iZ > MAXATOM ) THEN 
       CALL error('ERROR: iZ out of range!! ',WARNLev=5,i1=iZ)
    END IF
  END SUBROUTINE RangeCheck

  FUNCTION AtomicNumber(name) RESULT(iZ)
    USE Datatypes
    USE UtilitiesMod
    IMPLICIT NONE
    CHARACTER(LEN=2), INTENT(IN) :: name
    INTEGER(I4B) :: iZ
    iZ = IndexSearch(ElementNameData,name,case_sensitive=.FALSE.)
  END FUNCTION AtomicNumber

  FUNCTION ElementName_iZ(iZ) RESULT(name)
    USE Datatypes
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: iZ
    CHARACTER(LEN=2) :: name
    CALL RangeCheck(iZ)
    name = ElementNameData(iZ)
  END FUNCTION ElementName_iZ

END MODULE periodictablemod





MODULE molecularpropertiesmod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE





CONTAINS
  FUNCTION DipoleMoment(q,Crd,ctr) RESULT(DM)
    USE DataTypes
    Use UtilitiesMod
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: Crd(:,:), q(:)
    REAL(SP), optional, INTENT(IN) :: ctr(3)
    REAL(SP) :: DM(3)
    INTEGER(I4B) :: n,i
    REAL(SP) :: orig(3)
    n = assert_eq(SIZE(q),SIZE(Crd,2),'Dipoleoment')
    if(present(ctr)) THEN
       orig(:)=ctr(:)
    ELSE
       orig(1:3) = 0.0_SP
       DO i=1,n
          orig(1:3) = orig(1:3) + ABS(q(i))*Crd(1:3,i)
       END DO
       ! BG: n should be replaced by SUM(ABS(q(:))) ???
       orig(1:3) = orig(1:3)/n
    ENDIF
    DM(1:3) = 0.0_SP
    DO i=1,n
       DM(1:3) = DM(1:3) + q(i) * (Crd(1:3,i)-orig(1:3))
    END DO
  END FUNCTION DipoleMoment

END MODULE molecularpropertiesmod





MODULE coulombmod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

   REAL(SP), PRIVATE, PARAMETER :: TOL_HIGH = 5.0_SP
   REAL(SP), PRIVATE, PARAMETER :: TOL_LOW = 1.E-2_SP 
   REAL(SP), PRIVATE, PARAMETER :: TWO_D_SQRTPI = 2.0_SP/SQRT_PI

INTERFACE pointS_pointS
 Module Procedure pointS_pointS_pot,pointS_pointS_pgrad
END INTERFACE

INTERFACE pointS_gausS
 Module Procedure pointS_gausS_pot,pointS_gausS_pgrad
END INTERFACE

INTERFACE gausS_pointS
 Module Procedure gausS_pointS_pot,gausS_pointS_pgrad
END INTERFACE

INTERFACE gausS_gausS
 Module Procedure gausS_gausS_pot,gausS_gausS_pgrad
END INTERFACE

CONTAINS
subroutine pointS_gausS_pgrad(Crd1,q1,pot1,grad1,Crd2,q2,zeta2,pot2,grad2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:),zeta2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),grad1(:,:)
  REAL(SP), INTENT(INOUT) :: pot2(:),grad2(:,:)
  LOGICAL(LGD),intent(in) :: SameSet
  CALL gausS_pointS_pgrad(Crd2,q2,zeta2,pot2,grad2,Crd1,q1,pot1,grad1,sameset)
end subroutine pointS_gausS_pgrad

subroutine gausS_gausS_pot(Crd1,q1,zeta1,pot1,Crd2,q2,zeta2,pot2,sameset)
  USE DataTypes
  Use ConstantsMod
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:),zeta1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:),zeta2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),pot2(:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,j0,n1,n2
  REAL(SP) :: tphi,phi1,zeti,zetj,zetp,r
  REAL(SP) :: qi,xi(1:3),dx(1:3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  j0=1
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     zeti=zeta1(i)
     IF (SameSet) THEN
        j0=i+1
        zetp=(zeti/SQRT2)
        CALL GaussianInteraction(zetp,0.0_SP,tphi)
        phi1=qi*tphi
     END IF
     DO j=j0,n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j); zetj=zeta2(j)
        zetp=zeti*zetj/SQRT(zeti*zeti+zetj*zetj)
        r=SQRT(DOT_PRODUCT(dx,dx))
        Call GaussianInteraction(zetp,r,tphi)
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + q1(i)*tphi
     END DO
     pot1(i) = pot1(i) + phi1
  END DO
end subroutine gausS_gausS_pot

subroutine pointS_pointS_pot(Crd1,q1,pot1,Crd2,q2,pot2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),pot2(:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,js(2),n1,n2,jidx
  REAL(SP) :: tphi,phi1
  REAL(SP) :: qi,xi(1:3),dx(1:3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  js(1)=1
  IF(sameset) THEN
     jidx = 2
  ELSE
     jidx = 1
  ENDIF
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     js(2)=i+1
     DO j=js(jidx),n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j)
        tphi = SQRT(1.0_SP/DOT_PRODUCT(dx,dx))
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + qi*tphi
     END DO
     pot1(i) = pot1(i) + phi1
  END DO
end subroutine pointS_pointS_pot

  Subroutine dGaussianInteraction(zeta,r,phi,dphi)
    USE DataTypes
    USE MathMod, ONLY : ERF
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: zeta,r
    REAL(SP), INTENT(OUT) :: PHI,dphi
    REAL(SP) :: x, x2
    x = zeta*r
    ! all this junk is to prevent divide by zero, numerical stability issuse
    !  and calling erf
    IF ( x > TOL_HIGH ) THEN
       PHI  =  zeta/x
       ! remember... dphi below is 1/r * d/dr phi
       dphi = -zeta*zeta*PHI/(x*x)
    ELSE IF ( x > TOL_LOW ) THEN
       PHI = zeta*ERF(x)/x
       x2 = x*x
       ! remember... dphi below is 1/r * d/dr phi
       dphi = zeta*zeta*(TWO_D_SQRTPI * zeta * EXP(-x2) - PHI)/x2 
    ELSE
       ! for small x, we use a 3 term taylor expansion of erf(x)/x around 0
       x2 = x*x
       PHI = zeta*TWO_D_SQRTPI * ( 1.0_SP -x2/3.0_SP + x2*x2/10.0_SP )
       ! remember... dphi below is 1/r * d/dr phi 
       dphi = zeta*zeta*zeta*TWO_D_SQRTPI * ( -2.0_SP / 3.0_SP + 4.0_SP * x2 / 10.0_SP )
    END IF
  end Subroutine dGaussianInteraction

subroutine gausS_pointS_pgrad(Crd1,q1,zeta1,pot1,grad1,Crd2,q2,pot2,grad2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:),zeta1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),grad1(:,:)
  REAL(SP), INTENT(INOUT) :: pot2(:),grad2(:,:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,j0,n1,n2
  REAL(SP) :: tphi,phi1,tdphi,gx1(3),zetp,r
  REAL(SP) :: qi, xi(3),dx(3),tgrad(3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  j0=1
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP; gx1(1:3)=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     zetp=zeta1(i)
     IF (SameSet) THEN
        j0=i+1
        CALL GaussianInteraction(zetp,0.0_SP,tphi)
        phi1=qi*tphi
     END IF
     DO j=j0,n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j)
        r=SQRT(DOT_PRODUCT(dx,dx))
        CALL dGaussianInteraction(zetp,r,tphi,tdphi)
        tgrad(1:3) = tdphi * dx(1:3)
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + qi*tphi
        gx1(1:3) = gx1(1:3) + q2(j)*tgrad(1:3)
        grad2(1:3,j) = grad2(1:3,j) - qi*tgrad(1:3)
     END DO
     pot1(i) = pot1(i) + phi1
     grad1(1:3,i) = grad1(1:3,i) + gx1(1:3)
  END DO
end subroutine gausS_pointS_pgrad

subroutine pointS_pointS_pgrad(Crd1,q1,pot1,grad1,Crd2,q2,pot2,grad2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),grad1(:,:)
  REAL(SP), INTENT(INOUT) :: pot2(:),grad2(:,:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,js(2),jidx,n1,n2
  REAL(SP) :: tphi,phi1,gx1(3),rm2
  REAL(SP) :: qi, xi(3),dx(3),tgrad(3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  js(1)=1
  IF(sameset) THEN
     jidx = 2
  ELSE
     jidx = 1
  ENDIF
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP; gx1(1:3)=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     js(2)=i+1
     DO j=js(jidx),n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j)
        rm2 = 1.0_SP/DOT_PRODUCT(dx,dx)
        tphi = SQRT(rm2)
        tgrad(1:3) = -rm2 * tphi * dx(1:3)
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + qi*tphi
        gx1(1:3) = gx1(1:3) + q2(j) * tgrad(1:3)
        grad2(1:3,j) = grad2(1:3,j) - qi * tgrad(1:3)
     END DO
     pot1(i) = pot1(i) + phi1
     grad1(1:3,i) = grad1(1:3,i) + gx1(1:3)
  END DO
end subroutine pointS_pointS_pgrad

subroutine pointS_gausS_pot(Crd1,q1,pot1,Crd2,q2,zeta2,pot2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:),zeta2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),pot2(:)
  LOGICAL(LGD),intent(in) :: SameSet
  CALL gausS_pointS_pot(Crd2,q2,zeta2,pot2,Crd1,q1,pot1,sameset)
end subroutine pointS_gausS_pot

  Subroutine GaussianInteraction(zeta,r,phi)
    USE DataTypes
    USE MathMod, ONLY : ERF
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: zeta,r
    REAL(SP), INTENT(OUT) :: PHI
    REAL(SP) ::x, x2
    x = zeta*r
    ! all this junk is to prevent divide by zero, numerical stability issuse
    !  and calling erf
    IF ( x > TOL_HIGH ) THEN
       PHI = zeta/x
    ELSE IF ( x > TOL_LOW ) THEN
       PHI = zeta*ERF(x)/x
    ELSE
       ! for small x, we use a 3 term taylor expansion of erf(x)/x around 0
       x2 = x*x
       PHI = zeta * TWO_D_SQRTPI * ( 1.0_SP -x2/3.0_SP + x2*x2/10.0_SP )
    END IF
  END Subroutine GaussianInteraction

subroutine gausS_gausS_pgrad(Crd1,q1,zeta1,pot1,grad1,Crd2,q2,zeta2,pot2,grad2,sameset)
  USE DataTypes
  Use ConstantsMod
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:),zeta1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:),zeta2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),grad1(:,:)
  REAL(SP), INTENT(INOUT) :: pot2(:),grad2(:,:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,j0,n1,n2
  REAL(SP) :: tphi,phi1,tdphi,gx1(3),zetp,zeti,zetj,r
  REAL(SP) :: qi, xi(3),dx(3),tgrad(3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  j0=1
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP; gx1(1:3)=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     zeti=zeta1(i)
     IF (SameSet) THEN
        j0=i+1
        zetp=(zeti/SQRT2)
        CALL GaussianInteraction(zetp,0.0_SP,tphi)
        phi1=qi*tphi
     END IF
     DO j=j0,n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j); zetj=zeta2(j)
        zetp=zeti*zetj/SQRT(zeti*zeti+zetj*zetj)
        r=SQRT(DOT_PRODUCT(dx,dx))
        CALL dGaussianInteraction(zetp,r,tphi,tdphi)
        tgrad(1:3) = tdphi * dx(1:3)
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + qi*tphi
        gx1(1:3) = gx1(1:3) + q2(j)*tgrad(1:3)
        grad2(1:3,j) = grad2(1:3,j) - qi*tgrad(1:3)
     END DO
     pot1(i) = pot1(i) + phi1
     grad1(1:3,i) = grad1(1:3,i) + gx1(1:3)
  END DO
end subroutine gausS_gausS_pgrad

subroutine gausS_pointS_pot(Crd1,q1,zeta1,pot1,Crd2,q2,pot2,sameset)
  USE DataTypes
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: Crd1(:,:),q1(:),zeta1(:)
  REAL(SP), INTENT(IN) :: Crd2(:,:),q2(:)
  REAL(SP), INTENT(INOUT) :: pot1(:),pot2(:)
  LOGICAL(LGD),intent(in) :: SameSet
  INTEGER(I4B) :: i,j,j0,n1,n2
  REAL(SP) :: tphi,phi1,zetp,r
  REAL(SP) :: qi,xi(1:3),dx(1:3)
  n1 = SIZE(Crd1,2); n2 = SIZE(Crd2,2)
  j0=1
  DO i=1,n1
     IF(q1(i)==0.0_SP) CYCLE
     phi1=0.0_SP
     qi=q1(i); xi(1:3)=Crd1(1:3,i)
     zetp=zeta1(i)
     IF (SameSet) THEN
        j0=i+1
        Call GaussianInteraction(zetp,0.0_SP,tphi)
        phi1=qi*tphi
     END IF
     DO j=j0,n2
        dx(1:3)=xi(1:3)-Crd2(1:3,j)
        r=SQRT(DOT_PRODUCT(dx,dx))
        Call GaussianInteraction(zetp,r,tphi)
        phi1 = phi1 + q2(j)*tphi
        pot2(j) = pot2(j) + qi*tphi
     END DO
     pot1(i) = pot1(i) + phi1
  END DO
end subroutine gausS_pointS_pot

END MODULE coulombmod




MODULE SphericalHarmonicsMod
  USE DataTypes
  IMPLICIT NONE
  
  REAL(SP), pointer :: RT(:)=>NULL()

  PRIVATE 

  PUBLIC :: gen_roots, get_sh_root,&
       RegSolidHarmonic, RegSolidHarmonicDeriv,&
       IrregSolidHarmonic


contains
  Subroutine RegSolidHarmonic(Crd, lmax, Norm, Glm, Glmd )
    ! regular solid harmonics are returned:
    ! Glm = r^l * H_lm
    ! Norm == .true.
    !   H_lm = Ylm(theta,phi) s.t.:
    !      \int H_{lm}* H_{l'm'} \dddr = \delta{ll'}\delta{mm'}
    ! Norm == false
    !   H_lm = Clm{theta,phi) = SQRT(4*PI/(2*l+1)) * Ylm(theta,phi) s.t.:
    !      H_l0(0,0) = 1
    ! Return Values in order of
    !  Index     1   2   3    4    5   6    7    8    9   10  11
    !  G_{lm}   00  10  11c  11s  20  21c  21s  22c  21s  30  31c ...
    !    Where c and s stands for the cos and sin components of the complex harmonic.
    !    (The returned real values are linear combinations of the Complex harmonics)
    USE Datatypes
    USE ConstantsMod
    IMPLICIT NONE
    
    REAL(SP), INTENT(IN) :: Crd(:,:)              ! 1:3, 1:nq
    INTEGER(I4B), INTENT(IN) :: lmax
    LOGICAL(LGD), INTENT(IN) :: Norm
    REAL(SP), INTENT(OUT) :: Glm(:,:)             ! 1:ldim, 1:nq
    REAL(SP), INTENT(OUT),OPTIONAL :: Glmd(:,:,:) ! 1:3, 1:ldim, 1:nq
    
    REAL(SP) :: R2, RR, RF0, RF1
    REAL(SP) :: norm_factor(SIZE(Glm,1))
    INTEGER(I4B) :: i

    call gen_roots(lmax)
    
    do i=1,SIZE(Crd,2)
       R2=crd(1,i)**2+crd(2,i)**2+crd(3,i)**2
       !  Regular
       RR = R2
       RF0 = 1.0_SP
       RF1 = 1.0_SP
       Call GenHarmonic(Glm(:,i), Crd(:,i), lmax, R2, RR, RF0, RF1)
    enddo
    
    IF(PRESENT(Glmd))THEN
       CALL RegSolidHarmonicDeriv(SIZE(Crd,2), lmax,Glm,Glmd)
       
       IF(norm)THEN
          norm_factor = SHNorm(lmax)
          do i=1,SIZE(Crd,2)
             Glm(:,i)  = Glm(:,i) * norm_factor(:)
             Glmd(1:3,:,i) = Glmd(1:3,:,i)*spread(norm_factor(:),1,3)
          enddo
       ENDIF
    ELSEIF(norm)THEN
       norm_factor = SHNorm(lmax)
       do i=1,SIZE(Crd,2)
          Glm(:,i)=Glm(:,i) * norm_factor(:)
       enddo
    ENDIF
    
  end Subroutine RegSolidHarmonic



  Subroutine IrregSolidHarmonic(Crd, lmax, Norm,Glm, Glmd)
    ! irregular solid harmonics are returned:
    ! Glm = r^(-l+1) * H_lm
    ! Norm == .true.
    !   H_lm = Ylm(theta,phi) s.t.:
    !      \int H_{lm}* H_{l'm'} \dddr = \delta{ll'}\delta{mm'}
    ! Norm == false
    !   H_lm = Clm{theta,phi) = SQRT(4*PI/(2*l+1)) * Ylm(theta,phi) s.t.:
    !      H_l0(0,0) = 1
    ! Return Values in order of
    !  Index     1   2   3    4    5   6    7    8    9   10  11
    !  G_{lm}   00  10  11c  11s  20  21c  21s  22c  21s  30  31c ...
    !    Where c and s stands for the cos and sin components of the complex harmonic.
    !    (The returned real values are linear combinations of the Complex harmonics)
    USE Datatypes
    IMPLICIT NONE
 
    REAL(SP), INTENT(IN) :: Crd(:,:)              ! 1:3, 1:nq
    INTEGER(I4B), INTENT(IN) :: lmax
    LOGICAL(LGD), INTENT(IN) :: Norm
    REAL(SP), INTENT(OUT) :: Glm(:,:)             ! 1:ldim, 1:nq
    REAL(SP), INTENT(OUT),OPTIONAL :: Glmd(:,:,:) ! 1:3, 1:ldim, 1:nq
  
    REAL(SP) :: R2, RR, RF0, RF1
    REAL(SP) :: norm_factor(SIZE(Glm,1))
    REAL(SP), allocatable :: Glmtmp(:,:)
    INTEGER(I4B) :: i,ldim
    IF(PRESENT(Glmd))THEN
       ! Derivatives for Irregular harmonics require recursions to l+1
       call gen_roots(lmax+1)
       
       ldim = (lmax+2)**2
       allocate(Glmtmp(ldim,SIZE(Glm,2)))
       Glmtmp = 0.0_SP
       do i=1,SIZE(Crd,2)
          R2=crd(1,i)**2+crd(2,i)**2+crd(3,i)**2
          !  Irregular
          RR = 1.0_SP/R2
          RF0 = SQRT(RR)
          RF1 = RR
          Call GenHarmonic(Glmtmp(:,i), Crd(:,i), lmax+1, R2, RR, RF0, RF1)
       enddo
       
       CALL IrregSolidHarmonicDeriv(SIZE(Crd,2), lmax, Glmtmp, Glmd)
       
       ldim = (lmax+1)**2
       Glm(1:ldim,:) = Glmtmp(1:ldim,:)
       deallocate(Glmtmp)
       IF(norm)THEN
          norm_factor = SHNorm(lmax)
          do i=1,SIZE(Crd,2)
             Glm(:,i)  = Glm(:,i) * norm_factor(:)
             Glmd(1:3,:,i) = Glmd(1:3,:,i)*spread(norm_factor(:),1,3)
          enddo
       ENDIF
    ELSE
       call gen_roots(lmax)
       
       do i=1,SIZE(Crd,2)
          R2=crd(1,i)**2+crd(2,i)**2+crd(3,i)**2
          !  Irregular
          RR = 1.0_SP/R2
          RF0 = SQRT(RR)
          RF1 = RR
          Call GenHarmonic(Glm(:,i), Crd(:,i), lmax, R2, RR, RF0, RF1)
       enddo
       
       IF(norm)THEN
          norm_factor = SHNorm(lmax)
          do i=1,SIZE(Crd,2)
             Glm(:,i)=Glm(:,i) * norm_factor(:)
          enddo
       ENDIF
    ENDIF
    
  end Subroutine IrregSolidHarmonic

  Subroutine GenHarmonic(Glm, Crd, lmax, R2, RR, RF0, RF1)
    ! Generic routine for generating Regular and Irregular Solid Harmonics
    ! in addition to the Real spherical harmonics.
    USE Datatypes
    USE ConstantsMod
    IMPLICIT NONE
    
    REAL(SP), INTENT(OUT) :: Glm(:)
    REAL(SP), INTENT(IN) :: Crd(:)
    INTEGER(I4B), INTENT(IN) :: lmax
    REAL(SP), INTENT(IN) :: R2, RR, RF0, RF1
    
    
    INTEGER(I4B) :: ILC ! current L index vals
    INTEGER(I4B) :: V2LP1,VLp1,V2m ! variables: 2*L+1, L+1, 2*m
    INTEGER(I4B) :: L, M
    REAL(SP) :: Plm,plm0,prt,crt, tmp  ! previous lm
    REAL(SP) :: RFX,RFY,RFZ,S
    
    RFX = crd(1) * RF1
    RFY = crd(2) * RF1
    RFZ = crd(3) * RF1
    
    !  Obtain Glm(l+1,0) from Glm(l,0)*Glm(1,0) and Glm(l-1,0)
    plm=RF0; plm0=0.0_SP; ILC = 1
    LL: DO L=0, lmax-1
       Vlp1= L+1
       V2LP1=L+Vlp1
       Glm(ILC) = plm
       tmp=(V2LP1*RFZ*plm-L*RR*Plm0) / Vlp1
       plm0=plm; plm=tmp
       ILC= ILC + V2LP1
    enddo LL
    Glm(ILC) = plm
    ! For each value of m, accend through l space
    do m=1, lmax
       V2m = 2*m
       plm=RF0; plm0=0.0_SP; ILC = m*m+V2m; prt = RT(V2m)
       DO L=m,Lmax-1
          Vlp1 = L+1
          V2LP1=L+Vlp1
          !  Obtain Glm(l+1,m) from Glm(l,m)*Glm(1,0) and Glm(l-1,m)
          Glm(ILC) = plm * RFX
          Glm(ILC+1) = plm * RFY
          crt = RT(Vlp1+M)*RT(Vlp1-M)
          tmp=(V2LP1*RFZ*plm-prt*RR*Plm0)/ crt
          plm0=plm; plm=tmp; prt=crt
          ILC= ILC + V2LP1
       ENDDO
       Glm(ILC) = plm * RFX
       Glm(ILC+1) = plm * RFY
       
       S=RF1*RT(V2m+1)/RT(V2m+2)
       tmp = RFX
       RFX = S*(tmp*crd(1)-RFY*crd(2))
       RFY = S*(tmp*crd(2)+RFY*crd(1))
       
    ENDDO
    
  end Subroutine GenHarmonic
  subroutine RegSolidHarmonicDeriv(n,lmax,Nlm,Nlmd)
    ! X,Y and Z derivatives of the Regular Solid Harmonics
    USE Datatypes
    USE ConstantsMod
    Use MathMod
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) ::lmax,n
    REAL(SP),INTENT(INOUT) :: Nlm(:,:)    ! 1:ldim, 1:nq
    REAL(SP),INTENT(OUT) :: Nlmd(:,:,:) ! 1:3, 1:ldim,1:nq
    INTEGER(I4B) :: i,m1,l1,lind,lindp1,V2m,lp1,lp2,idx
    REAL(SP) :: a,b,factor
    
    IF(lmax==0)RETURN
    
    DO i=1,n
       a = Nlm(1,i)
       Nlmd(3,2,i)=Nlm(1,i)
       Nlmd(1,3,i)=a
       Nlmd(2,4,i)=a
       
       DO l1=1,Lmax-1
          lind=L1*L1 + 1
          lindp1=lind+2*l1+1
          lp1 = l1 + 1
          lp2 = l1 + 2
          
          ! terms for L'=L
          Nlmd(3,lindp1,i)=Nlm(lind,i) *(lp1)
          
          a = Nlm(lind,i) /SQRT2 *RT(lp1)*RT(lp2)
          Nlmd(1,lindp1+1,i)=a
          Nlmd(2,lindp1+2,i)=a
          ! terms for L'=L+1,M=+-1
          Nlmd(3,lindp1+1,i)=Nlm(lind+1,i) *RT(lp2)*RT(l1)
          Nlmd(3,lindp1+2,i)=Nlm(lind+2,i) *RT(lp2)*RT(l1)
          
          ! Terms for L'=L+1,M=+-2 from l,1
          a=Nlm(lind+1,i)
          b=Nlm(lind+2,i)
          
          ! terms for L'=L+1,M=0
          factor =  -RT(l1)*RT(lp1)/SQRT2
          Nlmd(1,lindp1,i) = a * factor
          Nlmd(2,lindp1,i) = b * factor
          
          factor =  0.5_SP * RT(l1+2)*RT(l1+3)
          a=a * factor
          b=b * factor
          idx=lindp1+3
          Nlmd(1,idx,i)=Nlmd(1,idx,i)+a
          Nlmd(2,idx,i)=Nlmd(2,idx,i)-b
          idx=lindp1+4
          Nlmd(1,idx,i)=Nlmd(1,idx,i)+b
          Nlmd(2,idx,i)=Nlmd(2,idx,i)+a
          
          
       END DO
       DO m1=2,Lmax-1
          V2m= 2*m1
          DO l1=m1,Lmax-1
             lind=L1*L1 + 1
             lindp1=lind+2*l1+1
             lp1 = l1 + 1
             lp2 = l1 + 2
             
             factor = RT(lp1+m1)*RT(lp1-m1)
             Nlmd(3,lindp1+V2m-1,i)=factor*Nlm(lind+V2m-1,i)
             Nlmd(3,lindp1+V2m,i)=factor*Nlm(lind+V2m,i)
             
             factor= 0.5_SP
             a=Nlm(lind+V2m-1,i)    *factor *RT(lp2-m1)*RT(lp1-m1)
             b=Nlm(lind+V2m,i)      *factor *RT(lp2-m1)*RT(lp1-m1)
             
             idx=lindp1+V2m-3
             Nlmd(1,idx,i)=Nlmd(1,idx,i)-a
             Nlmd(2,idx,i)=Nlmd(2,idx,i)-b
             idx=idx + 1
             Nlmd(2,idx,i)=Nlmd(2,idx,i)+a
             Nlmd(1,idx,i)=Nlmd(1,idx,i)-b
             
             factor = 0.5_SP*RT(l1+m1+1)*RT(l1+m1+2)
             a=Nlm(lind+V2m-1,i)*factor
             b=Nlm(lind+V2m,i) *factor
             idx=lindp1+V2m+1
             Nlmd(1,idx,i)=Nlmd(1,idx,i)+a
             Nlmd(2,idx,i)=Nlmd(2,idx,i)-b
             idx=idx + 1
             Nlmd(2,idx,i)=Nlmd(2,idx,i)+a
             Nlmd(1,idx,i)=Nlmd(1,idx,i)+b
          END DO
       END DO
    end DO
    
  end Subroutine RegSolidHarmonicDeriv

  subroutine IrregSolidHarmonicDeriv(n,lmax,Mlm,Mlmd)
    ! X,Y and Z derivatives of the Irregular Solid Harmonics
    USE Datatypes
    USE ConstantsMod
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) ::lmax,n
    REAL(SP),INTENT(INOUT) :: Mlm(:,:)    ! 1:ldim, 1:nq
    REAL(SP),INTENT(OUT) :: Mlmd(:,:,:) ! 1:3, 1:ldim,1:nq
    INTEGER(I4B) :: i,m1,l1,lind,lindp1,V2m,lp1,lp2
    REAL(SP) :: a,b,factor,factor2

    DO i=1,n
       Mlmd(1,1,i)=-Mlm(3,i)  
       Mlmd(2,1,i)=-Mlm(4,i)  
       Mlmd(3,1,i)=-Mlm(2,i) 

       DO l1=1,Lmax
          lind=L1*L1 + 1
          lindp1=lind+2*l1+1
          lp1 = l1 + 1
          lp2 = l1 + 2

          ! terms for L'=l+1 m=0 
          Mlmd(3,lind,i)=-(lp1) * Mlm(lindp1,i) 
          factor = -RT(lp1)*RT(lp2)/SQRT2
          Mlmd(1,lind,i)= factor * Mlm(lindp1+1,i)
          Mlmd(2,lind,i)= factor * Mlm(lindp1+2,i)

          ! Terms for L'=L+1,M=+-1 
          factor = -RT(lp2)*RT(l1)
          Mlmd(3,lind+1,i)= factor * Mlm(lindp1+1,i) 
          Mlmd(3,lind+2,i)= factor * Mlm(lindp1+2,i) 

          factor=0.5_SP*RT(lp2)*RT(l1+3)
          a = factor*Mlm(lindp1+3,i)
          b = 0.5_SP*RT(lp1)*RT(l1)*Mlm(lindp1,i)*SQRT2
          Mlmd(1,lind+1,i)=(-a+b)
          Mlmd(2,lind+2,i)=( a+b)

          a = -factor*Mlm(lindp1+4,i)
          Mlmd(1,lind+2,i)=a
          Mlmd(2,lind+1,i)=a

       END DO
       DO l1=2,Lmax
          lind=L1*L1 + 1
          lindp1=lind+2*l1+1
          lp1 = l1 + 1
          lp2 = l1 + 2
          DO m1=2,L1
             V2m= 2*m1
             factor = - RT(lp1+m1)*RT(lp1-m1)
             Mlmd(3,lind+V2m-1,i)= factor * Mlm(lindp1+V2m-1,i) 
             Mlmd(3,lind+V2m,i)  = factor * Mlm(lindp1+V2m,i) 

             factor  = 0.5_SP*RT(lp1+m1)*RT(lp2+m1)
             factor2 = 0.5_SP*RT(lp2-m1)*RT(lp1-m1)
             a = factor  * Mlm(lindp1+V2m+1,i)
             b = factor2 * Mlm(lindp1+V2m-3,i)
             Mlmd(1,lind+V2m-1,i) = -a + b 
             Mlmd(2,lind+V2m,i)   =  a + b

             a = -factor * Mlm(lindp1+V2m+2,i)
             b = factor2 * Mlm(lindp1+V2m-2,i)
             Mlmd(1,lind+V2m,i)   = a + b
             Mlmd(2,lind+V2m-1,i) = a - b
          END DO
       END DO
    end DO

  end Subroutine IrregSolidHarmonicDeriv


  subroutine gen_roots(lmax)
    USE DataTypes
    USE ConstantsMod
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: lmax
    INTEGER(I4B) :: i,istart,nval
    REAL(SP),pointer :: rt_new(:)=>NULL()
    
    ! space for entries from 0:2*lmax+2
    nval = max(2*(lmax+2),5)
    IF(associated(RT))THEN
       IF(size(RT) >= nval) RETURN
       istart = size(rt)
       allocate(rt_new(0:nval-1))
       rt_new(0:size(rt)-1) = rt(0:)
       deallocate(rt)
       rt => rt_new
       rt_new => NULL()
    ELSE
       allocate(rt(0:nval-1))
       RT(0)= 0.0_SP
       RT(1)= 1.0_SP
       RT(2)= SQRT2
       RT(3)= SQRT3
       RT(4)= 2.0_SP
       istart=5
    ENDIF
    DO i=istart,nval-1
       RT(i)=SQRT(REAL(i,SP))
    ENDDO
    
  end subroutine gen_roots
  
  function get_sh_root(i) result(root)
    USE DataTypes
    IMPLICIT NONE
      INTEGER(I4B),INTENT(IN) :: i
      REAL(SP) :: root
      root=RT(i)
  end function get_sh_root

  function SHNorm(lmax) result(nf)
    USE Datatypes
    USE ConstantsMod
    IMPLICIT NONE
    INTEGER(I4B),INTENT(IN):: lmax
    REAL(SP) :: nf(lmax**2+2*lmax+1),tmp
    REAL(SP),parameter :: SQRT_4PI_INV = 1.0_SP/(2.0_SP*SQRT_PI)
    INTEGER(I4B) :: L
    
    DO L=0, lmax
       tmp = RT(2*L+1)*SQRT_4PI_INV
       nf(l**2+1:l**2+2*l+1) = tmp
    enddo
    
  end function SHNorm


END MODULE SphericalHarmonicsMod


MODULE multipolemod
   USE DATATYPES
   USE CONSTANTSMOD
   IMPLICIT NONE

CONTAINS


  subroutine NearField(Crd,orig,lmax,Nlm, Nlmd)
    USE Datatypes
    Use SphericalHarmonicsMod, only : RegSolidHarmonic
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: lmax
    REAL(SP), INTENT(IN) :: Crd(:,:),orig(:)
    REAL(SP), INTENT(OUT) :: Nlm(:,:)
    REAL(SP), INTENT(OUT), optional :: Nlmd(:,:,:)
    
    CALL RegSolidHarmonic(crd(1:3,:)-spread(orig(1:3), DIM=2,NCOPIES=SIZE(Crd,2)),lmax,.FALSE.,Nlm, Nlmd)
    
  end subroutine NearField

  subroutine FarField(Crd,orig,lmax,Mlm, Mlmd)
    USE Datatypes
    Use SphericalHarmonicsMod, only : IrregSolidHarmonic
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: lmax
    REAL(SP), INTENT(IN) :: Crd(:,:),orig(:)
    REAL(SP), INTENT(OUT) :: Mlm(:,:)
    REAL(SP), INTENT(OUT), optional :: Mlmd(:,:,:)
    CALL IrregSolidHarmonic(crd(1:3,:)-spread(orig(1:3), DIM=2,NCOPIES=SIZE(Crd,2)),lmax,.FALSE.,Mlm, Mlmd)
  end subroutine FarField

  Subroutine qTimesNearField(q, Crd, orig, lmax, FMM,QNlm, QNlmd)
    ! Code adapted from spherical harmonic routines for
    ! L1 normalized regular solid harmonics.
    USE Datatypes
    USE SphericalHarmonicsMod , only : rt=>get_sh_root, gen_roots,RegSolidHarmonicDeriv
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: q(:),Crd(:,:),orig(3)              ! 1:3, 1:nq
    INTEGER(I4B), INTENT(IN) :: lmax
    LOGICAL(LGD),INTENT(IN) :: FMM
    REAL(SP), INTENT(OUT) :: QNlm(:)             ! 1:ldim
    REAL(SP), INTENT(OUT),OPTIONAL :: QNlmd(:,:) ! 1:3, 1:ldim
    INTEGER(I4B) :: ILC ! current L index vals
    INTEGER(I4B) :: V2LP1,VLp1,V2m ! variables: 2*L+1, L+1, 2*m
    INTEGER(I4B) :: i,L, M,ldim,lderiv
    REAL(SP) :: Plm,plm0,prt,crt, tmp
    REAL(SP) ::S
    REAL(SP) :: R2,RFX,RFY,RFZ,RF0
    REAL(SP) :: xyz(1:3)
    REAL(SP) :: QLM0(1:(lmax+2)**2,1)
    REAL(SP) :: QLM0d(1:3,1:(lmax+2)**2,1)    
    IF(FMM)THEN
       ldim=(lmax+2)**2
       CALL gen_roots(lmax+1)
       lderiv=lmax+1
    ELSE
       ldim=(lmax+1)**2
       CALL gen_roots(lmax)
       lderiv=lmax
    ENDIF
    QLM0 = 0.0_SP
    DO i=1,SIZE(crd,2)
       !  Some initialization ...
       xyz(1:3)=Crd(1:3,i)-orig(1:3)
       r2=DOT_PRODUCT(xyz,xyz)
       RF0 = q(i)
       RFX = xyz(1)
       RFY = xyz(2) 
       RFZ = xyz(3) 
       !  Obtain Qlm(l+1,0) from Qlm(l,0)*Qlm(1,0) and Qlm(l-1,0)
       plm=RF0; plm0=0.0_SP; ILC = 1
       LL: DO L=0, lmax-1
          Vlp1= L+1
          V2LP1=L+Vlp1
          Qlm0(ILC,1) = Qlm0(ILC,1)+plm 
          tmp=(V2LP1*RFZ*plm-L*R2*Plm0) / Vlp1
          plm0=plm; plm=tmp
          ILC= ILC + V2LP1
       enddo LL
       Qlm0(ILC,1) = Qlm0(ILC,1)+plm
       ! For each value of m, accend through l space
       do m=1, lmax
          V2m = 2*m
          plm=RF0; plm0=0.0_SP; ILC = m*m+V2m; prt = RT(V2m)
          DO L=m,Lmax-1
             Vlp1 = L+1
             V2LP1=L+Vlp1
             !  Obtain Qlm(l+1,m) from Qlm(l,m)*Qlm(1,0) and Qlm(l-1,m)
             Qlm0(ILC,1) = Qlm0(ILC,1)+plm * RFX
             Qlm0(ILC+1,1) = Qlm0(ILC+1,1)+plm * RFY
             crt = RT(Vlp1+M)*RT(Vlp1-M)
             tmp=(V2LP1*RFZ*plm-prt*R2*Plm0)/ crt
             plm0=plm; plm=tmp; prt=crt
             ILC= ILC + V2LP1
          ENDDO
          Qlm0(ILC,1) = Qlm0(ILC,1)+plm * RFX
          Qlm0(ILC+1,1) = Qlm0(ILC+1,1)+plm * RFY
          S=RT(V2m+1)/RT(V2m+2) 
          tmp = RFX
          RFX = S*(tmp*xyz(1)-RFY*xyz(2))
          RFY = S*(tmp*xyz(2)+RFY*xyz(1))
       ENDDO
    enddo
    QNlm(1:ldim) = QLM0(1:ldim,1)
    IF(PRESENT(QNlmd))THEN
       QLM0d= 0.0_SP
       CALL RegSolidHarmonicDeriv(1,lderiv,Qlm0,Qlm0d)
       QNlmd(1:3,1:ldim) = QLM0d(1:3,1:ldim,1)
    ENDIF
  end Subroutine qTimesNearField

  Subroutine qTimesFarField(q, Crd, orig, lmax, QMlm)
    USE Datatypes
    USE SphericalHarmonicsMod , only : rt=>get_sh_root, gen_roots
    IMPLICIT NONE
    ! Code adapted from spherical harmonic routines for
    ! L1 normalized irregular solid harmonics.
    REAL(SP), INTENT(IN) :: q(:),Crd(:,:),orig(3)              ! 1:3, 1:nq
    INTEGER(I4B), INTENT(IN) :: lmax
    REAL(SP), INTENT(OUT) :: QMlm(:)             ! 1:ldim
!    REAL(SP), INTENT(OUT),OPTIONAL :: QNlmd(:,:) ! 1:3, 1:ldim
    INTEGER(I4B) :: ILC ! current L index vals
    INTEGER(I4B) :: V2LP1,VLp1,V2m ! variables: 2*L+1, L+1, 2*m
    INTEGER(I4B) :: i,L, M, ldim
    REAL(SP) :: Plm,plm0,prt,crt, tmp
    REAL(SP) ::S
    REAL(SP) :: RR,R2,RFX,RFY,RFZ,RF0
    REAL(SP) :: xyz(1:3)
    REAL(SP) :: QLM0(1:(lmax+1)**2)
    CALL gen_roots(lmax)
    QLM0=0.0_SP
    ldim=(lmax+1)**2
    do i=1,SIZE(Crd,2)
       xyz(1:3)=Crd(1:3,i)-orig(1:3)
       r2=DOT_PRODUCT(xyz,xyz)
       !  Irregular
       RR = 1.0_SP/R2
       RF0 = q(i)*SQRT(RR)
       RFX = xyz(1) * RR
       RFY = xyz(2) * RR
       RFZ = xyz(3) * RR
       !  Obtain QMlm(l+1,0) from QMlm(l,0)*QMlm(1,0) and QMlm(l-1,0)
       plm=RF0; plm0=0.0_SP; ILC = 1
       LL: DO L=0, lmax-1
          Vlp1= L+1
          V2LP1=L+Vlp1
          Qlm0(ILC) = plm 
          tmp=(V2LP1*RFZ*plm-L*RR*Plm0) / Vlp1
          plm0=plm; plm=tmp
          ILC= ILC + V2LP1
       enddo LL
       Qlm0(ILC) = plm
       ! For each value of m, accend through l space
       do m=1, lmax
          V2m = 2*m
          plm=RF0; plm0=0.0_SP; ILC = m*m+V2m; prt = RT(V2m)
          DO L=m,Lmax-1
             Vlp1 = L+1
             V2LP1=L+Vlp1
             !  Obtain QMlm(l+1,m) from QMlm(l,m)*QMlm(1,0) and QMlm(l-1,m)
             Qlm0(ILC) = plm * RFX
             Qlm0(ILC+1) = plm * RFY
             crt = RT(Vlp1+M)*RT(Vlp1-M)
             tmp=(V2LP1*RFZ*plm-prt*RR*Plm0)/ crt
             plm0=plm; plm=tmp; prt=crt
             ILC= ILC + V2LP1
          ENDDO
          Qlm0(ILC) = plm * RFX
          Qlm0(ILC+1) = plm * RFY
          S=RR*RT(V2m+1)/RT(V2m+2) 
          tmp = RFX
          RFX = S*(tmp*xyz(1)-RFY*xyz(2))
          RFY = S*(tmp*xyz(2)+RFY*xyz(1))
       ENDDO
    enddo
    QMlm(1:ldim) = QLM0(1:ldim)
  end Subroutine qTimesFarField

END MODULE multipolemod


