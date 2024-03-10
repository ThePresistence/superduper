MODULE SurfaceMod
!!!****h* SurfaceMod
!!!
!!! NAME
!!!     SurfaceMod -- module for generating molecular surfaces
!!!
!!! COPYRIGHT
!!!     Prof. York's Group
!!!     Department of Chemistry
!!!     University of Minnesota
!!! AUTHOR
!!!     Darrin York
!!!     Kevin Range
!!!     Brent Gregersen
!!! CREATION DATE
!!!     1998
!!! DESCRIPTION
!!!     Module to generate molecular surfaces
!!!     It might be useful to read York_JPhysChemA_1999_v103_p11060
!!! USES
!!!     YorkLib: DataTypes, UtilitiesMod, AngularQuadratureMod, ConstantsMod
!!!
!!!***

  USE DataTypes
  USE UtilitiesMod

  IMPLICIT NONE
  ! everything private by default
  private

  ! used for determining whether to split set when clipping surface points
  integer(i4b), parameter :: nmax = 80 

  ! Crds and radii of atomic centers, around which surface is constructed
  TYPE ParticleSetT
     REAL(SP), POINTER :: Crd(:,:) =>NULL()
     REAL(SP), POINTER :: Rad(:) =>NULL()
     INTEGER(I4B) :: SPLIT
     INTEGER(I4B) :: NP
     REAL(SP) :: Orig(3)
     REAL(SP) :: EigV(3)
     INTEGER(I4B), POINTER :: Index(:) =>NULL()
  END TYPE ParticleSetT

  ! Stores the different angular quadrature rules to apply to 
  ! the different atom centers, for constant pps Ntype=1
  TYPE SphereTypes
     INTEGER(I4B) :: Ntype
     INTEGER(I4B), POINTER :: Npt(:) =>NULL()
     REAL(SP), POINTER :: RadMax(:) =>NULL()
     TYPE(SP2Dptr), POINTER :: Pt(:) =>NULL()
     TYPE(SP1Dptr), POINTER :: Wt(:) =>NULL()
  END TYPE SphereTypes

  ! These components are 1:npt
  TYPE AtomDiscritization 
     INTEGER(I4B) :: Npt      ! number of kept points around atom
     INTEGER(I4B) :: Type     ! type used (from SphereTypes)
     INTEGER(I4B), POINTER :: KeepPoint(:) => NULL()  ! Indexes of points to keep
     REAL(SP), POINTER :: Winv(:) => NULL()           ! Inverse Weights of points
     INTEGER(I4B), POINTER :: SwitchCount(:) =>NULL() ! number of times each point was switched
     real(sp) :: Rsw          ! switching radius for atom
     real(sp) :: res_eff      ! effective resolution
  END TYPE AtomDiscritization

  ! Contains 
  TYPE localSurfaceT
     integer(i4b) :: AngQuad_type
     integer(i4b) :: discretization_type
     integer(i4b) :: pts_per_sphere
     REAL(SP) :: Resolution
     REAL(SP) :: ProbeRadius
     REAL(SP) :: GammaSwitch
     REAL(SP) :: Pexp
     REAL(SP) :: AlphaShift
     LOGICAL(LGD) :: verbose = .FALSE.
     TYPE(SphereTypes) :: Sphere
     ! these are dimensioned 1:natom centers
     TYPE(AtomDiscritization), POINTER :: sphere_ctrs(:) =>NULL()
  END TYPE LocalSurfaceT

!!!****s* SurfaceMod/surface_t
!!!
!!! NAME
!!!     surface_t -- a data type for molecular surfaces
!!! SEE ALSO
!!!     SurfaceMod/SASurface for member definitions
!!! SOURCE
!!!
  type :: surface_t
     ! inputs
     real(SP) :: Resolution = 1.0_sp
     real(SP) :: ProbeRadius = 0.0_sp 
     real(SP) :: GammaSwitch = 0.0_sp
     real(SP) :: Pexp = 1.00_sp
     real(sp) :: r_thresh = 0.0_sp
     REAL(SP) :: AlphaShift = -1.0_SP ! shift parameter in Eq. 62 and E-6
     REAL(SP) :: rscale = 1.0_SP
     LOGICAL(LGD) :: verbose = .FALSE. ! print out info if true
     character(len=256) :: surf_file =''
     integer(i4b) :: AngQuad_type = 0
     integer(i4b) :: discretization_type = 0
     integer(i4b) :: pts_per_sphere = 110
     integer(i4b) :: iCDR = 0
     ! factors for CavitEnergy 
     REAL(SP) :: g_cav =0.10368_SP ! from 72 dyn/cm in kcal/mol/Ang^2
     REAL(SP) :: cav_rscale = 1.0_SP 
     REAL(SP) :: dr_rscale = 1.0_SP
     REAL(SP) :: pressure = 1.0_SP ! Atm

     ! outputs
     real(SP) :: SA                                 ! Molecular surface Area =SUM(ASA)
     real(SP) :: Vol                                ! Molecular volume
     real(SP), pointer :: SECRD(:,:) => Null()      ! Surface Element Crds
     real(SP), pointer :: SESCALE(:) => Null()      ! Surface Element Scale Factors
     real(SP), pointer :: COSMOZETA(:) => Null()    ! Cosmo Gaussian Exponents
     real(SP), pointer :: ASA(:) => Null()          ! Atomic Surface Areas
     integer(i4b), pointer :: KeptIndex(:)=> Null() ! Dim(1:# of kept atoms) in values in
     integer(i4b), pointer :: ASAINDEX(:) => Null() ! Index of Last surface pt for a given atom
     ! Wts stores weights including r^2
     REAL(SP), POINTER :: Wts(:) => Null()          ! Surface element weights for computing surface area
     ! GradSwf stores dS/dr_hat times gradr_hat without the Delta function
     REAL(SP), POINTER :: GradSwf(:,:,:) => Null()  ! gradient of switching function (second term in Eq 74)
  end type surface_t
!!!***

  ! parameter to use original definition of Rin and Rout for switching regions
  ! or to use revised definitions, which have proper limit in the case
  ! of an infinitly discritized sphere intersecting with a finite discritized sphere
  LOGICAL(LGD), PARAMETER :: org_Rin_Rout_def = .TRUE.

  ! parameter to use original Sik definition in construction of diag. elements of A or
  ! use sqrt(Sik) for diag. elements of A so that 'dimentionality of Sik' is 'consistent' 
  ! between use in A matrix and computing Surface Area.
  LOGICAL(LGD), PARAMETER,public :: org_Sik_in_A0_diag = .TRUE.
  
  ! sescale values smaller than this are set to this number
  ! this prevents problems in computation of the A diag elements
  ! when for all intents and purposes the point should be off
  REAL(SP), PARAMETER :: sescale_tol = 1.0E-37_SP

  PUBLIC :: SASurface, visualize_surface, switch_pot, &
       deallocate_surface, deallocateSurfaceArrays,cpt_A,cpt_Ainv, &
       write_surface, read_surface, surface_t, copy_surface_scalar_inputs, &
       read_surface_parms, print_surface_parms, &
       MolecularNonelecEnergy, AtomicNonelecEnergy, AnalytMolNonelecGradient
  
  interface SASurface
     module procedure SASurface_long, SASurface_short
  end interface
  
  INTERFACE SplitSets
     MODULE PROCEDURE SplitSets2, SplitSets3
  END INTERFACE
  
  INTERFACE cpt_A
     MODULE PROCEDURE cpt_A_surface1, cpt_A_surface2
  END INTERFACE
  
  INTERFACE cpt_Ainv
     MODULE PROCEDURE cpt_Ainv_switch, cpt_Ainv_noswitch
  END INTERFACE
  
CONTAINS
  
!!!****f* SurfaceMod/SASurface
!!!
!!! NAME
!!!     SASurface
!!! USAGE
!!!     SASurface( Crd, Rad, DT, PTS_PER_SPHERE, RESOLUTION, &
!!!       PROBERADIUS, GAMMASWITCH, SA, ASA, ASAINDEX, SECRD, &
!!!       SESCALE, PEXP, COSMOZETA, WTS,GCAVATOM, MAXPTS, MINPTS, AVGPTS )
!!!                                 or
!!!     SASurface( Crd, Rad, Surface, MAXPTS, MINPTS, AVGPTS)
!!!
!!! DESCRIPTION
!!!     Generates a molecular surface for a molecule of atoms at Crd with
!!!     radii Rad at a resolution of RESOLUTION points per square bohr of
!!!     surface area.
!!!
!!!     Crd(1:3,1:number of atoms) is the coordinates of the atoms
!!!     Rad(1:number of atoms) are the radii of the atoms
!!!
!!!     Surface is of type surface_t and contains the most commonly used input
!!!             and output variables in an easy to digest datatype.
!!!     
!!!     DT is the discretization type.  0 -> constant pts/sphere
!!!                                     1 -> constant resolution
!!!     PTS_PER_SPHERE specifies the minimum number of points per sphere
!!!     RESOLUTION is the number of points per square bohr on the surface
!!!                The default is 26 / FOUR_PI 
!!!     PROBERADIUS - probe radius in bohr
!!!                The default is 0.0 bohr
!!!     GAMMASWITCH - controls the degree of switching (\gamma_s in the paper)
!!!                The default is 0.0
!!!     SA - surface area of the molecule
!!!     ASA(1:number of kept atoms) - atomic surface area.  sum(asa)=sa
!!!     ASAINDEX(1:number of kept atoms) - index of the last surface point
!!!                of the a given atom
!!!     SECRD(1:3,1:number of surface elements) - coordinates of the
!!!                surface elements
!!!     SESCALE(1:number of surface elements) - scale factor for switching
!!!            function. (eqn. 70 in the paper) 
!!!     WTS(1:number of surface elements) - weights of the surface elements.
!!!     SEAREA(1:number of surface elements) - Area of each surface element. Can be
!!!          calculated by SESCALE * WTS
!!!     COSMOZETA(1:numer of suface elements) is the gaussian exponents of the
!!!              surface elements
!!!     GradSwf(1:3, 1:number of atoms, 1:number of surface elements) is the
!!!          gradient of switching function as the second term in Eq 74
!!!     PEXP - "overlap factor".  This parameter tries to account for 
!!!          "buried" atoms.  The default is 1.0
!!!     MAXPTS, MINPTS, and AVGPTS are the maximum, minimum, and average
!!!             number of points per sphere for the system.
!!!     SURF_FILE is the filename for atom/surface statistics
!!!*** 

  subroutine SASurface_short( Crd, Rad, MolSurf,MAXPTS, MINPTS,AVGPTS)
    USE Datatypes
    USE ConstantsMod
    USE AngularQuadratureMod

    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: Crd(:,:)
    REAL(SP), INTENT(IN) :: Rad(:)
    type(surface_t), intent(inout) :: MolSurf
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: MAXPTS, MINPTS
    REAL(SP), OPTIONAL, INTENT(OUT) :: AVGPTS

    ! Local variables
    TYPE(localSurfaceT) :: LocSurf
    TYPE(ParticleSetT) :: set
    INTEGER(I4B) :: i, n, nkeep, keep_index

    real(sp) :: radtmp

    n=SIZE(Crd,2)
    IF(n /= SIZE(Rad)) THEN
      WRITE(6,*)'SurfaceMod: Number of atom centers doesnt equal number of radii',n,SIZE(Rad)
      STOP 'Fatal Error in Call to SurfaceMod:SASurface'
    ENDIF

    ! Check Parameter Ranges
    MolSurf%GAMMASWITCH = MIN(1.0_SP,MAX(MolSurf%GAMMASWITCH,0.0_SP))

    ! only keep particles with radii greater than r_thresh
    nkeep = count( (Rad(:)*MolSurf%rscale) > MolSurf%r_thresh )
    ! Allocate local variables
    allocate(set%Crd(3,nkeep),set%Rad(nkeep),set%Index(nkeep)) 

    if(nkeep==0)THEN
       WRITE(6,*)'WARNING: After scaling radii by ',MolSurf%rscale,&
            ' no particles remain with r > ', MolSurf%r_thresh
    ENDIF
    keep_index = 0
    DO i = 1,n
       radtmp=Rad(i)*MolSurf%rscale
       if( radtmp > MolSurf%r_thresh) then
          keep_index = keep_index + 1
          set%Crd(:,keep_index) = Crd(1:3,i)
          set%rad(keep_index) = radtmp
          set%Index(keep_index) = i
       else
          IF (MolSurf%VERBOSE) write(6,100)'Throwing away particle ',i,' with radius ',Rad(i),' scaled to ',radtmp
       end if
    END DO
100 FORMAT (A,I9,2(A,F11.3))

    IF(associated(MolSurf%KeptIndex))deallocate(MolSurf%KeptIndex)
    allocate(MolSurf%KeptIndex(nkeep))
    MolSurf%KeptIndex(1:nkeep) = set%Index(1:nkeep)

    ! re-initalize index for use with surface construction routines
    set%Index(:)    = (/ (i, i=1,nkeep) /)
    set%NP = nkeep

    CALL SetupLocalSurface(MolSurf,set%Rad, LocSurf)

    IF (MolSurf%VERBOSE) THEN
       write(6,*)'Untrimmed surface has ',SUM( LocSurf%sphere_ctrs(:)%Npt ),' points'
       write(6,'(A23,F0.2,A)')'          ', &
            & real(SUM( LocSurf%sphere_ctrs(:)%Npt ), sp) / nkeep ,' points per atom'
    END IF


    CALL GetOrigin(set%NP,set%Crd,set%Orig)
    CALL ClipSpheres(set,set,LocSurf)
    CALL ReorderParticleIndeces(1,set%NP,set)

    CALL ComputeProperties( set, LocSurf, MolSurf,MAXPTS, MINPTS,AVGPTS)

!!!*** Compute dS/dr_hat times gradr_hat divided by Swf without Delta function as in Eq 75 and Eq 76
    IF (MolSurf%GammaSwitch>0) THEN
       IF (MolSurf%VERBOSE) write(6,*)'Computing switching func gradients'
       CALL CptSwfGradient(set, LocSurf, MolSurf%SECRD, MolSurf%GradSwf)
    END IF

    ! DEALLOCATE local variables from heap
    DO i=1,nkeep
       DEALLOCATE( LocSurf%sphere_ctrs(i)%KeepPoint)
       DEALLOCATE( LocSurf%sphere_ctrs(i)%Winv)
       DEALLOCATE( LocSurf%sphere_ctrs(i)%SwitchCount) 
    END DO
    deallocate(LocSurf%sphere_ctrs)

    DEALLOCATE( LocSurf%Sphere%Npt )
    DEALLOCATE( LocSurf%Sphere%RadMax )
    ! Need to deallocate pointer subcomponents first
    DO i=1,size(LocSurf%Sphere%Wt)
       deallocate(LocSurf%Sphere%Wt(i)%Sp1dptr)
       deallocate(LocSurf%Sphere%Pt(i)%Sp2dptr)
    ENDDO
    DEALLOCATE( LocSurf%Sphere%Pt )
    DEALLOCATE( LocSurf%Sphere%Wt )    
    deallocate(set%Crd,set%Rad,set%Index)

  end subroutine SASurface_short

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SASurface_long( Crd, Rad, AQTYPE, DT, PTS_PER_SPHERE, RESOLUTION, &
       PROBERADIUS, GAMMASWITCH, VOL, SA, ASA, ASAINDEX,KeptIndex, SECRD, SESCALE, &
       PEXP, COSMOZETA, MAXPTS, MINPTS, AVGPTS, WTS, R_THRESH,RSCALE, &
       GradSWF, SURF_FILE, VERBOSE, ALPHASHIFT)

    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: Crd(:,:)
    REAL(SP), INTENT(IN) :: Rad(:)
    
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: DT
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: AQTYPE
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: PTS_PER_SPHERE

    REAL(SP), OPTIONAL, INTENT(IN) :: RESOLUTION
    REAL(SP), OPTIONAL, INTENT(IN) :: PROBERADIUS
    REAL(SP), OPTIONAL, INTENT(IN) :: GAMMASWITCH
    REAL(SP), OPTIONAL, INTENT(IN) :: PEXP
    REAL(SP), OPTIONAL, INTENT(IN) :: ALPHASHIFT
    REAL(SP), OPTIONAL, INTENT(IN) :: R_THRESH
    REAL(SP), OPTIONAL, INTENT(IN) :: RSCALE
    CHARACTER(LEN=*), optional, intent(in) :: SURF_FILE
    LOGICAL, OPTIONAL, INTENT(IN) :: VERBOSE 

    REAL(SP), OPTIONAL, INTENT(OUT) :: VOL
    REAL(SP), OPTIONAL, INTENT(OUT) :: SA
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: MAXPTS, MINPTS
    REAL(SP), OPTIONAL, INTENT(OUT) :: AVGPTS

    INTEGER(I4B), OPTIONAL, POINTER :: KeptIndex(:)
    INTEGER(I4B), OPTIONAL, POINTER :: ASAINDEX(:)

    REAL(SP), OPTIONAL, POINTER :: ASA(:)
    REAL(SP), OPTIONAL, POINTER :: SECRD(:,:)
    REAL(SP), OPTIONAL, POINTER :: COSMOZETA(:)
    REAL(SP), OPTIONAL, POINTER :: SESCALE(:)
    REAL(SP), OPTIONAL, POINTER :: WTS(:)

    REAL(SP), OPTIONAL, POINTER :: GradSwf(:,:,:)

    type(surface_t) :: Surface

    ! Set intent IN variables
    IF ( PRESENT(DT) ) surface%discretization_type = DT
    IF ( PRESENT(AQTYPE) )surface%AngQuad_type=AQTYPE
    IF ( PRESENT(PTS_PER_SPHERE) ) surface%pts_per_sphere = PTS_PER_SPHERE
    IF ( PRESENT(RESOLUTION) ) Surface%Resolution = RESOLUTION
    IF ( PRESENT(PROBERADIUS) ) Surface%ProbeRadius = PROBERADIUS
    IF ( PRESENT(GAMMASWITCH) ) Surface%GammaSwitch = GAMMASWITCH
    IF ( PRESENT(PEXP) ) Surface%Pexp = PEXP
    IF ( PRESENT(ALPHASHIFT) ) Surface%AlphaShift = ALPHASHIFT
    IF ( PRESENT(R_THRESH) ) Surface%r_thresh = R_THRESH
    IF ( PRESENT(RSCALE) ) Surface%rscale = RSCALE
    IF ( PRESENT(VERBOSE) ) Surface%verbose = VERBOSE
    IF ( PRESENT(SURF_FILE) ) Surface%surf_file = SURF_FILE

    ! Call Short routine
    CALL SASurface_short(Crd, Rad, Surface,MAXPTS,MINPTS,AVGPTS)

    IF ( PRESENT(VOL))VOL = Surface%VOL
    IF ( PRESENT( SA))SA  = Surface%SA

    ! Set intent OUT variables
    IF ( PRESENT(KeptIndex))THEN
       IF(ASSOCIATED(KeptIndex))DEALLOCATE(KeptIndex)    
       KeptIndex => Surface%KeptIndex
       Surface%KeptIndex=>null()
    endif
    IF ( PRESENT(ASAINDEX))THEN
       IF(ASSOCIATED(ASAINDEX))DEALLOCATE(ASAINDEX)
       ASAINDEX=>Surface%ASAINDEX
       Surface%ASAINDEX=>NULL()
    ENDIF
    IF ( PRESENT(ASA))THEN
       IF(ASSOCIATED(ASA))DEALLOCATE(ASA)
       ASA=>Surface%ASA
       Surface%ASA=>NULL()
    ENDIF
    IF ( PRESENT(SECRD))THEN
       IF(ASSOCIATED(SECRD))DEALLOCATE(SECRD)
       SECRD=>Surface%SECRD
       Surface%SECRD=>NULL()
    ENDIF
    IF ( PRESENT(COSMOZETA))THEN
       IF(ASSOCIATED(COSMOZETA))DEALLOCATE(COSMOZETA)
       COSMOZETA=>Surface%COSMOZETA
       Surface%COSMOZETA=>NULL()
    ENDIF
    IF ( PRESENT(SESCALE))THEN
       IF(ASSOCIATED(SESCALE))DEALLOCATE(SESCALE)
       SESCALE=>Surface%SESCALE
       Surface%SESCALE=>NULL()
    ENDIF
    if ( PRESENT(WTS))THEN
       IF(ASSOCIATED(WTS))DEALLOCATE(WTS)
       WTS=>Surface%WTS
       Surface%WTS=>NULL()
    ENDIF
    IF ( PRESENT(GradSwf))THEN
       IF(ASSOCIATED(GradSwf))DEALLOCATE(GradSwf)
       GradSwf=>Surface%GradSwf
       Surface%GradSwf=>null()
    ENDIF
    CALL deallocateSurfaceArrays(Surface)

  END SUBROUTINE SASurface_long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ComputeProperties( set, LSurf, MSurf,MINPTS, MAXPTS, AVGPTS)
    use filemod

    TYPE(ParticleSetT), INTENT(IN) :: set
    TYPE(localSurfaceT), INTENT(IN) :: LSurf
    type(surface_t), intent(inout) :: MSurf
    INTEGER(I4B), OPTIONAL, INTENT(OUT) :: MAXPTS, MINPTS
    REAL(SP), OPTIONAL, INTENT(OUT) :: AVGPTS

    INTEGER(I4B) :: i,n,npt,ipt,TypeI,IndexI,KeepI, NSE
    REAL(SP) :: Sarea, sea, asea,wt, wtR2, scale, radi, zeta, percent,vtmp

    integer (i4b) :: switched, switch_count, sw_min, sw_max, sw_atom, unit

    switched = 0
    switch_count = 0
    sw_min = huge(sw_min)
    sw_max = 0

    n = set%NP

    NSE = SUM( LSurf%sphere_ctrs(:)%Npt )

    ! Allocate is necessary
    IF ( ASSOCIATED(MSurf%ASA) )       DEALLOCATE(MSurf%ASA)
    allocate(MSurf%ASA(n))
    IF ( ASSOCIATED(MSurf%ASAINDEX) )  DEALLOCATE(MSurf%ASAINDEX)
    allocate(MSurf%ASAINDEX(n))
    IF ( ASSOCIATED(MSurf%SECRD) )     DEALLOCATE(MSurf%SECRD)
    allocate(MSurf%SECRD(3,NSE))
    IF ( ASSOCIATED(MSurf%COSMOZETA) ) DEALLOCATE(MSurf%COSMOZETA)
    allocate(MSurf%COSMOZETA(NSE))
    if ( associated(MSurf%WTS) )       DEALLOCATE(MSurf%WTS)
    allocate(MSurf%WTS(NSE))
    IF ( ASSOCIATED(MSurf%SESCALE) )   DEALLOCATE(MSurf%SESCALE)
    allocate(MSurf%SESCALE(NSE))

    unit = -1
    IF(len_trim(MSurf%SURF_FILE)>0)THEN
       call OpenFile(trim(MSurf%SURF_FILE) // '.stat', unit)
       write(unit,'(A4,2A6,2A12,2A12)')&
            '# ','radius','Rsw','npt before','npt after','% switched','res_eff'
    ENDIF

    Sarea = 0.0_SP
    vtmp = 0.0_SP
    npt = 0
    ATOM_LOOP : DO i = 1, n
       IndexI = set%Index(i)
       radi = set%Rad(i)
       TypeI = LSurf%sphere_ctrs(IndexI)%Type

       sea = 0.0_SP
       asea=0.0_SP       

       CALL CosmoGaussianExponent(LSurf%AngQuad_type,LSurf%Sphere%Npt(TypeI),zeta)
       
       sw_atom = 0

       POINT_LOOP : DO ipt = 1, LSurf%sphere_ctrs(IndexI)%Npt
          KeepI = LSurf%sphere_ctrs(IndexI)%KeepPoint(ipt)

          npt = npt + 1
          MSurf%SECRD(1:3,npt) = set%Crd(:,i) &
               + radi * LSurf%sphere%pt(TypeI)%sp2dptr(1:3,KeepI) 

          IF ( LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI) > 0 ) THEN
             ! Collect switching stats
             switched = switched + 1
             sw_atom = sw_atom + 1
             switch_count = switch_count &
                  + LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI)
             if ( LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI) > sw_max ) &
                  then
                sw_max = LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI)
             elseif (LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI) < sw_min) &
                  then
                sw_min = LSurf%sphere_ctrs(IndexI)%SwitchCount(KeepI)
             end if

             scale = EXP(LSurf%Pexp * LSurf%sphere_ctrs(IndexI)%Winv(KeepI))
          ELSE
             scale = 1.0_SP
          END IF

          wt    = LSurf%Sphere%Wt(TypeI)%SP1Dptr(KeepI)
          wtR2  = wt*set%Rad(i)**2

          !$ AlphaShift is derived using Swf (not Swf^2) so scale must NOT be scale^2
          sea = wtR2 * scale                     ! individual surface element area
          
          asea = asea + sea                      ! running atomic surface area
          vtmp = vtmp + (set%Rad(i)*sea)/3.0_SP  ! running molecular volume
 
          MSurf%SESCALE(npt) = MAX(scale,sescale_tol) 
          MSurf%COSMOZETA(npt) = zeta / ( SQRT( wt ) * set%Rad(i) )
          MSurf%WTS(npt) = wtR2
       END DO POINT_LOOP

       ! sum up total molecular surface area
       Sarea = Sarea + asea

       MSurf%ASAINDEX(i) = npt
       MSurf%ASA(i) = asea

       if ( unit > -1) then
          ! write atom stats
          if ( sw_atom > 0 ) then
             percent = real(sw_atom) / real(LSurf%sphere_ctrs(IndexI)%Npt) * 100.0_sp
          else
             percent = 0.0_sp
          end if

          write(unit,10)IndexI,set%Rad(IndexI),LSurf%sphere_ctrs(IndexI)%Rsw,&
               LSurf%Sphere%Npt(TypeI), LSurf%sphere_ctrs(IndexI)%Npt, percent, &
               LSurf%sphere_ctrs(IndexI)%res_eff
       end if
    END DO ATOM_LOOP

    MSurf%SA = Sarea
    MSurf%VOL = vtmp

    IF ( PRESENT(MINPTS) ) MINPTS = MINVAL( LSurf%Sphere%Npt(:) )    
    IF ( PRESENT(MAXPTS) ) MAXPTS = MAXVAL( LSurf%Sphere%Npt(:) )
    IF ( PRESENT(AVGPTS) ) AVGPTS = &
         & REAL( SUM( LSurf%Sphere%Npt(LSurf%sphere_ctrs(1:n)%Type) ), SP ) / n

    IF (LSurf%verbose) THEN
       write(6,*)'Molecular Surface Area=',MSurf%sa
       write(6,*)'Molecular Volume=',MSurf%Vol
       write(6,*)'Final Nse=',nse 
       if (LSurf%gammaswitch > 0.0_sp ) then
          write(6,*)'Switching stats -'
          write(6,*)switched,'points are in a switching region'
          write(6,'(A19, F8.1)')'Percent switched:', &
               (real(switched,sp) / real(nse,sp)) * 100.0_sp
          if ( switched > 0 ) then
             write(6,'(A19, F8.1)')'          Average:' &
                  ,real(switch_count,sp) / switched
             write(6,*)'          Minimum:',sw_min
             write(6,*)'          Maximum:',sw_max
          end if
       end if
    END IF
    if (unit > -1 ) close(unit)

10  format(I4,2F6.2,2I12,2F12.1)

  END SUBROUTINE ComputeProperties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! set up local surface type
  SUBROUTINE SetupLocalSurface(MSurf, Rad, LSurf)

    USE ConstantsMod
    IMPLICIT NONE

    type(surface_t), intent(in) :: MSurf
    REAL(SP), INTENT(INOUT) :: Rad(:)
    TYPE(localSurfaceT),intent(inout) :: LSurf
    INTEGER(I4B) :: i, ipt, n, npt
    REAL(SP) :: rnse

    Lsurf%AngQuad_type = Msurf%AngQuad_type
    Lsurf%discretization_type = Msurf%discretization_type
    Lsurf%pts_per_sphere = MSurf%PTS_PER_SPHERE
    LSurf%Resolution = MSurf%RESOLUTION
    LSurf%ProbeRadius = MSurf%PROBERADIUS
    LSurf%GammaSwitch = MSurf%GammaSwitch 
    LSurf%Pexp = MSurf%PEXP
    LSurf%AlphaShift = MSurf%ALPHASHIFT
    LSurf%verbose = MSurf%VERBOSE

    n = SIZE(Rad)
    allocate(LSurf%sphere_ctrs(n))

    Rad = Rad + LSurf%ProbeRadius

    CALL SetupSpheres(Rad, LSurf%Sphere, LSurf%AngQuad_TYPE, &
         LSurf%discretization_type, Lsurf%pts_per_sphere, LSurf%Resolution, &
         LSurf%sphere_ctrs(:)%Type)

    DO i=1,n
       ! Get the number of points for this sphere
       npt  = LSurf%Sphere%Npt( LSurf%sphere_ctrs(i)%Type )
       rnse = real(npt,sp)
       LSurf%sphere_ctrs(i)%Npt = npt
       
       ! Compute effective resolution
       LSurf%sphere_ctrs(i)%res_eff = rnse / ( FOUR_PI * (Rad(i))**2)
       
       ! Compute switching radius for this atom
       LSurf%sphere_ctrs(i)%Rsw = LSurf%GammaSwitch * Rad(i) * sqrt(14.0_SP/rnse)
       
       allocate(LSurf%sphere_ctrs(i)%KeepPoint(npt))
       LSurf%sphere_ctrs(i)%KeepPoint(:) = (/ (ipt, ipt=1,npt) /)
       
       allocate(LSurf%sphere_ctrs(i)%Winv(npt))
       LSurf%sphere_ctrs(i)%Winv(:) = 0.0_SP
       
       allocate(LSurf%sphere_ctrs(i)%SwitchCount(npt))
       LSurf%sphere_ctrs(i)%SwitchCount(:)= 0 
    END DO
  END SUBROUTINE SetupLocalSurface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SetupSpheres(Rad,Sphere,AQ_TYPE,DT_TYPE,PTS_PER_SPHERE,resolution,ATOMTYPE)

    ! Fills DiscretizedSpheresT member of localSurfaceT
    USE ConstantsMod
    USE AngularQuadratureMod
    use errormod

    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: Rad(:)
    integer(i4b), intent(in) :: DT_type
    integer(i4b), intent(in) :: PTS_PER_SPHERE
    integer(i4b), intent(in) :: AQ_type
    REAL(SP), INTENT(IN) :: RESOLUTION ! points/area
    TYPE(SphereTypes) :: Sphere
    INTEGER(I4B), INTENT(OUT) :: ATOMTYPE(:)
    INTEGER(I4B), pointer :: npt_range(:)=>NULL()
    INTEGER(I4B) :: i, j, n
    INTEGER(I4B) :: npt_min, npt_max, npt

    npt=0
    ! Get the discretization type, 0 -> constant pts/sphere
    !                              1 -> constant res.
    select case(DT_type)

    case(0) ! constant pts/sphere
       npt_min = PTS_PER_SPHERE

       call GetNumAngularQuadPts(AQ_Type,npt,NMIN=npt_min,cosmo=.true.)

       ! How many rules do we have?
       Sphere%Ntype = 1

       ALLOCATE( Sphere%Npt(Sphere%Ntype) )
       ALLOCATE( Sphere%RadMax(Sphere%Ntype) )
       ALLOCATE( Sphere%Pt(Sphere%Ntype) )
       ALLOCATE( Sphere%Wt(Sphere%Ntype) )

       ! Fill in the different sphere types
       Sphere%Npt(1) = npt
       CALL AngularQuad(AQ_Type, npt, Sphere%Pt(1)%SP2Dptr, &
            Sphere%Wt(1)%SP1Dptr, ALLOC = .TRUE.)
       ATOMTYPE = 1

    case(1) ! constant resolution
       ! Get the smallest number of points there will be on a sphere and the
       ! corresponding angular quadrature rule number
       npt_min = CEILING( resolution * FOUR_PI * MINVAL( Rad )**2 )
       ! Get the largest number of points there will be on a sphere and the
       ! corresponding angular quadrature rule number
       npt_max = CEILING( resolution * FOUR_PI * MAXVAL( Rad )**2 )

       CALL GetAngularQuadPtRange(AQ_Type,npt_range,npt_min,npt_max,NPTS=.TRUE.,cosmo=.TRUE.)

       ! How many schemes do we have?
       Sphere%Ntype = SIZE(npt_range)

       ALLOCATE( Sphere%Npt(Sphere%Ntype) )
       ALLOCATE( Sphere%RadMax(Sphere%Ntype) )
       ALLOCATE( Sphere%Pt(Sphere%Ntype) )
       ALLOCATE( Sphere%Wt(Sphere%Ntype) )
       ! Fill in the different sphere types
       DO i=1,Sphere%Ntype
          Sphere%Npt(i) = npt_range(i)
          Sphere%RadMax(i) = SQRT( Sphere%Npt(i) / (resolution * FOUR_PI) )
          CALL AngularQuad(AQ_Type, Sphere%Npt(i), Sphere%Pt(i)%SP2Dptr, &
               Sphere%Wt(i)%SP1Dptr, ALLOC = .TRUE.)
       END DO
       deallocate(npt_range)

       ! Find the minimal rule number for each atom
       n = SIZE(Rad)
       IF(n /= SIZE(ATOMTYPE)) THEN
          CALL Error('SurfaceMod: n Radii and AtomTypes not equal',i1=n,i2=SIZE(ATOMTYPE))
          STOP 'Fatal error in call to SurfaceMod:SetupSpheres'
       ENDIF
       DO i=1,n
          ATOMTYPE(i) = -1
          DO j=1,Sphere%Ntype
             IF ( Rad(i) < Sphere%RadMax(j) ) THEN
                ATOMTYPE(i) = j
                EXIT
             END IF
          END DO
       END DO
       IF(ANY(ATOMTYPE <= 0))THEN
          CALL Error('SetupSpheres: Some atoms do not have sufficient resolution',WARNLEV=0)
          WHERE(ATOMTYPE <= 0)
             ATOMTYPE = Sphere%Ntype
          END WHERE
       ENDIF
    end select

  END SUBROUTINE SetupSpheres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE ClipSpheres(set1,set2,S)
    IMPLICIT NONE

    TYPE(localSurfaceT), INTENT(INOUT) :: S
    TYPE(ParticleSetT), INTENT(INOUT) :: set1, set2

    INTEGER(I4B) :: iDecision
    TYPE(ParticleSetT) :: set11, set12, set13, set21, set22, set23
    LOGICAL(LGD) :: SameSet

    SameSet = ASSOCIATED(set1%Crd,set2%Crd)
    CALL SplitDecision(set1,set2,iDecision)

    SELECT CASE(iDecision)

    CASE (0)
       CALL DirectTrim(set1,set2,S)
    CASE (2)
       IF ( set1%SPLIT == 1 ) THEN
          IF ( .NOT. SameSet ) THEN 
             IF ( set2%SPLIT == 1 ) THEN
                CALL ClipSpheres(set1,set2,S)
             ELSE IF ( set2%SPLIT == 2 ) THEN
                CALL SplitSets(set2,set21,set22)

                CALL ClipSpheres(set1,set21,S)
                CALL ClipSpheres(set1,set22,S)
             ELSE IF ( set2%SPLIT == 3 ) THEN
                CALL SplitSets(set2,set21,set22,set23)

                CALL ClipSpheres(set1,set22,S)
                CALL ClipSpheres(set1,set21,S)
                CALL ClipSpheres(set1,set23,S)
             END IF
          ELSE
             CALL ClipSpheres(set1,set1,S)
          END IF

       ELSE IF ( set1%SPLIT == 2 ) THEN
          CALL SplitSets(set1,set11,set12)
          IF ( .NOT. SameSet ) THEN
             IF ( set2%SPLIT == 1 ) THEN
                CALL ClipSpheres(set11,set2,S)
                CALL ClipSpheres(set12,set2,S)
             ELSE IF ( set2%SPLIT == 2 ) THEN
                CALL SplitSets(set2,set21,set22)

                CALL ClipSpheres(set11,set21,S)
                CALL ClipSpheres(set11,set22,S)
                CALL ClipSpheres(set12,set21,S)
                CALL ClipSpheres(set12,set22,S)
             ELSE IF ( set2%SPLIT == 3 ) THEN
                CALL SplitSets(set2,set21,set22,set23)

                CALL ClipSpheres(set11,set22,S)
                CALL ClipSpheres(set12,set22,S)
                CALL ClipSpheres(set11,set21,S)
                CALL ClipSpheres(set11,set23,S)
                CALL ClipSpheres(set12,set21,S)
                CALL ClipSpheres(set12,set23,S)
             END IF
          ELSE
             CALL ClipSpheres(set12,set12,S)
             CALL ClipSpheres(set11,set11,S)
             CALL ClipSpheres(set11,set12,S)
          END IF

       ELSE IF ( set1%SPLIT == 3 ) THEN
          CALL SplitSets(set1,set11,set12,set13)
          IF ( .NOT. SameSet ) THEN
             IF ( set2%SPLIT == 1 ) THEN
                CALL ClipSpheres(set12,set2,S)
                CALL ClipSpheres(set11,set2,S)
                CALL ClipSpheres(set13,set2,S)
             ELSE IF ( set2%SPLIT == 2 ) THEN
                CALL SplitSets(set2,set21,set22)

                CALL ClipSpheres(set12,set21,S)
                CALL ClipSpheres(set12,set22,S)
                CALL ClipSpheres(set11,set21,S)
                CALL ClipSpheres(set11,set22,S)
                CALL ClipSpheres(set13,set21,S)
                CALL ClipSpheres(set13,set22,S)
             ELSE IF ( set2%SPLIT == 3 ) THEN
                CALL SplitSets(set2,set21,set22,set23)

                CALL ClipSpheres(set12,set21,S)
                CALL ClipSpheres(set12,set22,S)
                CALL ClipSpheres(set12,set23,S)
                CALL ClipSpheres(set11,set21,S)
                CALL ClipSpheres(set11,set22,S)
                CALL ClipSpheres(set11,set23,S)
                CALL ClipSpheres(set13,set21,S)
                CALL ClipSpheres(set13,set22,S)
                CALL ClipSpheres(set13,set23,S)
             END IF
          ELSE
             CALL ClipSpheres(set12,set12,S)
             CALL ClipSpheres(set11,set11,S)
             CALL ClipSpheres(set13,set13,S)

             CALL ClipSpheres(set11,set12,S)
             CALL ClipSpheres(set12,set13,S)

             CALL ClipSpheres(set11,set13,S)
          END IF

       ELSE
          CALL error('ERROR in ClipSpheres: set1%SPLIT = set2%SPLIT = 0', &
               & WARNLev=5)
          STOP
       END IF

    END SELECT

  END SUBROUTINE ClipSpheres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SplitDecision(set1,set2,iDecision)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set1,set2

    INTEGER(I4B), INTENT(OUT) :: iDecision
    REAL(SP), PARAMETER :: ANI_TOL = 0.75_SP

    INTEGER(I4B) :: i
    REAL(SP) :: Rvec(3), r12, maxr1, maxr2, r, aniso
    LOGICAL(LGD) :: SameSet

    !    WRITE(*,*)'SplitDecision'
    CALL GetOrigin(set1%NP,set1%Crd,set1%Orig)
    CALL GetOrigin(set2%NP,set2%Crd,set2%Orig)

    SameSet = ASSOCIATED(set1%Crd,set2%Crd)

    IF ( .NOT. SameSet ) THEN

       Rvec = set1%Orig - set2%Orig
       r12 = SQRT( DOT_PRODUCT( Rvec, Rvec ) )
       Rvec = Rvec / r12

       ! Calculate maxr1, the maximum distance of an atomic sphere in set
       ! 1 from its origin in the direction of the origin of set 2
       maxr1 = 0.0_SP
       DO i=1,set1%NP
          r = DOT_PRODUCT( (set1%Orig-set1%Crd(:,i)), Rvec ) &
               + set1%Rad(i)
          maxr1 = MAX(r,maxr1)
       END DO

       ! Calculate maxr2, the maximum distance of an atomic sphere in set
       ! 2 from its origin in the direction of the origin of set 1
       maxr2 = 0.0_SP
       DO i=1,set2%NP
          r = DOT_PRODUCT( (set2%Crd(:,i)-set2%Orig), Rvec ) &
               + set2%Rad(i)
          maxr2 = MAX(r,maxr2)
       END DO

       IF ( (maxr1 + maxr2) < r12 ) THEN
          iDecision = 1 ! there is no overlap, do nothing
       ELSE IF (set1%NP*set2%NP <= Nmax*Nmax) THEN
          ! small set, direct trim
          iDecision = 0
       ELSE
          iDecision = 2 ! there is probably overlap, split

          IF ( set1%NP <= Nmax ) THEN
             set1%SPLIT = 1
          ELSE
             CALL GetSplitDirection(set1,ANISO=aniso)
             IF ( aniso < ANI_TOL .OR. set1%NP < 1.5_SP * Nmax ) THEN
                set1%SPLIT = 2
             ELSE
                set1%SPLIT = 3
             END IF
          END IF

          IF ( set2%NP <= Nmax ) THEN
             set2%SPLIT = 1
          ELSE
             CALL GetSplitDirection(set2,ANISO=aniso)
             IF ( aniso < ANI_TOL .OR. set2%NP < 1.5_SP * Nmax ) THEN
                set2%SPLIT = 2
             ELSE
                set2%SPLIT = 3
             END IF
          END IF
       END IF

    ELSE ! SameSet
       IF ( set1%NP*(set2%NP-1)/2 > Nmax*Nmax ) THEN
          iDecision = 2 ! there is probably overlap, split
          CALL GetSplitDirection(set2,ANISO=aniso)
          IF ( aniso < ANI_TOL .OR. set1%NP < 1.5_SP * Nmax ) THEN
             set1%SPLIT = 2
          ELSE
             set1%SPLIT = 3
          END IF
       ELSE
          iDecision = 0 ! who cares, just check it 
       END IF
    END IF

  END SUBROUTINE SplitDecision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SplitSets2(set,set1,set2)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set
    TYPE(ParticleSetT), INTENT(OUT) :: set1, set2


    CALL SplitSetAlongDirection(set,set1,set2,set%Orig,set%EigV)

  END SUBROUTINE SplitSets2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SplitSets3(set,set1,set2,set3)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set
    TYPE(ParticleSetT), INTENT(OUT) :: set1, set2, set3

    INTEGER(I4B) :: i, n1
    REAL(SP) :: Orig1(3),OD,dot
    TYPE(ParticleSetT) :: set12

    n1 = 0; Orig1 = 0.0_SP
    OD = DOT_PRODUCT( set%Orig, set%EigV )
    DO i=1,set%NP
       dot = DOT_PRODUCT( set%Crd(:,i), set%EigV )
       IF ( dot < OD ) THEN
          Orig1 = Orig1 + set%Crd(:,i); n1 = n1 + 1
       END IF
    END DO
    Orig1 = Orig1 / n1
    Orig1 = 0.33_SP * Orig1 + 0.67_SP * set%Orig      

    CALL SplitSetAlongDirection(set,set12,set3,Orig1,set%EigV)
    CALL SplitSetAlongDirection(set12,set1,set2,set12%Orig,set%EigV)

  END SUBROUTINE SplitSets3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SplitSetAlongDirection(set,set1,set2,Orig,Dir)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT),TARGET :: set
    TYPE(ParticleSetT), INTENT(OUT) :: set1, set2
    REAL(SP), INTENT(IN) :: Orig(3), Dir(3)

    INTEGER(I4B) :: i,j,np,np1
    REAL(SP), POINTER :: crd(:,:)=>NULL()
    REAL(SP) :: OD, dot

    i=0; np = set%NP; j = set%NP+1; np1=0; crd => set%Crd
    OD = DOT_PRODUCT(Orig,Dir)
    DO
       DO
          i=i+1
          dot = DOT_PRODUCT(crd(:,i),Dir)
          IF ( dot <= OD .OR. i == np ) EXIT
       END DO
       DO
          j=j-1
          dot = DOT_PRODUCT(crd(:,j),Dir)
          IF ( dot > OD ) np1 = j
          IF ( dot > OD .OR. j == 1 ) EXIT
       END DO
       IF ( i < j ) THEN
          CALL SwapParticles(set,i,j)
       ELSE
          EXIT
       END IF

    END DO

    set1%NP = np1
    set1%Crd      => set%Crd(1:3,1:np1)
    set1%Rad      => set%Rad(1:np1)
    set1%Index    => set%Index(1:np1)

    set2%NP = np - np1
    set2%Crd      => set%Crd(1:3,np1+1:)
    set2%Rad      => set%Rad(np1+1:)
    set2%Index    => set%Index(np1+1:)

    CALL GetOrigin(set1%NP,set1%Crd,set1%Orig)
    CALL GetOrigin(set2%NP,set2%Crd,set2%Orig)
  END SUBROUTINE SplitSetAlongDirection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DirectTrim(set1,set2,Surface)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set1, set2
    TYPE(localSurfaceT) :: Surface

    IF ( Surface%gammaswitch > 0 ) THEN
       CALL DirectTrimSwitch(set1,set2,Surface)
    ELSE
       CALL DirectTrimNoSwitch(set1,set2,Surface)
    END IF

  END SUBROUTINE DirectTrim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DirectTrimNoSwitch(seta,setb,Surface)
    use constantsmod
    IMPLICIT NONE

    TYPE(ParticleSetT), TARGET, INTENT(INOUT) :: seta, setb
    TYPE(localSurfaceT) :: Surface
    TYPE(ParticleSetT), pointer :: set1=>NULL(), set2=>NULL()

    INTEGER(I4B) :: i, j, ipt, Nkeep, IndexI, TypeI,nset
    LOGICAL(LGD) :: SameSet
    REAL(SP) :: radi, radi2, radj, radj2, sum_rad2, r2, rrv(3), rv(3)
    INTEGER(I4B), POINTER :: KeepI(:)
    REAL(SP) :: RinJ2,RinJ,RoutJ

    SameSet = ASSOCIATED(seta%Crd,setb%Crd)

    set1=>seta
    set2=>setb
    ! Loop through twice, 1st clip points in seta, second in setb
    do nset=1,2 
       ! trim surface elements around set 1
       SET1_LOOP_I : DO i=1,set1%NP
          IndexI = set1%Index(i); TypeI = Surface%sphere_ctrs(IndexI)%Type
          KeepI => Surface%sphere_ctrs(IndexI)%KeepPoint
          
          radi  = set1%Rad(i); radi2 = radi*radi
          IF ( Surface%sphere_ctrs(IndexI)%Npt == 0 ) CYCLE SET1_LOOP_I
          
          SET2_LOOP_J : DO j=1,set2%NP
             IF ( SameSet .AND. i == j ) CYCLE SET2_LOOP_J
             rv = set1%Crd(:,i) - set2%Crd(:,j)
             r2 = DOT_PRODUCT( rv, rv )
             
             RinJ = set2%Rad(j)
             RoutJ = RinJ
             radj  = RoutJ; radj2 = radj**2
             sum_rad2 = (radi + radj)**2
             
             IF ( r2 < sum_rad2 ) THEN
                Nkeep = 0
                rrv = set1%Crd(:,i) - set2%Crd(:,j)
                RinJ2 = RinJ**2
                DO ipt = 1, Surface%sphere_ctrs(IndexI)%Npt
                   rv = rrv + radi * Surface%Sphere%pt(TypeI)%sp2dptr(1:3,KeepI(ipt))
                   r2 = DOT_PRODUCT( rv, rv )
                   IF ( r2 > RinJ2 ) THEN
                      Nkeep = Nkeep + 1
                      KeepI(Nkeep) = KeepI(ipt)
                   END IF
                END DO
                Surface%sphere_ctrs(IndexI)%Npt = Nkeep
                IF ( Nkeep == 0 ) CYCLE SET1_LOOP_I       
             END IF
          END DO SET2_LOOP_J
       END DO SET1_LOOP_I
       IF (SameSet) EXIT
       set1=>setb
       set2=>seta
    enddo
    set1=>NULL()
    set2=>NULL()

  END SUBROUTINE DirectTrimNoSwitch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DirectTrimSwitch(seta,setb,Surface)
    use constantsmod
    IMPLICIT NONE

    TYPE(ParticleSetT), TARGET, INTENT(INOUT) :: seta, setb
    TYPE(localSurfaceT) :: Surface
    TYPE(ParticleSetT), pointer :: set1=>NULL(), set2=>NULL()

    INTEGER(I4B) :: i, j, ipt, IndexI, TypeI, Nkeep, nset, sel_rsw
    LOGICAL(LGD) :: SameSet
    REAL(SP) :: radi, radi2, radj, radj2, sum_rad2, r, r2, rrv(3), rv(3)
    REAL(SP) :: Swf, dSwf
    INTEGER(I4B), POINTER :: KeepI(:)
    REAL(SP), POINTER :: WinvI(:)
    REAL(SP) :: RinJ2, Rsw_vals(2), Rsw, AlphaJ, RinJ, RoutJ

    SameSet = ASSOCIATED(seta%Crd,setb%Crd)

    IF(.NOT. org_Rin_Rout_def)THEN
       sel_rsw=1
    ELSE
       sel_rsw=2
    ENDIF

    set1=>seta
    set2=>setb
    ! Loop through twice, 1st clip points in seta, second in setb
    do nset=1,2 
       ! trim surface elements around set 1
       SET1_LOOP_I : DO i=1,set1%NP
          ! Nkeep is keeping track of the number of points we are keeping.
          ! KeepI is reindexing them.
          IndexI = set1%Index(i); TypeI = Surface%sphere_ctrs(IndexI)%Type
          KeepI => Surface%sphere_ctrs(IndexI)%KeepPoint
          WinvI => Surface%sphere_ctrs(IndexI)%Winv
          
          Rsw_vals(1) = Surface%sphere_ctrs(IndexI)%Rsw
          
          radi  = set1%Rad(i); radi2 = radi*radi
          IF ( Surface%sphere_ctrs(IndexI)%Npt == 0 ) CYCLE SET1_LOOP_I
          
          SET2_LOOP_J : DO j=1,set2%NP
             IF ( SameSet .AND. i == j ) CYCLE SET2_LOOP_J
             Rsw_vals(2) = Surface%sphere_ctrs(set2%Index(j))%Rsw
             
             Rsw = Rsw_vals(sel_rsw)
             AlphaJ = AlphaF(set2%Rad(j), Rsw, Surface%AlphaShift)
             RinJ = set2%Rad(j) - AlphaJ * Rsw
             RoutJ = RinJ + Rsw
             
             rv = set1%Crd(:,i) - set2%Crd(:,j)
             r2 = DOT_PRODUCT( rv, rv )
             
             radj  = RoutJ; radj2 = radj**2
             sum_rad2 = (radi + radj)**2
             
             IF ( r2 < sum_rad2 ) THEN
                ! There is overlap
                Nkeep = 0
                RinJ2 = sign(RinJ**2,RinJ)
                IPT_LOOP: DO ipt = 1, Surface%sphere_ctrs(IndexI)%Npt
                   rrv = rv + radi * Surface%sphere%pt(typeI)%sp2dptr(1:3,KeepI(ipt))
                   
                   r2 = DOT_PRODUCT( rrv, rrv )
                   
                   IF ( r2 > radj2 ) THEN
                      Nkeep = Nkeep + 1
                      KeepI(Nkeep) = KeepI(ipt)
                   ELSE IF ( r2 > RinJ2 ) THEN
                      r = ( SQRT(r2) - RinJ ) / Rsw
                      Swf = SwitchF( r, dSwf )
                      if ( Swf <= 0.0_sp ) then
                         ! We really do NOT want this point
                         !write(6,*)'Skipped a zero weight'              
                         cycle IPT_LOOP
                      end if
                      
                      Nkeep = Nkeep + 1
                      KeepI(Nkeep) = KeepI(ipt)
                      
                      if ( Swf >= 1.0_sp ) then
                         ! don't switch, just keep it
                         !write(6,*)'did not switch a one weight', ipt
                         cycle ipt_loop
                      end if
                      
                      WinvI(KeepI(ipt)) = WinvI(KeepI(ipt)) + WinvF(Swf)
                      Surface%sphere_ctrs(IndexI)%SwitchCount(KeepI(ipt)) = &
                           & Surface%sphere_ctrs(IndexI)%SwitchCount(KeepI(ipt)) + 1
                      
                   END IF
                END DO IPT_LOOP
                Surface%sphere_ctrs(IndexI)%Npt = Nkeep
                IF ( Nkeep == 0 ) CYCLE SET1_LOOP_I
             END IF
          END DO SET2_LOOP_J
       END DO SET1_LOOP_I
       IF (SameSet) EXIT
       set1=>setb
       set2=>seta
    enddo
    set1=>NULL(); set2=>NULL()

  END SUBROUTINE DirectTrimSwitch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION SwitchF(r,dSwf) RESULT(Swf)
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: r
    REAL(SP) :: Swf, dSwf

    REAL(SP) :: r2

    IF ( r < 0 ) THEN
       Swf = 0.0_SP
       dSwf = 0.0_SP
    ELSE IF ( r > 1 ) THEN
       Swf = 1.0_SP
       dSwf = 0.0_SP
    ELSE
       r2 = r*r
       Swf = r * r2 * ( 10.0_sp - 15.0_sp * r + 6.0_sp * r2 )
       !Swf = r * r2 * ( 10.0_sp + r * (6.0_sp * r - 15.0_sp))
       !Swf = MIN(MAX(Swf,0.0_SP),1.0_SP)
       dSwf = 30.0_SP*r2 * (r-1) * (r-1)
    END IF

    Swf = MIN(MAX(Swf,0.0_SP),1.0_SP)

  END FUNCTION SwitchF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE ReorderParticleIndeces(nl,nu,set)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set
    INTEGER(I4B), INTENT(IN) :: nl,nu

    INTEGER(I4B) :: i,j,np,nmid

    !INTEGER(I4B), TARGET, SAVE :: PROC=0, TLEV=4, VLEV=4

    np = nu-nl+1

    IF ( np == 2 ) THEN
       IF ( set%index(nl) > set%index(nu) ) THEN
          CALL SwapParticles(set,nl,nu)
       END IF

    ELSE IF ( np > 2 ) THEN

       nmid = (nl+nu)/2
       i = nl-1; j=nu+1

       DO
          DO
             i=i+1
             IF ( set%index(i) >= nmid .OR. i == nu ) EXIT
          END DO
          DO
             j=j-1
             IF ( set%index(j) < nmid .OR. j == nl ) EXIT
          END DO
          IF ( i < j ) THEN
             CALL SwapParticles(set,i,j)
          ELSE
             EXIT
          END IF

       END DO

       CALL ReorderParticleIndeces(nl,nmid,set)
       CALL ReorderParticleIndeces(nmid,nu,set)

    END IF

  END SUBROUTINE ReorderParticleIndeces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GetOrigin(N,Crd,Orig)
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: N
    REAL(SP), INTENT(IN) :: Crd(:,:)
    REAL(SP), INTENT(OUT) :: Orig(3)

    INTEGER(I4B) :: i


    Orig = 0.0_SP
    if ( n == 0 ) RETURN
    DO i=1,N
       Orig = Orig + Crd(:,i)
    END DO
    Orig = Orig / N

  END SUBROUTINE GetOrigin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GetFluctuationMatrix(N,Crd,Orig,Fmat)
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: N
    REAL(SP), INTENT(IN) :: Crd(:,:)
    REAL(SP), INTENT(OUT) :: Orig(3), Fmat(3,3)

    INTEGER(I4B) :: i,k

    Orig = 0.0_SP; Fmat = 0.0_SP
    DO i=1,N
       Orig = Orig + Crd(:,i)
       DO k=1,3
          Fmat(1:3,k) = Fmat(1:3,k) + Crd(k,i)*Crd(:,i)
       END DO
    END DO
    Orig = Orig / N; Fmat = Fmat / N
    DO k=1,3
       Fmat(1:3,k) = Fmat(1:3,k) - Orig(k)*Orig(1:3)
    END DO

  END SUBROUTINE GetFluctuationMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GetSplitDirection(set,ANISO)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set
    REAL(SP), OPTIONAL, INTENT(OUT) :: ANISO

    REAL(SP) :: Orig(3),Fmat(3,3)

    CALL GetFluctuationMatrix(set%NP,set%Crd,Orig,Fmat)
    CALL GetLargestEigenvectorOf3x3(Fmat,set%EigV,ANISO=ANISO)

  END SUBROUTINE GetSplitDirection

  SUBROUTINE GetLargestEigenvectorOf3x3(M,EigMax,ANISO)
    USE ConstantsMod
    IMPLICIT NONE

    ! Subroutine gets the eigenvector (EigMax) corresponding to the largest
    ! eigenvalue of a 3x3 matrix (M).  The norm of EigMax gives the
    ! anisotropy of the distribution: 
    !    aniso = 1 - SQRT((Eig1 + Eig2)/(2 * EigMax))
    REAL(SP), INTENT(IN) :: M(3,3)
    REAL(SP), INTENT(OUT) :: EigMax(3)
    REAL(SP), OPTIONAL, INTENT(OUT) :: ANISO

    REAL(SP) :: a,b,c,q,r,theta,x1,x2,x3,t,sqrtq

    a = -(M(1,1)+M(2,2)+M(3,3))
    b = M(1,1)*(M(2,2)+M(3,3))+M(2,2)*M(3,3) &
         - M(1,2)*M(2,1) - M(1,3)*M(3,1) - M(2,3)*M(3,2)
    c = -M(1,1)*M(2,2)*M(3,3) - M(2,1)*M(3,2)*M(3,1) &
         -  M(2,1)*M(3,2)*M(3,1) + M(2,2)*M(3,1)*M(1,3) &
         +  M(3,3)*M(2,1)*M(1,2) + M(1,1)*M(3,2)*M(3,3)
    q = ABS(a*a-3*b)/9; sqrtq = SQRT(q)
    r = (2*a**3-9*a*b+27*c)/54
    t = r/sqrtq**3; IF (t > 1.0_SP) t = 1.0_SP; IF (t < -1.0_SP)&
         & t = -1.0_SP
    theta = ACOS(t)
    x1 = -2*sqrtq*COS(theta/3)-a/3
    x2 = -2*sqrtq*COS((theta+2*PI)/3)-a/3
    x3 = -2*sqrtq*COS((theta-2*PI)/3)-a/3
    r = MAX(x1,x2,x3)
    EigMax(1)=(M(2,2)-r)*(M(3,3)-r)-M(3,2)**2
    EigMax(2)=(-M(2,1)*(M(3,3)-r)+M(3,1)*M(3,2))
    EigMax(3)=(-M(3,1)*(M(2,2)-r)+M(2,1)*M(3,2))
    t = SQRT(EigMax(1)**2+EigMax(2)**2+EigMax(3)**2)
    EigMax = (EigMax / t)
    IF ( PRESENT(ANISO) ) THEN
       ANISO = 1.0_SP - SQRT( (x1 + x2 + x3 - r) / (2 * r) )
    END IF

  END SUBROUTINE GetLargestEigenvectorOf3x3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SwapParticles(set, i1, i2)
    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(INOUT) :: set
    INTEGER(I4B), INTENT(IN) :: i1, i2
    REAL(SP) :: tmp_crd(3), tmp_rad
    INTEGER(I4B)::tmp_idx
  
    tmp_idx = set%index(i1)
    set%index(i1) = set%index(i2)
    set%index(i2) = tmp_idx

    tmp_crd = set%crd(1:3,i1)
    set%crd(1:3,i1) = set%crd(1:3,i2)
    set%crd(1:3,i2) = tmp_crd

    tmp_rad = set%rad(i1)
    set%rad(i1) = set%rad(i2)
    set%rad(i2) = tmp_rad
 
 
  END SUBROUTINE SwapParticles



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION PowerF(a,y) RESULT( p )
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a
    REAL(SP), INTENT(IN) :: y

    REAL(SP) :: p

    REAL(SP), PARAMETER :: TOL=0.0_SP

    IF ( a <= TOL ) THEN
       p = 0.0_SP
    ELSE
       p = a**y
    END IF

  END FUNCTION PowerF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION WinvF(Swf) RESULT( p )
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: Swf

    REAL(SP) :: p

    REAL(SP) :: TOL=tiny(1.0_sp)
    real(sp), parameter :: margin = 16.0_sp

    if ( Swf > tol ) then
       p = LOG(Swf)
    else
       p = margin * log(tol)
    end if

  END FUNCTION WinvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CosmoGaussianExponent(type,npt,zeta)
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: type, npt
    REAL(SP), INTENT(OUT) :: zeta
    REAL(SP) :: zet
    REAL(SP), PARAMETER :: a0=0.219217_SP,a1=5.64287_SP,a2=4.72383_SP,third=1.0_SP/3.0_SP

    Select Case(type) 
    CASE(1)
       ! Gaussian Product Rules
       zet = a0/((npt+a1)**(third)) + a2
    CASE default
       ! Lebedev Grids
       Select Case(npt)
       CASE(   6)
          zet = 4.8456607786812596_SP
       CASE(  14)
          zet = 4.8645871433400361_SP
       CASE(  26)
          zet = 4.8547822621992047_SP
       CASE(  38)
          zet = 4.9010581268545996_SP
       CASE(  50)
          zet = 4.8925067329562086_SP
       ! CASE(  74)  ! this has negative weights
       CASE(  86)
          zet = 4.8974137258096730_SP
       CASE( 110)
          zet = 4.9010106098795623_SP
       CASE( 146)
          zet = 4.8982518739257808_SP
       CASE( 170)
          zet = 4.9068551772536679_SP
       CASE( 194)
          zet = 4.9033764424897530_SP
       ! CASE( 230) ! this has negative weights
       ! CASE( 266) ! this has negative weights
       CASE( 302)
          zet = 4.9049808816981910_SP
       CASE( 350)
          zet = 4.8687947483240288_SP
       CASE( 434)
          zet = 4.9056734908003916_SP
       CASE( 590)
          zet = 4.9062407135926875_SP
       CASE( 770)
          zet = 4.9065643577907592_SP
       CASE( 974)
          zet = 4.9068516799891162_SP
       CASE(1202)
          zet = 4.9070409821651158_SP
       CASE(1454)
          zet = 4.9072102386914338_SP
       CASE(1730)
          zet = 4.9073327069138495_SP
       CASE(2030)
          zet = 4.9074449914280249_SP
       CASE(2354)
          zet = 4.9075308282568502_SP
       CASE(2702)
          zet = 4.9076097276694908_SP
       CASE(3074)
          zet = 4.9076728239484062_SP
       CASE(3470)
          zet = 4.9077314137132584_SP
       CASE(3890)
          zet = 4.9077796598116947_SP
       CASE(4334)
          zet = 4.9078246952619979_SP
       CASE(4802)
          zet = 4.9074912555330155_SP
       CASE(5294)
          zet = 4.9076207345293437_SP
       CASE(5810)
          zet = 4.9079290252272934_SP
       CASE DEFAULT
          write(6,*)"Lebedev CosmoGaussianExponent: I don't know what to do with npt=",npt
          zet = 0.0_SP
       end Select
    end select
    
    zeta = zet

  END SUBROUTINE CosmoGaussianExponent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine visualize_surface(surface,filename,in_ang)

    ! Write an xyz file for visualization in gopenmol to filename

    use FileMod
    use PeriodicTableMod

    type(surface_t), intent(in) :: surface
    character(len=*), intent(in) :: filename
    logical(LGD), optional,intent(in) :: in_ang

    integer(i4b) :: unit, nse, i, j
    real(sp) :: pos,scale
    integer(i4b), parameter :: mult=5, min=5
    logical(lgd) :: ang


    if(present(in_ang))THEN
      ang=in_ang
    else
      ang=.FALSE.
    endif

    if(ang)THEN
      scale= 1.0_SP/ANGSTROM
    else
      scale=1.0_SP
    endif

    unit = -1
    call OpenFile(FILE=filename,UNIT=unit)

    nse = size(surface%cosmozeta)

    if (nse <= 0) then
       write(6,*)'No surface to visualize'
       return
    end if

    write(unit,*)nse

    do i = 1, nse
       write(unit,'(1X,A5,3F10.3)')map(surface%sescale(i)),surface%secrd(1:3,i) * scale
    end do

    call CloseFile(unit)

    unit = -1
    call OpenFile(FILE='key.xyz',UNIT=unit)

    write(unit,*)mult+1

    pos = 0.0_sp
    do i = min, mult+min
       j = i
       if ( j == 10 ) j = 1
       write(unit,'(1X,A5,3F10.3)')ElementName(j),0.0,0.0,pos
       pos = pos + 1
    end do

    call CloseFile(unit)

  contains

    function map(x) result(a)

      ! maps x to an element

      use PeriodicTableMod

      real(sp), intent(in) :: x

      character(len=2) :: a

      integer(i4b) :: i

      ! 0<x<1 mapped to 5-10 with 10 = 1
      i = int(x*mult + min)
      if ( i == 10 ) i = 1
      a = ElementName(i)

    end function map

  end subroutine visualize_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine copy_surface_scalar_inputs(surf,surf_new)
     type(surface_t), intent(in) :: surf
     type(surface_t), intent(out) :: surf_new

     surf_new%Resolution  = surf%Resolution
     surf_new%ProbeRadius = surf%ProbeRadius
     surf_new%GammaSwitch = surf%GammaSwitch
     surf_new%Pexp        = surf%Pexp
     surf_new%r_thresh    = surf%r_thresh
     surf_new%AlphaShift  = surf%AlphaShift
     surf_new%rscale      = surf%rscale
     surf_new%verbose     = surf%verbose
     surf_new%surf_file   = surf%surf_file
     surf_new%AngQuad_type = surf%AngQuad_type
     surf_new%discretization_type = surf%discretization_type
     surf_new%pts_per_sphere = surf%pts_per_sphere
     surf_new%iCDR        = surf%iCDR
     surf_new%g_cav       = surf%g_cav
     surf_new%cav_rscale  = surf%cav_rscale
     surf_new%dr_rscale   = surf%dr_rscale
  end subroutine copy_surface_scalar_inputs



  subroutine deallocate_surface(surface)

    ! Deallocates surafce type and sets all scalars back to defaults

    type(surface_t), intent(inout) :: surface
    type(surface_t):: my_surf

    CALL deallocateSurfaceArrays(surface)
    ! default initalize all scalar variables
    surface=my_surf

  end subroutine deallocate_surface


  subroutine deallocateSurfaceArrays(surface)
    type(surface_t), intent(inout) :: surface
 
    if ( associated(surface%secrd) ) then
       deallocate(surface%secrd)
    end if
    if ( associated(surface%sescale) ) then
       deallocate(surface%sescale)
    end if
    if ( associated(surface%cosmozeta) ) then
       deallocate(surface%cosmozeta)
    end if
    IF ( ASSOCIATED(surface%ASA )) THEN
       DEALLOCATE(surface%ASA)
    ENDIF
    IF ( ASSOCIATED(surface%ASAINDEX) ) then
       DEALLOCATE(surface%ASAINDEX)
    END IF
    IF ( ASSOCIATED(surface%Wts) ) then
       DEALLOCATE(surface%Wts)
    END IF
    IF ( ASSOCIATED(surface%GradSwf) ) then
       DEALLOCATE(surface%GradSwf)
    END IF
    IF ( ASSOCIATED(surface%KeptIndex) ) then 
       DEALLOCATE(surface%KeptIndex)
    ENDIF

  end subroutine deallocateSurfaceArrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_surface(surface, FILE, UNIT, FORMATTED)

    use FileMod
    use StringMod
    use errormod
    use array_io
    type(surface_t), intent(in) :: surface
    character(len=*), optional :: FILE
    integer(i4b), optional :: UNIT
    logical, optional :: FORMATTED

    integer(i4b) :: unit_l
    character(len=8) :: form

    unit_l = -1

    if (present(FILE)) then
       if (optional_flag(FORMATTED,DEFAULT=.true.)) then
          call OpenFile(FILE,unit_l)
       else
          call OpenFile(FILE,UNIT=unit_l,FORM='unformatted')
       end if
    elseif (present(UNIT)) then
       unit_l = UNIT
    else
       unit_l = 6
    end if

    inquire(UNIT=unit_l,FORM=form)

    ! scalars
    if ( CompareString(form,'form') ) then
       ! Formatted I/O

       write(unit_l,*)surface%discretization_type, surface%pts_per_sphere
       write(unit_l,*)surface%resolution, surface%proberadius
       write(unit_l,*)surface%gammaswitch, surface%pexp, surface%r_thresh
       write(unit_l,*)surface%sa

    elseif ( CompareString(form,'unform') ) then
       ! Unformatted I/O

       write(unit_l)surface%discretization_type, surface%pts_per_sphere
       write(unit_l)surface%resolution, surface%proberadius
       write(unit_l)surface%gammaswitch, surface%pexp, surface%r_thresh
       write(unit_l)surface%sa

    else
       call error('unable to determine file format in write_surface')
    end if

    ! arrays
    call write_array(surface%secrd,unit_l)
    call write_array(surface%sescale,unit_l)
    call write_array(surface%cosmozeta,unit_l)

  end subroutine write_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_surface(surface, FILE, UNIT, FORMATTED)


    use FileMod
    use StringMod
    use errormod
    use array_io

    type(surface_t), intent(out) :: surface
    character(len=*), optional :: FILE
    integer(i4b), optional :: UNIT
    logical, optional :: FORMATTED

    integer(i4b) :: unit_l, iostat
    character(len=8) :: form

    unit_l = -1

    if (present(FILE)) then
       if (optional_flag(FORMATTED,DEFAULT=.true.)) then
          call OpenFile(FILE,UNIT=unit_l,STAT='old',IOSTAT=iostat)
       else
          call OpenFile(&
               FILE,UNIT=unit_l,FORM='unformatted',STAT='old',IOSTAT=iostat)
       end if

       if (iostat > 0) then
          call error('Error opening file in read_surface')
       end if

    elseif (present(UNIT)) then
       unit_l = UNIT
    else
       unit_l = 5
    end if

    inquire(UNIT=unit_l,FORM=form)

    ! scalars
    if ( CompareString(form,'form') ) then
       ! Formatted I/O
       read(unit_l,*)surface%discretization_type, surface%pts_per_sphere
       read(unit_l,*)surface%resolution, surface%proberadius
       read(unit_l,*)surface%gammaswitch, surface%pexp, surface%r_thresh
       read(unit_l,*)surface%sa
    elseif ( CompareString(form,'unform') ) then
       ! Unformatted I/O

       read(unit_l)surface%discretization_type, surface%pts_per_sphere
       read(unit_l)surface%resolution, surface%proberadius
       read(unit_l)surface%gammaswitch, surface%pexp, surface%r_thresh
       read(unit_l)surface%sa

    else
       call error('unable to determine file format in read_surface')
    end if

    ! arrays
    call read_array(surface%secrd,unit_l)
    call read_array(surface%sescale,unit_l)
    call read_array(surface%cosmozeta,unit_l)

  end subroutine read_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CptSwfGradient(set, Surface, SECRD, GradSwf)
!!!****f* SurfaceMod/CptSwfGradient
!!!
!!! Description
!!!     computes dS/dr_hat times gradr_hat divided by Swf as in 
!!!     Eq 75 and Eq 76 without Delta functions (leave out dimension m).
!!!
!!! Output
!!!     GradSwf - dS/dr_hat * gradr_hat / Swf
!!! Scratch
!!!     dSwf - dS/dr_hat
!!!***

    USE ConstantsMod

    IMPLICIT NONE

    TYPE(ParticleSetT), INTENT(IN) :: set
    TYPE(localSurfaceT), INTENT(IN) :: Surface
    REAL(SP), INTENT(IN):: SECRD(:,:)
    REAL(SP), POINTER :: GradSwf(:,:,:)

    INTEGER(I4B) :: n, i, j, ipt, IndexI, npt, nse, keepI,sel_rsw
    REAL(SP) :: radi, radi2, radj, radj2, sum_rad2, r2, r2new, rrv(3), rv(3)
    REAL(SP) :: Swf, rhat, dSwf, tmp
    REAL(SP) :: RinJ2, Rsw_vals(2), Rsw, AlphaJ, RinJ, RoutJ

    n = set%NP
    nse = SUM( Surface%sphere_ctrs(:)%Npt )
    npt = 0

    IF(.NOT. org_Rin_Rout_def)THEN
       sel_rsw=1
    ELSE
       sel_rsw=2
    ENDIF

    IF( ASSOCIATED(GradSwf)) Deallocate(GradSwf)
    ALLOCATE(GradSwf(1:nse,1:n,1:3))
    GradSwf(1:nse,1:n,1:3) = 0.0_SP

    AtomLoopI: DO i = 1, n
       IndexI = set%Index(i)
       Rsw_vals(1) = Surface%sphere_ctrs(IndexI)%Rsw

!!!**** loop over all not trimmed points on atom i
       PointLoop: DO ipt = 1, Surface%sphere_ctrs(IndexI)%Npt

          KeepI = Surface%sphere_ctrs(IndexI)%KeepPoint(ipt)
          npt = npt + 1

!!!**** Compute GradSwf if the point is partially switched else GradSwf=0.0
          IF ( Surface%sphere_ctrs(IndexI)%SwitchCount(KeepI) > 0 ) THEN
             radi  = set%Rad(i); radi2 = radi*radi
             
             AtomLoopJ: DO j = 1, n
                
                IF ( i == j ) CYCLE AtomLoopJ
                
                Rsw_vals(2) = Surface%sphere_ctrs(set%Index(j))%Rsw
                Rsw = Rsw_vals(sel_rsw)
                
                AlphaJ = AlphaF(set%Rad(j), Rsw, Surface%AlphaShift)
                RinJ = set%Rad(j) - AlphaJ * Rsw
                RoutJ = RinJ + Rsw
                
                rv = set%Crd(:,i) - set%Crd(:,j)
                r2 = DOT_PRODUCT( rv, rv )
                
                radj  = RoutJ; radj2 = radj**2
                sum_rad2 = (radi + radj)**2
                
!!!**** only if there is overlap
                IF ( r2 < sum_rad2 ) THEN

                   ! Sign function to make sure r2new stays bigger than
                   ! negative RinJ 
                   RinJ2 = sign(RinJ**2,RinJ)

                   rrv = SECRD(:,npt)-set%Crd(:,j)
                   r2new = DOT_PRODUCT( rrv, rrv )

                   IF ( (r2new < radj2) .AND. (r2new > RinJ2) ) THEN

                      rhat = ( SQRT(r2new) - RinJ ) / Rsw
                      Swf = SwitchF(rhat,dSwf)
                      IF ( Swf <= 0.0_sp ) THEN
                         !GradSwf(npt,j, :) = 0.0_SP
                         CYCLE PointLoop
                      END IF
                      tmp = SQRT(r2new) * Rsw * Swf

                      GradSwf(npt,j,:) = dSwf* rrv(:)/tmp

                   END IF !r2new
                END IF ! r2

             END DO AtomLoopJ
          END IF ! point switched
       END DO PointLoop

    END DO AtomLoopI
 
  END SUBROUTINE CptSwfGradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION AlphaF(rad, Rsw, AlphaShift) Result(AlphaJ)

    USE ConstantsMod

    REAL(SP), INTENT(IN) :: rad, Rsw, AlphaShift
    REAL(SP) :: AlphaJ
    REAL(SP) :: Rd
    REAL(SP), PARAMETER :: ONE_28TH = 1.0_SP/28.0_SP

    IF (AlphaShift > 1.0_SP .OR. AlphaShift < 0.0_SP ) THEN
       IF(rsw > 0.0_SP) THEN
          Rd = Rad / Rsw
          AlphaJ = 0.5_SP + Rd - SQRT( Rd*Rd - ONE_28TH)
       ELSE
          AlphaJ = 0.0_SP
       ENDIF
    ELSE
       AlphaJ = AlphaShift

    END IF
  END FUNCTION AlphaF

  SUBROUTINE AtomicNonelecEnergy(Surface,AST,ABC,Ecdr)
    USE ConstantsMod
    !
    ! Noneletrostatic Energy = SUM ( atomic surface tension * 
    !            atomic exposed surface area + atomic buried atom constant)
    ! 
    ! in atomic units

    TYPE(surface_t) :: Surface
    ! Ecdr = Nonelectrostatic energy from cavitation, dispersion, and repulsion
    REAL(SP), INTENT(OUT) :: Ecdr
    ! Atomic Surface Tensions for each atom
    REAL(SP), INTENT(IN) :: AST(:) 
    ! Atomic Buried Contributions
    REAL(SP), INTENT(IN) :: ABC(:)
    INTEGER(I4B) :: i,nea

    nea = SIZE(SURFACE%KeptIndex(:)) !number of exposed atoms

    Ecdr = 0.0_SP

    DO i=1,nea
       Ecdr = Ecdr + AST(SURFACE%KeptIndex(i)) * SURFACE%ASA(i) * (AU_DISTANCE_IN_ANGSTROM) ** 2
    END DO

    Ecdr = Ecdr + SUM(ABC(:))
    Ecdr = Ecdr * KCAL_PER_MOL

  END SUBROUTINE AtomicNonelecEnergy

  SUBROUTINE MolecularNonelecEnergy(Surface,Ecdr)
!!!****f* SurfaceMod/Cavitation
!!!
!!! Description
!!!     computes cavitation, dispersion and repulsion energy 
!!!       from surface and or volume
!!!    cavitation Energy = P*V*f^3 + gamma*S*f^2
!!!    Dispersion Energy =
!!!    Repulsion  Energy = 
!!! Parameters
!!!     Surface%g_cav are in kcal/(mol*ang*ang). 
!!!     Surface%cav_factor factor to scale 
!!!     Surface%dr_factor factor to scale
!!!     Need to be converted to a.u.first
!!! Output
!!!     Ecav in atomic units
!!!     EDisp in atomic units
!!!     ERep in atomic units
!!!     Ecdr in atomic units 
!!!
!!!***
    USE ConstantsMod

    TYPE(surface_t) :: Surface
    ! Parameter to convert volume dependent energy term to a.u.
    REAL(SP), PARAMETER :: vol_factor = 0.0_SP ! 
    REAL(SP) :: Surface_in_Ang2
    REAL(SP) :: Ecav,Edis,Erep
    REAL(SP), INTENT(OUT) :: Ecdr
    REAL(SP) :: cscale2
    REAL(SP) :: drscale2

    cscale2= Surface%cav_rscale * Surface%cav_rscale
    drscale2 = Surface%dr_rscale * Surface%dr_rscale

    Surface_in_Ang2 = Surface%sa * (AU_DISTANCE_IN_ANGSTROM) ** 2

    ! surface term
    Ecav = Surface%g_cav * Surface_in_Ang2 * cscale2 * KCAL_PER_MOL
    ! volume Term (too small to include)
    !Ecav = Ecav + vol_factor*Surface%Pressure*Surface%Vol*(cscale2*Surface%cav_rscale) 

    ! EDis and Erep all in the same equation
    !
    ! Tomasi J. Comp. Chem. 1991 v12 p784
    ! Linear fit to alkane p 787
    Edis = ( -0.03208_SP - 0.0767_SP * Surface_in_Ang2 * drscale2 ) * KCAL_PER_MOL
    Erep = 0.0_SP

    Ecdr = Ecav + Edis + Erep

  END SUBROUTINE MolecularNonelecEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AnalytMolNonelecGradient(Surface, CDRGrad)
!!!****f* SurfaceMod/CavitGradient
!!!
!!! Description
!!!     computes gradients of the cavitation energy
!!!     This only happens when switching is on
!!!     (If switching is off, pts are either on or off. IE
!!!     cavitation energy is constant for a given number 
!!!     of surface pts, ie NO GRADIENT!!)
!!! Output
!!!     CDRGrad
!!!***
    USE ConstantsMod

    TYPE(surface_t) :: Surface

    INTEGER  :: nkatm, nse, NseAtm, first, last, m, ik, k,index
    REAL(SP) :: tmp
    REAL(SP) :: cavtmp, distmp, reptmp
    REAL(SP), PARAMETER :: factor =  KCAL_PER_MOL * (AU_DISTANCE_IN_ANGSTROM **2)
    REAL(SP), INTENT(OUT) :: CDRGrad(:,:)
    REAL(SP) :: TMPGrad(1:3)

    nse = SIZE(surface%cosmozeta)
    nkatm =  SIZE(Surface%keptIndex)
 
    cavtmp = Surface%g_cav * Surface%PEXP * Surface%cav_rscale * Surface%cav_rscale * factor
    distmp = -0.0767_SP * Surface%PEXP * Surface%dr_rscale * Surface%dr_rscale * factor
    reptmp = 0.0_SP
    CDRGrad(:,:)=0.0_SP

!!!**** compute gradients of the surface area. This only applies when switching is on
    first = 1 
    GradAtomLoop: DO m = 1, nkatm
       index = Surface%keptIndex(m)
       last = surface%asaindex(m)
       NseAtm = last - first + 1
       TMPGrad(:) = 0.0_SP
       GradComponentLoop: DO k = 1, 3
!!!**** first term in my notes
          IF (NseAtm /=0) THEN
             DO ik = first, last
                tmp = SUM(Surface%GradSwf(ik,:,k))
                !$ Must use non squared Sik!
                TMPGrad(k)  = TMPGrad(k) + &
                     Surface%WTS(ik) * Surface%SESCALE(ik) * tmp
             END DO !ik
          END IF
!!!**** second term in my notes
          DO ik = 1, nse
             tmp = Surface%GradSwf(ik,m,k)
             TMPGrad(k)  = TMPGrad(k) - &
                  Surface%WTS(ik)* Surface%SESCALE(ik) * tmp
          END DO
       END DO GradComponentLoop
!!!**** compute cavitation dispersion repulsion energy gradients
       CDRGrad(1:3,index) = TMPGrad(1:3) * ( cavtmp + distmp + reptmp )
       first = last + 1
    END DO GradAtomLoop
  END SUBROUTINE AnalytMolNonelecGradient

  function Ainv_ok(A,Ainv,low_memory,verbose)
    use errormod

    REAL(SP),pointer :: A(:,:), Ainv(:,:)
    LOGICAL(LGD), INTENT(IN) :: low_memory,verbose
    LOGICAL(LGD) :: Ainv_ok
    INTEGER(I4B) :: nse

    Ainv_ok = .TRUE.
    if ( (.not. associated(A) ) ) then
       if (.not. associated(Ainv)) then
          call error("I have to compute Ainv but I don't have A!")
       else
          ! I don't have A, but I do have Ainv, so Ainv better be right
          return
       end if
    end if
    
    ! If we get here, we have A, and maybe Ainv.
    nse = assert_eq(size(A,1),size(A,2),'A is not square!')
    if ( associated(Ainv) ) then
       ! We might be able to use it
       if ( (size(Ainv,1) /= nse) .or. &
            (size(Ainv,2) /= nse) ) then
          ! AINV is not the right size/shape
          IF (verbose) write(6,*)&
               "WARNING: AINV had incorrect size/shape.  I'll fix that."
          deallocate(Ainv)
       else
          IF (verbose) write(6,*)'Ainv appears valid'
          return
       end if
    endif
    
    ! If we get here then we dont have Ainv
    Ainv_ok = .FALSE.
    if (low_memory) then
       Ainv => A
       A => null()
    else
       ! We need to allocate Ainv
       allocate(Ainv(1:nse,1:nse))
       Ainv = A
    end if
  end function Ainv_ok

  subroutine cpt_Ainv_switch(gammaswitch,sescale,A,Ainv,low_memory,verbose)
    ! Computes Ainv
    use UtilitiesMod, only : assert_eq
    use cosmo_util

    REAL(SP),pointer :: A(:,:), Ainv(:,:)
    LOGICAL(LGD), INTENT(IN) :: low_memory, verbose
    REAL(SP), intent(in) :: gammaswitch,sescale(:)
    integer(i4b) :: nse, i
    REAL(SP),allocatable :: said(:)

    IF(Ainv_ok(A,Ainv,low_memory,verbose))return
    nse=SIZE(Ainv,1)
    if ( gammaswitch > 0.0_sp ) then
       ! make sure sescale is there
       nse = assert_eq(nse,size(sescale),&
            'sescale vector is messed up in cpt_Ainv')

       ! handle switching
       allocate(said(1:nse))
       ! said will have the diagonal elements of A and Ainv will be A_off
       forall (i = 1:nse)
          said(i) =  Ainv(i,i)
          Ainv(i,i) = 0.0_sp
       end forall

       ! said is the square root of the scaled A inverse diagonal elements
       IF(org_Sik_in_A0_diag)THEN
          said = sqrt(sescale/said)
       ELSE
          said = sqrt(sqrt(sescale)/said)
       ENDIF
       ! Construct the matrix we need to invert (the middle of eqn. 81)
       IF (verbose) write(6,*)'Constructing middle of A...'
       Ainv = diagadd(diagmul(said,diagmul(Ainv, said)),1.0_sp)
       ! invert that puppy
       CALL do_syminv(nse,Ainv,verbose)
       ! multiply on both sides by said... (the rest of eqn 81)
       Ainv = diagmul(said,diagmul(Ainv,said))
       deallocate(said)
    else
       CALL do_syminv(nse,Ainv,verbose)
    end if

  end subroutine cpt_Ainv_switch

  subroutine cpt_Ainv_noswitch(A,Ainv,low_memory,verbose)
    ! Computes Ainv
    use cosmo_util

    REAL(SP),pointer :: A(:,:), Ainv(:,:)
    LOGICAL(LGD), INTENT(IN) :: low_memory, verbose
    integer(i4b) :: nse

    IF(Ainv_ok(A,Ainv,low_memory,verbose))return
    nse=SIZE(Ainv,1)
    CALL do_syminv(nse,Ainv,verbose)

  end subroutine cpt_Ainv_noswitch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_A_surface1(s1_crd,s1_zeta,A)
    USE CoulombMod, only : GaussianInteraction
    REAL(SP), INTENT(IN) :: s1_crd(:,:),s1_zeta(:)
    REAL(SP),pointer :: A(:,:)
    integer(i4b) :: nse, i, j
    real(sp) :: zetai, zetaj, zeta, r2, r

    nse = size(s1_zeta)

    ! allocate space for and construct A matrix
    if ( associated(A) ) DEALLOCATE(A)
    allocate(A(1:nse,1:nse))

    ! Construct A matrix (self-interaction matrix)
    ! Taking some advantage of the fact that it is symmetric
    !$OMP PARALLEL DO
    do i = 1, nse
       zetai = s1_zeta(i)
       do j = i, nse
          zetaj = s1_zeta(j)
          zeta = zetai*zetaj/sqrt(zetai*zetai+zetaj*zetaj)
          r2 = (s1_crd(1,i) - s1_crd(1,j))**2 + &
               (s1_crd(2,i) - s1_crd(2,j))**2 + &
               (s1_crd(3,i) - s1_crd(3,j))**2
          r = sqrt(r2)
          CALL GaussianInteraction(zeta,r,A(j,i))
          A(i,j) = A(j,i)
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine cpt_A_surface1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cpt_A_surface2(s1_crd,s1_zeta,s2_crd,s2_zeta,A_12)
    USE CoulombMod, only : GaussianInteraction
    REAL(SP), INTENT(IN) :: s1_crd(:,:),s1_zeta(:),s2_crd(:,:),s2_zeta(:)
    REAL(SP),pointer :: A_12(:,:)

    integer(i4b) :: nse_1,nse_2, i, j 
    real(sp) :: zeta1, zeta2, zeta, r2, r

    nse_1 = size(s1_zeta)
    nse_2 = size(s2_zeta)
    call estimate_real_size(nse_1,nse_2,THRESH=0.0_sp)

    if ( associated(A_12) )deallocate(A_12)
    allocate(A_12(1:nse_1,1:nse_2))

    write(6,*)'Computing A_12...'
    ! Construct A_12 matrix (interaction matrix between 2 surfaces)
    do i = 1, nse_2
       zeta2 = s2_zeta(i)
       do j = 1, nse_1
          zeta1 = s1_zeta(j)
          zeta = zeta1*zeta2/sqrt(zeta1*zeta1+zeta2*zeta2)
          r2 = (s2_crd(1,i) - s1_crd(1,j))**2 + &
               (s2_crd(2,i) - s1_crd(2,j))**2 + &
               (s2_crd(3,i) - s1_crd(3,j))**2
          r = sqrt(r2)
          CALL GaussianInteraction(zeta,r,A_12(j,i))
       end do
    end do

  end subroutine cpt_A_surface2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine switch_pot(pot, rho, sescale, zeta)
    use ConstantsMod
    use UtilitiesMod, only: assert_eq
    IMPLICIT NONE
    real(SP), intent(inout) :: pot(:)
    real(sp), intent(in) :: rho(:), sescale(:), zeta(:)
    real(SP) :: tmp
    integer(i4b) :: i, n

    n = assert_eq( (/ size(pot), size(sescale), size(zeta), size(rho) /), &
         'Size mismatch in switch_pot')

    IF(org_Sik_in_A0_diag) THEN
       do i=1,n
          pot(i) = pot(i) + 2.0_sp * (rho(i) * zeta(i) / (SQRT2 * SQRT_PI)) &
               * ((1.0_sp - sescale(i)) / sescale(i))
       end do
    ELSE
       do i=1,n
          tmp = SQRT(sescale(i))
          pot(i) = pot(i) + 2.0_sp * (rho(i) * zeta(i) / (SQRT2 * SQRT_PI)) &
               * ((1.0_sp - tmp) / tmp)
       end do       
    ENDIF
  end subroutine switch_pot


  subroutine read_surface_parms(newsurface,unit)
    USE Datatypes
    Use ConstantsMod
    IMPLICIT NONE
    type(surface_t), intent(inout) :: newsurface
    INTEGER(I4B),optional,intent(IN) :: unit
    TYPE(surface_t) :: surf

    INTEGER(I4B) :: iunit
    real(sp) :: res, probe, gs, pexp, rscale,r_thresh,AlphaShift, &
         g_cav, cav_rscale,dr_rscale,pressure
    integer(i4b) :: dt,pts_per_sphere,iCDR
    character(len=256) :: surface_file
    logical :: Verbose

    namelist /surface/ surface_file,rscale,r_thresh,dt,pts_per_sphere,res,probe,gs,pexp, &
         AlphaShift,verbose, iCDR !, g_cav,cav_rscale,dr_rscale,pressure

    surface_file = surf%surf_file
    rscale = surf%rscale
    r_thresh = surf%r_thresh
    dt = surf%discretization_type
    pts_per_sphere = surf%pts_per_sphere
    res = surf%resolution
    probe = surf%proberadius / ANGSTROM
    gs = surf%gammaswitch
    pexp = surf%pexp
    AlphaShift = surf%AlphaShift
    verbose = surf%verbose
    iCDR = surf%iCDR
    g_cav = surf%g_cav
    cav_rscale = surf%cav_rscale
    dr_rscale = surf%dr_rscale
    pressure = surf%pressure

    iunit = 5
    IF(present(unit))THEN
       IF(unit>=0)iunit=unit
    ENDIF

    read(iunit,nml=surface)

    newsurface%surf_file=surface_file
    newsurface%rscale = rscale
    newsurface%r_thresh = r_thresh
    newsurface%discretization_type = dt
    newsurface%pts_per_sphere = pts_per_sphere
    newsurface%resolution = res
    newsurface%proberadius = probe * ANGSTROM
    newsurface%gammaswitch = gs
    newsurface%pexp = pexp
    newsurface%AlphaShift = AlphaShift
    newsurface%verbose = verbose
    newsurface%iCDR = iCDR
    newsurface%g_cav = g_cav
    newsurface%cav_rscale = cav_rscale
    newsurface%dr_rscale = dr_rscale
    newsurface%pressure = pressure
    
  end subroutine read_surface_parms

  subroutine print_surface_parms(surface,unit)
    USE Datatypes
    Use ConstantsMod
    IMPLICIT NONE
    type(surface_t), intent(in) :: surface
    INTEGER(I4B),optional,intent(in) :: unit
    INTEGER(I4B) :: iunit

    iunit = 6
    IF(present(unit))THEN
       IF(unit>=0)iunit=unit
    ENDIF

    write(iunit,*)'surface_file= ',trim(surface%surf_file)
    if (  surface%discretization_type == 1 ) then
       write(iunit,10)' resolution=', surface%Resolution
    else
       write(iunit,13)' points per sphere=', surface%pts_per_sphere
    end if
    write(iunit,11)' probe radius(Ang)=', surface%ProbeRadius / ANGSTROM
    write(iunit,10)' radii scale=', surface%rscale
    write(iunit,10)' radii threshold=', surface%r_thresh
    write(iunit,10)' gamma switch=', surface%GammaSwitch
    write(iunit,10)' pexp=', surface%pexp
    if(surface%AlphaShift <0.0_SP .OR. surface%AlphaShift >1.0_SP)THEN
       WRITE(iunit,*)'AlphaShift=        Variable'
    ELSE
       WRITE(iunit,10)'AlphaShift=',surface%AlphaShift
    ENDIF

   IF(surface%iCDR == 1) THEN
      WRITE(iunit,*)'Cavitation Dispersion Repulsion gradients with g_cav = ', surface%g_cav
      WRITE(iunit,*)'                  cav_rscale = ', surface%cav_rscale
      WRITE(iunit,*)'                   dr_rscale = ', surface%dr_rscale
   ELSE IF(surface%iCDR == 2) THEN
      WRITE(iunit,*)'Nonelectrostatics calculated using atomic surface tensions'
   ELSE
      WRITE(iunit,*)'No Cavitation Dispersion Repulsion Requested'
   ENDIF


10  format(A,T24,F5.2)
11  format(A,T24,F5.2)
13  format(A,T24,I5)
  end subroutine print_surface_parms




END MODULE SurfaceMod

