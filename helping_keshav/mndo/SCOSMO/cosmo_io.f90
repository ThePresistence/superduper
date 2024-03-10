Module cosmo_io
  use Datatypes
  use constantsMod
  use cosmomod
  use surfaceMod
  USE MinimizeMod
  implicit none
private
  public ::read_cosmo_parms, &
           print_cosmo_parms

CONTAINS 

  Subroutine read_cosmo_parms(cosmo_data,infile,energyf,RadType,GaussType,unit)
    type(cosmo_t), intent(inout) :: cosmo_data
    CHARACTER(LEN=256),optional,intent(OUT)::infile,energyf
    ! Two extensions for mndo... If they are useful in regular cosmo,
    ! they should be rolled into cosmo_t
    INTEGER(I4B),optional,intent(INOUT):: radType
    LOGICAL(LGD),optional,intent(INOUT):: GaussType
    INTEGER(I4B),optional,intent(IN) :: unit
!!! The cosmo namelist 
!!! (to be replaced by a more user friendly interface in the future) 
!!! 
!!! filename         The name of the pdb file containing your molecule with
!!!                       radii and charges in the occupation number and Bvalue
!!!                       columns respectively.
!!! 
!!! surface_file     The root name of the surface data files.  If null, nothing
!!!                       happens.  If it is not null, you will get a .stat,
!!!                       .xyz, and .pdb file.  The .stat contains surface 
!!!                       construction stats.  The xyz and pdb are for
!!!                       visualizing the surface in gOpenMol.
!!!
!!! res              The resolution in points per square bohr.  (1.5)
!!!
!!! probe            Probe radius in Angstroms (1.4)
!!!
!!! gs               Gamma switch (0.0)
!!!
!!! pexp             Correction factor for overswitching (1.0)
!!!
!!! e1               epsilon inside (1.0)
!!!
!!! e2               epsilon outside (80.0)
!!!
!!! dt               Discretization type (0)
!!!                       0 -> constant pts/sphere
!!!                       1 -> constant resolution
!!!
!!! pts_per_sphere   Points per sphere (110)
!!!
!!! rscale           Scale radii by this factor (1.0)
!!!
!!! r_thresh         only keep particles with radii greater than r_thresh (0.0)
!!!
!!! l                Constrain COSMO to this multipole on the surface (-1)
!!!
!!! variational      variational constrained COSMO (Appendix B)
!!!
!!! low_memory       discard A matrix and other "unnecessary" data to save
!!!                      memory (false)
!!!
!!! restart          Is this a restart?
!!!
!!! tol_for_minimize tolerance for direct energy minimization (1.0E-5)
!!!
!!! itmax            maximum number of direct energy 
!!!                       minimization iterations (100)
!!! 
!!! itol             tolerance mode (1)
!!!
!!! ipreconditioner   Choose a ipreconditioner (Details, Jana?) (1)
!!!
!!! imultipole        use recursive fast imultipole gets 1, else 0 (0)
!!!
!!! WS               Well separatedness parameter for LS (2.0)
!!!
!!! L_rb             maximum l in recursive bisection-fast multipole (8)
!!! 
!!! IfAmatrix        Compute A matrix in Minimization (false)
!!! 
!!! iCosmo           inversion COSMO or minimization (0)
!!!                       0 = invert
!!!                       1 = minimization
!!!
!!! energy_file      Write energy and gradients to this file (energy)
!!!
!!! igrad            Numerical (2), analytic (1) gradients, both (3), 
!!!                     or none (default) (0)
!!!
!!! iRadii           Type of COSMO radii (0)
!!!                  -1 = Emsley; 0 = default; 1 = input; 
!!!                   2 = variable radii: R = R0+a*Q+b*Q^2
!!!                   3 = variable radii: 1/R = 1/R0 +a*Q+b*Q^2
!!!                   4 = variable radii: R = R0 * exp[-(a*Q+b*Q^2)]
!!!
!!! FDstep           finite difference step size for numerical gradients (1.0E-6)
!!!
!!! AlphaShift       Alpha parameter used in Eq. 62 and 63
!!!                   < 0.0_SP  or > 1.0_SP use the value as determined in Eq.E-6
!!!                   >= 0.0_SP and <= 1.0_SP use input value.
!!!                   0.0_SP (default)
!!!
!!! iCDR             Calculate Nonelectrostatics: cavitation dispersion repulsion terms 
!!!                   (0) = no nonelectrostatics (default)
!!!                   (1) = molecular noneletrostatics
!!!                   (2) = atomic nonelectrostatics
!!!                   (3) = constrained atom nonelectrostatics (better ask ATM)
!!!
!!! g_cav           Surface tension parameter for molecular cavitation in kcal/mol/Ang^2
!!!                   = 0.10368 (default) = 72 dyn/cm
!!!
!!! cav_rscale      Scale surface area by cav_rscale^2 for cavitation
!!! 
!!! dr_rscale       Scale surface area by dr_rscale^2 for dispersion/repulsion
!!!


  ! This will contain the defaults we use for cosmo/surface
  TYPE(cosmo_t) :: cd
  character(len=256) :: filename
  CHARACTER(LEN=256) :: energy_file, g_file
  INTEGER(I4B):: iunit
  real(sp) :: e1, e2, tol_for_minimize,WS, FDstep
  integer(i4b) :: l,itmax,itol,ipreconditioner,imultipole,iCosmo, &
       L_rb,igrad,iRadii
  logical :: variational,IfAmatrix,low_memory,Verbose,lGaus

  namelist /cosmo/ filename, e1, e2, l, variational,  &
       tol_for_minimize, itmax, itol,ipreconditioner,imultipole, WS, L_rb, &
       IfAmatrix,iCosmo,energy_file,igrad,FDstep, low_memory, &
       g_file, verbose, iRadii,lgaus

  filename=''
  energy_file='energy'
  g_file=''
!!! **** Initialization

  e1 = cd%e1
  e2 = cd%e2
  l = cd%l
  variational = cd%variational
  low_memory = cd%low_memory
  verbose = cd%verbose
  iCosmo = cd%icosmo
  iGrad = cd%iGrad

  tol_for_minimize = cd%minpara%tol 
  WS = cd%minpara%WS
  itmax = cd%minpara%itmax
  itol = cd%minpara%itol
  ipreconditioner = cd%minpara%precond 
  imultipole =  cd%minpara%multipole
  L_rb = cd%minpara%L 
  IfAmatrix = cd%minpara%IfAmatrix 
  FDstep = cd%minpara%FDStep

  IF(PRESENT(GaussType)) lgaus = GaussType
  iunit=5
  IF(present(unit))THEN
     IF(unit>=0)iunit=unit
  ENDIF

read(iunit,nml=cosmo)

!!! **** Initialization

  cosmo_data%e1 = e1
  cosmo_data%e2 = e2
  cosmo_data%l = l
  cosmo_data%variational = variational
  cosmo_data%low_memory = low_memory
  cosmo_data%verbose = verbose
  cosmo_data%icosmo = iCosmo
  cosmo_data%iGrad = iGrad

  cosmo_data%minpara%tol = tol_for_minimize 
  cosmo_data%minpara%WS = WS
  cosmo_data%minpara%itmax = itmax
  cosmo_data%minpara%itol = itol
  cosmo_data%minpara%precond = ipreconditioner 
  cosmo_data%minpara%multipole = imultipole
  cosmo_data%minpara%L = L_rb
  cosmo_data%minpara%IfAmatrix = IfAmatrix 
  cosmo_data%minpara%FDStep=FDstep
  IF(present(infile))infile=filename
  IF(present(energyf))energyf=energy_file
  IF(present(RadType))radtype=iRadii
  IF(present(GaussType))GaussType=lgaus
  end Subroutine read_cosmo_parms


  subroutine print_cosmo_parms(cosmo_data,unit,infile,energyf,radtype,gausstype)
    type(cosmo_t), intent(in) :: cosmo_data
    INTEGER(I4B),optional,intent(in) :: unit
    CHARACTER(LEN=256),optional,intent(in) :: infile, energyf
    INTEGER(I4B),optional,intent(IN):: radType
    LOGICAL(LGD),optional,intent(IN):: GaussType

    INTEGER(I4B) :: iunit

    iunit = 6
    IF(present(unit))THEN
       IF(unit>=0)iunit=unit
    ENDIF

    write(iunit,*)'****SMOOTH COSMO INPUT PARAMETERS****'
    IF(present(infile)) write(iunit,*)'filename= ',trim(infile)
    IF(present(energyf)) write(iunit,*)'energy and gradient_file= ',trim(energyf)
    WRITE(iunit,12)'Type of COSMO= ', cosmo_data%iCosmo
    ! Specific extensions for mndo... If they become popular
    ! they should be rolled into cosmo_t
    IF(present(radType))WRITE(iunit,*)'Type of radii= ', radtype
    IF(present(GaussType))WRITE(iunit,*)'Gaussians used for density-srf interaction=', GaussType

    if (cosmo_data%e1 > 1.0_sp) then
       write(iunit,10)' e1=', cosmo_data%e1
    end if
    write(iunit,10)' e2=', cosmo_data%e2

    if (  cosmo_data%l > -1 ) then
       write(iunit,13)' constrained to l=', cosmo_data%l
       if (  cosmo_data%variational ) then
          write(iunit,*)'using variational method'
       else
          write(iunit,*)'using non-variational method'
       end if
    else
       write(iunit,*)'unconstrained'
    end if
    if ( cosmo_data%low_memory) then
       write(iunit,*)'using low memory algorithm'
    end if

    IF (cosmo_data%igrad /= 0) THEN
       WRITE(iunit,*)'Gradient calculation'
       IF (cosmo_data%igrad ==1)THEN
            WRITE(iunit,*)'compute analytic gradients'
       ELSEIF (cosmo_data%igrad ==2) THEN
          WRITE(iunit,*) 'compute numerical gradients'
          WRITE(iunit,'(A,E8.2)') 'NumGrad Finite difference step size=', cosmo_data%minpara%FDstep
       ELSEIF (cosmo_data%igrad ==3) THEN
          WRITE(iunit,*)'compute both analytic and numerical gradients'
          WRITE(iunit,'(A,E8.2)') 'NumGrad Finite difference step size=', cosmo_data%minpara%FDstep
       ENDIF
   END IF

    WRITE(iunit,*)'Energy minimization input parameters:'
    WRITE(iunit,*)'Use A matrix=', cosmo_data%minpara%IfAmatrix
    WRITE(iunit,12)'Option for tolerance=', cosmo_data%minpara%itol
    WRITE(iunit,14)'Tolerance=',cosmo_data%minpara%tol
    WRITE(iunit,12)'Max iteration number=',cosmo_data%minpara%itmax
    WRITE(iunit,12)'Type of preconditioner=',cosmo_data%minpara%precond
    IF (cosmo_data%minpara%multipole == 0) THEN
       WRITE(iunit,*) 'Use direct Coulomb for electrostatics'
    ELSE
       WRITE(iunit,*)'Use recursive bisection for electrostatics'
       WRITE(iunit,'(A,F5.2,I5)')'WS,L_rb=', cosmo_data%minpara%WS, cosmo_data%minpara%L
    END IF


10  format(A,T24,F5.2)
11  format(A,T24,F5.2)
12  format(A,T24,I5)
13  format(A,T24,I5)
14  format(A,T24,E8.2)
  end subroutine print_cosmo_parms

end Module cosmo_io
