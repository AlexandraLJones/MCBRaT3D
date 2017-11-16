! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 31 $, $Date: 2010-06-14 23:12:00 +0100 (Mon, 14 Jun 2010) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/opticalProperties.f95 $
! --------------------------------------------
module opticalProperties
  ! Provide a representation of the three dimensional optical properties 
  !   of the atmosphere with each component (e.g. cloud, aerosol) 
  !   represented separately. 
  ! There are two objects: the domain in which the optical properties are specified, 
  !   and a set of one or more optical components. 
  use characterUtils
  use ErrorMessages
  use scatteringPhaseFunctions
  use netcdf
  use inversePhaseFunctions
  use numericUtilities

  implicit none
  private
  
  integer, parameter :: maxNameLength = 256
  real,    parameter :: Pi = 3.14159265358979312 
  real(8),    parameter :: light_spd = 2.99792458E8 ![m/s]
  real(8), parameter :: Na = 6.02214129E23 ![mol^-1]
  real(8), parameter :: Rstar = 8.3144621 ![J K^-1 mol^-1]
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  ! Matrix and it's functions are also used by the integrator 
  !   One of the components of the public type (Domain) is of type matrix.

  type matrix
    integer                        :: numX = 0, numY  = 0
    real, dimension(:, :), pointer :: values => null()
  end type matrix


  type opticalComponent
    ! Each component has a name. The optical properties may be defined for a subset of the 
    !   z levels in the domain, and may be uniform or variable  horiontally. 
    !   the values of the optical properties (extinction, single scattering 
    !   albedo, phase function) of each component at each point. 
    ! The vectors defining the domain's cell boundaries should be one element
    !   larger than the corresponding array dimension. TIn other words, a field (extinction, say)
    !   defined everywhere is an an array of size nx, ny, nz, while the xPosition, xPosition, 
    !   and zPostion vectors are of length nx+1, ny+1, and nz+1.
    !
    private
    character (len = maxNameLength)      :: name       = ""
    integer                              :: zLevelBase = 0
    logical                              :: horizontallyUniform = .false.
    real(8),    dimension(:, :, :), pointer :: extinction             => null()
    real(8),    dimension(:, :, :), pointer :: singleScatteringAlbedo => null()
    integer, dimension(:, :, :), pointer :: phaseFunctionIndex     => null()
    type(phaseFunctionTable)             :: table
  end type opticalComponent

  type commonDomain
      ! The domain contains the x, y, and z cell boundaries and temperatures, these things are common to all domains in the simulation regardless of wavelength
    
    real(8), allocatable          :: xPosition(:) 
    real(8), allocatable          :: yPosition(:) 
    real(8), allocatable          :: zPosition(:) 
    real(8), allocatable          :: temps(:,:,:) 
    real(8), allocatable          :: rho(:,:,:)
    real(8), allocatable          :: numConc(:,:,:) !  total molecular numberconcentration(x, y, z)
    real(8), allocatable          :: massConc(:,:,:,:)! component mass concentration (component,x,y,z)
    real(8), allocatable          :: Reff(:,:,:,:)! component effective radius (component,x,y,z)

  end type commonDomain

  type domain
!
    private
    real(8), pointer, dimension(:)          :: xPosition => null()
    real(8), pointer, dimension(:)          :: yPosition => null()
    real(8), pointer, dimension(:)          :: zPosition => null()
    real(8), pointer, dimension(:,:,:)      :: temps => null()   
    real(8)				    :: lambda
    integer				    :: lambdaI, nlambda
    real(8)				    :: surfaceAlbedo
    logical                              :: xyRegularlySpaced = .false., zRegularlySpaced = .false. 
    type(opticalComponent), &
                   dimension(:), pointer :: components  => null()
    real(8),    dimension(:, :, :),    pointer :: totalExt           => null()
    real(8),    dimension(:, :, :, :), pointer :: cumulativeExt      => null()
    real(8),    dimension(:, :, :, :), pointer :: ssa                => null()
    integer, dimension(:, :, :, :), pointer :: phaseFunctionIndex => null()
    !
    ! We store the original forward phase function tables even though calculations
    !   inside the module use the matrix representations below. The originals don't
    !   take much rooms (a few Mb at most) and this allows us to recompute them to
    !   arbitrary accuracy at any time.
    !
    type(phaseFunctionTable), &
             dimension(:),    pointer :: forwardTables => null()

    ! We store tabulated phase function and inverse (cumulative) phase functions in
    !   two-D arrays, but these arrays can be different sizes for different components.
    !   We define a derived type to represent each 2D array and make a vector of these
    !   (one matrix for each component).
    !
    type(matrix), &
          dimension(:),       pointer :: tabulatedPhaseFunctions => null()
    type(matrix), &
          dimension(:),       pointer :: tabulatedOrigPhaseFunctions => null()
    type(matrix), &
          dimension(:),       pointer :: inversePhaseFunctions => null()
    ! -------------------------------------------------------------------------=    
  end type domain
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  
  interface addOpticalComponent
    module procedure addOpticalComponent3D, addOpticalComponent1D
  end interface ! addOpticalComponent

  interface replaceOpticalComponent
    module procedure replaceOpticalComponent3D, replaceOpticalComponent1D
  end interface ! replaceOpticalComponent

 interface new_Domain
    module procedure new_DomainMono, new_DomainBB
  end interface
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------

  ! The types...
  public :: domain, matrix, commonDomain 
  ! ... and these procedures
  public :: new_Domain, getInfo_Domain, write_Domain, read_Domain, finalize_Domain, &
            addOpticalComponent, deleteOpticalComponent, replaceOpticalComponent,   &
            getOpticalPropertiesByComponent, finalize_Matrix, new_Matrix,           &
            accumulateExtinctionAlongPath, tabulateInversePhaseFunctions,           &
            tabulateForwardPhaseFunctions, read_Common, read_SSPTable, calc_RayleighScattering !,  getAverageOpticalProperties
            

contains

   subroutine read_SSPTable(ncFileIds, lambdaIndex, commonD, thisDomain, setup, calcRayl, status)
!    character(len = *), intent(in   ) :: fileName
    integer, dimension(4), intent(in)    :: ncFileIds
    integer,            intent(in)    :: lambdaIndex
    type(commonDomain), intent(in)    :: commonD
    type(Domain), intent(out)         :: thisDomain
    logical, intent(in)               :: setup, calcRayl
    type(ErrorMessage), intent(inout) :: status

    integer, dimension(30)            :: ncStatus
    type(phaseFunctionTable)          :: table
    logical                           :: horizontallyUniform, fillsVerticalDomain
    character(len = maxNameLength)    :: name, extType
    integer                           :: nLambda, nComponents, zLevelBase, i, nZGrid, nXEdges, nYEdges, nZEdges, j, length, n
    real(8)                           :: lambda, albedo, freq, f
    integer                           :: nDims, ncVarID, ncDimId, zGridDimId, err, ncFileId, dimId, &
					nReff, il, comp, gasComp, ix, iy, iz
    integer, dimension(3)             :: dimIds
    real(8), allocatable              :: extinction(:,:,:), singleScatteringAlbedo(:,:,:)
    real(8), allocatable              :: xsec(:), extinctionT(:), singleScatteringAlbedoT(:)
	real, allocatable                 :: key(:)
    integer, allocatable              :: phaseFunctionIndex(:,:,:)
	real, dimension(2)                                  :: LG
   type(phaseFunction), dimension(1)                   :: phaseFunc

!    ierr = nf_set_log_level(3)
!    PRINT *, trim(nf90_strerror(err))

    nXEdges = size(commonD%xPosition)
    nYEdges = size(commonD%yPosition)
    nZEdges = size(commonD%zPosition)

    nZGrid = nZEdges-1

 n=1 ! SSP file counter
 comp = 1 ! total number of components counter
 gasComp = 0 ! total number of gaseous components
 DO    ! loop over SSPTable files
    ncFileId = ncFileIds(n)
    ncStatus(:) = nf90_NoErr
    ncStatus(:) = nf90_NoErr
!    err = nf90_open(trim(fileName), nf90_NoWrite, ncFileID)
!    if(err /= nf90_NoErr) then
!      call setStateToFailure(status, "read_SSPTable: Can't open file " // trim(fileName))
!      PRINT *, "read_SSPTable: error opening file ", err
!      PRINT *, trim(nf90_strerror(err))
!    end if

    ncStatus( 1) = nf90_inq_dimid(ncFileId, "f_grid_nelem", ncDimId)
    ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nLambda)
    ncStatus( 3) = nf90_inq_varid(ncFileId, "f_grid", ncVarId)
    ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, freq, start = (/lambdaIndex/))
    lambda = (light_spd * (10**6))/freq  ![microns]
    ncStatus( 5) = nf90_inq_varid(ncFileId, "surfaceAlbedo", ncVarId)
    ncStatus( 6) = nf90_get_var(ncFileID, ncVarId, albedo, start = (/lambdaIndex/))

    if(any(ncStatus(:) .ne. nf90_NoErr)) then
        call setStateToFailure(status, "read_SSPTable: doesn't look an optical properties file.")
        PRINT *, "read_SSPTable: ncStatus properties file errors ", ncStatus(1:8)
    end if
    if(n .eq. 1) thisDomain = new_Domain(commonD, lambda, lambdaIndex, nlambda, albedo, status)

    ncStatus( 9) = nf90_get_att(ncFileID, nf90_Global, "numberOfComponents", nComponents)
    do i = 1, nComponents ! loop over components in file n
        ncStatus(10) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "Name", name)
!PRINT *, "readSSP_Table: name= ", trim(makePrefix(i)), "Name"
        ncStatus(11) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "zLevelBase", zLevelBase)
!PRINT *, "readSSP_Table: name= ", trim(makePrefix(i)), "zLevelBase"
	ncStatus(12) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "extType", extType)
		
        if(extType .eq. "absXsec")then! Read in the profile of absorption cross section; we assume that there is no associated scattering for this component
		gasComp = gasComp + 1
		allocate(xsec(1:nZGrid))
		ncStatus(13) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "xsec", ncVarId)
		ncStatus(14) = nf90_get_var(ncFileId, ncVarId, xsec(:), start = (/1,lambdaIndex/), count = (/nZGrid,1/))
		allocate(extinction(1,1,nZGrid),singleScatteringAlbedo(1, 1, nZGrid), phaseFunctionIndex(1, 1, nZGrid))
           	extinction(1,1,:) = xsec * commonD%numConc(1,1,:) * 1000.0 ! xsec should have units of m^2 per molecule and the number concentrations should be in molecules per meter cubed, which means the resulting volume extinction coefficient is in m^-1, however, physical distances are in km and the units need to cancel with volume extinction coefficient, so that means converting to units of km^-1, thus the factor of 1000 multiplication
            	deallocate(xsec)
		singleScatteringAlbedo = 0.0_8
		phaseFunctionIndex = 1
		LG = (/0.0, 0.0/)
		phaseFunc(1) = new_PhaseFunction(LG,status=status)
		table = new_PhaseFunctionTable(phaseFunc(1:1),key=(/0.0/),tableDescription="Molecular Absorption", status=status)
		call finalize_PhaseFunction(phaseFunc(1))
		if(any(ncStatus(:) .ne. nf90_NoErr)) then
			PRINT *, "read_SSPTable: Error reading scalar fields from file ", ncStatus(9:), "lambdaIndex= ", lambdaIndex
			call setStateToFailure(status, "read_SSPTable: Error reading scalar fields from file")
		end if
	elseif(extType .eq. "volExt")then ! Read in volume extinction
		ncStatus(13) = nf90_inq_dimid(ncFileId, trim(makePrefix(i)) // "phaseFunctionNumber", dimId)
		ncStatus(14) = nf90_Inquire_Dimension(ncFileId, dimId, len = nReff)
		allocate(extinctionT(nReff),singleScatteringAlbedoT(nReff), key(nReff))
		ncStatus(15) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "ExtinctionT", ncVarId)
		ncStatus(16) = nf90_get_var(ncFileId, ncVarId, extinctionT(:), start = (/1,lambdaIndex/), count = (/nReff,1/))
		ncStatus(17) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedoT", ncVarId)
		ncStatus(18) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedoT(:), start = (/1,lambdaIndex/), count = (/nReff,1/))
		ncStatus(19) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "phaseFunctionKeyT", ncVarId)
		ncStatus(20) = nf90_get_var(ncFileId, ncVarId, key(:))
		if(any(ncStatus(:) .ne. nf90_NoErr)) then
			PRINT *, "read_SSPTable: Error reading scalar fields from file ", ncStatus(9:), "lambdaIndex= ", lambdaIndex
			call setStateToFailure(status, "read_SSPTable: Error reading scalar fields from file")
		end if
		allocate(phaseFunctionIndex(nXEdges-1, nYEdges-1, nZGrid), &
		extinction(nXEdges-1, nYEdges-1,nZGrid),singleScatteringAlbedo(nXEdges-1, nYEdges-1, nZGrid))
		phaseFunctionIndex = 1
                extinction = 0.0_8
                singleScatteringAlbedo = 0.0_8
		if(.not. stateIsFailure(status) .and. .not. setup) then ! We don't need the time consuming phase function calculations for the setup step
			call read_PhaseFunctionTable(fileId = ncFileId, spectIndex = lambdaIndex, table = table,  &
                                       prefix = "Component" // trim(IntToChar(i)) // "_", status = status)
		else
			LG = (/0.0, 0.0/)
                	phaseFunc(1) = new_PhaseFunction(LG,status=status)
			table = new_PhaseFunctionTable(phaseFunc(1:1),key=(/0.0/),tableDescription="dummy table", status=status)
                	call finalize_PhaseFunction(phaseFunc(1))
		endif
		! Fill Domain vars according to physical properties and phase function key nearest neighbor	
		do iz = 1, nZGrid
			do iy = 1, nYEdges-1
				do ix = 1, nXEdges-1
					if(commonD%MassConc(comp-gasComp,ix,iy,iz) .gt. 0.0_8 .and. &
					commonD%Reff(comp-gasComp,ix,iy,iz) .lt. MAXVAL(key) .and. &
					commonD%Reff(comp-gasComp,ix,iy,iz) .ge. MINVAL(key))then
						! Binary search to find effective radius entry in table
						il = findIndex(commonD%Reff(comp-gasComp,ix,iy,iz), REAL(key(:),8))

						! Interpolate optical properties linearly in Reff for this component
						f = (commonD%Reff(comp-gasComp,ix,iy,iz)-key(il)) / (key(il+1)-key(il))
						extinction(ix,iy,iz) = commonD%MassConc(comp-gasComp,ix,iy,iz) * &
                                    		((1-f)*extinctionT(il) + f*extinctionT(il+1))
						singleScatteringAlbedo(ix,iy,iz) = &
						     (1-f)*singleScatteringAlbedoT(il) + &
						     f*singleScatteringAlbedoT(il+1)
						if(.not. setup)then  ! the presense of raylExt indicates that we are not in the setup step and therefore information about the phase function is needed
						   ! Chose the closest phase function
						   if (f < 0.5) then
							phaseFunctionIndex(ix,iy,iz) = il
						   else
							phaseFunctionIndex(ix,iy,iz) = il+1
						   endif
						endif
					else if(commonD%MassConc(comp-gasComp,ix,iy,iz) .gt. 0.0_8)then
						call setStateToFailure(status, "read_SSPTable: Effective radius outside of table range")
					end if
				end do
			end do
		end do	
		DEALLOCATE(extinctionT, singleScatteringAlbedoT, key)
			
			
	else ! unrecognizable format
		call setStateToFailure(status, "read_SSPTable: unrecognizable extType")
	end if
        
             
        !
        ! Add the new component to the domain.
        !
        if(.not. stateIsFailure(status)) then

!PRINT *, "read_Domain: extinction size ", size(extinction, 1), size(extinction, 2), size(extinction, 3)
		call addOpticalComponent(thisDomain, name, extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, table, zLevelBase = zLevelBase,   &
                                   status = status)
        else
		call setStateToFailure(status, "read_SSPTable: Error reading phase function table.")
        end if
        deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex)
	CALL finalize_PhaseFunctionTable(table)
	comp = comp+1
      end do
      n=n+1
      if(ncFileIds(n).eq.-9999)EXIT
 END DO


 ! do we need to add rayleigh scattering?
 if(calcRayl .and. .not. setup)then
    zLevelBase=1
    name = "Rayleigh Scattering"
    allocate(extinction(1,1,nZGrid), singleScatteringAlbedo(1,1,nZGrid), phaseFunctionIndex(1,1,nZGrid))
    CALL calc_RayleighScattering(lambda,commonD%rho(1,1,:),commonD%numConc(1,1,:), &
                                  ext=extinction(1,1,:), ssa=singleScatteringAlbedo(1,1,:), &
                                  phaseInd=phaseFunctionIndex(1,1,:),table=table,status=status)
    if(.not. stateIsFailure(status)) then
!PRINT *, "read_Domain: extinction size ", size(extinction, 1), size(extinction, 2), size(extinction, 3)
          call addOpticalComponent(thisDomain, name, extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, table, zLevelBase = zLevelBase,   &
                                   status = status)
    else
          call setStateToFailure(status, "read_SSPTable: Error calculating rayleigh scattering.")
    end if
    deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex)
	CALL finalize_PhaseFunctionTable(table)
 end if

 if(.not. stateIsFailure(status))  call getOpticalPropertiesByComponent(thisDomain, status)
 if(.not. stateIsFailure(status)) call setStateToSuccess(status)
   end subroutine read_SSPTable

   subroutine read_Common(filename, commonD, status)
    character(len = *), intent(in   ) :: fileName
    type(commonDomain), intent(out)   :: commonD
    type(ErrorMessage), intent(inout) :: status

    integer, dimension(16)             :: ncStatus
    integer                           :: ncFileID, ncDimID, zGridDimID, nXEdges, nYEdges, &
					 nZEdges, nZGrid, ncVarID, nComponents, i, nDims
    real(8), allocatable              :: prssr(:,:,:)

    ncStatus(:) = nf90_NoErr
    if(nf90_open(trim(fileName), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      call setStateToFailure(status, "read_Common: Can't open file " // trim(fileName))
    end if

    if(.not. stateIsFailure(status)) then
      ncStatus( 1) = nf90_inq_dimid(ncFileId, "x-edges", ncDimId)
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nXEdges)
      ncStatus( 3) = nf90_inq_dimid(ncFileId, "y-edges", ncDimId)
      ncStatus( 4) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nYEdges)
      ncStatus( 5) = nf90_inq_dimid(ncFileId, "z-edges", ncDimId)
      ncStatus( 6) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nZEdges)
      ncStatus( 7) = nf90_inq_dimid(ncFileId, "z-grid", zGridDimId)
      ncStatus( 8) = nf90_Inquire_Dimension(ncFileId, zGridDimId, len = nZGrid)
      if(any(ncStatus(:) /= nf90_NoErr)) then
        call setStateToFailure(status, "read_Common: " // trim(fileName) // &
                               " problem reading dimensions.")
        PRINT *, "read_Common: ncStatus dimension errors ", ncStatus(1:9) 
      end if
      allocate(commonD%xPosition(nXEdges), commonD%yPosition(nYEdges), commonD%zPosition(nZEdges), &
		 commonD%temps(nXEdges-1,nYEdges-1,nZEdges-1))
      ncStatus( 1) = nf90_inq_varid(ncFileId, "x-edges", ncVarId)
      ncStatus( 2) = nf90_get_var(ncFileId, ncVarId, commonD%xPosition)
      ncStatus( 3) = nf90_inq_varid(ncFileId, "y-edges", ncVarId)
      ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, commonD%yPosition)
      ncStatus( 5) = nf90_inq_varid(ncFileId, "z-edges", ncVarId)
      ncStatus( 6) = nf90_get_var(ncFileId, ncVarId, commonD%zPosition)
      ncStatus( 7) = nf90_inq_varid(ncFileId, "Temperatures", ncVarId)
      ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, commonD%temps)
!      if( COUNT(temps .le. 0.0_8) .gt. 0)PRINT *, 'readDomain: there are temps at or below 0.0 K'
!PRINT *, 'read_Common: status before Nc', ncStatus(:)
      if(any(ncStatus(:) /= nf90_NoErr)) then
        call setStateToFailure(status, "read_Common: " // trim(fileName) // &
                               " doesn't look an optical properties file.")
        PRINT *, "read_Common: ncStatus optical property file errors ", ncStatus(1:8)
      else
        ncStatus( 1) = nf90_inq_varid(ncFileId,"Pressures", ncVarID)

        if(ncStatus(1) .eq.  nf90_NoErr)then
	  allocate(prssr(1:nXEdges-1,1:nYEdges-1,1:nZEdges-1))
	  allocate(commonD%numConc(nXEdges-1,nYEdges-1,nZEdges-1))
	  ncStatus(2)= nf90_Inquire_Variable(ncFileId, ncVarId, ndims = nDims)
	  if (nDims .eq. 3)then
 	     ncStatus(3) = nf90_get_var(ncFileId, ncVarId, prssr(:,:,:))
 	  else if(nDims .eq. 1)then
	     ncStatus(3) = nf90_get_var(ncFileId, ncVarId, prssr(1,1,:))
	     prssr = spread(spread(prssr(1,1,:),1, nCopies=nXEdges-1), 2, nCopies=nYEdges-1)
	  else
	     call setStateToFailure(status, "read_Common: " // trim(fileName) // &
                               " strange number of dimensions for pressure")
	  end if  
!	  do i= 1, nComponents THIS BLOCK IS A PLACEHOLDER FOR BRINGING IN MIXING RATIOS OF SOME SORT IF WANTED/NEEDED TO PARTITION THE TOTAL NUMBER CONCENTRATION
!	     ncStatus(2*i) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "numberConcen", ncVarId)
!             ncStatus((2*i)+1) = nf90_get_var(ncFileId, ncVarId, commonD%numConc(1,1,:)) 
!	  end do     
!PRINT *, 'read_Common: status after Nc', ncStatus(:)
 	  commonD%numConc(:,:,:)=(prssr*100.0*Na)/(Rstar*commonD%temps) ! factor of 100 converts from hPa in domain to the Pa needed for proper unit cancelation
!PRINT *, "number concentration: ", commonD%numConc(1,1,:)
          if(any(ncStatus(:) /= nf90_NoErr)) then
          	call setStateToFailure(status, "read_Common: " // trim(fileName) // &
                               " problem reading PRESSURE.") 
	  else
	    ncStatus( 1) = nf90_inq_dimid(ncFileId, "nonGasComps", ncDimId)
      	    ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nComponents)
	    if (ncStatus(1) .eq. nf90_NoErr .and. nComponents .gt. 0) then
		allocate(commonD%massConc(nComponents,nXEdges-1,nYEdges-1,nZEdges-1), commonD%Reff(nComponents,nXEdges-1,nYEdges-1,nZEdges-1))
		ncStatus( 3) = nf90_inq_varid(ncFileId,"massConc", ncVarID)
		ncStatus( 4) = nf90_get_var(ncFileId, ncVarId,commonD%massConc)
		ncStatus( 5) = nf90_inq_varid(ncFileId,"Reff", ncVarID)
                ncStatus( 6) = nf90_get_var(ncFileId, ncVarId,commonD%Reff)
	    end if
	  end if
        end if

	ncStatus( 1) = nf90_inq_varid(ncFileId,"Density", ncVarID)
	if(ncStatus(1) .eq.  nf90_NoErr)then
	   ncStatus( 2) = nf90_Inquire_Variable(ncFileId, ncVarId, ndims = nDims)
	   allocate(commonD%rho(nXEdges-1,nYEdges-1,nZEdges-1))
	   if(ndims.eq.1)then
		ncStatus(3) = nf90_get_var(ncFileId, ncVarId,commonD%rho(1,1,:))
		commonD%rho=spread(spread(commonD%rho(1,1,:),1, nCopies=nXEdges-1), 2, nCopies=nYEdges-1)
	   else
		ncStatus(3) = nf90_get_var(ncFileId, ncVarId,commonD%rho)
	   end if
	end if

	if(any(ncStatus(:) /= nf90_NoErr))then
	   call setStateToFailure(status, "read_Common: " // trim(fileName) // &
                               " problem reading DENSITY, massConc, or Reff.")
	else
	   call setStateToSuccess(status)
	end if
      end if
    end if
   end subroutine
  !------------------------------------------------------------------------------------------
  ! Initialization: Routine to create new domains 
  !------------------------------------------------------------------------------------------
   function new_DomainBB(commonD, lambda, lambdaI, nlambda, albedo, status)
!  function new_Domain(xPosition, yPosition, zPosition, lambda, lambdaI, nlambda, albedo, temps, status)
!    real(8),  dimension(:), intent(in   ) :: xPosition, yPosition, zPosition
!    real(8),  dimension(:,:,:), intent(in) :: temps
    type(commonDomain), TARGET,  intent(in)       :: commonD
    real(8), intent(in)			 :: lambda, albedo
    integer, intent(in)			 :: lambdaI, nlambda
    type(ErrorMessage),    intent(inout) :: status
    type(domain)                         :: new_DomainBB
    
    ! Local variables
    integer :: numX, numY, numZ

      
    ! -------------------------
    ! Checks: always increasing, within limits; 
    !   other checks performed by addOpticalComponent
!    numX = size(xPosition); numY = size(yPosition); numZ = size(zPosition)
    numX = size(commonD%xPosition); numY = size(commonD%yPosition); numZ = size(commonD%zPosition)
    if(any(commonD%xPosition(2:) - commonD%xPosition(:numX-1) <= 0.) .or. &
       any(commonD%yPosition(2:) - commonD%yPosition(:numY-1) <= 0.) .or. &
       any(commonD%zPosition(2:) - commonD%zPosition(:numZ-1) <= 0.))     &
      call setStateToFailure(status, "new_Domain: Positions must be increasing, unique.")
     
    ! -------------------------
     if(.not. stateIsFailure(status)) then
!       allocate(new_Domain%xPosition(numX), new_Domain%yPosition(numY), new_Domain%zPosition(numZ), new_Domain%temps(numX-1,numY-1,numZ-1))
       new_DomainBB%xPosition(1:) => commonD%xPosition(:)
       new_DomainBB%yPosition(1:) => commonD%yPosition(:)
       new_DomainBB%zPosition(1:) => commonD%zPosition(:)
       new_DomainBB%temps(1:,1:,1:) => commonD%temps(:,:,:)       
       new_DomainBB%lambda       = lambda
       new_DomainBB%lambdaI      = lambdaI
       new_DomainBB%nlambda      = nlambda
!PRINT *, 'new_Domain: nlambdas match= ', new_Domain%nlambda, nlambda
       new_DomainBB%surfaceAlbedo = albedo

       ! Are the grids regularly spaced? Compare the distance between each pair 
       !   of array elements to the distance between the first two. 
       ! The default value is false. 
       !
       if(all(abs( (commonD%xPosition(2:) - commonD%xPosition(:numX-1)) -                                &
                   (commonD%xPosition(2)  - commonD%xPosition(1)) ) <= 2 * spacing(commonD%xPosition(2:))) .and. &
          all(abs( (commonD%yPosition(2:) - commonD%yPosition(:numY-1)) -                                &
                   (commonD%yPosition(2)  - commonD%yPosition(1)) ) <= 2 * spacing(commonD%yPosition(2:))))      &
         new_DomainBB%xyRegularlySpaced = .true.
       if(all(abs( (commonD%zPosition(2:) - commonD%zPosition(:numZ-1)) -                           &
                   (commonD%zPosition(2)  - commonD%zPosition(1)) ) <= 2 * spacing(commonD%zPosition(2:)))) &
         new_DomainBB%zRegularlySpaced = .true.
       call setStateToSuccess(status)
     end if

  end function new_DomainBB
  !------------------------------------------------------------------------------------------
  function new_DomainMono(xPosition, yPosition, zPosition, temps, status)
    real(8),  dimension(:), intent(in   ) :: xPosition, yPosition, zPosition
    real(8),  dimension(:,:,:), intent(in) :: temps
    type(ErrorMessage),    intent(inout) :: status
    type(domain)                         :: new_DomainMono

    ! Local variables
    integer :: numX, numY, numZ


    ! -------------------------
    ! Checks: always increasing, within limits;
    !   other checks performed by addOpticalComponent
    numX = size(xPosition); numY = size(yPosition); numZ = size(zPosition)
    if(any(xPosition(2:) - xPosition(:numX-1) <= 0.) .or. &
       any(yPosition(2:) - yPosition(:numY-1) <= 0.) .or. &
       any(zPosition(2:) - zPosition(:numZ-1) <= 0.))     &
      call setStateToFailure(status, "new_Domain: Positions must be increasing, unique.")

    ! -------------------------
     if(.not. stateIsFailure(status)) then
       allocate(new_DomainMono%xPosition(numX), new_DomainMono%yPosition(numY), &
	new_DomainMono%zPosition(numZ), new_DomainMono%temps(numX-1,numY-1,numZ-1))
       new_DomainMono%xPosition(1:) = xPosition(:)
       new_DomainMono%yPosition(1:) = yPosition(:)
       new_DomainMono%zPosition(1:) = zPosition(:)
       new_DomainMono%temps(1:,1:,1:) = temps(:,:,:)

       ! Are the grids regularly spaced? Compare the distance between each pair
       !   of array elements to the distance between the first two.
       ! The default value is false.
       !
       if(all(abs( (xPosition(2:) - xPosition(:numX-1)) -                                &
                   (xPosition(2)  - xPosition(1)) ) <= 2 * spacing(xPosition(2:))) .and. &
          all(abs( (yPosition(2:) - yPosition(:numY-1)) -                               &
                   (yPosition(2)  - yPosition(1)) ) <= 2 * spacing(yPosition(2:))))      &
         new_DomainMono%xyRegularlySpaced = .true.
       if(all(abs( (zPosition(2:) - zPosition(:numZ-1)) -                           &
                   (zPosition(2)  - zPosition(1)) ) <= 2 * spacing(zPosition(2:)))) &
         new_DomainMono%zRegularlySpaced = .true.
       call setStateToSuccess(status)
     end if

  end function new_DomainMono



  !------------------------------------------------------------------------------------------
  subroutine addOpticalComponent3D(thisDomain, componentName,          &
                                   extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, phaseFunctions, &
                                   zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real(8),    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer, optional,           intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    ! 
    ! Add a new optical component to the domain. We check the optical 
    !   properties to make sure they make sense, and also that they can be 
    !   associated with the domain. 
    ! Implementation note: the domain type has an array containing the optical components; 
    !   each new component means allocating a slightly bigger array and copying the 
    !   optical component objects.  

    ! Local variables
    integer                               :: nComponents, i
    type(opticalComponent), &
                    dimension(:), allocatable :: tempComponents
    integer                               :: baseLevel 

    ! -----------------
    ! Checks are done in validateOpticalComponent; zBaseLevel is assumed to be 1
    !   if it isn't supplied, if this doesn't make sense it'll get caught in validateOpticalComponent
    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if
!PRINT *, 'addOpticalComponent3D: size of extinction array', size(extinction,1), size(extinction,2), size(extinction,3)
    call validateOpticalComponent(thisDomain, componentName,          &
                                  extinction, singleScatteringAlbedo, &
                                  phaseFunctionIndex, phaseFunctions, &
                                  baseLevel, status)
    ! -----------------
    if(.not. stateIsFailure(status)) then 
      ! How many components will we have when we're done adding this one? 
      !   Allocate new memory, then copy the data from the pre-existing array. 
      ! 
      if(containsComponents(thisDomain)) then
        nComponents = size(thisDomain%components)
        allocate(tempComponents(nComponents + 1))
        tempComponents(:nComponents) = thisDomain%components(:)
      else
        nComponents = 0
        allocate(tempComponents(1))
      end if 
      
      tempComponents(nComponents + 1) =                                          &
          newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                              phaseFunctionIndex, baseLevel, phaseFunctions)
      
      ! Free the memory associated with the old array
      if(associated(thisDomain%components)) then
	Do i = 1, nComponents
	   call finalizeComponent(thisDomain%components(i))
	End Do
	deallocate(thisDomain%components)       
      end if
      allocate(thisDomain%components(nComponents + 1))
      thisDomain%components = tempComponents 
      Do i = 1, nComponents+1
           call finalizeComponent(tempComponents(i))
      End Do
      deallocate(tempComponents)
      
      call setStateToSuccess(status)
    else
      call setStateToFailure(status, "addOpticalComponent: optical properties aren't valid.")
    end if
  end subroutine addOpticalComponent3D
  !------------------------------------------------------------------------------------------
  subroutine addOpticalComponent1D(thisDomain, componentName,          &
                                   extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, phaseFunctions, &
                                   zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real(8),    dimension(:),       intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:),       intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer,           optional, intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Add a new one-dimensional component to the domain. 
    !   We do this by pretending it's a 3D component with dimension size 1 in the 
    !   x and y directions.

    ! -----------------
    ! Local variables
    integer :: baseLevel

    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if

    call addOpticalComponent3D(thisDomain, componentName,                                                 &
                               reshape(extinction,             (/ 1, 1, size(extinction) /)),             &
                               reshape(singleScatteringAlbedo, (/ 1, 1, size(singleScatteringAlbedo) /)), &
                               reshape(phaseFunctionIndex,     (/ 1, 1, size(phaseFunctionIndex) /)),     &
                               phaseFunctions, zLevelBase = baseLevel, status = status)

  end subroutine addOpticalComponent1D
  !------------------------------------------------------------------------------------------
  subroutine replaceOpticalComponent3D(thisDomain, componentNumber, componentName,                &
                                       extinction, singleScatteringAlbedo,                        &
                                       phaseFunctionIndex, phaseFunctions, &
                                       zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    character (len = *),         intent(in   ) :: componentName
    real(8),    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer, optional,           intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    ! 
    ! Replace one of the domains optical components. This is like adding except we don't have 
    !   to allocate memory.  
    integer                               :: baseLevel 

    ! -----------------
    if(.not. containsComponents(thisDomain)) then
      call setStateToFailure(status, "replaceOpticalComponent: no components to replace" )
    else if(componentNumber < 1 .or. componentNumber > size(thisDomain%components)) then
      call setStateToFailure(status, "replaceOpticalComponent: no components to replace" )
    end if
    ! Checks are done in validateOpticalComponent; zBaseLevel is assumed to be 1
    !   if it isn't supplied, if this doesn't make sense it'll get caught in validateOpticalComponent
    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if
    call validateOpticalComponent(thisDomain, componentName,          &
                                  extinction, singleScatteringAlbedo, &
                                  phaseFunctionIndex, phaseFunctions, &
                                  baseLevel, status)

    ! -----------------
    if(.not. stateIsFailure(status)) then 
      call finalizeComponent(thisDomain%components(componentNumber))
      thisDomain%components(componentNumber) =                                   &
          newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                             phaseFunctionIndex, baseLevel, phaseFunctions)
      call setStateToSuccess(status)
    end if
  end subroutine replaceOpticalComponent3D
  !------------------------------------------------------------------------------------------
  subroutine replaceOpticalComponent1D(thisDomain, componentNumber, componentName, &
                                       extinction, singleScatteringAlbedo,         &
                                       phaseFunctionIndex, phaseFunctions,         &
                                       zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    character (len = *),         intent(in   ) :: componentName
    real(8),    dimension(:),       intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:),       intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer,           optional, intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Replace a new component in the domain with something one-dimensional. It's like adding but with no 
    !   memory headaches. 

    ! -----------------
    ! Local variables
    integer :: baseLevel

    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if

    call replaceOpticalComponent3D(thisDomain, componentNumber, componentName,                            &
                               reshape(extinction,             (/ 1, 1, size(extinction) /)),             &
                               reshape(singleScatteringAlbedo, (/ 1, 1, size(singleScatteringAlbedo) /)), &
                               reshape(phaseFunctionIndex,     (/ 1, 1, size(phaseFunctionIndex) /)),     &
                               phaseFunctions, zLevelBase = baseLevel, status = status)

  end subroutine replaceOpticalComponent1D
  !------------------------------------------------------------------------------------------
  subroutine deleteOpticalComponent(thisDomain, componentNumber, status)
    ! Delete a component from the domain
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Delete a component from the domain. By analogy to the process when 
    !   we add a new component, we create a shorter array in the domain 
    !   object, then copy the remaining optical components to the new array 
    !   and loose the old one. 
    ! 

    ! Local variables
    type(opticalComponent), &
      dimension(size(thisDomain%components) - 1), target :: tempComponents
    
    ! --------------------
    if(.not. isValid(thisDomain)) &
      call setStateToFailure(status, "deleteOpticalComponent: domain hasn't been initialized.")
    if(containsComponents(thisDomain)) then
      if(componentNumber < 1 .or. componentNumber > size(thisDomain%components)) &
        call setStateToFailure(status, "deleteOpticalComponent: non-existent component.")
    else
      call setStateToFailure(status, "deleteOpticalComponent: no components to delete.")
    end if 
    
    ! --------------------
    if(.not. stateIsFailure(status)) then
      if(size(thisDomain%components) == 1) then 
        ! There is only one component, so we deallocate (and nullify) the 
        !   array that holds the components. 
        !
        call finalizeComponent(thisDomain%components(1))
        deallocate(thisDomain%components)
      else 
        tempComponents(:) = (/ thisDomain%components(:(componentNumber - 1)), &
                               thisDomain%components((componentNumber + 1):) /)
        ! This deallocates the pointer to the array that holds the optical components
        !   but not the underlying data. 
        call finalizeComponent(thisDomain%components(componentNumber))
        deallocate(thisDomain%components); allocate(thisDomain%components(size(tempComponents)))
        thisDomain%components(:) = tempComponents(:)
      end if 
      
      call setStateToSuccess(status)
    end if
  end subroutine deleteOpticalComponent
  !------------------------------------------------------------------------------------------
  !  Getting information back from the object
  !------------------------------------------------------------------------------------------
  subroutine getInfo_Domain(thisDomain, numX, numY, numZ, albedo,   &
                            lambda, lambdaIndex, namelistNumLambda,         &
                            domainNumLambda, xPosition, yPosition, zPosition, temps, &
                            numberOfComponents, componentNames,     &
                            totalExt, cumExt, ssa, ext, phaseFuncI, &
 			    inversePhaseFuncs, tabPhase, tabOrigPhase, status)
 
    type(domain),                    intent(in   ) :: thisDomain
    integer,               optional, intent(  out) :: numX, numY, numZ, lambdaIndex
    real(8), optional, intent(out)                    :: lambda, albedo
    real(8),    dimension(:), optional, intent(  out) :: xPosition, yPosition, zPosition
    real(8), dimension(:,:,:), optional, intent( out) :: temps, totalExt
    integer,               optional, intent(  out) :: numberOfComponents
    integer,		optional, intent(in)	   :: namelistNumLambda
    integer,            optional, intent(out)       :: domainNumLambda
    character(len = *), &
             dimension(:), optional, intent(  out) :: componentNames
    real(8), dimension(:,:,:,:), optional, intent(  out) :: cumExt, ssa, ext
    integer, dimension(:,:,:,:), optional, intent(  out) :: phaseFuncI
    type(matrix), dimension(:), optional, intent(out)  :: inversePhaseFuncs, tabPhase, tabOrigPhase
    type(ErrorMessage),              intent(inout) :: status
    
    integer                                        :: j
    ! What can you get back from the domain? The number of cells in the arrays, 
    !   and the locations of the x, y, z, boundaries (which is one bigger). 
    ! Also the number and names of the various optical components.
    ! Also the temperature field-added by Alexandra Jones UIUC Fall 2011 
    
    ! --------------------
    ! Checks: domain is valid, position arrays are big enough (the right size?)
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "getInfo_Domain: domain hasn't been initialized.")
    else
      if(present(namelistNumLambda)) then
	!PRINT *, 'namelistNumLambda= ', namelistNumLambda, 'thisDomain%nlambda= ', thisDomain%nlambda
	if(namelistNumLambda .ne. thisDomain%nlambda) call setStateToFailure(status, &
		"getInfo_Domain: number of wavelengths from namelist does not match number of wavelengths from domain file")
      end if
      if(present(domainNumLambda)) domainNumLambda = thisDomain%nlambda
      if(present(lambdaIndex)) lambdaIndex = thisDomain%lambdaI
      if(present(lambda)) lambda = thisDomain%lambda

      ! Number of positions in each dimension for the property arrays. It's one shorter 
      !   than the position arrays, which describe the boundaries. 
      if(present(numX)) numX = size(thisDomain%xPosition) - 1
      if(present(numY)) numY = size(thisDomain%yPosition) - 1 
      if(present(numZ)) numZ = size(thisDomain%zPosition) - 1
      
      if(present(albedo)) albedo = thisDomain%surfaceAlbedo

      if(present(totalExt)) then
	if(size(totalExt,1) .ne. size(thisDomain%xPosition)-1 .or. size(totalExt,2) .ne. size(thisDomain%yPosition)-1 .or. &
           size(totalExt,3) .ne. size(thisDomain%zPosition)-1) then
	  call setStateToFailure(status, "getInfo_Domain: array for totalExt is wrong dimensions.")
	else
	  totalExt = thisDomain%totalExt
	end if
      endif
      if(present(cumExt)) then
        if(size(cumExt,1) .ne. size(thisDomain%xPosition)-1 .or. size(cumExt,2) .ne. size(thisDomain%yPosition)-1 .or. &
           size(cumExt,3) .ne. size(thisDomain%zPosition)-1 .or. size(cumExt,4) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for cumExt is wrong dimensions.")
        else
	  cumExt = thisDomain%cumulativeExt
        end if
      endif
      if(present(ssa)) then
        if(size(ssa,1) .ne. size(thisDomain%xPosition)-1 .or. size(ssa,2) .ne. size(thisDomain%yPosition)-1 .or. &
           size(ssa,3) .ne. size(thisDomain%zPosition)-1 .or. size(ssa,4) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for ssa is wrong dimensions.")
        else
	  ssa = thisDomain%ssa
        end if
      endif
      if(present(ext)) then
        if(size(ext,1) .ne. size(thisDomain%xPosition)-1 .or. size(ext,2) .ne. size(thisDomain%yPosition)-1 .or. &
           size(ext,3) .ne. size(thisDomain%zPosition)-1 .or. size(ext,4) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for ext is wrong dimensions.")
        else
          ext(:,:,:,1) = thisDomain%totalExt(:,:,:) * thisDomain%cumulativeExt(:,:,:,1)
	  if (size(ext,4) .gt. 1) then
	    forall (j=2:size(ext,4))
		ext(:,:,:,j) = thisDomain%totalExt(:,:,:) * (thisDomain%cumulativeExt(:,:,:,j)-thisDomain%cumulativeExt(:,:,:,j-1))
	    end forall
	  end if
        end if
      endif
      if(present(phaseFuncI)) then
        if(size(phaseFuncI,1) .ne. size(thisDomain%xPosition)-1 .or. size(phaseFuncI,2) .ne. size(thisDomain%yPosition)-1 .or. &
           size(phaseFuncI,3) .ne. size(thisDomain%zPosition)-1 .or. size(phaseFuncI,4) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for phaseFuncIndex is wrong dimensions.")
        else
	  phaseFuncI = thisDomain%phaseFunctionIndex
        end if
      endif
      if(present(inversePhaseFuncs)) then
        if(size(inversePhaseFuncs) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for inversePhaseFuncs is wrong dimensions.")
        else
	  inversePhaseFuncs = thisDomain%inversePhaseFunctions
        end if
      endif
      if(present(tabPhase)) then
        if(size(tabPhase) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for tabulatedPhaseFuncs is wrong dimensions.")
        else
          tabPhase = thisDomain%tabulatedPhaseFunctions
        end if
      endif
      if(present(tabOrigPhase)) then
        if(size(tabOrigPhase) .ne. size(thisDomain%components)) then
          call setStateToFailure(status, "getInfo_Domain: array for tabulatedOriginalPhaseFuncs is wrong dimensions.")
        else
          tabOrigPhase = thisDomain%tabulatedOrigPhaseFunctions
        end if
      endif

      ! Location of boundaries in each dimension
      if(present(xPosition)) then 
        if(size(xPosition) /= size(thisDomain%xPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for x positions is wrong length.")
        else
          xPosition(:) = thisDomain%xPosition(:)
        end if
      end if
      if(present(yPosition)) then 
        if(size(yPosition) /= size(thisDomain%yPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for y positions is wrong length.")
        else
          yPosition(:) = thisDomain%yPosition(:)
        end if
      end if
      if(present(zPosition)) then 
        if(size(zPosition) /= size(thisDomain%zPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for z positions is wrong length.")
        else
          zPosition(:) = thisDomain%zPosition(:)
        end if
      end if
      if(present(temps))then
        if(size(temps) /= size(thisDomain%temps))then
!PRINT *, 'getInfo_Domain: array for temps not the right size'
          call setStateToFailure(status, "getInfo_Domain: array for temps not the right size.")
!	elseif (COUNT(thisDomain%temps .le. 0.0_8) .gt. 0) then
!PRINT *, 'getInfo_Domain: there are temperatures at or below 0.0 K. MIN(thisDomain%temps)= ', MINVAL(thisDomain%temps)
!	  call setStateToFailure(status, "getInfo_Domain: there are temperatures at or below 0.0 K.")
        else
          temps(:,:,:) = thisDomain%temps(:,:,:)
!PRINT *, 'getINfoDomain: min(thisDomain%temps), min(temps) ', MINVAL(thisDomain%temps), MINVAL(temps)
        end if
      end if
      if(containsComponents(thisDomain)) then 
        if(present(numberOfComponents)) numberOfComponents = size(thisDomain%components)
        if(present(componentNames)) then 
          if(size(componentNames) < size(thisDomain%components)) then
            call setStateToFailure(status, "getInfo_Domain: component names array is wrong length")
          else 
            componentNames(:size(thisDomain%components)) = thisDomain%components(:)%name
          end if
        end if
      else 
        if(present(numberOfComponents)) numberOfComponents = 0
        if(present(componentNames))     componentNames(:) = ""
      end if
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
  end subroutine getInfo_Domain
  !------------------------------------------------------------------------------------------
  ! Getting the optical properties of the domain
  !------------------------------------------------------------------------------------------
  subroutine getOpticalPropertiesByComponent(thisDomain, status)
    type(domain),                    intent(inout) :: thisDomain
    type(ErrorMessage),              intent(inout) :: status
    !real(8),    allocatable, dimension(:, :, :)		 :: totalExtinction
    !real(8),    allocatable, dimension(:, :, :, :)	 :: cumulativeExtinction, singleScatteringAlbedo
    !integer, allocatable, dimension(:, :, :, :)		 :: phaseFunctionIndex
    !type(phaseFunctionTable), allocatable, dimension(:)	 :: phaseFunctions
    
    ! 
    ! Fill the optical properties of the domain, component by component. 
    !   The properties returned are defined at all points in the domain. 
    !   If the component is horizontally homogeneous the properties are expanded into 
    !   3D; if it exists only on some subset of vertical layers the properties at 
    !   other layers are set to 0. 
    ! Arrays are dimensioned x, y, z, component. 
    ! Cumulative extinction for the i'th component is the proportion of the total extinction due 
    !   to components 1 to i. It's intended to be used to pick the component doing the extinction
    !   at each scattering/absorption event in a Monte Carlo scheme.  
    ! Phase function tables are returned in arrays. 
    !   
    
    ! Local variable
    integer :: numX, numY, numZ, numComponents
    integer :: i, minZ, maxZ
    
    !------------
    ! Sanity checks
    !   Domain is in use and contains components
    !   x, y, z, dimensions, number of components correct for all arrays
    !   number of entries in each table is correct
    !
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "getOpticalPropertiesByComponent: domain is not initialized.")
    else if (.not. containsComponents(thisDomain)) then
      call setStateToFailure(status, "getOpticalPropertiesByComponent: domain contains no optical components.")
    else 
      !
      ! Checks for x, y, z, sizes
      numX = size(thisDomain%xPosition) - 1; numY = size(thisDomain%yPosition) - 1
      numZ = size(thisDomain%zPosition) - 1; numComponents = size(thisDomain%components)
      allocate(thisDomain%totalExt(1:numX,1:numY,1:numZ))
      allocate(thisDomain%cumulativeExt(1:numX,1:numY,1:numZ,1:numComponents))
      allocate(thisDomain%ssa(1:numX,1:numY,1:numZ,1:numComponents))
      allocate(thisDomain%phaseFunctionIndex(1:numX,1:numY,1:numZ,1:numComponents))
      allocate(thisDomain%forwardTables(1:numComponents))
    end if 

    ! -----------------------
    !   Everything looks ok for now. 
    !
    if(.not. stateIsFailure(status)) then 
!      totalExtinction       (:, :, :)    = 0.
!      cumulativeExtinction  (:, :, :, :) = 0.
!      singleScatteringAlbedo(:, :, :, :) = 0.
!      phaseFunctionIndex    (:, :, :, :) = 0
      
      do i = 1, numComponents
        minZ = thisDomain%components(i)%zLevelBase
        maxZ = minZ + size(thisDomain%components(i)%extinction, 3) - 1
        
        !
        ! In this loop cumulativeExtinction is the extinction by component 
        !
        if(.not. thisDomain%components(i)%horizontallyUniform) then 
          thisDomain%cumulativeExt  (:, :, minZ:maxZ, i) = thisDomain%components(i)%extinction(:, :, :)
          thisDomain%ssa(:, :, minZ:maxZ, i) = thisDomain%components(i)%singleScatteringAlbedo(:, :, :)
          thisDomain%phaseFunctionIndex    (:, :, minZ:maxZ, i) = thisDomain%components(i)%phaseFunctionIndex(:, :, :)
        else
          thisDomain%cumulativeExt(:, :, minZ:maxZ, i) =                                                      &
                                             spread(spread(thisDomain%components(i)%extinction(1, 1, :),  &
                                                           1, nCopies = numX), 2, nCopies = numY)
          thisDomain%ssa(:, :, minZ:maxZ, i) =                                                               &
                                             spread(spread(thisDomain%components(i)%singleScatteringAlbedo(1, 1, :), &
                                                           1, nCopies = numX), 2, nCopies = numY)
          thisDomain%phaseFunctionIndex    (:, :, minZ:maxZ, i) =                                                           &
                                             spread(spread(thisDomain%components(i)%phaseFunctionIndex(1, 1, :), &
                                                           1, nCopies = numX), 2, nCopies = numY)
        end if
        ! We don't finalize the array element phaseFunctions(i) in case something else is pointing to the 
        !   underlying memory. Users should finalize before they pass the array in. 
        ! 
         thisDomain%forwardTables(i) = copy_PhaseFunctionTable(thisDomain%components(i)%table)
         call finalizeComponent(thisDomain%components(i))  
      end do 
	!!! But don't deallocate thisDomain%Components !!!
      !
      ! Compute total extinction and convert cumulative extinction to fractional cumulative extinction. 
      !   Every element of cumulativeExtinction(:, :, :, numComponents) should equal 1 by definition
      !
      do i = 2, numComponents
        thisDomain%cumulativeExt(:, :, :, i) = thisDomain%cumulativeExt(:, :, :, i)  + thisDomain%cumulativeExt(:, :, :, i-1) 
      end do
      thisDomain%totalExt(:, :, :) = thisDomain%cumulativeExt(:, :, :, numComponents)
      where(spread(thisDomain%totalExt, 4, nCopies = numComponents) > tiny(thisDomain%totalExt)) &
        thisDomain%cumulativeExt(:, :, :, :) = thisDomain%cumulativeExt(:, :, :, :)  / &
						spread(thisDomain%totalExt, 4, nCopies = numComponents)
    end if

    !thisDomain%totalExt = totalExtinction
    !thisDomain%cumulativeExt = cumulativeExtinction
    !thisDomain%ssa = singleScatteringAlbedo
    !thisDomain%phaseFunctionIndex = phaseFunctionIndex
    !thisDomain%forwardTables = phaseFunctions
  
  !!!!! CAN'T FIGURE OUT HOW TO FREE THE MEMORY ASSOCIATED WITH THISDOMAIN%COMPONENTS SINCE IT'S A POINTER 
    
  end subroutine getOpticalPropertiesByComponent
  !------------------------------------------------------------------------------------------
!   subroutine getAverageOpticalProperties(thisDomain, extinction, singleScatteringAlbedo, &
!                 phaseFunctionIndex, phaseFunctions, status)
!     type(domain),                 intent( in) :: thisDomain
!     real,    dimension(:, :, :),  intent(out) :: extinction, singleScatteringAlbedo
!     integer, dimension(:, :, :),  intent(out) :: phaseFunctionIndex
!     type(phaseFunctionTable),     intent(out) :: phaseFunctions
!     type(ErrorMessage),           intent(inout) :: status
!     
!     
!   end subroutine getAverageOpticalProperties
  !------------------------------------------------------------------------------------------
  ! Storage and retrieval
  !------------------------------------------------------------------------------------------
  subroutine write_Domain(thisDomain, fileName, status)
   ! has been edited to write out the new temperature field that was added to the domain-AJ UIUC Fall 2011
    type(domain),       intent(in   ) :: thisDomain
    character(len = *), intent(in   ) :: fileName
    type(ErrorMessage), intent(inout) :: status
    
    ! Local variables
    integer                        :: i, j
    logical                        :: fillsDomainInVertical
    
    ! Netcdf-related local variables
    integer, dimension(20) :: ncStatus
    integer                :: ncFileId, xEdgeDimID, yEdgeDimId, zEdgeDimId,  &
                                        xGridDimId, yGridDimId, zGridDimId,  &
                                        extinctionVarId, ssaVarId, indexVarId, ncVarId
    
    ! Checks: domain is valid, contains components(?)
    ncStatus(:) = nf90_NoErr
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "write_Domain: domain hasn't been initialized.") 
    else
      ncStatus( 1) = nf90_create(trim(fileName), nf90_Clobber, ncFileId)
      !
      ! Domain dimensions 
      !
      ncStatus( 2) = nf90_def_dim(ncFileId, "x-Edges", size(thisDomain%xPosition),    xEdgeDimId) 
      ncStatus( 3) = nf90_def_dim(ncFileId, "y-Edges", size(thisDomain%yPosition),    yEdgeDimId) 
      ncStatus( 4) = nf90_def_dim(ncFileId, "z-Edges", size(thisDomain%zPosition),    zEdgeDimId) 
      ncStatus( 5) = nf90_def_dim(ncFileId, "x-Grid",  size(thisDomain%xPosition) - 1, xGridDimId) 
      ncStatus( 6) = nf90_def_dim(ncFileId, "y-Grid",  size(thisDomain%yPosition) - 1, yGridDimId) 
      ncStatus( 7) = nf90_def_dim(ncFileId, "z-Grid",  size(thisDomain%zPosition) - 1, zGridDimId)
      !
      ! Domain variables
      ! 
      ncStatus( 8) = nf90_def_var(ncFileId, "x-Edges", nf90_double, xEdgeDimId, ncVarId) 
      ncStatus( 9) = nf90_def_var(ncFileId, "y-Edges", nf90_double, yEdgeDimId, ncVarId) 
      ncStatus(10) = nf90_def_var(ncFileId, "z-Edges", nf90_double, zEdgeDimId, ncVarId) 
      ncStatus(11) = nf90_def_var(ncFileId, "Temperatures", nf90_double, (/ xGridDimId, yGridDimId, zGridDimId /), ncVarId)
      !
      ! Domain attributes
      !
      ncStatus(12) = nf90_put_att(ncFileID, nf90_Global, "xyRegularlySpaced", asInt(thisDomain%xyRegularlySpaced)) 
      ncStatus(13) = nf90_put_att(ncFileID, nf90_Global,  "zRegularlySpaced", asInt(thisDomain%zRegularlySpaced))
      ncStatus(14) = nf90_put_att(ncFileID, nf90_Global, "lambda", thisDomain%lambda)
      ncStatus(15) = nf90_put_att(ncFileID, nf90_Global,  "lambdaIndex", thisDomain%lambdaI)
      ncStatus(16) = nf90_put_att(ncFileID, nf90_Global,  "numberOfLambdas", thisDomain%nlambda)
      ncStatus(17) = nf90_put_att(ncFileID, nf90_Global, "surfaceAlbedo", thisDomain%surfaceAlbedo)
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "write_Domain: error writing domain information") 
        
      if(.not. stateIsFailure(status) .and. containsComponents(thisDomain)) then
        ncStatus( 1) = nf90_put_att(ncFileID, nf90_Global, "numberOfComponents", size(thisDomain%components))
        do i = 1, size(thisDomain%components)
          ! 
          ! For each component: Attributes 
          !
          ncStatus( 1) = nf90_put_att(ncFileId, nf90_global, trim(makePrefix(i)) // "Name",       &
                                      trim(thisDomain%components(i)%name))
          ncStatus( 2) = nf90_put_att(ncFileId, nf90_global, trim(makePrefix(i)) // "zLevelBase", &
                                      thisDomain%components(i)%zLevelBase)
          !
          ! Dimension definition -  the length of the z dimension can be different than the domains
          ! 
          fillsDomainInVertical = thisDomain%components(i)%zLevelBase == 1 .and. &
                                  size(thisDomain%components(i)%extinction, 3) == (size(thisDomain%zPosition) -1)
          if(fillsDomainInVertical) then
            ncStatus( 6) = nf90_inq_dimid(ncFileId, "z-Grid", zGridDimId)
          else
            ncStatus( 6) = nf90_def_dim(ncFileId, trim(makePrefix(i)) // "z-Grid", &
                                        size(thisDomain%components(i)%extinction, 3), zGridDimId)
          end if 
          !
          ! Variables
          !   Doesn't seem like there will ever be more than 2^15 possible phase functions, 
          !   so we can use a two-byte integer to store the phaseFunctionIndex. 
          ! 
          if(thisDomain%components(i)%horizontallyUniform) then
            ! Variables are 1D
            !
            ncStatus( 7) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "Extinction",             nf90_double, &
                                                       zGridDimId, extinctionVarId)
            ncStatus( 8) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", nf90_double, &
                                                       zGridDimId, ssaVarId)
            ncStatus( 9) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",     nf90_short, &
                                                       zGridDimId, indexVarId)
          else
            ! Variables are 3D
            !
            ncStatus( 7) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "Extinction",             nf90_double, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , extinctionVarId)
            ncStatus( 8) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", nf90_double, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , ssaVarId)
            ! Doesn't seem like there will ever be more than 2^15 possible phase functions...
            ncStatus( 9) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",     nf90_short, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , indexVarId)
          end if
          if(any(ncStatus(:) /= nf90_NoErr)) &
            call setStateToFailure(status,   &
                                   "write_Domain: Error creating definitions for component" // trim(IntToChar(i)))
        end do 
      end if
      !
      ! All the attributes have been written, and the dimensions and variables defined. Now we'll take the
      !   file out of define mode and write the domain data, then the data for each component.
      !
      ncStatus( 1) = nf90_EndDef(ncFileId)
      ncStatus( 2) = nf90_inq_varid(ncFileId, "x-Edges", ncVarID)
      ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, thisDomain%xPosition)
      ncStatus( 4) = nf90_inq_varid(ncFileId, "y-Edges", ncVarID)
      ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, thisDomain%yPosition)
      ncStatus( 6) = nf90_inq_varid(ncFileId, "z-Edges", ncVarID)
      ncStatus( 7) = nf90_put_var(ncFileId, ncVarId, thisDomain%zPosition)
      ncStatus( 8) = nf90_inq_varid(ncFileId, "Temperatures", ncVarID)
      ncStatus( 9) = nf90_put_var(ncFileId, ncVarId, thisDomain%temps(:,:,:))
      if(any(ncStatus(:) /= nf90_NoErr)) then
!        print *, ncStatus(1:9)
        call setStateToFailure(status, "write_Domain: error writing domain data") 
      end if
  
      if(.not. stateIsFailure(status) .and. containsComponents(thisDomain)) then
        do i = 1, size(thisDomain%components)
          ncstatus( 1) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "Extinction",  extinctionVarId)
          ncstatus( 2) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo",  ssaVarId)
          ncstatus( 3) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",  indexVarId)
          if(thisDomain%components(i)%horizontallyUniform) then
            ncStatus( 4) = nf90_put_var(ncFileId, extinctionVarId, &
                                        thisDomain%components(i)%extinction(1, 1, :))
            ncStatus( 5) = nf90_put_var(ncFileId, ssaVarId,        &
                                        thisDomain%components(i)%singleScatteringAlbedo(1, 1, :))
            ncStatus( 6) = nf90_put_var(ncFileId, indexVarId,      &
                                        thisDomain%components(i)%phaseFunctionIndex(1, 1, :))
          else
            ncStatus( 4) = nf90_put_var(ncFileId, extinctionVarId, &
                                        thisDomain%components(i)%extinction(:, :, :))
            ncStatus( 5) = nf90_put_var(ncFileId, ssaVarId,        &
                                        thisDomain%components(i)%singleScatteringAlbedo(:, :, :))
            ncStatus( 6) = nf90_put_var(ncFileId, indexVarId,      &
                                        thisDomain%components(i)%phaseFunctionIndex(:, :, :))
          end if
          call add_PhaseFunctionTable(thisDomain%components(i)%table,            &
                                      fileId = ncFileId,                         &
                                      prefix = "Component" // trim(IntToChar(i)) // "_", &
                                      status = status)
          do j = 1, 6
            if(ncStatus(j) /= nf90_NoErr) &
              call setStateToFailure(status, "write_Domain: " // trim(nf90_StrError(ncStatus(j))))
          end do
          if(stateIsFailure(status)) then 
            !
            ! File write has failed somehow - delete file and exit
            ! 
            ncstatus(1) = nf90_close(ncFileId) 
            open(20, file = trim(fileName))
            close(20, status = "delete")
            exit
          end if
        end do
      end if
      ncStatus(1) = nf90_close(ncFileId)
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
    
  end subroutine write_Domain
  !------------------------------------------------------------------------------------------
  subroutine read_Domain(fileName, commonD, thisDomain, status)
    character(len = *), intent(in   ) :: fileName
    type(commonDomain), intent(in) :: commonD
    type(domain), intent(  out) :: thisDomain
    type(ErrorMessage), intent(inout) :: status
    
    ! Local variables
    integer(kind = selected_int_kind(2)) &
                       :: oneByte
    integer            :: i
    integer            :: nXEdges, nYEdges, nZEdges, nComponents, nZGrid
!    real(8), dimension(:), &
!           allocatable :: xEdges, yEdges, zEdges
!    real(8), dimension(:,:,:), allocatable :: temps

    ! Variable for each component
    real(8),    dimension(:, :, :), allocatable :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), allocatable :: phaseFunctionIndex
    type(phaseFunctionTable)                 :: table
    character(len = maxNameLength)           :: name
    integer                                  :: zLevelBase, lambdaI, nlambda
    logical                                  :: fillsVerticalDomain, horizontallyUniform
    real(8)				     :: lambda, albedo
    
    ! Netcdf-related local variables
    integer, dimension(16) :: ncStatus
    integer                :: ncFileId, ncDimId, ncVarId, zGridDimId, nDims, err
    integer, dimension(3)  :: dimIds
   
    ! There is nothing to check a priori
    ncStatus(:) = nf90_NoErr
    err = nf90_open(trim(fileName), nf90_NoWrite, ncFileID)
    if(err /= nf90_NoErr) then
      call setStateToFailure(status, "read_Domain: Can't open file " // trim(fileName)) 
      PRINT *, "read_Domain: file open error ", err
    end if
    
    if(.not. stateIsFailure(status)) then 
      ncStatus( 1) = nf90_inq_dimid(ncFileId, "x-Edges", ncDimId) 
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nXEdges)
      ncStatus( 3) = nf90_inq_dimid(ncFileId, "y-Edges", ncDimId) 
      ncStatus( 4) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nYEdges)
      ncStatus( 5) = nf90_inq_dimid(ncFileId, "z-Edges", ncDimId) 
      ncStatus( 6) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nZEdges)
      ncStatus( 7) = nf90_inq_dimid(ncFileId, "z-Grid", zGridDimId)
      ncStatus( 8) = nf90_Inquire_Dimension(ncFileId, zGridDimId, len = nZGrid)
!PRINT *, 'n_Edges', nXEdges, nYEdges, nZEdges, nZGrid
!      allocate(commonD%xPosition(nXEdges), commonD%yPosition(nYEdges), commonD%zPosition(nZEdges), commonD%temps(nXEdges-1,nYEdges-1,nZEdges-1))
!      ncStatus( 9) = nf90_inq_varid(ncFileId, "x-Edges", ncVarId)
!      ncStatus(10) = nf90_get_var(ncFileId, ncVarId, commonD%xPosition)
!      ncStatus(11) = nf90_inq_varid(ncFileId, "y-Edges", ncVarId)
!      ncStatus(12) = nf90_get_var(ncFileId, ncVarId, commonD%yPosition)
!      ncStatus(13) = nf90_inq_varid(ncFileId, "z-Edges", ncVarId)
!      ncStatus(14) = nf90_get_var(ncFileId, ncVarId, commonD%zPosition)
      
!      ncStatus(15) = nf90_inq_varid(ncFileId, "Temperatures", ncVarId)
!      ncStatus(16) = nf90_get_var(ncFileId, ncVarId, commonD%temps)
!      if( COUNT(temps .le. 0.0_8) .gt. 0)PRINT *, 'readDomain: there are temps at or below 0.0 K'
!PRINT *, 'readDomain: min temp = ', MINVAL(temps)
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "read_Domain: " // trim(fileName) // &
                               " doesn't look an optical properties file.") 
    end if

!    ncStatus( 7) = nf90_inq_dimid(ncFileId, "z-Grid", zGridDimId) ! need to know zGridDimId to know if later components fill the domain completely 
!    nXEdges = size(commonD%xPosition)
!    nYEdges = size(commonD%yPosition)
!    nZEdges = size(commonD%zPosition) ! need these dimensions for loater allocations even if not first read

    if(.not. stateIsFailure(status)) then 
!PRINT *, 'read_Domain: getting attributes'
      !
      ! Create a new domain using the initialization procedure. 
      !   The domain may have been regularly spaced when it was written to the file, but 
      !   this isn't guaranteed any more. 
      !
      ncStatus( 1) = nf90_get_att(ncFileID, nf90_Global,  "lambda", lambda)
      ncStatus( 2) = nf90_get_att(ncFileID, nf90_Global,  "lambdaIndex", lambdaI)
      ncStatus( 3) = nf90_get_att(ncFileID, nf90_Global,  "numberOfLambdas", nlambda)
      ncStatus( 4) = nf90_get_att(ncFileID, nf90_Global,  "surfaceAlbedo", albedo)

!      call finalize_Domain(thisDomain)
!      thisDomain = new_Domain(xEdges, yEdges, zEdges, lambda, lambdaI, nlambda, albedo, temps, status)
!PRINT *, 'read_Domain: commonD max/mins ', MAXVAL(commonD%xPosition), MINVAL(commonD%xPosition), MAXVAL(commonD%yPosition), MINVAL(commonD%yPosition), MAXVAL(commonD%zPosition), MINVAL(commonD%zPosition), MAXVAL(commonD%temps), MINVAL(commonD%temps)
      thisDomain = new_Domain(commonD, lambda, lambdaI, nlambda, albedo, status)
      ncStatus( 5) = nf90_get_att(ncFileID, nf90_Global, "xyRegularlySpaced", oneByte) 
      if(asLogical(oneByte) .neqv. thisDomain%xyRegularlySpaced) &
        call setStateToWarning(status, "read_Domain: file and new domain don't agree on regularity of x-y spacing.") 
      ncStatus( 6) = nf90_get_att(ncFileID, nf90_Global,  "zRegularlySpaced", oneByte)
      if(asLogical(oneByte) .neqv. thisDomain%zRegularlySpaced) &
        call setStateToWarning(status, "read_Domain: file and new domain don't agree on regularity of z spacing.") 
!      deallocate(xEdges, yEdges, zEdges, temps)
    end if 

    if(.not. stateIsFailure(status)) then 
      !
      ! Read in the data for each component, and then add the component to the domain using the 
      !   same function we'd use from outside the module. 
      !
      ncStatus( 1) = nf90_get_att(ncFileID, nf90_Global, "numberOfComponents", nComponents)
!PRINT *, 'read_Domain: reading ', nComponents, ' components'
      do i = 1, nComponents
        ncStatus( 2) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "Name", name)
        ncStatus( 3) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "zLevelBase", zLevelBase)
        !
        ! Is the component horizontally homogeneous? Does it fill the entire domain vertically? 
        ! 
        ncStatus( 4) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "Extinction", ncVarId)
        ncStatus( 5) = nf90_Inquire_Variable(ncFileId, ncVarId, ndims = nDims, dimids = dimIds)
        horizontallyUniform = (ndims == 1)
        fillsVerticalDomain = (dimIds(ndims) == zGridDimId)
!PRINT *, dimIds, zGridDimId
        if(fillsVerticalDomain) then
          nZGrid = nZEdges - 1
!PRINT *, 'fillsVerticalDomain=TRUE', nZgrid, nZEdges
        else
          ncStatus( 6) = nf90_inq_dimId(ncFileId, trim(makePrefix(i)) // "z-Grid", ncDimId)
          ncStatus( 7) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nZGrid)
!PRINT *, 'fillsVerticalDomain=FALSE', ncDimId, nZgrid, nZEdges
        end if
        !
        ! Read in the scalar 3D or 1D fields
        !
        if(horizontallyUniform) then
          allocate(extinction(1, 1, nZGrid), singleScatteringAlbedo(1, 1, nZGrid), phaseFunctionIndex(1, 1, nZGrid)) 
!PRINT *, 'read_DOmain: horizontally uniform. extinction size: ', size(extinction,1), size(extinction,2), size(extinction,3)
          ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, extinction(1, 1, :))
          ncStatus( 9) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", ncVarId)
          ncStatus(10) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedo(1, 1, :))
          ncStatus(11) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex", ncVarId)
          ncStatus(12) = nf90_get_var(ncFileId, ncVarId, phaseFunctionIndex(1, 1, :))
        else 
          allocate(            extinction(nXEdges - 1, nYEdges - 1, nZGrid), &
                   singleScatteringAlbedo(nXEdges - 1, nYEdges - 1, nZGrid), &
                       phaseFunctionIndex(nXEdges - 1, nYEdges - 1, nZGrid))
!PRINT *, 'read_Domain: horizontally nonuniform. extinction size: ', size(extinction,1), size(extinction,2), size(extinction,3)
          ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, extinction(:, :, :))
          ncStatus( 9) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", ncVarId)
          ncStatus(10) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedo(:, :, :))
          ncStatus(11) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex", ncVarId)
          ncStatus(12) = nf90_get_var(ncFileId, ncVarId, phaseFunctionIndex(:, :, :))
        end if
        if(any(ncStatus(:) /= nf90_NoErr)) then
          call setStateToFailure(status, "read_Domain: Error reading scalar fields from file " // trim(fileName))
	  !PRINT *, ncStatus(:)
	end if
        !
        ! Read in the phase function table(s) 
        !
        if(.not. stateIsFailure(status)) &
          call read_PhaseFunctionTable(fileId = ncFileId, table = table,                  &
                                       prefix = "Component" // trim(IntToChar(i)) // "_", &
                                       status = status)
          
        !
        ! Add the new component to the domain. 
        !
        if(.not. stateIsFailure(status)) then

!PRINT *, "read_Domain: extinction size ", size(extinction, 1), size(extinction, 2), size(extinction, 3)
          call addOpticalComponent(thisDomain, name, extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, table, zLevelBase = zLevelBase,   &
                                   status = status)
        else
          call setStateToFailure(status, "read_Domain: Error reading phase function table.") 
        end if
        deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex)
      end do
      if(.not. stateIsFailure(status))  call getOpticalPropertiesByComponent(thisDomain, status)
!      if(.not. stateIsFailure(status)) then
!         do i = 1, 4
!	    call deleteOpticalComponent(thisDomain, i, status)
!         end do
!      end if
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if
  end subroutine read_Domain
  !------------------------------------------------------------------------------------------
  ! Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_Domain(thisDomain) 
    ! Return the variable to it uninitialized state
    type(domain), intent(inout) :: thisDomain
    
    ! Loca variable
    integer :: i
    
    thisDomain%xyRegularlySpaced = .false.
    thisDomain%zRegularlySpaced  = .false.
    
    ! Position vectors
    if(associated(thisDomain%xPosition)) NULLIFY(thisDomain%xPosition)
    if(associated(thisDomain%yPosition)) NULLIFY(thisDomain%yPosition)
    if(associated(thisDomain%zPosition)) NULLIFY(thisDomain%zPosition)
    if(associated(thisDomain%temps)) NULLIFY(thisDomain%temps)
	
	! Cumulative arrays
	if(associated(thisDomain%totalExt)) DEALLOCATE(thisDomain%totalExt)
	if(associated(thisDomain%cumulativeExt)) DEALLOCATE(thisDomain%cumulativeExt)
	if(associated(thisDomain%ssa)) DEALLOCATE(thisDomain%ssa)
	if(associated(thisDomain%phaseFunctionIndex)) DEALLOCATE(thisDomain%phaseFunctionIndex)
	
    !  Optical components
    if(containsComponents(thisDomain)) then
      do i = 1, size(thisDomain%components)
        call finalizeComponent(thisDomain%components(i))
      end do
      if(associated(thisDomain%components)) deallocate(thisDomain%components)
    end if
    ! phaseTables
     do i = size(thisDomain%components), 1, -1
	if(associated(thisDomain%forwardTables)) CALL finalize_PhaseFunctionTable(thisDomain%forwardTables(i))
	if(associated(thisDomain%tabulatedPhaseFunctions)) CALL finalize_Matrix(thisDomain%tabulatedPhaseFunctions(i))
	if(associated(thisDomain%tabulatedOrigPhaseFunctions)) CALL finalize_Matrix(thisDomain%tabulatedOrigPhaseFunctions(i))
	if(associated(thisDomain%inversePhaseFunctions)) CALL finalize_Matrix(thisDomain%inversePhaseFunctions(i))
     end do
      if(associated(thisDomain%tabulatedPhaseFunctions)) DEALLOCATE(thisDomain%tabulatedPhaseFunctions)
      if(associated(thisDomain%forwardTables)) DEALLOCATE(thisDomain%forwardTables)
      if(associated(thisDomain%tabulatedOrigPhaseFunctions)) DEALLOCATE(thisDomain%tabulatedOrigPhaseFunctions)
      if(associated(thisDomain%inversePhaseFunctions)) DEALLOCATE(thisDomain%inversePhaseFunctions)
    
  end subroutine finalize_Domain
  !------------------------------------------------------------------------------------------
  subroutine finalizeComponent(component) 
    ! Return the variable to it uninitialized state
    type(opticalComponent), intent(out) :: component
    
    if(associated(component%extinction))             deallocate(component%extinction)
    if(associated(component%singleScatteringAlbedo)) deallocate(component%singleScatteringAlbedo)
    if(associated(component%phaseFunctionIndex))     deallocate(component%phaseFunctionIndex)
    
    call finalize_PhaseFunctionTable(component%table)

    component%horizontallyUniform     = .false. 
    component%zLevelBase              = 0
    component%name                    = ""
    
  end subroutine finalizeComponent
  !------------------------------------------------------------------------------------------
  ! Utilities, used inside the module
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  ! Optical component objects aren't visible outside the module. so this creation routine isn't
  !   visible either.
  !------------------------------------------------------------------------------------------
  function newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                               phaseFunctionIndex, zLevelBase, phaseFunctions) result(newOne)
    character (len = *),         intent( in) :: componentName
    real(8),    dimension(:, :, :), intent( in) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent( in) :: phaseFunctionIndex
    integer,           optional, intent( in) :: zLevelBase
    type(phaseFunctionTable),    intent( in) :: phaseFunctions
    type(opticalComponent)                   :: newOne
        
    ! Local variables
    integer :: nx, ny, nz
    
    ! Optical components are always associated with a domain. 
    !   So all error checking (grid sizes, etc.) that might be done here is instead done
    !   in validateOpticalComponent, which knows about the grid). 
    nx = size(extinction, 1); ny = size(extinction, 2); nz = size(extinction, 3)
    newOne%name = trim(componentName)
    newOne%zLevelBase = zLevelBase
    allocate(newOne%extinction(nx, ny, nz), newOne%singleScatteringAlbedo(nx, ny, nz), &
             newOne%phaseFunctionIndex(nx, ny, nz))
    newOne%extinction(:, :, :)             = extinction(:, :, :)
    newOne%singleScatteringAlbedo(:, :, :) = singleScatteringAlbedo(:, :, :)
    newOne%phaseFunctionIndex(:, :, :)     = phaseFunctionIndex(:, :, :)
    newOne%table = copy_PhaseFunctionTable(phaseFunctions)
    if(nx == 1 .and. ny == 1) newOne%horizontallyUniform = .true.
    
    if(present(zLevelBase)) then
      newOne%zLevelBase = zLevelBase
    else 
      newOne%zLevelBase = 1
    end if  
    
  end function newOpticalComponent
  !------------------------------------------------------------------------------------------
  subroutine validateOpticalComponent(thisDomain, componentName,          &
                                      extinction, singleScatteringAlbedo, &
                                      phaseFunctionIndex, table,          &
                                      zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real(8),    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: table
    integer,                     intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Check to be sure that data being used to add or replace components is compatible 
    !   with the existing domain. This is folded into a separate subroutine because the logic
    !   is complicated enough that we don't want to replicate it everywhere. 
    

    ! Local variables
    integer :: numX, numY, numZ, numPhaseFunctions
        
    if(isValid(thisDomain)) then 
      numX = size(extinction, 1); numY = size(extinction, 2); numZ = size(extinction, 3)
      
      ! Are the arrays all the same size? 
      if(any((/ size(singleScatteringAlbedo, 1), size(phaseFunctionIndex, 1) /) /= numX) .or. &
         any((/ size(singleScatteringAlbedo, 2), size(phaseFunctionIndex, 2) /) /= numY) .or. &
         any((/ size(singleScatteringAlbedo, 3), size(phaseFunctionIndex, 3) /) /= numZ)) &
        call setStateToFailure(status, "validateOpticalComponent: optical property grids must be the same size.") 
        
      ! Do the arrays conform to the grid in the domain? 
      if(.not. any(numX == (/ 1, size(thisDomain%xPosition) - 1 /)) .or. &
         .not. any(numY == (/ 1, size(thisDomain%yPosition) - 1/)))  then !   &
!PRINT *, 'numX, size(thisDomain%xPosition), numY, size(thisDomain%yPosition) ', numX, size(thisDomain%xPosition), numY, size(thisDomain%yPosition)          
        call setStateToFailure(status, "validateOpticalComponent: arrays don't conform to horizontal extent of domain.")
      end if
      if(zLevelBase + numZ - 1 > size(thisDomain%zPosition) .or. zLevelBase < 1) &
        call setStateToFailure(status, "validateOpticalComponent: arrays don't conform to vertical extent of domain.")
     
      ! Resonable values for the properties
      if(any(extinction(:, :, :) < 0.)) &
        call setStateToFailure(status, "validateOpticalComponent: LambdaI= " // &
		IntToChar(thisDomain%lambdaI) // "extinction must be >= 0.")
      if(any(singleScatteringAlbedo(:, :, :) < 0.) .or. any(singleScatteringAlbedo(:, :, :) > 1. )) &
        call setStateToFailure(status, "validateOpticalComponent: singleScatteringAlbedo must be between 0 and 1")
      
      ! Check the phase function table
      if(.not. stateIsFailure(status)) then 
        call getInfo_PhaseFunctionTable(table, nEntries = numPhaseFunctions, status = status)
        if(any(phaseFunctionIndex(:, :, :) < 0 .or. phaseFunctionIndex(:, :, :) > numPhaseFunctions)) &
          call setStateToFailure(status, "validateOpticalComponent: phase function index is out of bounds")
        ! Are the phase functions ready to go ? 
        if(.not. isReady_PhaseFunctionTable(table)) &
          call setStateToFailure(status, "validateOpticalComponent: phase function table is not ready.")
      end if 
      
      ! We could check to see if the component names overlap. 
    else
      call setStateToFailure(status, "validateOpticalComponent: domain hasn't been initialized.")
    end if

    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine validateOpticalComponent
  !------------------------------------------------------------------------------------------
  logical function isValid(thisDomain) 
    type(domain), intent(in) :: thisDomain
    ! Checks to see if domain is associated with an initialized
    !   object. 
    
    isValid = associated(thisDomain%xPosition) .and. &
              associated(thisDomain%yPosition) .and. &
              associated(thisDomain%zPosition) .and. &
              associated(thisDomain%temps)
  end function isValid
  !------------------------------------------------------------------------------------------
  logical function containsComponents(thisDomain) 
    type(domain), intent(in) :: thisDomain
    ! Checks to see if this domain contains any optical components. 
    
    containsComponents = associated(thisDomain%components)
  end function containsComponents
  !------------------------------------------------------------------------------------------
  function makePrefix(i) 
    integer, intent( in) :: i
    character(len = 32)  :: makePrefix
    ! Constructs a unique prefix for each component
    !   We put this in a function so reading and writing subroutine are consistent
    !   Uses utility IntToChar from module characterUtils 
    
    character(len = 16), parameter :: prefixBase = "Component"
    
    makePrefix = trim(prefixBase) // trim(IntToChar(i)) // "_"
  end function makePrefix
  !------------------------------------------------------------------------------------------
  elemental function asInt(inValue)
    logical,                 intent( in) :: inValue
    integer(kind = selected_Int_Kind(2)) :: asInt
    
    asInt = 0; if(inValue) asInt = 1
    
  end function asInt
  !------------------------------------------------------------------------------------------
  elemental function asLogical(inValue)
    integer(kind = selected_Int_Kind(2)), intent( in) :: inValue
    logical :: asLogical
    
    asLogical = .not. (inValue == 0)
    
  end function asLogical
  !------------------------------------------------------------------------------------------
  function new_Matrix(array)
    real, dimension(:, :), intent(in) :: array
    type(matrix)          :: new_Matrix

    new_Matrix%numX         = size(array, 1)
    new_Matrix%numY         = size(array, 2)
    allocate(new_Matrix%values(size(array, 1), size(array, 2)))
    new_Matrix%values(:, :) = array(:, :)
  end function new_Matrix
  !------------------------------------------------------------------------------------------
  subroutine finalize_Matrix(thisMatrix)
    type(matrix), intent(inout) :: thisMatrix

    thisMatrix%numX = 0; thisMatrix%numY = 0
    if(associated(thisMatrix%values)) deallocate(thisMatrix%values)
  end subroutine finalize_Matrix
  !------------------------------------------------------------------------------------------
  subroutine accumulateExtinctionAlongPath(thisDomain, directionCosines,         &
                                                xPos, yPos, zPos, xIndex, yIndex, zIndex, &
                                                extAccumulated, extToAccumulate)
!  pure subroutine accumulateExtinctionAlongPath(thisIntegrator, directionCosines,         &
!                                                xPos, yPos, zPos, xIndex, yIndex, zIndex, &
!                                                extAccumulated, extToAccumulate)
    !
    ! Trace through the medium in a given direction accumulating extinction
    !   along the path. Tracing stops either when the boundary is reached or
    !   when the accumulated extinction reaches the (optional) value extToAccumulate.
    !   Reports the final position and optical path accumulated.
    !
    
    type(Domain),       intent(in   ) :: thisDomain
    real, dimension(3), intent(in   ) :: directionCosines
    real(8),               intent(inout) :: xPos, yPos, zPos
    integer,            intent(inout) :: xIndex, yIndex, zIndex
    real,               intent(  out) :: extAccumulated
    real, optional,     intent(in   ) :: extToAccumulate

    ! Local variables
    integer               :: nXcells, nYcells, nZcells
    real(8)                  :: thisStep, totalPath
    real(8)               :: z0, zMax,thisCellExt
    real(8),    dimension(3) :: step
    integer, dimension(3) :: SideIncrement, CellIncrement

    extAccumulated = 0.; totalPath = 0.

     ! Make some useful parameters
    nXcells = size(thisDomain%xPosition) - 1
    nYcells = size(thisDomain%yPosition) - 1
    nZcells = size(thisDomain%zPosition) - 1
    ! Side increment is 1 where directionCosines is > 0; 0 otherwise
    sideIncrement(:) = merge(1, 0, directionCosines(:) >= 0.)
    ! Cell increment is +1 where direction cosines > 0; -1 otherwise
    CellIncrement(:) = merge(1, -1, directionCosines(:) >= 0.)

    z0   = thisDomain%zPosition(1)
    zMax = thisDomain%zPosition(nZCells + 1)

    accumulationLoop: do
      !
      ! step - how far away is the closest cell boundary along the direction of travel
      !
      ! How big a step must we take to reach the next edge?  Go the right direction (+ or -)
      !    and find out for each dimension (but don't divide by zero).
      !
!if(zIndex+SideIncrement(3) .gt. 37) PRINT *, zIndex, SideIncrement(3)
      where(abs(directionCosines) >= 2. * tiny(directionCosines))
        step(:) = (/ thisDomain%xPosition(xIndex + SideIncrement(1)) - xPos,     &
                     thisDomain%yPosition(yIndex + SideIncrement(2)) - yPos,     &
                     thisDomain%zPosition(zIndex + SideIncrement(3)) - zPos /) / &
                  directionCosines(:)
      elsewhere
        step(:) = huge(step)
      end where

       ! The step length across the cell is the smallest of the three directions
       !   We guard against step being negative, which can happen if the
       !   direction cosine or the distance to the boundary is very small
       !
      thisStep = minval(step(:))
      if (thisStep <= 0.0_8) then
        extAccumulated = -2.0  ! Error flag
        exit accumulationLoop
      end if

       ! If this cell pushes the optical path past the desired extToAccumulate,
       !  then find how large a step is needed, take it, and exit

      thisCellExt = thisDomain%totalExt(xIndex, yIndex, zIndex)

      if (present(extToAccumulate)) then
        if(extAccumulated + thisStep * thisCellExt > extToAccumulate) then
          thisStep = dble(extToAccumulate - extAccumulated) / thisCellExt
          xPos = xPos + dble(thisStep * directionCosines(1))
          yPos = yPos + dble(thisStep * directionCosines(2))
          zPos = zPos + dble(thisStep * directionCosines(3))
          totalPath = totalPath + thisStep
          extAccumulated = extToAccumulate
          exit accumulationLoop
        end if
      end if

       ! Add this cell crossing to the accumulated optical path and distance

      extAccumulated = extAccumulated + thisStep*thisCellExt
      totalPath = totalPath + thisStep

       ! Determine which side of the cell we're going to hit first (the smallest step).
       !   Set the position to that side of the cell and increment the cell indices to
       !   the next cell, which means extra code for periodic boundaries in X and Y.
       ! If we wind up within spacing() of the coming boundary just say the position is
       !   in the next cell. (This slows things down but protects againt rounding errors.)

      if(step(1) <= thisStep) then
        xPos = thisDomain%xPosition(xIndex + SideIncrement(1))
        xIndex = xIndex + CellIncrement(1)
      else
        xPos = xPos + dble(thisStep * directionCosines(1))
        if(abs(thisDomain%xPosition(xIndex + sideIncrement(1)) - xPos) <= 2 * spacing(xPos)) &
             xIndex = xIndex + cellIncrement(1)      !
      end if

      if(step(2) <= thisStep) then
        yPos = thisDomain%yPosition(yIndex+SideIncrement(2))
        yIndex = yIndex + CellIncrement(2)
      else
        yPos = yPos + dble(thisStep * directionCosines(2))
        if(abs(thisDomain%yPosition(yIndex + sideIncrement(2)) - yPos) <= 2 * spacing(yPos)) &
                yIndex = yIndex + cellIncrement(2)
      end if

      if(step(3) <= thisStep) then
        zPos = thisDomain%zPosition(zIndex+SideIncrement(3))
        zIndex = zIndex + CellIncrement(3)
      else
        zPos = zPos + dble(thisStep * directionCosines(3))
        if(abs(thisDomain%zPosition(zIndex + sideIncrement(3)) - zPos) <= 2 * spacing(zPos)) &
          zIndex = zIndex + cellIncrement(3)
      end if

      !
      ! Enforce periodicity
      !
      if (xIndex <= 0) then
        xIndex = nXcells
        xPos = thisDomain%xPosition(xIndex+1) + cellIncrement(1) * 2 * spacing(xPos)
      else if (xIndex >= nXcells+1) then
        xIndex = 1
        xPos = thisDomain%xPosition(xIndex) + cellIncrement(1) * 2 * spacing(xPos)
      end if

      if (yIndex <= 0) then
        yIndex = nYcells
        yPos = thisDomain%yPosition(yIndex+1) + cellIncrement(1) * 2 * spacing(yPos)
      else if (yIndex >= nYcells+1) then
        yIndex = 1
        yPos = thisDomain%yPosition(yIndex) + cellIncrement(1) * 2 * spacing(yPos)
      end if

      !
      ! Out the top?
      !
      if(zIndex > nZCells) then
        zPos = zMax + 2 * spacing(zMax)
        exit accumulationLoop
      end if

      !
      ! Hit the bottom?
      !
      if(zIndex < 1) then
        zPos = z0
        exit accumulationLoop
      end if

    end do accumulationLoop
  end subroutine accumulateExtinctionAlongPath
  !------------------------------------------------------------------------------------------                  
  subroutine tabulateInversePhaseFunctions(thisDomain, tableSize, status)
    !
    ! Tabulate the inverse (cummulative) phase functions (i.e. scattering angle
    !   as a function of the position in the cumulative distribution).
    !
    type(domain),  intent(inout) :: thisDomain
    integer, intent(in)		:: tableSize
    type(ErrorMessage),intent(inout) :: status

    ! Local variables
    integer                            :: i, numComponents, nEntries, nSteps
    logical                            :: computeThisTable
    real, dimension(:, :), allocatable :: tempMatrix

    numComponents = size(thisDomain%forwardTables)
    if(.not. associated(thisDomain%inversePhaseFunctions)) &
        allocate(thisDomain%inversePhaseFunctions(numComponents))

    componentLoop: do i = 1, numComponents
      !
      ! Does the table already exist at high enough resolution?
      !
      computeThisTable = .true.
      if(associated(thisDomain%inversePhaseFunctions(i)%values)) then
        if(thisDomain%inversePhaseFunctions(i)%numX >= tableSize) computeThisTable = .false.
      end if

      if(computeThisTable) then
        !
        ! Compute at the minimum desired "resolution" (number of intervals bet. 0 and 1)
        !   computeInversePhaseFuncTable expects a simple 2D array, dimensioned nSteps, nEntries
        !   We need to copy this into our "matrix" type
        !
        call getInfo_phaseFunctionTable(thisDomain%forwardTables(i), nEntries = nEntries, status = status)
        if(stateIsFailure(status)) exit componentLoop
        nSteps = tableSize
        allocate(tempMatrix(nSteps, nEntries))
        !
        ! Use the code in the inversePhaseFunctions module to compute the inverse table
        !
        call computeInversePhaseFuncTable(thisDomain%forwardTables(i), tempMatrix, status)
        if(stateIsFailure(status)) exit componentLoop
        thisDomain%inversePhaseFunctions(i) = new_Matrix(tempMatrix)
        deallocate(tempMatrix)
      end if
    end do componentLoop

    if(stateIsFailure(status)) then
      call setStateToFailure(status, "tabulateInversePhaseFunctions: failed on component" // trim(intToChar(i)) )
    else
      call setStateToSuccess(status)
    end if

  end subroutine tabulateInversePhaseFunctions
  !------------------------------------------------------------------------------------------
  subroutine tabulateForwardPhaseFunctions(thisDomain, tableSize, hybrid, hybridWidth, status)
    type(Domain),  intent(inout) :: thisDomain
    integer, intent(in)			:: tableSize
    logical, intent(in)			:: hybrid
    real, intent(in)			:: hybridWidth
    type(ErrorMessage),intent(inout) :: status

    ! Local variables
    integer                            :: i, j, numComponents, nEntries, nSteps
    logical                            :: computeThisTable
    real, dimension(:),    allocatable :: angles
    real, dimension(:, :), allocatable :: tempMatrix

    !
    ! We store the forward tables as matrices evenly spaced in angle to make interpolation simple.
    !
    numComponents = size(thisDomain%forwardTables)
    if(.not. associated(thisDomain%tabulatedPhaseFunctions)) &
        allocate(thisDomain%tabulatedPhaseFunctions(numComponents))
    if(.not. associated(thisDomain%tabulatedOrigPhaseFunctions)) &
        allocate(thisDomain%tabulatedOrigPhaseFunctions(numComponents))

    componentLoop: do i = 1, numComponents
      !
      ! Does the matrix form of the table already exist high enough resolution?
      !   The tabulated phase functions and the original versions are always computed at the same resolution
      !
      computeThisTable = .true.
      if(associated(thisDomain%tabulatedPhaseFunctions(i)%values)) then
        if(thisDomain%tabulatedPhaseFunctions(i)%numX >= tableSize) computeThisTable = .false.
      end if

      if(computeThisTable) then
        !
        ! Compute the phase function values at nSteps points equally spaced in angle from 0 to pi radians
        !
        nSteps = tableSize
        call getInfo_PhaseFunctionTable(thisDomain%forwardTables(i), nEntries = nEntries, status = status)
        if(stateIsFailure(status)) exit componentLoop
        allocate(angles(nSteps), tempMatrix(nSteps, nEntries))
        angles(:) = (/ (j, j = 0, nSteps - 1) /) / real(nSteps - 1) * Pi
        call getPhaseFunctionValues(thisDomain%forwardTables(i), angles(:), tempMatrix(:, :), status)

        thisDomain%tabulatedOrigPhaseFunctions(i) = new_Matrix(tempMatrix)

        if(hybrid) then
          if(hybridWidth > 0.)                                   &
            tempMatrix(:, :) = computeHybridPhaseFunctions(angles(:), tempMatrix(:, :), hybridWidth)
        end if
        !
        ! Copy the tabulated phase functions into a matrix
        !
        thisDomain%tabulatedPhaseFunctions(i) = new_Matrix(tempMatrix)
        deallocate(angles, tempMatrix)
      end if
    end do componentLoop

    if(stateIsFailure(status)) then
      call setStateToFailure(status, "tabulatePhaseFunctions: failed on component" // trim(intToChar(i)) )
    else
      call setStateToSuccess(status)
    end if
  end subroutine tabulateForwardPhaseFunctions
  !------------------------------------------------------------------------------------------
  pure function computeHybridPhaseFunctions(angles, values,  GaussianWidth) result(newValues)
    real, dimension(:),    intent( in) :: angles         ! phase function angles in radians
    real, dimension(:, :), intent( in) :: values         ! Phase  function, nEntries by nAngles
    real,                  intent( in) :: GaussianWidth  ! In degrees
    real, dimension(size(values, 1), size(values, 2) ) :: newValues
    !
    ! Creates a phase function that's a hybrid of the original and a Gaussian of
    !   specified width. The Gaussian part replaces the forward  peak, and is
    !   continuous with the original phase function

    ! Local variables.
    real, dimension(size(angles)) :: gaussianValues, angleCosines
    real                          :: P0, lowDiff, upDiff, midDiff
    integer                       :: nAngles, transitionIndex
    integer                       :: i, lowerBound, upperBound,  midPoint, increment

    nAngles = size(angles)
    angleCosines(:) = cos(angles(:))
    !
    ! Gaussian phase function values - we won't need most of these
    !   but it's easier to compute them all than keep adding
    !
    gaussianValues(:) = exp(-( angles(:)/(GaussianWidth * Pi/180) )**2)

    ! First set the output phase function in input one in case there is no root
    newValues(:, :) = values(:, :)
    entryLoop: do i = 1, size(values, 2)

      ! Set the lower transition bound according to width of the Gaussian
      lowerBound = findIndex(GaussianWidth * Pi/180., angles(:)) + 1
      if(lowerBound >= nAngles - 2) exit entryLoop
      ! We haven't found the position of the Gaussian width in the table

      ! We want the transition angle at which the two phase functions are the same
      !   (i.e. the difference between them is 0).
      !
      ! First we "hunt", expanding the range in which we're trying to bracket the value
      !
      lowDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), lowerBound)
      increment = 1
      huntingLoop: do
        upperBound = min(lowerBound + increment, nAngles - 1)
        upDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), upperBound)

        if (lowerBound == nAngles - 1) cycle entryLoop  ! There's no root, so use the original phase function
        if (lowDiff * upDiff < 0) exit huntingLoop

        lowerBound = upperBound
        lowDiff = upDiff
        increment = increment * 2
      end do huntingLoop

      ! Bisection: figure out which half of the remaining interval holds the
      !   root, discard the other half, and repeat
      bisectionLoop: do
        if (upperBound <= lowerBound + 1) exit bisectionLoop
        midPoint = (lowerBound + upperBound)/2
        midDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), midPoint)
        if (midDiff * upDiff < 0) then
          lowerBound = midPoint
          lowDiff = midDiff
        else
          upperBound = midPoint
          upDiff = midDiff
        end if
      end do bisectionLoop

      transitionIndex = lowerBound
      P0 = computeNormalization(angleCosines(:), values(:, i),  gaussianValues(:), transitionIndex)
      newValues(:transitionIndex,   i) = P0 * gaussianValues(:transitionIndex)
      newValues(transitionIndex+1:, i) = values(transitionIndex+1:, i)
    end do entryLoop

  end function computeHybridPhaseFunctions
  !------------------------------------------------------------------------------------------
  pure function phaseFuncDiff(angleCosines, values, gaussianValues,  transitionIndex) &
                                        result(d)
    !
    ! Compute the difference between the normalized Gaussian phase function and the
    !   original phase function at the transition index.
    !
    real, dimension(:), intent(in) :: angleCosines, values,  gaussianValues
    integer,            intent(in) :: transitionIndex
    real                           :: d

    real :: P0

    P0 = computeNormalization(angleCosines, values, gaussianValues,  transitionIndex)
    d = P0 * gaussianValues(transitionIndex) - values(transitionIndex)
  end function phaseFuncDiff
  !------------------------------------------------------------------------------------------
  pure function computeNormalization(angleCosines, values,  gaussianValues, transitionIndex) result(P0)
    real, dimension(:), intent(in) :: angleCosines, values,  gaussianValues
    integer,            intent(in) :: transitionIndex
    real                           :: P0

    integer :: nAngles
    real    :: IntegralGaus, IntegralOrig

    ! Normalization for the Gaussian part of the phase function, computed by
    !   forcing the complete phase function to be normalized.
    !
    nAngles = size(angleCosines)
    IntegralGaus = dot_product( &
                   0.5*(gaussianValues(1:transitionIndex-1) + gaussianValues(2:transitionIndex)), &
                   angleCosines(1:transitionIndex-1) - angleCosines(2:transitionIndex) )
    IntegralOrig = dot_product( &
                   0.5*(values(transitionIndex:nAngles-1) + values(transitionIndex+1:nAngles)), &
                   angleCosines(transitionIndex:nAngles-1) - angleCosines(transitionIndex+1:nAngles) )
    if (IntegralOrig >= 2.0) then
      P0 = 1.0/IntegralGaus
    else
      P0 = (2. - IntegralOrig) / IntegralGaus
    end if
  end function computeNormalization
  !------------------------------------------------------------------------------------------
  subroutine calc_RayleighScattering(lambda,rho,N,ext,ssa,phaseInd,table,status)
   real(8),               intent(in)                   :: N(:), rho(:)
   real(8),               intent(in)                   :: lambda
   real(8), dimension(:), intent(inout)                :: ext, ssa
   integer, dimension(:), intent(inout)                :: phaseInd
   type(ErrorMessage),    intent(inout)                :: status
   type(phaseFunctionTable), intent(out)               :: table
   

   real(8)                                             :: mr1, f=1.060816681_8, rho0=1.275_8
   real, dimension(2)                                  :: LG
   type(phaseFunction), dimension(1)                   :: phaseFunc
   integer                                             :: err, ncid, dimid, varid, nz
   
!   err=nf90_open("/mnt/a/u/sciteam/aljones4/I3RC/Tools/density.nc",nf90_NoWrite, ncid)
!seFunc(1) = new_PhaseFunction(LG,status=status)V   err=nf90_inq_dimid(ncid, "z-grid", dimid)
!   err=nf90_Inquire_Dimension(ncid, dimid, len=nz)
!!   allocate(rho(1:nz))
!   err=nf90_inq_varid(ncid, "density", varid)
!   err=nf90_get_var(ncid, varid, rho)
   
!   mrsq = (1 + 6.4328E-5 + (2.94981E-2/(146-(lambda**(-2)))) + (2.554E-4/(41-(lambda**(-2)))))**2
!   ext = ((8.0E24_8)*f*(Pi**3)*(mrsq-1)**2)/(3.0_8*(lambda**4)*(N**2)) 
   mr1 = 6.4328E-5 + (2.94981E-2/(146-(lambda**(-2)))) + (2.554E-4/(41-(lambda**(-2))))
   ext = (32.0E27_8)*f*(Pi**3)*(rho**2)*(mr1**2)/(3.0_8*N*(rho0**2)*(lambda**4))
   ssa = 1.0_8 
!   ext = ext*N*1000.0_8
   phaseInd = 1
   LG = (/0.0, 0.5/)/(/2.0*1.0+1.0, 2.0*2.0+1.0/)
   phaseFunc(1) = new_PhaseFunction(LG,status=status)
   table = new_PhaseFunctionTable(phaseFunc(1:1),key=(/0.0/),tableDescription="Rayleigh Scattering", status=status)
   call finalize_PhaseFunction(phaseFunc(1))
   if(.not. stateIsFailure(status)) call setStateToSuccess(status)

  end subroutine calc_RayleighScattering
  !-----------------------------------------------------------------------------------------
end module opticalProperties   
