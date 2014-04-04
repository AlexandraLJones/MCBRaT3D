! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

!! Extra Modules added to store, write, and read a file containing surfaceProperty information
!! Adapted from routines read_Domain and write_Domain from opticalProperties.mod
!! Maxwell Smith (maxwell.a.smith@gmail.com)
!! Research Assistant
!! University of Illinois at Urbana-Champaign, Department of Atmospheric Sciences
!! 10 August 2012, REVISION: 0 LAST EDIT: 10 AUG 2012

! $Revision: 2.1 $, $Date: 2012/09/25 00:19:45 $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/surfaceProperties.f95 $
  ! This module represent surface  models in which the reflectance 
  !   (the "bi-directional reflectance distribution function," or BRDF) can be modeled with
  !    just a few parameters. These parameters can vary with horizontal position.
  ! The BRDF is a function of the incoming and outgoing polar cosines and azimuthal 
  !   angles.
  ! This particular module implements Lambertian surface reflectance but 
  !    can provide a template for more complicated surface BRDFs. 
  ! The horizontal coordinate system is local to the surface reflection object. 
  !   In particular, it's up to the user to guarantee that the coordinate system
  !   used to define the surface model is the same as the one used to define 
  !   the atmosphere. 
! --------------------------------------------
module surfaceProperties
  use ErrorMessages
  use numericUtilities
  use netcdf
  implicit none
  private
  !------------------------------------------------------------------------------------------
  ! Module constants
  !------------------------------------------------------------------------------------------
  integer, parameter :: numberOfParameters = 1
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type surfaceDescription
    private
    real(8), dimension(:),       pointer :: xPosition, yPosition
    real(8), dimension(:, :, :), pointer :: BRDFParameters
  end type surfaceDescription
  
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  !  Surface properties can vary with location or be constant. 
  ! 
  interface new_SurfaceDescription
    module procedure newSurfaceDescriptionXY, newSurfaceUniform
  end interface ! new_SurfaceDescription
  
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: surfaceDescription
  public :: new_SurfaceDescription, copy_surfaceDescription, finalize_SurfaceDescription, &
            isReady_surfaceDescription, &
            computeSurfaceReflectance
  !! addition of extra I/O modules, see header for more info
  public :: getInfo_SurfaceDescription, write_SurfaceDescription, read_SurfaceDescription

contains  
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create new surface description objects
  !------------------------------------------------------------------------------------------
  function newSurfaceDescriptionXY(surfaceParameters, xPosition, yPosition, status) &
    result(thisSurfaceDescription)
    real(8), dimension(:, :, :), intent(in   ) :: surfaceParameters
    real(8), dimension(:),       intent(in   ) :: xPosition, yPosition
    type(ErrorMessage),       intent(inout) :: status
    type(surfaceDescription)                :: thisSurfaceDescription 

    integer :: numX, numY
    
    ! Checks : array sizes 
    numX = size(xPosition); numY = size(yPosition)
    if(size(surfaceParameters, 1) /= numberOfParameters) &
      call setStateToFailure(status, "new_SurfaceDescription: Wrong number of parameters supplied for surface BRDF.") 
    if(size(surfaceParameters, 2) /= numX - 1 .or. &
       size(surfaceParameters, 3) /= numY - 1)     &
       call setStateToFailure(status, "new_SurfaceDescription: position vector(s) are incorrect length.")
    if(any(xPosition(2:) - xPosition(:numX-1) <= 0.0_8) .or. &
       any(yPosition(2:) - yPosition(:numY-1) <= 0.0_8))     &
      call setStateToFailure(status, "new_SurfaceDescription: positions must be unique, increasing.")

    ! Check to ensure that surface property parameter values make sense
    !   For a Lambertian surface (our example) the reflectance must be between 0 and 1, incl. 
    if(any(surfaceParameters(1, :, :) < 0.) .or. &
       any(surfaceParameters(1, :, :) > 1.)) &
      call setStateToFailure(status, "new_SurfaceDescription: surface reflectance must be between 0 and 1")
       
       
    if(.not. stateIsFailure(status)) then
      allocate(thisSurfaceDescription%xPosition(numX), &
               thisSurfaceDescription%yPosition(numY), &
               thisSurfaceDescription%BRDFParameters(numberOfParameters, numX - 1, numY - 1)) 
      thisSurfaceDescription%xPosition(:) = xPosition(:)
      thisSurfaceDescription%yPosition(:) = yPosition(:)  
      thisSurfaceDescription%BRDFParameters(:, :, :) = surfaceParameters(:, :, :)
      call setStateToSuccess(status) 
    end if 
  end function newSurfaceDescriptionXY
  !------------------------------------------------------------------------------------------
  function newSurfaceUniform(surfaceParameters, status) result(thisSurfaceDescription)
    real, dimension(:), intent(in   ) :: surfaceParameters
    type(ErrorMessage), intent(inout) :: status
    type(surfaceDescription)          :: thisSurfaceDescription
    !
    ! Define a horizontally uniform surface
    !
    real(8), dimension(2), parameter :: xPosition = (/ 0., huge(surfaceParameters) /), &
                                     yPosition = (/ 0., huge(surfaceParameters) /)
    real(8), dimension(numberOfParameters, 1, 1) :: surfaceParams
 
    if(size(surfaceParameters) /= numberOfParameters) then
      call setStateToFailure(status, "new_SurfaceDescription: Wrong number of parameters supplied for surface BRDF.") 
    else 
      surfaceParams(:, 1, 1) = surfaceParameters(:)
      thisSurfaceDescription = &
        newSurfaceDescriptionXY(surfaceParams, xPosition, yPosition, status)
    end if 
    
  end function newSurfaceUniform
  !------------------------------------------------------------------------------------------
  ! Compute surface reflectance at a given x-y location
  !------------------------------------------------------------------------------------------
  pure function computeSurfaceReflectance(thisSurfaceDescription, xPos, yPos, &
                                     incomingMu, outgoingMu, incomingPhi, outgoingPhi) result(surfaceReflectance)
    type(surfaceDescription), intent(in) :: thisSurfaceDescription
    real(8),		      intent(in) :: xPos, yPos
    real,                     intent(in) :: incomingMu, outgoingMu, incomingPhi, outgoingPhi
    real                                 :: surfaceReflectance

    real(8)    :: x0, xMax, y0, yMax
    integer :: xIndex, yIndex
    
    !
    ! Find the square on the surface that's doing the reflecting 
    !
    x0 = thisSurfaceDescription%xPosition(1)
    y0 = thisSurfaceDescription%yPosition(1)
    xMax = thisSurfaceDescription%xPosition(size(thisSurfaceDescription%xPosition))
    yMax = thisSurfaceDescription%yPosition(size(thisSurfaceDescription%yPosition))
    xIndex = findIndex(makePeriodic(xPos, x0, xMax), thisSurfaceDescription%xPosition)
    yIndex = findIndex(makePeriodic(yPos, y0, yMax), thisSurfaceDescription%yPosition)
    
    !
    ! Compute the reflectance from the BRDF parameters
    !   Developers should make no assumptions about the sign of value 
    !   of incomingMu and outgoingMu and should be prepared to treat 
    !   values of phi between 360 and -360. 
    !
    surfaceReflectance = R(thisSurfaceDescription%BRDFParameters(:, xIndex, yIndex))
    
  end function computeSurfaceReflectance
  !------------------------------------------------------------------------------------------
  !   Compute reflectance given a set of BRDF parameters 
  !     This is where the work is done and it's the main section of the code
  !     that needs changing when implementing a new BRDF. 
  !------------------------------------------------------------------------------------------
  pure function R(surfaceParameters) 
    real(8), dimension(numberOfParameters), intent(in) :: surfaceParameters
    real                                            :: R
    !
    ! This example implements a Lambertian surface albedo
    !
    R = surfaceParameters(1)
    
  end function R 
  !------------------------------------------------------------------------------------------
  !   Readiness
  !------------------------------------------------------------------------------------------
  !elemental 
    function isReady_surfaceDescription(thisSurface) 
    type(surfaceDescription), intent(in) :: thisSurface
    logical                              :: isReady_surfaceDescription
    
    isReady_surfaceDescription = associated(thisSurface%xPosition) .and. &
                                 associated(thisSurface%yPosition) .and. &
                                 associated(thisSurface%BRDFParameters)
  end function isReady_surfaceDescription
  !------------------------------------------------------------------------------------------
  !   Copy
  !------------------------------------------------------------------------------------------
  function copy_surfaceDescription(original) result(copy)
    type(surfaceDescription), intent( in) :: original
    type(surfaceDescription)              :: copy
    
    integer :: numX, numY 
    
    if(isReady_surfaceDescription(original)) then 
      numX = size(original%xPosition); numY = size(original%yPosition)
      allocate(copy%xPosition(numX), &
               copy%yPosition(numY), &
               copy%BRDFParameters(numberOfParameters, numX - 1, numY - 1)) 
      copy%xPosition(:) = original%xPosition(:)
      copy%yPosition(:) = original%yPosition(:)  
      copy%BRDFParameters(:, :, :) = original%BRDFParameters(:, :, :)
    end if
  end function copy_surfaceDescription

  !------------------------------------------------------------------------------------------
  !   Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_surfaceDescription(thisSurface)
    type(surfaceDescription), intent(out) :: thisSurface
    
    if(associated(thisSurface%xPosition)) deallocate(thisSurface%xPosition)
    if(associated(thisSurface%yPosition)) deallocate(thisSurface%yPosition)
    if(associated(thisSurface%BRDFParameters)) &
                                          deallocate(thisSurface%BRDFParameters)

  end subroutine finalize_surfaceDescription
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Utility procedures private to the module 
  !------------------------------------------------------------------------------------------
  elemental function makePeriodic(a, aMin, aMax)
    real(8), intent(in) :: a, aMin, aMax
    real             :: makePeriodic
    !
    ! Ensure that a position is within domain when the boundary conditions are periodic
    !
    ! makePeriodic = aMin + mod(a - aMin, aMax - aMin)
    ! if(makePeriodic < aMin) makePeriodic = aMax - abs(makePeriodic - aMin)
    makePeriodic = a
    do 
      if(makePeriodic <= aMax .and. makePeriodic > aMin) exit
      if(makePeriodic > aMax) then
        makePeriodic = makePeriodic - (aMax - aMin)
      else if (makePeriodic == aMin) then 
        makePeriodic = aMax 
      else
        makePeriodic = makePeriodic + (aMax - aMin)
      end if 
    end do
  end function makePeriodic
  !------------------------------------------------------------------------------------------



  !! Extra Modules added to store, write, and read a file containing surfaceProperty information
  !! see header for more info
   !------------------------------------------------------------------------------------------
  !  Getting information back from the object after it is loaded
  !------------------------------------------------------------------------------------------
  subroutine getInfo_SurfaceDescription(thisSurfaceDescription, numX, numY, numberOfParameters, &
                                xPosition, yPosition, BRDFParameters, status) 
    type(surfaceDescription),            intent(in   ) :: thisSurfaceDescription
    integer,                   optional, intent(  out) :: numX, numY, numberOfParameters
    real(8),    dimension(:),     optional, intent(  out) :: xPosition, yPosition
    real(8),    dimension(:,:,:), optional, intent(  out) :: BRDFParameters
    type(ErrorMessage),                  intent(inout) :: status
    
    ! Available information: The number of cells in the arrays, 
    !   and the locations of the x, y boundaries (which is one bigger). 
    ! The number of BRDF parameters is available (==1 for lambertian) and the
    ! array of parameters themselves at dimensionality(numberOfParameterss, xPosition, yPosition)
    
    ! --------------------
    ! Checks: description is valid, position arrays are big enough (the right size?)
    if(.not. isReady_SurfaceDescription(thisSurfaceDescription)) then
      call setStateToFailure(status, "getInfo_SurfaceDescription: surfaceDescription hasn't been initialized.")
    else
      ! Number of positions in each dimension for the property arrays. It's one shorter 
      !   than the position arrays, which describe the boundaries. 
      if(present(numX)) numX = size(thisSurfaceDescription%xPosition) -1
      if(present(numY)) numY = size(thisSurfaceDescription%yPosition) -1 
      if(present(numberOfParameters)) &
                 numberOfParameters = size(thisSurfaceDescription%BRDFParameters, 1)
      ! Location of boundaries in each dimension
      if(present(xPosition)) then 
        if(size(xPosition) /= size(thisSurfaceDescription%xPosition)) then
          call setStateToFailure(status, "getInfo_SurfaceDescription: vector for x positions is wrong length.")
        else
          xPosition(:) = thisSurfaceDescription%xPosition(:)
        end if
      end if
      if(present(yPosition)) then 
        if(size(yPosition) /= size(thisSurfaceDescription%yPosition)) then
          call setStateToFailure(status, "getInfo_SurfaceDescription: vector for y positions is wrong length.")
        else
          yPosition(:) = thisSurfaceDescription%yPosition(:)
        end if
      end if
      if(present(BRDFParameters)) then 
        if(size(BRDFParameters, 1) /= numberOfParameters) &
                  call setStateToFailure(status, "getInfo_SurfaceDescription: Wrong number of parameters supplied for surface BRDF.") 
        if(size(BRDFParameters, 2) /= size(thisSurfaceDescription%xPosition) - 1 .or. &
           size(BRDFParameters, 3) /= size(thisSurfaceDescription%yPosition) - 1)     &
           call setStateToFailure(status, "getInfo_SurfaceDescription: position vector(s) are incorrect length.")
        if( .not. stateIsFailure(status) ) &
          BRDFParameters(:,:,:) = thisSurfaceDescription%BRDFParameters(:,:,:)
      end if
            
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
  end subroutine getInfo_SurfaceDescription
  
    !------------------------------------------------------------------------------------------
  ! Storage and retrieval
  !------------------------------------------------------------------------------------------
  subroutine write_SurfaceDescription(thisSurfaceDescription, fileName, status)
    type(surfaceDescription),       intent(in   ) :: thisSurfaceDescription
    character(len = *),             intent(in   ) :: fileName
    type(ErrorMessage),             intent(inout) :: status
    
    ! Local variables
    integer                        :: i, j
        
    ! Netcdf-related local variables
    integer, dimension(9) :: ncStatus
    integer                :: ncFileId, xEdgeDimID, yEdgeDimId, nParamDimId,  &
                                        xGridDimId, yGridDimId, &
                                        ncVarId
    
    ! Checks: surface description is valid
    ncStatus(:) = nf90_NoErr
    if(.not. isReady_surfaceDescription(thisSurfaceDescription)) then
      call setStateToFailure(status, "write_SurfaceDescription: surface hasn't been initialized.") 
    else
      ncStatus( 1) = nf90_create(trim(fileName), nf90_Clobber, ncFileId)
      !
      ! Domain dimensions 
      !
      ncStatus( 2) = nf90_def_dim(ncFileId, "x-Edges", size(thisSurfaceDescription%xPosition),    xEdgeDimId) 
      ncStatus( 3) = nf90_def_dim(ncFileId, "y-Edges", size(thisSurfaceDescription%yPosition),    yEdgeDimId) 
      ncStatus( 4) = nf90_def_dim(ncFileId, "x-Grid",  size(thisSurfaceDescription%xPosition) - 1, xGridDimId) 
      ncStatus( 5) = nf90_def_dim(ncFileId, "y-Grid",  size(thisSurfaceDescription%yPosition) - 1, yGridDimId) 
      ncStatus( 6) = nf90_def_dim(ncFileId, "Parameters",size(thisSurfaceDescription%BRDFParameters, 1) , nParamDimId)
      !
      ! Domain variables
      ! 
      ncStatus( 7) = nf90_def_var(ncFileId, "x-Edges", nf90_double, xEdgeDimId, ncVarId) 
      ncStatus( 8) = nf90_def_var(ncFileId, "y-Edges", nf90_double, yEdgeDimId, ncVarId) 
      ncStatus( 9) = nf90_def_var(ncFileId, "BRDF", nf90_double, (/ nParamDimId, xGridDimId, yGridDimId /), ncVarId)
       
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "write_SurfaceDescription: error defining surface information") 
        
      if(.not. stateIsFailure(status)) then
        ncStatus( 1) = nf90_put_att(ncFileId, nf90_Global, "numberOf-BRDF-Parameters", size(thisSurfaceDescription%BRDFParameters, 1))
        ncStatus( 2) = nf90_put_att(ncFileId, nf90_Global, "numXcells", size(thisSurfaceDescription%xPosition) - 1)
        ncstatus( 3) = nf90_put_att(ncFileId, nf90_Global, "numYcells", size(thisSurfaceDescription%yPosition) - 1)
        
        if(any(ncStatus(1:3) /= nf90_NoErr)) &
          call setStateToFailure(status, "write_surfaceDescription: could not define global attributes")
      end if
      !
      ! All the attributes have been written, and the dimensions and variables defined. Now we'll take the
      !   file out of define mode and write the BRDF data.
      !
      ncStatus( 1) = nf90_EndDef(ncFileId)
      ncStatus( 2) = nf90_inq_varid(ncFileId, "x-Edges", ncVarID)
      ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, thisSurfaceDescription%xPosition)
      ncStatus( 4) = nf90_inq_varid(ncFileId, "y-Edges", ncVarID)
      ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, thisSurfaceDescription%yPosition)
      ncStatus( 6) = nf90_inq_varid(ncFileId, "BRDF", ncVarID)
      ncStatus( 7) = nf90_put_var(ncFileId, ncVarId, thisSurfaceDescription%BRDFParameters)
      
      if(any(ncStatus(1:7) /= nf90_NoErr)) &
        call setStateToFailure(status, "write_SurfaceProperties: error writing surface data") 
        
      if(stateIsFailure(status)) then 
          !
          ! File write has failed somehow - delete file and exit
          ! 
          ncstatus(1) = nf90_close(ncFileId) 
          open(20, file = trim(fileName))
          close(20, status = "delete")
      end if
        
      
      ncStatus(1) = nf90_close(ncFileId)
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
  end subroutine write_SurfaceDescription
  
  !------------------------------------------------------------------------------------------
  subroutine read_SurfaceDescription(fileName, thisSurfaceDescription, status)
    character(len = *),       intent(in   ) :: fileName
    type(surfaceDescription), intent(  out) :: thisSurfaceDescription
    type(ErrorMessage),       intent(inout) :: status
    
    ! Local variables
    integer            :: nXEdges, nYEdges, nParams
    real(8), dimension(:), &
           allocatable :: xEdges, yEdges
    real(8), dimension(:,:,:), &
           allocatable :: BRDFParameters
    
    
    ! Netcdf-related local variables
    integer, dimension(12) :: ncStatus
    integer                :: ncFileId, ncDimId, ncVarId
    
   
    ! There is nothing to check a priori
    ncStatus(:) = nf90_NoErr
    if(nf90_open(trim(fileName), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      call setStateToFailure(status, "read_SurfaceDescription: Can't open file " // trim(fileName)) 
    end if
    
    if(.not. stateIsFailure(status)) then 
      ncStatus( 1) = nf90_inq_dimid(ncFileId, "x-Edges", ncDimId) 
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nXEdges)
      ncStatus( 3) = nf90_inq_dimid(ncFileId, "y-Edges", ncDimId) 
      ncStatus( 4) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nYEdges)
      ncStatus( 5) = nf90_inq_dimid(ncFileId, "Parameters", ncDimId) 
      ncStatus( 6) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nParams)
      allocate(xEdges(nXEdges), yEdges(nYEdges))
      allocate(BRDFParameters(nParams,nxEdges-1, nyEdges-1))
      ncStatus( 7) = nf90_inq_varid(ncFileId, "x-Edges", ncVarId)
      ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, xEdges)
      ncStatus( 9) = nf90_inq_varid(ncFileId, "y-Edges", ncVarId)
      ncStatus(10) = nf90_get_var(ncFileId, ncVarId, yEdges)
      ncStatus(11) = nf90_inq_varid(ncFileId, "BRDF", ncVarId)
      ncStatus(12) = nf90_get_var(ncFileId, ncVarId, BRDFParameters)
       
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "read_SurfaceDescription: " // trim(fileName) // &
                               ". Error reading NetCDF file. Is it a surfaceDescription file?") 
    end if 
    
    if(.not. stateIsFailure(status)) then 
      !
      ! Create a new surfaceDescription using the initialization procedure. 
      !   
      !
      call finalize_SurfaceDescription(thisSurfaceDescription)
      thisSurfaceDescription = new_SurfaceDescription(BRDFParameters, xEdges, yEdges, status)
      
    end if 

    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine read_SurfaceDescription  
  !! End additional modules for I/O
  !-------------------------------------------------------------------------------
  
end  module surfaceProperties
