! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 25 $, $Date: 2010-01-04 03:28:27 +0000 (Mon, 04 Jan 2010) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Integrators/monteCarloRadiativeTransfer.f95 $
module monteCarloRadiativeTransfer
  use CharacterUtils
  use ErrorMessages
  use RandomNumbers
  use numericUtilities
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use surfaceProperties
  use monteCarloIllumination
  implicit none
  private
  
  !------------------------------------------------------------------------------------------
  ! Constants
  !------------------------------------------------------------------------------------------
  integer, parameter :: defaultMinForwardTableSize = 9001, &
                        defaultMinInverseTableSize = 9001
  real,    parameter :: defaultHybridPhaseFunWidth =  7., & 
                        maxHybridPhaseFunWidth     = 30.
  integer, parameter :: defOrdersOrigPhaseFunIntenCalcs = 0
  real,    parameter :: defaultZetaMin             =  0.3
  real,    parameter :: defaultMaxIntensityContrib =  huge(1.) 
  real,    parameter :: Pi = 3.14159265358979312
  integer, parameter :: defNumRecScatOrd           = -1 !default number of scattering orders to
                                                       !record for output
  
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  ! The public type
  !
  type integrator
    private
    ! Shortcuts - the status of the integrator could be obtained by looking to see if 
    !   arrays have been allocated, etc., but we also keep track with these flags
    logical                                 :: readyToCompute   = .false., &
                                               computeIntensity = .false. 
    ! Integrator specific parameters
    integer                                 :: minForwardTableSize = defaultMinForwardTableSize, &
                                               minInverseTableSize = defaultMinInverseTableSize
    ! -------------------------------------------------------------------------=                                           
    ! Algorithmic choices - these can be modified through specifyParameters()
    !
    ! Use ray-tracing? The alternative is max cross-section
    logical                                 :: useRayTracing = .true. 
    ! Use Russian roulette? 
    logical                                 :: useRussianRoulette = .true. 
    real                                    :: RussianRouletteW = 1. 
        
    ! -------------------------------------------------------------------------=                                           
    ! The atmosphere and surface 
    logical                                 :: xyRegularlySpaced = .false., &
                                                zRegularlySpaced = .false. 
    real(8)                                    :: deltaX = 0.0_8, deltaY = 0.0_8 , deltaZ = 0.0_8, &
                                               x0 = 0.0_8, y0= 0.0_8, z0 = 0.0_8
    real(8),    allocatable,dimension(:)           :: xPosition  
    real(8),    allocatable,dimension(:)           :: yPosition  
    real(8),    allocatable,dimension(:)           :: zPosition 
    real(8), allocatable, dimension(:,:,:)      :: temps 
    real                                    :: LW_flag = -1.
    ! Surface reflection BDRF
    logical                          :: useSurfaceBDRF = .false.
    type(surfaceDescription)         :: surfaceBDRF

    ! -------------------------------------------------------------------------=                                           
    !
    ! Direction cosines at which to compute intensity
    !
    real, dimension(:, :),    allocatable :: intensityDirections 
    
    ! -------------------------------------------------------------------------= 
    ! Variables for variance reduction for intensity calculations

    ! Build hybrid phase functions for intensity calculations for at least 
    !   one component

    logical :: useHybridPhaseFunsForIntenCalcs = .false. 
    real    :: hybridPhaseFunWidth             = defaultHybridPhaseFunWidth
    integer :: numOrdersOrigPhaseFunIntenCalcs = defOrdersOrigPhaseFunIntenCalcs

    ! Russian roulette for intensity 
    logical :: useRussianRouletteForIntensity = .False. 
    real    :: zetaMin                        = defaultZetaMin
    
    ! Redistribute large intenstity contributions in space
    logical                              :: limitIntensityContributions = .false. 
    real                                 :: maxIntensityContribution    = &
                                                        defaultMaxIntensityContrib
    real, dimension(:, :),       allocatable:: intensityExcess 

    ! -------------------------------------------------------------------------=                                           
    ! Output arrays
    
    real, dimension(:, :),    allocatable :: fluxUp , fluxDown , &
                                         fluxAbsorbed 
    real, dimension(:, :, :), allocatable :: volumeAbsorption 
    real, dimension(:, :, :), allocatable :: intensity 
    real, dimension(:, :, :, :), &
                              allocatable :: intensityByComponent 

    !output arrays for upward flux, downward flux, intensiies by number of scattering events
    real, dimension(:, :, :), allocatable :: fluxUpByScatOrd , fluxDownByScatOrd 
    real, dimension(:,:,:,:), allocatable :: intensityByScatOrd 

    !logical values to decide whether or not to record the number of scatterings
    logical                           :: recScatOrd = .false.
    integer                           :: numRecScatOrd = defNumRecScatOrd
    
  end type integrator
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: integrator
  public :: new_Integrator, copy_Integrator, isReady_Integrator, finalize_Integrator, &
            specifyParameters, computeRadiativeTransfer, reportResults !, getinfo_integrator
contains
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create new integrators, specifying 
  !   integrator-specific parameters, making copies. 
  !------------------------------------------------------------------------------------------
  function new_Integrator(atmosphere, status) result(new)
    type(domain),                       intent(in   ) :: atmosphere
    type(ErrorMessage),                 intent(inout) :: status
    type(integrator)                                  :: new
    ! 
    ! Create an object of type integrator. This holds all the information
    !   about the atmosphere and surface that's needed to do the 
    !   problem at hand. 
    
    ! Local variables
    integer :: numX, numY, numZ, numComponents
    real    :: deltaX, deltaY, deltaZ
      
    ! Sanity checks - be sure that the atmosphere is initialized. 
    !--------------------------------------------------------------
    ! Atmosphere
    !   How big is the atmospheric domain, and what are the cell boundaries? 
    !
    call getInfo_Domain(atmosphere, numX, numY, numZ, &
                        numberOfComponents = numComponents,  status = status)  
!print *, 'newIntegrator: got info 1'
    if(.not. stateIsFailure(status)) then   
      allocate(new%xPosition(numX + 1), new%yPosition(numY + 1), &
               new%zPosition(numZ + 1), new%temps(numX, numY, numZ)) 
      call getInfo_Domain(atmosphere,                       &
                          xPosition = new%xPosition, &
                          yPosition = new%yPosition, &
                          zPosition = new%zPosition, &
                          temps = new%temps,  status = status)
!print *, 'newIntegrator: got info 2'
      !
      ! Are the cells all the same size in the X and Y or in the Z directions? 
      !   This enables us to take some shortcuts 
      !
      new%x0 = new%xPosition(1)
      new%y0 = new%yPosition(1)
      new%z0 = new%zPosition(1)
      deltaX = new%xPosition(2)  - new%xPosition(1)
      deltaY = new%yPosition(2)  - new%yPosition(1)
      deltaZ = new%zPosition(2)  - new%zPosition(1)
      if(all(abs( (new%xPosition(2:) - new%xPosition(:numX)) -  deltaX ) <= &
             2. * spacing(new%xPosition(2:))) .and.                         &
         all(abs( (new%yPosition(2:) - new%yPosition(:numY)) -  deltaY ) <= &
             2. * spacing(new%yPosition(2:)))) then
        new%xyRegularlySpaced = .true. 
        new%deltaX            = deltaX
        new%deltaY            = deltaY
      end if
      if(all(abs( (new%zPosition(2:) - new%zPosition(:numZ)) -  deltaZ ) <= &
             spacing(new%zPosition(2:)))) then
        new%zRegularlySpaced = .true. 
        new%deltaZ           = deltaZ
      end if
    end if

    if(stateIsFailure(status)) then
      call setStateToFailure(status, "new_Integrator: Problems reading domain.")
      new%readyToCompute = .false. 
    else
      !--------------------------------------------------------------
      ! Create and zero the output arrays
      !
      allocate(new%fluxUp(numX, numY), new%fluxDown(numX, numY), new%fluxAbsorbed(numX, numY))
      new%fluxUp(:, :) = 0.; new%fluxDown(:, :) = 0.; new%fluxAbsorbed(:, :) = 0.
      allocate(new%volumeAbsorption(numX, numY, numZ))
      new%volumeAbsorption(:, :, :) = 0.
      !--------------------------------------------------------------
      ! Ready to go? 
      !
      new%readyToCompute = .true. 
      call setStateToSuccess(status)
    end if
  end function new_Integrator
  ! ------------------------------------------------------- 

  !------------------------------------------------------------------------------------------
  ! Computation: Routines to trace the photons and compute the radiative
  !   fluxes. There actual work might be done in other subroutines (to reflect 
  !   algorithmic choices, say). 
  !------------------------------------------------------------------------------------------
  subroutine computeRadiativeTransfer(thisIntegrator,thisDomain, randomNumbers, incomingPhotons,& 
		numPhotonsPerBatch, numPhotonsProcessed, status, option2)
    type(integrator),           intent(inout) :: thisIntegrator
    type(domain),		intent(inout)    :: thisDomain
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(photonStream),         intent(inout) :: incomingPhotons  
    integer(8),                    intent(in   ) :: numPhotonsPerBatch
    integer(8),                    intent(  out) :: numPhotonsProcessed
    type(ErrorMessage),         intent(inout) :: status
    integer, dimension(:,:,:), optional, intent(out) :: option2

    !
    ! Monte Carlo "integrator" to compute flux up at the top boundary, flux down at the
    !   bottom boundary, and colum absorption. (Absorption is calculated separately, 
    !   and should agree with the difference in boundary fluxes.) 
    ! Intensity at a set of zenith angle cosine and azimuth angles at the top and bottom boundaries
    !   may be computed using local estimation. 
    ! The bottom boundary is a Lambertian surface with default albedo 0. That 
    !   value may be changed with a call to specifyParameters(). 
    ! This routine calls one of two solvers: computeRT_MaxCrossSection for maximum
    !   cross-section and computeRT_PhotonTracing for simple ray tracing. 
    
    ! Local variables
    integer :: numIntensityDirections, numX, numY, numZ, numComponents, j, k, d,p
    real, dimension(:, :), allocatable &
            :: numPhotonsPerColumn

    ! Sanity checks
    if(.not. isReady_integrator(thisIntegrator)) then
      call setStateToFailure(status, "computeRadiativeTransfer: problem not completely specified.")
    else
      numX = size(thisIntegrator%xPosition) - 1
      numY = size(thisIntegrator%yPosition) - 1
      numZ = size(thisIntegrator%zPosition) - 1
      call getInfo_Domain(thisDomain,  numberOfComponents=numComponents, status=status)
      if(thisIntegrator%computeIntensity) & 
        numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
      
      ! Zero output arrays
      thisIntegrator%fluxUp(:, :)       = 0.
      thisIntegrator%fluxDown(:, :)     = 0.
      thisIntegrator%fluxAbsorbed(:, :) = 0.
      thisIntegrator%volumeAbsorption(:, :, :) &
                                        = 0.
      !
      ! We want the intensity arrays to be zero'd even if we're not computing new 
      !   intensity values
      !
      if(allocated(thisIntegrator%intensity)) thisIntegrator%intensity(:, :, :) = 0. 
      if(allocated(thisIntegrator%intensityByComponent)) &
                                               thisIntegrator%intensityByComponent(:, :, :, :) = 0. 
      if(allocated(thisIntegrator%intensityExcess)) &
                                               thisIntegrator%intensityExcess(:, :) = 0. 

      !
      ! We want the arrays by scattering order to be zero'd even if we're not computing
      !    new intensity values
      !
      if(allocated(thisIntegrator%fluxUpByScatOrd)) thisIntegrator%fluxUpByScatOrd(:, :, :)   = 0.
      if(allocated(thisIntegrator%fluxDownByScatOrd)) &
                                               thisIntegrator%fluxDownByScatOrd(:, :, :) = 0.
        
      if(allocated(thisIntegrator%intensityByScatOrd)) &
                                               thisIntegrator%intensityByScatOrd(:, :, :,  :) = 0.
   

      
      !
      ! Compute tablulated forward and inverse phase functions
      !   Forward phase functions are only needed if we're going to compute intensity
      !
      call tabulateInversePhaseFunctions(thisDomain, thisIntegrator%minInverseTableSize, status)
      if(thisIntegrator%computeIntensity .and. .not. stateIsFailure(status)) then
        call tabulateForwardPhaseFunctions(thisDomain, thisIntegrator%minForwardTableSize, &
                                           thisIntegrator%useHybridPhaseFunsForIntenCalcs,   &
                                           thisIntegrator%hybridPhaseFunWidth, status)
      end if
            
      !------------------------------------------------------------------------------
      !
      ! Compute radiative transfer for this photon batch 
      !
      if(.not. stateIsFailure(status)) &
         call computeRT(thisIntegrator, thisDomain, randomNumbers, incomingPhotons, numPhotonsPerBatch, numPhotonsProcessed, status)

      if(thisIntegrator%computeIntensity .and. &
         thisIntegrator%limitIntensityContributions) then 
        !
        ! Redistribute large intensity contributions in space
        !
        do j  = 0, numComponents
          do d = 1, numIntensityDirections 
            if(thisIntegrator%intensityExcess(d, j) > 0.) then 
              
              thisIntegrator%intensity(:, :, d) =  thisIntegrator%intensity(:, :, d) + &
                (    thisIntegrator%intensityByComponent(:, :, d, j) /                 &
                 sum(thisIntegrator%intensityByComponent(:, :, d, j)) ) * thisIntegrator%intensityExcess(d, j)
              
!              if(thisIntegrator%recScatOrd) then
!                do p = 0, thisIntegrator%numRecScatOrd
!                  thisIntegrator%intensityByScatOrd(:,:,d, p) =   &
!                    thisIntegrator%intensityByScatOrd(:,:,d,p) +  &
!                    (    thisIntegrator%intensityByComponent(:, :, d, j) /              &
!                     sum(thisIntegrator%intensityByComponent(:, :, d, j)) ) * thisIntegrator%intensityExcess(d,j)
!                end do
!              end if 
              thisIntegrator%intensityByComponent(:, :, d, j) =   & 
                thisIntegrator%intensityByComponent(:, :, d, j) + &
                (    thisIntegrator%intensityByComponent(:, :, d, j) /                 &
                 sum(thisIntegrator%intensityByComponent(:, :, d, j)) ) * thisIntegrator%intensityExcess(d, j)
            end if
          end do 
        end do
      end if
      
      !------------------------------------------------------------------------------
      !
      ! Normalization - compute the average number of photons incident on each column
      !
      allocate(numPhotonsPerColumn(numX, numY))
      
      if(thisIntegrator%xyRegularlySpaced) then
        numPhotonsPerColumn(:, :) = numPhotonsProcessed / real(numX * numY)
!PRINT *, "average numPhotonsPerColumn =", sum(numPhotonsPerColumn)/real (numX * numY)
      else
        forall(j = 1:numY)
          ! Relative area of each column 
          numPhotonsPerColumn(:, j) = ( (thisIntegrator%yPosition(j+1) - thisIntegrator%yPosition(j)) *       &
                                        (thisIntegrator%xPosition(2:) - thisIntegrator%xPosition(1:numX)) ) / & 
                                      ( (thisIntegrator%xPosition(numX+1) - thisIntegrator%xPosition(1)) *    &
                                        (thisIntegrator%yPosition(numY+1) - thisIntegrator%yPosition(1)) ) 
        end forall
        ! Now the number of photons incident per column
        numPhotonsPerColumn(:, :) = numPhotonsPerColumn(:, :) * numPhotonsProcessed  
      end if 
      
      ! 
      ! Normalize fluxes by average number of photons per column
      !
      thisIntegrator%fluxUp      (:, :) = thisIntegrator%fluxUp      (:, :) / numPhotonsPerColumn(:, :)
      thisIntegrator%fluxDown    (:, :) = thisIntegrator%fluxDown    (:, :) / numPhotonsPerColumn(:, :)
      thisIntegrator%fluxAbsorbed(:, :) = thisIntegrator%fluxAbsorbed(:, :) / numPhotonsPerColumn(:, :) 

!      if(thisIntegrator%recScatOrd) then
!        forall (p = 0:thisIntegrator%numRecScatOrd)
!          thisIntegrator%fluxUpByScatOrd (:, :, p) = thisIntegrator%fluxUpByScatOrd(:, :, p) / numPhotonsPerColumn(:,:)
!          thisIntegrator%fluxDownByScatOrd(:,:,p) = thisIntegrator%fluxDownByScatOrd(:,:,p)  / numPhotonsPerColumn(:,:)
!        end forall
!      end if 
      !
      ! Normalize absorption profile by cell depth
      !
      forall(k = 1:numZ)  &
        thisIntegrator%volumeAbsorption(:, :, k) =      &
          thisIntegrator%volumeAbsorption(:, :, k) /    &
               (numPhotonsPerColumn(:, :) * (thisIntegrator%zPosition(k+1) - thisIntegrator%zPosition(k))*1000.0 ) ! since physical positions are provided in km and totalFlux is calculatted or provided in Wm^-2 we need to use a factor of 1000 to get the units to make sense for flux divergence. result will be Wm^-3
  
      !
      ! Intensity is also normalized by the average number of photons per column. 
      !
      if(thisIntegrator%computeIntensity) then
        forall(d = 1:numIntensityDirections)
          thisIntegrator%intensity(:, :, d) =  thisIntegrator%intensity(:, :, d) / numPhotonsPerColumn(:, :)
        end forall
!PRINT *, SIZE(thisIntegrator%intensityByComponent,1), SIZE(thisIntegrator%intensityByComponent,2), SIZE(thisIntegrator%intensityByComponent,3), SIZE(thisIntegrator%intensityByComponent,4)
!PRINT *, SIZE(numPhotonsPerColumn,1), SIZE(numPhotonsPerColumn,2)
        forall(d = 1:numIntensityDirections, j = 1:numComponents)
          thisIntegrator%intensityByComponent(:, :, d, j) =  &
                                               thisIntegrator%intensityByComponent(:, :, d, j) / &
                                                                                   numPhotonsPerColumn(:, :)
        end forall
!PRINT *, "average column intensity normalized by photonsPerCol =", sum(thisIntegrator%intensity(:,:,1))/real (numX * numY)
!        if(thisIntegrator%recScatOrd) then
!          forall(d = 1:numIntensityDirections, p = 0:thisIntegrator%numRecScatOrd)
!            thisIntegrator%intensityByScatOrd(:, :, d, p) =  &
!                                                 thisIntegrator%intensityByScatOrd(:, :, d, p) / &
!                                                                                   numPhotonsPerColumn(:, :)
!          end forall
!        end if
      end if
      deallocate(numPhotonsPerColumn)  
    end if
  end subroutine computeRadiativeTransfer
  !------------------------------------------------------------------------------------------
  subroutine computeRT(thisIntegrator, thisDomain, randomNumbers, incomingPhotons, &
                                numPhotonsPerBatch, numPhotonsProcessed, status, option2)
    type(integrator),           intent(inout) :: thisIntegrator
    type(domain),               intent(in)    :: thisDomain
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(photonStream),         intent(inout) :: incomingPhotons
    integer(8),                    intent(in   ) :: numPhotonsPerBatch
    integer(8),                    intent(  out) :: numPhotonsProcessed
    type(ErrorMessage),         intent(inout) :: status
    integer, dimension(:,:,:), optional, intent(out) :: option2
    !
    ! Implements a standard ray-tracing Monte Carlo algorthm or 
    !   the Marchuk (1980) maximum cross-section algorithm. The extinction
    !   along each path segment is scaled by the maximum extinction within the domain. 
    !   Scattering events are then classified as "mathemtical" or "physical" by comparing
    !   a random number to the ratio of local to maximum extinction. 

        
    ! Local variables
    real :: mu, phi
    real(8) :: x0, y0, z0, xMax, yMax, zMax,xPos, yPos, zPos, xNew, yNew, zNew, remainder, albedo
    real :: tauToTravel, photonWeight, scatteringAngle, tauAccumulated, ssa, maxExtinction
    real :: initialMu, initialPhi
    logical :: useRayTracing, useMaxCrossSection, scatterThisEvent
    integer :: xIndex, yIndex, zIndex, numX, numY, numZ, numComps, component, phaseFunctionIndex
    integer(8) :: nPhotons
    integer :: i, p, scatteringOrder 
    integer :: nBad
    real, dimension(3) :: directionCosines
    !
    ! Variables related to intensity calculations
    !
    integer            :: numIntensityDirections, current
    real, dimension(:), allocatable ::  contributions
    integer, dimension(:), allocatable :: xIndexF, yIndexF
    real(8), allocatable, dimension(:,:,:) :: totalExt
    real(8), allocatable, dimension(:,:,:,:) :: cumExt, singleScattAlbedo
    integer, allocatable, dimension(:,:,:,:) :: phaseFuncI
    type(matrix), allocatable, dimension(:)  :: inversePhaseFuncs
    
    ! ---------------------------------------------------------------------------------------
     call getInfo_Domain(thisDomain, numberOfComponents=numComps, numX=numX, numY=numY, numZ=numZ, status=status)
    allocate(totalExt(1:numX,1:numY,1:numZ))
    allocate(cumExt(1:numX,1:numY,1:numZ,1:numComps))
    allocate(singleScattAlbedo(1:numX,1:numY,1:numZ,1:numComps))
    allocate(phaseFuncI(1:numX,1:numY,1:numZ,1:numComps))
    allocate(inversePhaseFuncs(1:numComps))

    call getInfo_Domain(thisDomain,albedo=albedo,                                   &
                        totalExt=totalExt, cumExt=cumExt, ssa=singleScattAlbedo,           &
                        phaseFuncI=phaseFuncI, inversePhaseFuncs=inversePhaseFuncs, status=status)
!    option2=0
    useRayTracing = thisIntegrator%useRayTracing; useMaxCrossSection = .not. useRayTracing
    scatterThisEvent = .true. 
    if(useMaxCrossSection) &
      maxExtinction = maxval(totalExt(:, :, :))
    x0 = thisIntegrator%x0; xMax = thisIntegrator%xPosition(size(thisIntegrator%xPosition))
    y0 = thisIntegrator%y0; yMax = thisIntegrator%yPosition(size(thisIntegrator%yPosition))
    z0 = thisIntegrator%z0; zMax = thisIntegrator%zPosition(size(thisIntegrator%zPosition))
    if(thisIntegrator%computeIntensity) then
      numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
      allocate(contributions(numIntensityDirections), &
                xIndexF(numIntensityDirections), yIndexF(numIntensityDirections))
    end if 
    !
    ! Begin loop over photons
    !
    nPhotons = 0; nBad = 0
CALL getNextPhoton(incomingPhotons, xPos, yPos, zPos, mu, phi, status, current)
!PRINT *, 'Starting with Photon', current
    photonLoop: DO WHILE (nPhotons .lt. numPhotonsPerBatch)
!PRINT *, nPhotons
      if(.not. morePhotonsExist(incomingPhotons)) EXIT photonLoop! This means we've used all the photons
      call getNextPhoton(incomingPhotons, xPos, yPos, zPos, mu, phi, status)
!      if(zPos .ge. 36)PRINT *, zPos
      if(stateIsFailure(status)) EXIT photonLoop
      scatteringOrder = 0
!PRINT *, 'mu=', mu, ' phi=', phi
      directionCosines(:) = makeDirectionCosines(mu, phi)
      photonWeight = 1. 
      nPhotons = nPhotons + 1
      !
      ! Incoming xPos, yPos are between 0 and 1. 
      !   Translate these positions to the local domain
      !
      XIndex = 1; yIndex = 1; zIndex = 1
!PRINT *, 'xPos before=', xPos
      xPos = x0 + xPos * (xMax - x0) 
!PRINT *, 'xPos after=', xPos
      yPos = y0 + yPos * (yMax - y0)
      call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
      if(thisIntegrator%zRegularlySpaced)then
        zPos = z0 + zPos * (zMax - z0)  
        call findZIndex(thisIntegrator, zPos, zIndex)    
      else
!        zIndex = 1-(zPos-FLOOR(zPos))+((SIZE(thisIntegrator%zPosition)-1)*zPos)
!        zPos = thisIntegrator%zPosition(zIndex) + (zPos-FLOOR(zPos))* (thisIntegrator%zPosition(zIndex+1)-thisIntegrator%zPosition(zIndex)) ! convert the zPos to one that works for an irregularly spaced grid. This line must follow and not preceed the zPos= line above.
         
         remainder = (zPos-z0)*numZ - FLOOR((zPos-z0)*numZ) !These lines added 11/4/2013 to replace the above calculations that seemed inaccuarate for irregularly spaced vertical levels
         zIndex = MIN(FLOOR((zPos-z0)*numZ)+1, numZ)
         zPos = thisIntegrator%zPosition(zIndex) + remainder*(thisIntegrator%zPosition(zIndex+1)-thisIntegrator%zPosition(zIndex))
      end if

!if (zPos > 0.0_8)then
!  option2(xIndex,yIndex,Zindex)=option2(xIndex,yIndex,Zindex)+1
!end if

!write(14,"(I7, 2X, 3I5)") nPhotons, xIndex, yIndex, zIndex
!PRINT *, xIndex, yIndex, zIndex
!if(zIndex .ge. 36)PRINT *, zIndex, zPos, thisIntegrator%zRegularlySpaced

      if(thisIntegrator%LW_flag > 0)then ! if we are doing a LW simulation we want to calculate the emmission contribution to the radiance and the absorbed flux
        if(zPos > 0.0_8)then
	   thisIntegrator%fluxAbsorbed(xIndex, yIndex) = thisIntegrator%fluxAbsorbed(xIndex, yIndex) - 1.0
           thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) = thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) - 1.0   !These 2 lines added 1/17/13 to account for photon emission by the medium in the accurate accounting of flux divergence, ultimately
        endif

         if (thisIntegrator%computeIntensity)then
           if (zPos .eq. 0.0_8)then
            call computeIntensityContribution(thisIntegrator, thisDomain, photonWeight, &
                                                xPos,   yPos,   zPos,         &
                                                xIndex, yIndex, zIndex,       &
                                                directionCosines, 0,  &
                                               randomNumbers, scatteringOrder, &
                                                contributions, xIndexF(:), yIndexf(:))
           else 
            call computeIntensityContribution(thisIntegrator, thisDomain, photonWeight, &
                                                xPos,   yPos,   zPos,         &
                                                xIndex, yIndex, zIndex,       &
                                                directionCosines, -1,  &
                                               randomNumbers, scatteringOrder, &
                                                contributions, xIndexF(:), yIndexf(:))
           end if
!PRINT *, "COntribution to radiance from photon ", nPhotons, " is ", contributions
!           if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then
!            !only record scattering order if it is within bounds
!              forall(i = 1:numIntensityDirections)
!                thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder)   = &
!                thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder) + &
!                                                                                    contributions(i)
!              end forall
!           end if
           forall(i = 1:numIntensityDirections)
                thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) =  &
                  thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) + contributions(i)
                thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) = &
                  thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) + contributions(i)  ! as of right now this last line attributes these contributions to the 0th component aka the surafce, incorrectly-not sure how to resolve this, can we give it a value of -1? should we just drop these lines? should I add in this functionality for th atmosphere? should we randomly choose one of the components to ahve done the scattering?
           end forall
         end if
        end if


      !
      ! Loop over orders of scattering
      !
      scatteringLoop: do
        !
        ! The optical distance we need to travel. 
        !   It's possible for the random number generator to produce exactly 0; 
        !   we set a lower bound. 
        !
        tauToTravel = -log(max(tiny(tauToTravel), getRandomReal(randomNumbers)))
        if(useRayTracing) then 
          !
          ! Ray tracing  - travel until we have accumulated enough extinction
          !
          call accumulateExtinctionAlongPath(thisDomain, directionCosines, &
                                             xPos, yPos, zPos, xIndex, yIndex, zIndex, &
                                             tauAccumulated, tauToTravel)            
          if(tauAccumulated < 0.) nBad = nBad + 1 
          if(tauAccumulated < 0.) cycle photonLoop 
        else 
          !
          ! Max cross-section: move the photon to the new location according to the maximum extinction
          !
          xPos = makePeriodic(xPos + directionCosines(1) * tauToTravel/maxExtinction, x0, xMax)
          yPos = makePeriodic(yPos + directionCosines(2) * tauToTravel/maxExtinction, y0, yMax)
          zPos =              zPos + directionCosines(3) * tauToTravel/maxExtinction
        end if 
        
        if(zPos >= zMax) then
          !
          ! The photon has gone out the top of the domain.  
          !   Add to the upward flux at that point, then start on a new photon
          !
          if(useMaxCrossSection) then 
            !
            ! Trace backwards to domain top if using max cross-section 
            ! 
            xPos = makePeriodic(xPos - directionCosines(1) * abs((zPos - zMax)/directionCosines(3)), x0, xMax)
            yPos = makePeriodic(yPos - directionCosines(2) * abs((zPos - zMax)/directionCosines(3)), y0, yMax)
            call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
          end if 
          
          thisIntegrator%fluxUp(xIndex, yIndex) = thisIntegrator%fluxUp(xIndex, yIndex) + photonWeight
!          if (scatteringOrder .eq. 0 .and. thisIntegrator%computeIntensity)then
!            call computeIntensityContribution(thisIntegrator, photonWeight, &
!                                                xPos,   yPos,   zPos,         &
!                                                xIndex, yIndex, zIndex,       &
!                                                directionCosines, -1,  &
!                                                randomNumbers, scatteringOrder, &
!                                                contributions, xIndexF(:), yIndexf(:))
!            if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then
!            !only record scattering order if it is within bounds              
!              forall(i = 1:numIntensityDirections)
!                thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder)   = &
!                  thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder) + &
!                                                                                    contributions(i)
!              end forall
!            end if
!              forall(i = 1:numIntensityDirections)
!                thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) =  &
!                  thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) + contributions(i)
!                thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) = &
!                  thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) + contributions(i)  ! as of right now this last line attributes these contributions to the 0th component aka the surafce, incorrectly-not sure how to resolve this, can we give it a value of -1? should we just drop these lines? should I add in this functionality for th atmosphere? should we randomly choose one of the components to ahve done the scattering?
!              end forall 
!          end if

!          if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then 
!          !only record scattering order if it is within bounds
!              thisIntegrator%fluxUpByScatOrd(xIndex, yIndex, scatteringOrder)   = &
!                        thisIntegrator%fluxUpByScatOrd(xIndex, yIndex, scatteringOrder) + photonWeight
!            
!          end if
          cycle photonLoop

        else if(zPos <= z0 + spacing(z0)) then
          !
          ! The photon is at the surface. Add it to the surface flux. 
          !   directionCosines(3) will always be non-zero since the photon was traveling vertically. 
          !
          if(useMaxCrossSection) then 
            !
            ! Trace backwards to domain base if using max cross-section 
            ! 
            xPos = makePeriodic(xPos - directionCosines(1) * abs((zPos - z0)/directionCosines(3)), x0, xMax)
            yPos = makePeriodic(yPos - directionCosines(2) * abs((zPos - z0)/directionCosines(3)), y0, yMax)
            call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
          end if 
          zIndex = 1
          zPos = z0 + spacing(z0)
          thisIntegrator%fluxDown(xIndex, yIndex) = thisIntegrator%fluxDown(xIndex, yIndex) + photonWeight
!          if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then 
!          !only record scattering order if it is within bounds
!              thisIntegrator%fluxDownByScatOrd(xIndex, yIndex, scatteringOrder) = &
!                        thisIntegrator%fluxDownByScatOrd(xIndex, yIndex, scatteringOrder) + photonWeight
!          end if
          
          !the photon scatters off ground
          !increment scattering order
          scatteringOrder = scatteringOrder + 1

          !
          ! Compute new photon weight and a new direction to travel. 
          !   Save the old directions in case we're using a BDRF
          !
          initialMu  = directionCosines(3)
          if(initialMu <= 1.) then 
            initialPhi = acos(directionCosines(1)/sqrt(1. - initialMu**2))
          else
            initialPhi = 0. 
          end if 
          do 
            !
            ! Ensure that new trajectory has at least some vertical component - otherwise 
            !   the trajectory can get stuck forever if the lowest layer has no extinction
            !
            mu  = sqrt(getRandomReal(randomNumbers))
            if(abs(mu) > 2 * tiny(mu)) exit
          end do 
          phi = 2 * Pi * getRandomReal(randomNumbers)
          !
          ! New weight from the surface reflectance. 
          !
          if(thisIntegrator%useSurfaceBDRF) then 
            photonWeight = photonWeight * &
                           computeSurfaceReflectance(thisIntegrator%surfaceBDRF, xPos, yPos, &
                                                     initialMu, mu, initialPhi, phi)
          else
            !   Special case: Lambertian surface
            photonWeight = photonWeight * albedo
          end if
          if(photonWeight <= tiny(photonWeight)) cycle photonLoop
          directionCosines(:) = makeDirectionCosines(mu, phi)
          !
          ! Add contribution of surface reflection to intensity
          !
          if(thisIntegrator%computeIntensity) then
            call computeIntensityContribution(thisIntegrator, thisDomain, photonWeight, &
                                              xPos,   yPos,   zPos,         & 
                                              xIndex, yIndex, zIndex,       &
                                              directionCosines, 0,          &
                                              randomNumbers, scatteringOrder, &
                                              contributions, xIndexF(:), yIndexF(:))     
!PRINT *, "COntribution to radiance from photon ", nPhotons, " is ", contributions, " for scattering order ", scatteringorder, " off the surface."         
!            if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then
!            !only record scattering order if it is within bounds              
!              forall(i = 1:numIntensityDirections)
!                thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder)   = &
!                  thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder) + &
!                                                                                    contributions(i)
!              end forall
!            end if
              forall(i = 1:numIntensityDirections)
                thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) =  &
                  thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) + contributions(i)
                thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) = &
                  thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) + contributions(i)
              end forall 
          end if 
        else 
          !
          ! Scattering event. 
          !
          
          ! Max cross-section - test for "Physical scattering event" 
          if(useMaxCrossSection) scatterThisEvent = & 
            getRandomReal(randomNumbers) < totalExt(xIndex, yIndex, zIndex)/maxExtinction

          if(useRayTracing .or. scatterThisEvent ) then
            scatteringOrder = scatteringOrder + 1
            !
            ! Time for the next scattering event.
            ! The photon might have accumulated the right amount of extinction 
            !   just on the boundary of a cell. If it's traveling in the positive
            !   direction, the upper boundary is considered in the next cell, and 
            !   that cell might not have any extinction in it. 
            ! In this (very rare) case we take the smallest possible step backwards, 
            !   since the optical properties aren't defined in cells with no extinction. 
            !  Don't need to worry about max cross-section, since we can't have a physical 
            !    scattering event if extinction is 0. 
            ! 
            ! We need to enforce periodicity here 
            !
      
            if(totalExt(xIndex, yIndex, zIndex) <= 0.) then
              if(xPos - thisIntegrator%xPosition(xIndex) <= 0.0_8 .and. directionCosines(1) > 0. ) then
                xPos = xPos - spacing(xPos) 
                xIndex = xIndex - 1
                if (xIndex <= 0) then
                  xIndex = size(thisIntegrator%xPosition) - 1
                  xPos =  thisIntegrator%xPosition(xIndex) 
                  xPos = xpos - 2. * spacing(xPos)
                end if
              end if 
              if(yPos - thisIntegrator%yPosition(yIndex) <= 0.0_8 .and. directionCosines(2) > 0. ) then
                yPos = yPos - spacing(yPos) 
                yIndex = yIndex - 1
                if (yIndex <= 0) then
                  yIndex = size(thisIntegrator%yPosition) - 1
                  yPos =  thisIntegrator%xPosition(yIndex) 
                  yPos = ypos - 2. * spacing(yPos)  ! There was a typo here. I changed it from xPos to yPos on 11/3/13
                end if 
              end if 
              if(zPos - thisIntegrator%zPosition(zIndex) <= 0. .and. directionCosines(3) > 0. ) then
                zPos = zPos - spacing(zPos) 
                zIndex = zIndex - 1
              end if 
              !
              ! No need to worry about edge case - photons won't have come from below the domain
              !
            end if 
            !
            !   Figure out which component does the extinction, 
            !   and compute the new direction and weight of the photon. 
            !
            component = findIndex(getRandomReal(randomNumbers), &
                                  (/ 0.0_8, cumExt(xIndex, yIndex, zIndex, :) /))    
            !
            ! Absorption 
            !
            ssa = singleScattAlbedo(xIndex, yIndex, zIndex, component)
            if(ssa < 1.0_8) then 
              thisIntegrator%fluxAbsorbed(xIndex, yIndex) =   &
                thisIntegrator%fluxAbsorbed(xIndex, yIndex)             + photonWeight * (1.0_8 - ssa)
              thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) =   &
                thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) + photonWeight * (1.0_8 - ssa)
              photonWeight = photonWeight * ssa
            end if 
  
            !
            ! Compute the contribution to intensity from this scattering event if need be
            !
            if(thisIntegrator%computeIntensity) then
              call computeIntensityContribution(thisIntegrator, thisDomain, photonWeight, &
                                                xPos,   yPos,   zPos,         & 
                                                xIndex, yIndex, zIndex,       &
                                                directionCosines, component,  &
                                                randomNumbers, scatteringOrder, &
                                                contributions, xIndexF(:), yIndexf(:))   
!PRINT *, "COntribution to radiance from photon ", nPhotons, " is ", contributions, " for scattering order ", scatteringorder, " off the atmosphere."           
    
              forall (i = 1:numIntensityDirections)
                thisIntegrator%intensity(xIndexF(i), yIndexf(i), i) =  &
                  thisIntegrator%intensity(xIndexF(i), yIndexf(i), i) + contributions(i)
                thisIntegrator%intensityByComponent(xIndexF(i), yIndexf(i), i, component) = &
                  thisIntegrator%intensityByComponent(xIndexF(i), yIndexf(i), i, component) + contributions(i)
              end forall

!              if(thisIntegrator%recScatOrd .and. scatteringOrder <= thisIntegrator%numRecScatOrd) then
!              !only record scattering order if it is within bounds
!                forall(i = 1:numIntensityDirections)
!                  thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder)   = &
!                    thisIntegrator%intensityByScatOrd(xIndexF(i), yIndexF(i), i, scatteringOrder) + &
!                                                                                    contributions(i)
!                end forall
!              end if
            end if  ! end of local estimation for radiance contribution
    
            !
            ! "Russian roulette" 
            !
            if(thisIntegrator%useRussianRoulette .and. photonWeight < thisIntegrator%RussianRouletteW/2. ) then 
              if(getRandomReal(randomNumbers) >= photonWeight/thisIntegrator%RussianRouletteW) then 
                photonWeight = 0. 
              else
                photonWeight = thisIntegrator%RussianRouletteW
              end if 
            end if
            if(photonWeight <= tiny(photonWeight)) cycle photonLoop
            !
            ! Scattering - look up the scattering angle
            !
            phaseFunctionIndex = phaseFuncI(xIndex, yIndex, zIndex, component)
            scatteringAngle = computeScatteringAngle(getRandomReal(randomNumbers), &
                                      inversePhaseFuncs(component)%values(:, phaseFunctionIndex)) 
            call next_direct(randomNumbers, cos(scatteringAngle), directionCosines)
          end if
        end if
      end do scatteringLoop
    end do photonLoop

    
    if(thisIntegrator%computeIntensity) deallocate(contributions, xIndexF, yIndexF)
    deallocate(totalExt,cumExt,singleScattAlbedo,phaseFuncI,inversePhaseFuncs)
    !
    ! Status is only set when getting photons, so if that hasn't failed we know we're ok. 
    !
    if(.not. stateIsFailure(status) .and. nPhotons > 0) then
      call setStateToCompleteSuccess(status, "computeRadiativeTransfer: finished with photons")
      numPhotonsProcessed = nPhotons

    else if(nPhotons == 0) then 
        call setStateToFailure(status, "computeRadiativeTransfer: Didn't process any photons.")
    else 
        call setStateToFailure(status, "computeRadiativeTransfer: Error.")
    end if
!PRINT *, "Numphotons Processed=", nPhotons, " Numphotons not processed=", nbad, " sum of intensity across domain = ", sum(thisIntegrator%intensity)    
  end subroutine computeRT
  !------------------------------------------------------------------------------------------
  ! Reporting 
  !------------------------------------------------------------------------------------------
  subroutine reportResults(thisIntegrator,                             &
                           meanFluxUp, meanFluxDown, meanFluxAbsorbed, &
                               fluxUp,     fluxDown,     fluxAbsorbed, &
                           absorbedProfile, volumeAbsorption,          &
                           meanIntensity, intensity,                   &
                           meanFluxUpByScatOrd, meanFluxDownByScatOrd, &
                           fluxUpByScatOrd, fluxDownByScatOrd,         &
                           meanIntensityByScatOrd,intensityByScatOrd,  &
                           status)
    type(integrator),                      intent(in   ) :: thisIntegrator
    real,                        optional, intent(  out) :: meanFluxUp,          meanFluxDown, meanFluxAbsorbed
    real, dimension(:, :),       optional, intent(  out) ::     fluxUp,              fluxDown,     fluxAbsorbed
    real, dimension(:),          optional, intent(  out) :: meanFluxUpByScatOrd, meanFluxDownByScatOrd
    real, dimension(:, :, :),    optional, intent(  out) ::     fluxUpByScatOrd,     fluxDownByScatOrd
    real, dimension(:),          optional, intent(  out) :: absorbedProfile
    real, dimension(:, :, :),    optional, intent(  out) :: volumeAbsorption
    real, dimension(:),          optional, intent(  out) :: meanIntensity
    real, dimension(:, :, :),    optional, intent(  out) ::     intensity
    real, dimension(:, :),       optional, intent(  out) :: meanIntensityByScatOrd
    real, dimension(:, :, :, :), optional, intent(  out) ::     intensitybyScatOrd
    type(ErrorMessage),                    intent(inout) :: status
    
    !
    ! Local variables
    !
    integer :: direction, numDirections, numColumns, p
    
    !
    ! This integrator computes fluxes at the top and bottom boundaries, the total 
    !   column absorption, and the absorption profile.  Users can ask for any of 
    !   the pixel level fluxes as a domain mean and/or on a column-by-column basis. 
    !
    numColumns = size(thisIntegrator%fluxUp) 
    
    ! Domain averaged fluxes
    !
    if(present(meanFluxUp))   meanFluxUp       = sum(thisIntegrator%fluxUp)       / numColumns
    if(present(meanFluxDown)) meanFluxDown     = sum(thisIntegrator%fluxDown)     / numColumns
    if(present(meanFluxAbsorbed)) &
                              meanFluxAbsorbed = sum(thisIntegrator%fluxAbsorbed) / numColumns

    !domain averaged fluxes by scattering order
!    if(thisIntegrator%recScatOrd) then
!      if(present(meanFluxUpByScatOrd)) then
!        if(size(meanFluxUpByScatOrd) /= thisIntegrator%numRecScatOrd+1) then
!          call setStateToFailure(status, "reportResults: meanFluxUpByScatOrd is the wrong size")
!        else
!          forall(p = 0:thisIntegrator%numRecScatOrd)
!            meanFluxUpByScatOrd(p)   = sum(thisIntegrator%fluxUpByScatOrd(:,:,p))   / numColumns
!          end forall
!        end if
!      end if
!      if(present(meanFluxDownByScatOrd)) then
!        if(size(meanFluxDownByScatOrd) /= thisIntegrator%numRecScatOrd+1) then
!          call setStateToFailure(status, "reportResults: meanFluxDownByScatOrd is the wrong size")
!        else
!          forall (p = 0:thisIntegrator%numRecScatOrd)
!            meanFluxDownByScatOrd(p) = sum(thisIntegrator%fluxDownByScatOrd(:,:,p)) / numColumns
!          end forall
!        end if
!      end if
!    end if
    !
    ! Pixel-by-pixel fluxes
    !
    if(present(fluxUp)) then
      if(any((/ size(fluxUp, 1), size(fluxUp, 2) /) /= &
             (/ size(thisIntegrator%fluxUp, 1), size(thisIntegrator%fluxUp, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxUp array is the wrong size")
      else
        fluxUp(:, :) = thisIntegrator%fluxUp(:, :)
      end if 
    end if 
    
    if(present(fluxDown)) then
      if(any((/ size(fluxDown, 1), size(fluxDown, 2) /) /= &
             (/ size(thisIntegrator%fluxDown, 1), size(thisIntegrator%fluxDown, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxDown array is the wrong size")
      else
        fluxDown(:, :) = thisIntegrator%fluxDown(:, :)
      end if 
    end if 
    
    if(present(fluxAbsorbed)) then
      if(any((/ size(fluxAbsorbed, 1), size(fluxAbsorbed, 2) /) /= &
             (/ size(thisIntegrator%fluxAbsorbed, 1), size(thisIntegrator%fluxAbsorbed, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxAbsorbed array is the wrong size")
      else
        fluxAbsorbed(:, :) = thisIntegrator%fluxAbsorbed(:, :)
      end if 
    end if 
    
!    if(thisIntegrator%recScatOrd) then 
!      if(present(fluxUpByScatOrd)) then
!        if(any((/ size(fluxUpByScatOrd, 1), size(fluxUpByScatOrd, 2), size(fluxUpByScatOrd, 3)/) /= &
!               (/ size(thisIntegrator%fluxUpByScatOrd, 1), size(thisIntegrator%fluxUpByScatOrd, 2),&
!                  size(thisIntegrator%fluxUpbyScatOrd, 3) /))) then
!          call setStateToFailure(status, "reportResults: fluxUpByScatOrder array is the wrong size")
!        else
!          fluxUpByScatOrd(:, :, :) = thisIntegrator%fluxUpByScatOrd(:, :, :)
!        end if
!      end if

!      if(present(fluxDownByScatOrd)) then
!        if(any((/ size(fluxDownByScatOrd, 1), size(fluxDownByScatOrd, 2), size(fluxDownByScatOrd, 3)/) /= &
!               (/ size(thisIntegrator%fluxDownByScatOrd, 1), size(thisIntegrator%fluxDownByScatOrd, 2), & 
!                  size(thisIntegrator%fluxDownByScatOrd, 3) /))) then
!          call setStateToFailure(status, "reportResults: fluxDownByScatOrd array is the wrong size")
!        else
!          fluxDownByScatOrd(:, :, :) = thisIntegrator%fluxDownByScatOrd(:, :, :)
!        end if
!      end if
!    end if
    
    !
    ! Absorption - heating rate profile and volume absorption
    !
    if(present(absorbedProfile)) then
      if (size(absorbedProfile) /= size(thisIntegrator%volumeAbsorption, 3)) then
        call setStateToFailure(status, "reportResults: absorbedProfile array is the wrong size")
      else
        absorbedProfile(:) = sum(sum(thisIntegrator%volumeAbsorption(:, :, :), dim = 1), dim = 1) / numColumns
      end if 
    end if 
     
    if(present(volumeAbsorption)) then
      if(any((/ size(volumeAbsorption, 1), size(volumeAbsorption, 2), size(volumeAbsorption, 3) /) /= &
             (/ size(thisIntegrator%volumeAbsorption, 1), size(thisIntegrator%volumeAbsorption, 2),   &
                size(thisIntegrator%volumeAbsorption, 3)  /))) then
        call setStateToFailure(status, "reportResults: volumeAbsorption array is the wrong size")
      else
        volumeAbsorption(:, :, :) = thisIntegrator%volumeAbsorption(:, :, :)
      end if 
    end if 
    !
    ! Domain-averaged intensity
    !
    if(present(meanIntensity)) then
      if(.not. allocated(thisIntegrator%intensity)) then 
        call setStateToFailure(status, "reportResults: intensity information not available") 
      else if (size(thisIntegrator%intensity, 3) /= size(meanIntensity)) then 
        call setStateToFailure(status, "reportResults: requesting mean intensity in the wrong number of directions.") 
      else
        numDirections = size(thisIntegrator%intensity, 3)
        forall(direction = 1:numDirections)
          meanIntensity(direction) = sum(thisIntegrator%intensity(:, :, direction)) / numColumns
        end forall
      end if 
    end if
    
    !
    ! Pixel-by-pixel intensity
    !
    if(present(intensity)) then
      if(.not. allocated(thisIntegrator%intensity)) then 
        call setStateToFailure(status, "reportResults: intensity information not available") 
      else if (any( (/ size(               intensity, 1), size(               intensity, 2), &
                       size(               intensity, 3) /) /=                               &
                    (/ size(thisIntegrator%intensity, 1), size(thisIntegrator%intensity, 2), &
                       size(thisIntegrator%intensity, 3) /) ))  then 
        call setStateToFailure(status, "reportResults: intensity array has wrong dimensions.") 
      else
        intensity(:, :, :) = thisIntegrator%intensity(:, :, :)
!PRINT *, "mean intensity in reportResults =", sum(intensity(:,:,1))/numColumns
      end if
    end if

!    if(present(meanIntensityByScatOrd)) then
!      if(.not. associated(thisIntegrator%intensityByScatOrd)) then
!        call setStateToFailure(status, "reportResults: intensityByScatOrd information not available")
!      else if (size(thisIntegrator%intensityByScatOrd, 3) /= size(meanIntensity) .or. &
!              size(thisIntegrator%intensityByScatOrd, 4) /= thisIntegrator%numRecScatOrd) then
!        call setStateToFailure(status, "reportResults: requesting mean intensityByScatOrd in the wrong number of directions.")
!      else
!        numDirections = size(thisIntegrator%intensity, 3)
!        forall(direction = 1:numDirections, p=0:thisIntegrator%numRecScatOrd)
!          meanIntensityByScatOrd(direction, p) = sum(thisIntegrator%intensityByScatOrd(:, :, direction, p)) / numColumns
!        end forall
!      end if
!    end if

!    if(present(intensityByScatOrd)) then
!      if(.not. associated(thisIntegrator%intensityByScatOrd)) then
!        call setStateToFailure(status, "reportResults: intensityByScatOrd information not available")
!      else if (any( (/ size(intensityByScatOrd, 1), size(intensityByScatOrd, 2), &
!                       size(intensityByScatOrd, 3), size(intensityByScatOrd, 4)  /) /=       &
!                    (/ size(thisIntegrator%intensityByScatOrd, 1), size(thisIntegrator%intensityByScatOrd, 2), &
!                       size(thisIntegrator%intensityByScatOrd, 3), size(thisIntegrator%intensityByScatOrd, 4) /)))  then
!        call setStateToFailure(status, "reportResults: intensity array has wrong dimensions.")
!      else
!        intensityByScatOrd(:, :, :, :) = thisIntegrator%intensityByScatOrd(:, :, :, :)
!      end if
!    end if


    
    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine reportResults 
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  subroutine specifyParameters(thisIntegrator, surfaceBDRF,           &
                               minForwardTableSize, minInverseTableSize,             &
                               intensityMus, intensityPhis, computeIntensity,        &
                               useRayTracing,  useRussianRoulette,                   &
                               useRussianRouletteForIntensity, zetaMin,              &
                               useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                               numOrdersOrigPhaseFunIntenCalcs,                      &
                               limitIntensityContributions, maxIntensityContribution,&
                               recScatOrd, numRecScatOrd, LW_flag, numComps,                   &
                               status)
    type(integrator),   intent(inout) :: thisIntegrator
    type(surfaceDescription), &
              optional, intent(in   ) :: surfaceBDRF
    integer,  optional, intent(in   ) :: minForwardTableSize, minInverseTableSize
    real, dimension(:), &
              optional, intent(in   ) :: intensityMus, intensityPhis
    logical,  optional, intent(in   ) :: computeIntensity, useRayTracing, useRussianRoulette
    logical,  optional, intent(in   ) :: useRussianRouletteForIntensity
    real,     optional, intent(in   ) :: zetaMin
    logical,  optional, intent(in   ) :: useHybridPhaseFunsForIntenCalcs
    real,     optional, intent(in   ) :: hybridPhaseFunWidth
    integer,  optional, intent(in   ) :: numOrdersOrigPhaseFunIntenCalcs
    logical,  optional, intent(in   ) :: limitIntensityContributions
    real,     optional, intent(in   ) :: maxIntensityContribution
    logical,  optional, intent(in   ) :: recScatOrd
    integer,  optional, intent(in   ) :: numRecScatOrd, numComps
    real,     optional, intent(in   ) :: LW_flag
    type(ErrorMessage), intent(inout) :: status
    
    !
    !  Set integrator-specific parameter. The parameters might or might not need to be persistent. 
    !    If the medium is known to be highly variable within some vertical subset
    !    of the domain, for example, we might store the upper and lower boundaries 
    !    of this region, then do maximum cross section within it and photon tracing
    !    outside. 
    !  This might also be useful outside the core integrators, i.e. those that work 
    !    on a single batch of photons.
    
    ! Local variables
    integer :: i, p

    ! ------------------------------------------------------------------------          
    ! 
    ! Sanity checks for input variables
    !
    if(present(surfaceBDRF)) then
      if(.not. isReady_surfaceDescription(surfaceBDRF)) & 
        call setStateToFailure(status, "specifyParameters: surface description isn't valid.") 
    end if 
    
    !
    ! Table sizes must be bigger than 0; could put a minimum bound
    !
    if( present(minForwardTableSize)) then
      if(minForwardTableSize < defaultMinForwardTableSize) &
        call setStateToWarning(status, "specifyParameters: minForwardTableSize less than default. Value ignored.") 
    end if
    if( present(minInverseTableSize)) then
      if(minInverseTableSize < defaultMinInverseTableSize) &
        call setStateToWarning(status, "specifyParameters: minInverseTableSize less than default. Value ignored.") 
    end if
    
    if(present(hybridPhaseFunWidth)) then
      if(hybridPhaseFunWidth > maxhybridPhaseFunWidth .or. hybridPhaseFunWidth < 0.)        &
        call setStateToWarning(status,                                                         &
                               "specifyParameters: hybridPhaseFunWidth out of range (0 to " // &
                               trim(intToChar(int(maxhybridPhaseFunWidth))) // "degrees)."  // &
                               "Using default (" // trim(intToChar(int(defaultHybridPhaseFunWidth))) // ")")
    end if
    if(present(numOrdersOrigPhaseFunIntenCalcs)) then
      if(numOrdersOrigPhaseFunIntenCalcs < 0)                                                        &
        call setStateToWarning(status,                                                               &
                               "specifyParameters: numOrdersExactPhaseFunIntenCalcs less than 0." // &
                               "Using default (" // trim(intToChar(defOrdersOrigPhaseFunIntenCalcs)) // ")")
    end if
    
    if(present(maxIntensityContribution)) then
      if(maxIntensityContribution <= 0.)                                                      &
        call setStateToWarning(status,                                                       &
                               "specifyParameters: maxIntensityContribution <= 0. Value is unchanged.")
        
    end if 


    ! intensity direction arrays should be the same length
    !  both or neither supplied
    ! intensityMus must be between -1 and 1; can't be identically 0
    ! intensityPhis must be between 0 and 360 
    !
    if(present(intensityMus) .neqv. present(intensityPhis)) &
      call setStateToFailure(status, "specifyParameters: Both or neither of intensityMus and intensityPhis must be supplied") 
    if(present(intensityMus)) then
      if(size(intensityMus) /= size(intensityPhis)) &
        call setStateToFailure(status, "specifyParameters: intensityMus, intensityPhis must be the same length.") 
      if(any(intensityMus < -1.) .or. any(intensityMus > 1.)) &
        call setStateToFailure(status, "specifyParameters: intensityMus must be between -1 and 1") 
      if(any(abs(intensityMus) < tiny(intensityMus))) &
        call setStateToFailure(status, "specifyParameters: intensityMus can't be 0 (directly sideways)") 
      if(any(intensityPhis < 0.) .or. any(intensityPhis > 360.)) &
        call setStateToFailure(status, "specifyParameters: intensityPhis must be between 0 and 360") 
    end if
    
    !
    ! If someone specifies intensityDirections but trys to set computeIntensity to false we ignore them
    !
    if(present(computeIntensity)) then
      if(.not. computeIntensity .and. present(intensityMus))                                                             &
        call setStateToWarning(status, "specifyParameters: intensity directions *and* computeIntensity set to false." // &
                                       "Will compute intensity at given angles.") 
      if(computeIntensity .and. .not. present(intensityMus) .and. .not. allocated(thisIntegrator%intensityDirections)) &
        call setStateToFailure(status, "specifyParameters: Can't compute intensity without specifying directions.") 
    end if
    
    if(present(recScatOrd)) then
      if(.not. recScatOrd .and. present(numRecScatOrd)) then
        call setStateToWarning(status, "specifyParameters: set recScatOrd to false, but specified numRecScatOrders will override")
      else if (recScatOrd .and. .not. present(numRecScatOrd)) then
        call setStateToFailure(status, "specifyParameters: set recScatOrd to true, but did not provide number of orders to track.")
      end if
    end if
    if(present(numRecScatOrd)) then 
      if (numRecScatOrd < 0) then 
        call setStateToWarning(status, "specifyParameters: recScatOrd < 0. Will not record scattering orders.")
      end if
    end if
    ! ------------------------------------------------------------------------          
    if( .not. StateIsFailure(status)) then 
        if (present(surfaceBDRF)) then
        thisIntegrator%surfaceBDRF = copy_SurfaceDescription(surfaceBDRF)
        thisIntegrator%useSurfaceBDRF = .true. 
      end if 
      !
      ! Algorithmic choices
      !
      if(present(useRayTracing)) thisIntegrator%useRayTracing = useRayTracing
      if(present(minForwardTableSize)) &
        thisIntegrator%minForwardTableSize    = max(minForwardTableSize, defaultMinForwardTableSize)
      if(present(minInverseTableSize)) &
        thisIntegrator%minInverseTableSize    = max(minInverseTableSize, defaultMinInverseTableSize)

      if(present(useRussianRoulette)) &
        thisIntegrator%useRussianRoulette = useRussianRoulette
        
      !
      ! Russian roulette for intensity
      !
      if(present(useRussianRouletteForIntensity)) &
         thisIntegrator%useRussianRouletteForIntensity = useRussianRouletteForIntensity
      if(present(zetaMin)) then 
        if(zetaMin < 0.) then 
          call setStateToWarning(status, "specifyParameters: zetaMin must be >= 0. Value is unchanged.")
        else 
          thisIntegrator%zetaMin = zetaMin
          if(zetaMin > 1.) &
            call setStateToWarning(status, "specifyParameters: zetaMin > 1. That's kind of large.")
        end if 
      end if 

      !
      ! Hyrbid phase function for local estimation
      !
      if(present(useHybridPhaseFunsForIntenCalcs)) &
        thisIntegrator%useHybridPhaseFunsForIntenCalcs = useHybridPhaseFunsForIntenCalcs
      if(present(hybridPhaseFunWidth)) then
        if(hybridPhaseFunWidth > 0 .and. hybridPhaseFunWidth < maxhybridPhaseFunWidth) then
          thisIntegrator%hybridPhaseFunWidth = hybridPhaseFunWidth
        else
           thisIntegrator%hybridPhaseFunWidth = defaultHybridPhaseFunWidth
        end if
        !
        ! Need to re-tabulate the phase functions
        !
!        if(associated(thisIntegrator%tabulatedPhaseFunctions)) then 
!          do i = 1, size(thisIntegrator%tabulatedPhaseFunctions)
!            call finalize_Matrix(thisIntegrator%tabulatedPhaseFunctions(i))
!            deallocate(thisIntegrator%tabulatedPhaseFunctions)
!          end do
!        end if
      end if 
      if(present(numOrdersOrigPhaseFunIntenCalcs)) then
        if(numOrdersOrigPhaseFunIntenCalcs >= 0 ) then
          thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = numOrdersOrigPhaseFunIntenCalcs
        else
          thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = defOrdersOrigPhaseFunIntenCalcs
        end if
      end if 

      !
      ! Limited maximum local estimate
      !
      if(present(limitIntensityContributions)) & 
        thisIntegrator%limitIntensityContributions = limitIntensityContributions
      if(present(maxIntensityContribution)) then 
        if(maxIntensityContribution > 0.) &
          thisIntegrator%maxIntensityContribution = maxIntensityContribution
      end if
      !
      ! Intensity 
      !
      if(present(intensityMus)) then 
        if(allocated(thisIntegrator%intensityDirections))  deallocate(thisIntegrator%intensityDirections)
        allocate(thisIntegrator%intensityDirections(3, size(intensityMus)))

        if(allocated(thisIntegrator%intensity))            deallocate(thisIntegrator%intensity)
        allocate(thisIntegrator%intensity(size(thisIntegrator%xPosition)-1, &
                                          size(thisIntegrator%yPosition)-1, &
                                          size(intensityMus)))

        if(allocated(thisIntegrator%intensityByComponent)) deallocate(thisIntegrator%intensityByComponent) 
        allocate(thisIntegrator%intensityByComponent(size(thisIntegrator%xPosition)-1, &
                                                     size(thisIntegrator%yPosition)-1, &
                                                     size(intensityMus),               &
                                                     0:numComps))
        if(allocated(thisIntegrator%intensityByScatOrd)) deallocate(thisIntegrator%intensityByScatOrd)
        if(thisIntegrator%recScatOrd) then 
          allocate(thisIntegrator%intensityByScatOrd(size(thisIntegrator%xPosition)-1, &
                                                     size(thisIntegrator%yPosition)-1, &
                                                     size(intensityMus),&
                                                     0:numRecScatOrd                )   )
        endif

        forall(i = 1:size(intensityMus)) &
          thisIntegrator%intensitydirections(:, i) = makeDirectionCosines(intensityMus(i), &
                                                                          intensityPhis(i) * Pi/180.)
        thisIntegrator%computeIntensity = .true.
      end if 
      
      if(present(computeIntensity)) then 
        !
        ! If computeintensity is true we've already assured that the directions have been supplied at some point, 
        !    so there's nothing to do. 
        !
        if(.not. computeIntensity .and. .not. present(intensityMus)) then
          if(allocated(thisIntegrator%intensityDirections))  deallocate(thisIntegrator%intensityDirections)
          if(allocated(thisIntegrator%intensity))            deallocate(thisIntegrator%intensity)
          if(allocated(thisIntegrator%intensityByComponent)) deallocate(thisIntegrator%intensityByComponent)
          if(allocated(thisIntegrator%intensityByScatOrd))   deallocate(thisIntegrator%intensityByScatOrd)
          thisIntegrator%computeIntensity = .false.
        end if 
      end if 
            
      if(thisIntegrator%computeIntensity .and. thisIntegrator%limitIntensityContributions) then
        if(allocated(thisIntegrator%intensityExcess))      deallocate(thisIntegrator%intensityExcess)
        allocate(thisIntegrator%intensityExcess(size(thisIntegrator%intensityDirections, 2), &
                                                0:numComps))
      end if

      !allow recScatOrd to be set to true only if numRecScatOrd is set concurrently
      if(present(numRecScatOrd)) then 
        if(allocated(thisIntegrator%fluxUpByScatOrd))    deallocate(thisIntegrator%fluxUpByScatOrd)
        if(allocated(thisIntegrator%fluxDownByScatOrd))  deallocate(thisIntegrator%fluxDownByScatOrd)
        if(allocated(thisIntegrator%intensityByScatOrd)) deallocate(thisIntegrator%intensityByScatOrd)

        if(numRecScatOrd >= 0)  then
          thisIntegrator%recScatOrd=.true.
          thisIntegrator%numRecScatOrd = numRecScatOrd
              
          allocate(thisIntegrator%fluxUpByScatOrd   (size(thisIntegrator%xPosition)-1, &
                                                     size(thisIntegrator%yPosition)-1, &
                                                     0:numRecScatOrd                )   )
          allocate(thisIntegrator%fluxDownByScatOrd (size(thisIntegrator%xPosition)-1, &
                                                     size(thisIntegrator%yPosition)-1, &
                                                     0:numrecScatOrd                )   )
          if (thisIntegrator%computeIntensity) then
             allocate(thisIntegrator%intensityByScatOrd(size(thisIntegrator%intensity, 1), &
                                                        size(thisIntegrator%intensity, 2), &
                                                        size(thisIntegrator%intensity, 3), & 
                                                        0:numRecScatOrd                )   )
          end if
        else 
            thisIntegrator%recScatOrd=.false.
        end if
      end if

      if(present(recScatOrd)) then
        if(.not. recScatOrd .and. .not.present(numRecScatOrd)) then
          if(allocated(thisIntegrator%fluxUpByScatOrd))    deallocate(thisIntegrator%fluxUpByScatOrd)
          if(allocated(thisIntegrator%fluxDownByScatOrd))  deallocate(thisIntegrator%fluxDownByScatOrd)
          if(allocated(thisIntegrator%intensityByScatOrd)) deallocate(thisIntegrator%intensityByScatOrd)       
          thisIntegrator%recScatOrd = .false.
          thisIntegrator%numRecScatOrd = -1
        end if
      end if

      if(present(LW_flag))thisIntegrator%LW_flag=LW_flag  ! added by ALexandra Jones 4/4/12 to aid in calculation of emitted radiance

      call setStateToSuccess(status)
    end if 
      
       
    
  end subroutine specifyParameters 
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  function isReady_Integrator(thisIntegrator)
    type(integrator), intent( in) :: thisIntegrator
    logical                       :: isReady_Integrator
    
    isReady_Integrator = thisIntegrator%readyToCompute
  end function isReady_Integrator
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  function copy_Integrator(original) result(copy)
    type(integrator), intent(in) :: original
    type(integrator)             :: copy
    !
    ! Copy the state of one integrator to another
    !
    
    ! Local variables
    integer :: i
    
    ! Scalars
    copy%readyToCompute    = original%readyToCompute
    copy%computeIntensity  = original%computeIntensity
    
    
    if(isReady_surfaceDescription(original%surfaceBDRF)) &
      copy%surfaceBDRF     = copy_surfaceDescription(original%surfaceBDRF)
    copy%useSurfaceBDRF    = original%useSurfaceBDRF
    
    copy%xyRegularlySpaced = original%xyRegularlySpaced
    copy%zRegularlySpaced  = original%zRegularlySpaced
    
    copy%deltaX = original%deltaX; copy%deltaY = original%deltaY; copy%deltaZ = original%deltaZ
    copy%x0     = original%x0;     copy%y0     = original%y0;     copy%z0     = original%z0 
    
    copy%minForwardTableSize = original%minForwardTableSize
    copy%minInverseTableSize = original%minInverseTableSize
    !
    ! Algorithmic choices
    !
    copy%useRayTracing = original%useRayTracing
    
    copy%useRussianRoulette = original%useRussianRoulette
    copy%RussianRouletteW   = original%RussianRouletteW

    copy%useRussianRouletteForIntensity   = original%useRussianRouletteForIntensity
    copy%zetaMin                          = original%zetaMin
    
    !
    ! Evans phase function truncation
    !
    copy%useHybridPhaseFunsForIntenCalcs = original%useHybridPhaseFunsForIntenCalcs
    copy%hybridPhaseFunWidth             = original%hybridPhaseFunWidth
    copy%numOrdersOrigPhaseFunIntenCalcs = original%numOrdersOrigPhaseFunIntenCalcs
    

    copy%recScatOrd        = original%recScatOrd
    copy%numRecScatOrd     = original%numRecScatOrd
    !
    ! Location vectors
    !
    if(allocated(original%xPosition)) then
      allocate(copy%xPosition(size(original%xPosition)))
      copy%xPosition(:) = original%xPosition(:)
    end if

    if(allocated(original%yPosition)) then
      allocate(copy%yPosition(size(original%yPosition)))
      copy%yPosition(:) = original%yPosition(:)
    end if

    if(allocated(original%zPosition)) then
      allocate(copy%zPosition(size(original%zPosition)))
      copy%zPosition(:) = original%zPosition(:)
    end if
    
  
    !
    ! Intensity directions
    !
    if(allocated(original%intensityDirections)) then
      allocate(copy%intensityDirections(size(original%intensityDirections, 1), size(original%intensityDirections, 2)))
      copy%intensityDirections(:, :)  = original%intensityDirections(:, :)
    end if 
    
    !
    ! Output arrays
    !
    if(allocated(original%fluxUp)) then
      allocate(copy%fluxUp(size(original%fluxUp, 1), &
                           size(original%fluxUp, 2)))
      copy%fluxUp(:, :) = original%fluxUp(:, :)
    end if

    if(allocated(original%fluxDown)) then
      allocate(copy%fluxDown(size(original%fluxDown, 1), &
                             size(original%fluxDown, 2)))
      copy%fluxDown(:, :) = original%fluxDown(:, :)
    end if

    if(allocated(original%fluxAbsorbed)) then
      allocate(copy%fluxAbsorbed(size(original%fluxAbsorbed, 1), &
                                 size(original%fluxAbsorbed, 2)))
      copy%fluxAbsorbed(:, :) = original%fluxAbsorbed(:, :)
    end if
    
    if(allocated(original%volumeAbsorption)) then
      allocate(copy%volumeAbsorption(size(original%volumeAbsorption, 1), &
                                     size(original%volumeAbsorption, 2), &
                                     size(original%volumeAbsorption, 3)))
      copy%volumeAbsorption(:, :, :) = original%volumeAbsorption(:, :, :)
    end if
    
    if(allocated(original%intensity)) then
      allocate(copy%intensity(size(original%intensity, 1), &
                              size(original%intensity, 2), &
                              size(original%intensity, 3)))
      copy%intensity(:, :, :) = original%intensity(:, :, :)
    end if

    if(allocated(original%fluxUpByScatOrd)) then
       allocate(copy%fluxUpByScatOrd   (size(original%fluxUpByScatOrd, 1), &
                                        size(original%fluxUpByScatOrd, 2), &
                                        size(original%fluxUpByScatOrd, 3)))
       copy%fluxUpByScatOrd(:, :, :)       = original%fluxUpByScatOrd(:, :, :)
    end if
    if(allocated(original%fluxDownByScatOrd)) then   
       allocate(copy%fluxDownByScatOrd  (size(original%fluxDownByScatOrd, 1), &
                                         size(original%fluxDownByScatOrd, 2), &
                                         size(original%fluxDownByScatOrd, 3)))
       copy%fluxDownByScatOrd(:, :, :)     = original%fluxDownByScatOrd(:, :, :)
    end if
    if(allocated(original%intensityByScatOrd)) then   
       allocate(copy%intensityByScatOrd (size(original%intensityByScatOrd, 1), &
                                         size(original%intensityByScatOrd, 2), &
                                         size(original%intensityByScatOrd, 3), &
                                         size(original%intensityByScatOrd, 4)))
       copy%intensityByScatOrd(:, :, :, :) = original%intensityByScatOrd(:, :, :, :)
    end if


  end function copy_Integrator

  !------------------------------------------------------------------------------------------
  ! Finalization 
  !------------------------------------------------------------------------------------------
  subroutine finalize_Integrator(thisIntegrator)
    type(integrator), intent(out) :: thisIntegrator
    
    !
    ! Finalize by copying component values from a variable that's never been used
    !
    type(integrator)               :: pristineI
    
    ! Local variable
    integer :: i
    
    thisIntegrator%readyToCompute    =  pristineI%readyToCompute
    thisIntegrator%computeIntensity  =  pristineI%computeIntensity 
    
    thisIntegrator%xyRegularlySpaced = pristineI%xyRegularlySpaced 
    thisIntegrator%zRegularlySpaced  = pristineI%zRegularlySpaced
    
    thisIntegrator%useRayTracing              = pristineI%useRayTracing     
    thisIntegrator%useRussianRoulette         = pristineI%useRussianRoulette 
    thisIntegrator%RussianRouletteW           = pristineI%RussianRouletteW
    
    thisIntegrator%useRussianRouletteForIntensity = pristineI%useRussianRouletteForIntensity 
    thisIntegrator%zetaMin                        = pristineI%zetaMin 

    thisIntegrator%useHybridPhaseFunsForIntenCalcs = pristineI%useHybridPhaseFunsForIntenCalcs 
    thisIntegrator%hybridPhaseFunWidth             = pristineI%hybridPhaseFunWidth
    thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = pristineI%numOrdersOrigPhaseFunIntenCalcs 

    thisIntegrator%deltaX = 0.; thisIntegrator%deltaY = 0.; thisIntegrator%deltaZ = 0. 
    thisIntegrator%x0     = 0.; thisIntegrator%y0     = 0.; thisIntegrator%z0     = 0. 
    
    thisIntegrator%minForwardTableSize = pristineI%minForwardTableSize
    thisIntegrator%minInverseTableSize = pristineI%minInverseTableSize
    
    
    thisIntegrator%useSurfaceBDRF = pristineI%useSurfaceBDRF 
    call finalize_surfaceDescription(thisIntegrator%surfaceBDRF)

    thisIntegrator%recScatOrd   = pristineI%recScatOrd
    thisIntegrator%numrecScatOrd= pristineI%numRecScatOrd
    
    if(allocated(thisIntegrator%xPosition))          deallocate(thisIntegrator%xPosition)
    if(allocated(thisIntegrator%yPosition))          deallocate(thisIntegrator%yPosition)
    if(allocated(thisIntegrator%zPosition))          deallocate(thisIntegrator%zPosition)

    
    if(allocated(thisIntegrator%intensityDirections))    deallocate(thisIntegrator%intensityDirections)
    
    !
    ! Output arrays
    !
    if(allocated(thisIntegrator%fluxUp))             deallocate(thisIntegrator%fluxUp)
    if(allocated(thisIntegrator%fluxDown))           deallocate(thisIntegrator%fluxDown)
    if(allocated(thisIntegrator%fluxAbsorbed))       deallocate(thisIntegrator%fluxAbsorbed)
    if(allocated(thisIntegrator%volumeAbsorption))   deallocate(thisIntegrator%volumeAbsorption)
    if(allocated(thisIntegrator%intensity))          deallocate(thisIntegrator%intensity)
    if(allocated(thisIntegrator%fluxUpByScatOrd))    deallocate(thisIntegrator%fluxUpByScatOrd)
    if(allocated(thisIntegrator%fluxDownByScatOrd))  deallocate(thisIntegrator%fluxDownByScatOrd)
    if(allocated(thisIntegrator%intensityByScatOrd)) deallocate(thisIntegrator%intensityByScatOrd)
    if(allocated(thisIntegrator%intensityByComponent)) deallocate(thisIntegrator%intensityByComponent)

  end subroutine finalize_Integrator
  !------------------------------------------------------------------------------------------
  ! Functions for use inside the module 
  !------------------------------------------------------------------------------------------
 subroutine findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
! pure subroutine findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
    type(integrator), intent(in ) :: thisIntegrator
    real(8),             intent(in ) :: xPos, yPos
    integer,          intent(inout) :: xIndex, yIndex
    
    
    if(thisIntegrator%xyRegularlySpaced) then
      xIndex = min(int((xPos - thisIntegrator%x0)/thisIntegrator%deltaX) + 1, & 
                   size(thisIntegrator%xPosition)-1)
      yIndex = min(int((yPos - thisIntegrator%y0)/thisIntegrator%deltaY) + 1, &
                   size(thisIntegrator%yPosition)-1)
      !
      ! Insure against rounding errors
!PRINT *, 'findXindex xPos', xPos
      if(abs(thisIntegrator%xPosition(xIndex+1) - xPos) < spacing(xPos)) xIndex = xIndex + 1
      if(abs(thisIntegrator%yPosition(yIndex+1) - yPos) < spacing(yPos)) yIndex = yIndex + 1
      if(xIndex == size(thisIntegrator%xPosition)) xIndex = 1
      if(yIndex == size(thisIntegrator%yPosition)) yIndex = 1
    else
      xIndex = findIndex(xPos, thisIntegrator%xPosition, xIndex)
      yIndex = findIndex(yPos, thisIntegrator%yPosition, yIndex)
      if(abs(thisIntegrator%xPosition(xIndex) - xPos) < spacing(xPos)) xIndex = xIndex + 1
      if(abs(thisIntegrator%yPosition(yIndex) - yPos) < spacing(yPos)) yIndex = yIndex + 1
      if(xIndex .ge. size(thisIntegrator%xPosition)) xIndex = 1
      if(yIndex .ge. size(thisIntegrator%yPosition)) yIndex = 1
    end if
  end subroutine findXYIndicies
  !------------------------------------------------------------------------------------------
  pure subroutine findZIndex(thisIntegrator, zPos, zIndex)
    type(integrator), intent(in ) :: thisIntegrator
    real(8),             intent(in ) :: zPos
    integer,          intent(out) :: zIndex
    if(thisIntegrator%zRegularlySpaced) then
      zIndex = min(int((zPos - thisIntegrator%z0)/thisIntegrator%deltaZ) + 1, &
                   size(thisIntegrator%zPosition)-1)
      ! Insure against rounding errors
      if(abs(thisIntegrator%zPosition(zIndex+1) - zPos) < spacing(zPos)) zIndex = zIndex + 1
    else
      zIndex = findIndex(zPos, thisIntegrator%zPosition, zIndex)
    end if
  end subroutine findZIndex
  !------------------------------------------------------------------------------------------
  pure function computeScatteringAngle(randomDeviate, inversePhaseFunctionTable) result(scatteringAngle)
    !
    ! Linearly interpolate the scattering angle from a vector containing the 
    !   the angle as a function of the cumulative distribution (inverse phase function). 
    !   Recall that the first entry in the table is for CDF = 0 (hence angleIndex - 1) 
    ! 
    real,              intent (in ) :: randomDeviate
    real, dimension(:), intent(in ) :: inversePhaseFunctionTable
    real                            :: scatteringAngle
  
    ! Local variables
    integer :: angleIndex, numIntervals
    real    :: leftOver
    
    numIntervals = size(inversePhaseFunctionTable)
    angleIndex = int(randomDeviate * numIntervals) + 1
    if(angleIndex < numIntervals) then
      !
      ! Interpolate between entries in the inverse phase function table. 
      !   The first entry in the table is for CDF = 0 (hence angleIndex - 1) 
      !
      leftOver           = randomDeviate - real(angleIndex - 1)/real(numIntervals)
      scatteringAngle = (1. - leftOver) * inversePhaseFunctionTable(angleIndex) + &
                              leftOver  * inversePhaseFunctionTable(angleIndex + 1)
    else
      scatteringAngle = inversePhaseFunctionTable(numIntervals) 
    end if 
  end function computeScatteringAngle
  !------------------------------------------------------------------------------------------
  subroutine computeIntensityContribution(thisIntegrator, thisDomain, photonWeight,    &
                                          xPos,   yPos,   zPos,         & 
                                          xIndex, yIndex, zIndex,       &
                                          directionCosines, component,  &
                                          randomNumbers, scatteringOrder, &
                                          contributions, xIndexF, yIndexF)
    !
    ! Compute the contribution to the intensity in each direction 
    !   from a scattering event at xPos, xPos, zPos from a specified
    !   component. 
    !
    type(integrator),      intent(inout) :: thisIntegrator
    type(domain),             intent(in ) :: thisDomain
    real(8),                  intent(in   ) :: xPos,   yPos,   zPos
    real,                  intent(in   ) :: photonWeight
    integer,               intent(in   ) ::               xIndex, yIndex, zIndex
    real,    dimension(:), intent(in   ) :: directionCosines
    integer,               intent(in   ) :: component
    type(randomNumberSequence), &
                           intent(inout) :: randomNumbers ! Needed for Russian roulette, if used
    integer,               intent(in   ) :: scatteringOrder
    real,    dimension(:), intent(  out) :: contributions
    integer, dimension(:), intent(  out) :: xIndexF, yIndexF
    
    ! Local variables
    integer :: phaseFunctionIndex, numIntensityDirections, i
    real, dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: projections, scatteringAngles, phaseFunctionVals, tausToBoundary
    integer,  dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: zIndexF
    integer, allocatable, dimension(:,:,:,:) 	 :: phaseFuncI        
    type(ErrorMessage)				 :: status
    type(matrix), allocatable, dimension(:)	 :: tabulatedPhaseFunctions, tabulatedOrigPhaseFunctions
    !
    ! Variables for Russian Roulette as applied to intensity calculations
    !   Notation follows H. Iwabuchi, JAS 2006, Vol 63, pp 2324-2339
    !   contributions(:) corresponds to w_n zeta_n/Pi in his Eq 9
    !   We omit the cancelling factors of Pi that appear in his Eqs 9 and 10
    !
    real    :: tauFree, tauMax
    integer :: zIndexMax, numX, numY, numZ, numComps
    real(8)    :: xTemp, yTemp, zTemp, xPosI, yPosI, zPosI
    real, dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: normalizedPhaseFunc
   
    call getInfo_Domain(thisDomain, numX=numX, numY=numY, numZ=numZ, numberOfComponents=numComps, status=status)
    allocate(phaseFuncI(1:numX,1:numY,1:numZ,1:numComps))
    allocate(tabulatedPhaseFunctions(1:numComps))
    allocate(tabulatedOrigPhaseFunctions(1:numComps))
    call getInfo_Domain(thisDomain, phaseFuncI=phaseFuncI,tabPhase=tabulatedPhaseFunctions, &
		tabOrigPhase=tabulatedOrigPhaseFunctions, status=status)
 
    ! -------------------------------------
    numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
    zIndexMax = size(thisIntegrator%zPosition)
    xPosI = xPos; yPosI = yPos; zPosI = zPos
    
    xIndexF(:) = xIndex; yIndexF(:) = yIndex; zIndexF(:) = zIndex
    
    ! The photon weight as supplied has already been adjusted to account for the 
    !    albedo of the surface or scattering component
    !
    ! We rely on the calling routines setting component to 0 to indicate surface reflection. 
    !   There might be better ways to trap for that. 
    !
    if(component .eq. 0) then ! We're reflecting from the surface
      !
      ! This normalization may only work for Lambertian surfaces
      !   It's supposed to be the ratio of the BRDF to the albedo at 
      !   the incident solar zenith angle
      !
      normalizedPhaseFunc(:) = 1/Pi
    else if(component < 0)then  ! isotropic emission in the atmosphere
      normalizedPhaseFunc(:) = 1/(4.*Pi*abs(thisIntegrator%intensityDirections(3,:))) !not sure if this value is right
    else
      !
      ! Determine the scattering angle (in radians) from the photon direction to each intensity direction
      !   as the acos of dot product of the photon's direction cosine with the intensity direction's cosine
      !
      ! Numerics are such that once in a great while the projection can be slightly larger than 1. 
      !
      projections(:) = matmul(directionCosines, thisIntegrator%intensityDirections(:, :))
      where(abs(projections(:)) > 1.) projections(:) = sign(1., projections(:))
      scatteringAngles(:) = acos(projections)
      
      !
      ! Look up the phase function values in each of the scattering directions 
      !   from the tabulated phase functions
      !
      phaseFunctionIndex = phaseFuncI(xIndex, yIndex, zIndex, component)
      if(thisIntegrator%useHybridPhaseFunsForIntenCalcs .and. & 
         scatteringOrder <= thisIntegrator%numOrdersOrigPhaseFunIntenCalcs) then 
        phaseFunctionVals(:) =                                                                                          &
          lookUpPhaseFuncValsFromTable(tabulatedOrigPhaseFunctions(component)%values(:, phaseFunctionIndex), &
                                       scatteringAngles)
      else 
        !
        ! If useHybridPhaseFunsForIntenCalcs is false then the tabulated phase functions are the orginals
        ! 
        phaseFunctionVals(:) =                                                                                          &
          lookUpPhaseFuncValsFromTable(tabulatedPhaseFunctions(component)%values(:, phaseFunctionIndex), &
                                       scatteringAngles)
      end if 
      normalizedPhaseFunc(:) = phaseFunctionVals(:) / (4 * Pi * abs(thisIntegrator%intensityDirections(3, :)))
    end if 
    
    if(.not. thisIntegrator%useRussianRouletteForIntensity) then 
      !
      ! Find the integrated extinction from the current location to the boundary along each
      !   direction at which intensity is desired. 
      !
      do i = 1, numIntensityDirections
        xTemp = xPosI; yTemp = yPosI; zTemp = zPosI
        call accumulateExtinctionAlongPath(thisDomain, thisIntegrator%intensityDirections(:, i), &
                                           xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                           tausToBoundary(i))            
      end do
      !
      ! The contribution to intensity in each direction is the product of 
      !   photonWeight * phase function value * 1/(4 pi) * 1/abs(mu) * transmission
      !   The photon weight already contains the single scattering albedo for 
      !   this collision
      ! This contribution will be added to the intensity at the x-y location
      where(tausToBoundary(:) >= 0.) 
        contributions(:) = photonWeight * normalizedPhaseFunc * exp(-tausToBoundary(:)) 
      elsewhere
        ! The extinction subroutine signals errors by producing extinction values < 0 
        !   We could trap this but for the moment we'll just safeguard against disaster
        contributions(:) = 0. 
      end where 
    else
      ! 
      ! Russian roulette for intensity calculations (H. Iwabuchi, JAS 2006, Vol 63, pp 2324-2339) 
      !
      do i = 1, numIntensityDirections
        xTemp = xPosI; yTemp = yPosI; zTemp = zPosI
        tauFree = -log(max(tiny(tauFree), getRandomReal(randomNumbers)))
        !
        !   Small phase function contributions (Iwabuchi Eq 13)
        ! 
        if(Pi * normalizedPhaseFunc(i) <= thisIntegrator%zetaMin) then 
          call accumulateExtinctionAlongPath(thisDomain, thisIntegrator%intensityDirections(:, i), &
                                             xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                             tausToBoundary(i), tauFree)            
          !
          ! accumulateExtinctionAlongPath stops tracing once tausToBoundary(i) == tauFree
          !    so we use zF to see if the photon has escaped (i.e. if tau <= tauFree)
          !
          if(getRandomReal(randomNumbers) <= Pi * normalizedPhaseFunc(i)/thisIntegrator%zetaMin .and. &
             zIndexF(i) >= zIndexMax) then
            contributions(i) = photonWeight * thisIntegrator%zetaMin/Pi
          else
            contributions(i) = 0. 
          end if 
        else 
          !
          ! Use the full contribution if the optical depth is small and 
          !   play Russian roulette with contributions that undergo large extinction. 
          !   (Iwabuchi Eq 14). 
          !
          tauMax = -log(thisIntegrator%zetaMin/max(tiny(normalizedPhaseFunc), Pi * normalizedPhaseFunc(i)))
          call accumulateExtinctionAlongPath(thisDomain, thisIntegrator%intensityDirections(:, i), &
                                             xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                             tausToBoundary(i), tauMax)            
          if(zIndexF(i) >= zIndexMax .and. tausToBoundary(i) >= 0.) then 
            !
            !  This means tau <= tauMax
            ! 
            contributions(i) =  photonWeight * normalizedPhaseFunc(i) * exp(-tausToBoundary(i)) 
          else if (tausToBoundary(i) >= 0.) then 
            call accumulateExtinctionAlongPath(thisDomain, thisIntegrator%intensityDirections(:, i), &
                                               xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                               tausToBoundary(i), tauFree)            
            !
            ! accumulateExtinctionAlongPath stops tracing once tausToBoundary(i) == tauFree
            !    so we use zF to see if the photon has escaped (i.e. if tau <= tauMax + tauFree)
            !
            if(zIndexF(i) >= zIndexMax) then 
              contributions(i) = photonWeight * thisIntegrator%zetaMin/Pi
            else 
              contributions(i) = 0. 
            end if
          else
            !
            ! This branch means the ray tracing to tauMax failed
            !
            contributions(i) = 0.
          end if
        end if
      end do 
    end if 
    
    if(thisIntegrator%limitIntensityContributions) then 
      !
      ! Limit local estimate contribution to intensity; keep track of excess 
      !   so it can be redistributed after the fact
      !
      where(contributions(:) > thisIntegrator%maxIntensityContribution) 
        thisIntegrator%intensityExcess(:, component) = &
          thisIntegrator%intensityExcess(:, component) + & 
          contributions(:) - thisIntegrator%maxIntensityContribution
        contributions(:) = thisIntegrator%maxIntensityContribution
      end where
    end if
    DO i = 1, numComps
      call finalize_Matrix(tabulatedPhaseFunctions(i))
      call finalize_Matrix(tabulatedOrigPhaseFunctions(i)) 
    END DO
    deallocate(phaseFuncI, tabulatedPhaseFunctions,tabulatedOrigPhaseFunctions) 
  end subroutine computeIntensityContribution
  !------------------------------------------------------------------------------------------
  pure function lookUpPhaseFuncValsFromTable(tablulatedPhaseFunction, scatteringAngles) &
                    result(phaseFunctionVals)
    real, dimension(:),         intent(in ) :: tablulatedPhaseFunction, scatteringAngles
    real, dimension(size(scatteringAngles)) :: phaseFunctionVals
    !
    ! Find the value of phase function for the component doing the scattering at the 
    !   angle from the photon direction to each intensity direction.  
    ! Interpolate linearly in angle between the two closest tabulated values - 
    !   here we have passed in only the single relevant tabulated phase function
    !   (i.e. we have picked the correct element from the table corresponding to the 
    !    right component.) 
    ! The phase functions are tabulated equally in angle, with the first element corresponding 
    !   to scattering angle 0 and the last to scattering angle of pi. 
    !
    
    ! Local variables
    integer                                    :: nAngleSteps
    real                                       :: deltaTheta
    integer, dimension(size(scatteringAngles)) :: angleIndicies
    real,    dimension(size(scatteringAngles)) :: weights

    nAngleSteps = size(tablulatedPhaseFunction)
    deltaTheta  = Pi / (nAngleSteps - 1) 
    angleIndicies(:) = int(scatteringAngles(:) / deltaTheta) + 1

    !
    ! Since angleIndices could be >= nAngleSteps...
    !
    where(angleIndicies(:) < nAngleSteps)  
      ! Angle at corresponding index 
      weights(:) = 1. - (scatteringAngles(:) - (angleIndicies(:) - 1) * deltaTheta)/ deltaTheta
      phaseFunctionVals(:) =                                            &
              weights(:)  * tablulatedPhaseFunction(angleIndicies(:)) + &
        (1. - weights(:)) * tablulatedPhaseFunction(angleIndicies(:) + 1)
    elsewhere
      phaseFunctionVals(:) = tablulatedPhaseFunction(nAngleSteps) 
    end where
    
    
  end function lookUpPhaseFuncValsFromTable
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  pure function makeDirectionCosines(mu, phi)
    real,   intent(in) :: mu, phi
    real, dimension(3) :: makeDirectionCosines
    !
    ! Direction cosines S are defined so that
    !  S(1) = sin(theta) * cos(phi), projection in X direction 
    !  S(2) = sin(theta) * sin(phi), projection in Y direction
    !  S(3) = cos(theta),            projection in Z direction
    ! 
    ! Input is mu = cos(theta) and phi
    !
    
    real :: sinTheta, cosPhi, sinPhi
    
    sinTheta = sqrt(1. - mu**2)
    cosPhi   = cos(Phi) 
    sinPhi   = sin(Phi) ! sqrt(1 - cosPhi) is ambiguous at 90, 270 degrees. 
    makeDirectionCosines(:) = (/ sinTheta * cosPhi, sinTheta * sinPhi, mu /)
  end function makeDirectionCosines
  !------------------------------------------------------------------------------------------

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

  !------------------------------------------------------------------------------------------
  SUBROUTINE NEXT_DIRECT (randomNumbers, scatteringCosine, S)
    !
    !   Finds new set of direction cosines S from the original set and 
    !   the cosine of the scattering angle. 
    !   Algorithm from Marchuk et al (1980); implementation by Frank Evans
    !
    type(randomNumberSequence), intent(inout) :: randomNumbers
    real,                       intent(in   ) :: scatteringCosine
    real, dimension(:),         intent(inout) :: S
     
    REAL :: D, AX, AY, B
  
    D = 2.0
    DO WHILE (D .GT. 1.0)
      AX = 1.0 - 2.0*getRandomReal(randomNumbers)
      AY = 1.0 - 2.0*getRandomReal(randomNumbers)
      D = AX**2 + AY**2
    ENDDO
    B = SQRT((1.0 - scatteringCosine**2)/D)
    AX = AX*B
    AY = AY*B
    B = S(1)*AX-S(2)*AY
    D = scatteringCosine - B/(1.0+ABS(S(3)))
    
    S(1) = S(1)*D + AX
    S(2) = S(2)*D - AY
    S(3) = S(3)*scatteringCosine - sign(b, s(3) * b)
  END SUBROUTINE NEXT_DIRECT
  !------------------------------------------------------------------------------------------
!  subroutine getInfo_Integrator(thisIntegrator, ssa, cumExt)
!    type(Integrator), intent(in)                  :: thisIntegrator
!    real(8), dimension(:,:,:,:), intent(out)         :: ssa
!    real(8), dimension(:,:,:),   intent(out)           :: cumExt
!
!    ssa(:,:,:,:) = thisIntegrator%ssa(:,:,:,:)
!    cumExt(:,:,:) = thisIntegrator%totalExt(:,:,:)
!
!  end subroutine  getInfo_Integrator
!-------------------------------------------------------------------------------------------
end module monteCarloRadiativeTransfer
