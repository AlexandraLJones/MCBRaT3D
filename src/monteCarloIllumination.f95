! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/monteCarloIllumination.f95 $
module monteCarloIllumination
  ! Provides an object representing a series of photons. The object specifies the 
  !   initial x, y position and direction. 
  ! In principle any arbitrary illumination condition can be specified. 
  ! Initial x, y positions are between 0 and 1, and must be scaled by the Monte Carlo
  !   integrator. 
  ! On input, azimuth is in degrees and zenith angle is specified as the cosine. 
  ! On output, azimuth is in radians (0, 2 pi) and solar mu is negative (down-going). 
  
  use ErrorMessages
  use RandomNumbers
  use emissionAndBBWeights
  use numericUtilities
 
  implicit none
  private

  !------------------------------------------------------------------------------------------
  ! Constants
  !------------------------------------------------------------------------------------------
  logical, parameter :: useFiniteSolarWidth = .false. 
  real,    parameter :: halfAngleDiamaterOfSun = 0.25 ! degrees
  
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type photonStream
    integer                     :: currentPhoton = 0
    real(8), dimension(:), pointer :: xPosition    => null()
    real(8), dimension(:), pointer :: yPosition    => null()
    real(8), dimension(:), pointer :: zPosition    => null()
    real, dimension(:), pointer :: solarMu      => null()
    real, dimension(:), pointer :: solarAzimuth => null()
  end type photonStream
  
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  interface new_PhotonStream
    module procedure newPhotonStream_Directional, newPhotonStream_RandomAzimuth, newPhotonStream_BBEmission, &
                     newPhotonStream_Flux, newPhotonStream_Spotlight, newPhotonStream_LWemission
  end interface new_PhotonStream
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: photonStream 
  public :: new_PhotonStream, finalize_PhotonStream, morePhotonsExist, getNextPhoton, emission_weighting
contains
  !------------------------------------------------------------------------------------------
  ! Code
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create streams of incoming photons
  !------------------------------------------------------------------------------------------
  function newPhotonStream_Directional(solarMu, solarAzimuth, &
                                       numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine and  
    !   azimuth. 
    real,                       intent(in   ) :: solarMu, solarAzimuth
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
        ! Local variables
    integer :: i
    
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(solarAzimuth < 0. .or. solarAzimuth > 360.) &
      call setStateToFailure(status, "setIllumination: solarAzimuth out of bounds")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
                        
      do i = 1, numberOfPhotons
        ! Random initial positions 
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      ! Specified inital directions
      photons%solarMu( :) = -abs(solarMu)
      photons%solarAzimuth(:) = solarAzimuth * acos(-1.) / 180. 
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
   end if   
  end function newPhotonStream_Directional
  ! ------------------------------------------------------
  function newPhotonStream_RandomAzimuth(solarMu, numberOfPhotons, randomNumbers, status) &
           result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine but
    !  random initial azimuth. 
    real,                       intent(in   ) :: solarMu
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
    ! Local variables
    integer :: i
    
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
      do i = 1, numberOfPhotons
        ! Random initial positions 
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial azimuth
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.) 
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      ! but specified inital mu
      photons%solarMu( :) = -abs(solarMu)
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
    end if   
 end function newPhotonStream_RandomAzimuth
 ! ------------------------------------------------------
 function newPhotonStream_Flux(numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with random initial azimuth and initial
    !  mus constructed so the solar flux on the horizontal is equally weighted
    !  in mu (daytime average is 1/2 solar constant; this is "global" average weighting)
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
    
    ! Local variables
    integer :: i
    
    ! Checks
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
     
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons),  photons%solarAzimuth(numberOfPhotons))
               
      do i = 1, numberOfPhotons
        ! Random initial positions
        photons%xPosition(   i) = getRandomReal(randomNumbers)
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial directions
        photons%solarMu(     i) = -sqrt(getRandomReal(randomNumbers)) 
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
      end do
      photons%zPosition(:) = 1. - spacing(1.) 
      
      photons%currentPhoton = 1
      call setStateToSuccess(status)  
    end if     
  end function newPhotonStream_Flux
  !------------------------------------------------------------------------------------------
  function newPhotonStream_Spotlight(solarMu, solarAzimuth, solarX, solarY, &
                                     numberOfPhotons, randomNumbers, status) result(photons)
    ! Create a set of incoming photons with specified initial zenith angle cosine and  
    !   azimuth. 
    real,                       intent(in   ) :: solarMu, solarAzimuth, solarX, solarY
    integer                                   :: numberOfPhotons
    type(randomNumberSequence), optional, &
                                intent(inout) :: randomNumbers
    type(ErrorMessage),         intent(inout) :: status
    type(photonStream)                        :: photons
        
    ! Checks: are input parameters specified correctly? 
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
    if(solarAzimuth < 0. .or. solarAzimuth > 360.) &
      call setStateToFailure(status, "setIllumination: solarAzimuth out of bounds")
    if(abs(solarMu) > 1. .or. abs(solarMu) <= tiny(solarMu)) &
      call setStateToFailure(status, "setIllumination: solarMu out of bounds")
    if(abs(solarX) > 1. .or. abs(solarX) <= 0. .or. &
       abs(solarY) > 1. .or. abs(solarY) <= 0. )    &
      call setStateToFailure(status, "setIllumination: x and y positions must be between 0 and 1")
    
    if(.not. stateIsFailure(status)) then
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons), photons%solarAzimuth(numberOfPhotons))
                        
      ! Specified inital directions and position
      photons%solarMu( :) = -abs(solarMu)
      photons%solarAzimuth(:) = solarAzimuth * acos(-1.) / 180. 
      photons%xPosition(:) = solarX
      photons%yPosition(:) = solarY
      photons%zPosition(:) = 1. - spacing(1.) 
      photons%currentPhoton = 1
      
      call setStateToSuccess(status)
   end if   
  end function newPhotonStream_Spotlight

!-------------------------------------------------------------------------------------------
  function newPhotonStream_LWemission(numberOfPhotons, atms_photons, voxel_weights, col_weights, level_weights, nx, ny, nz, randomNumbers, status, option1) result(photons)
    ! Create a set of emitted photons with random initial azimuth, random 
    ! mus, and random x,y,z location within the domain. This is the LW source from the atmosphere
    ! and surface. The x,y,z locations are weighted based on the power emitted 
    ! from each  voxel, pixel, in both the atmosphere and surface.
    ! Written by Alexandra Jones at the University of Illinois, Urbana-Champaign. Fall 2011
    ! Updated Fall 2012 to remove predetermination of number of photons emitted per column
    implicit none
    
    integer, intent(in)                             :: numberOfPhotons, nx, ny, nz
    type(randomNumberSequence), intent(inout)       :: randomNumbers
    type(ErrorMessage),         intent(inout)       :: status
    type(photonStream)                              :: photons
    integer,                       intent(in)       :: atms_photons
    real(8), dimension(nx,ny,nz), intent(in)         :: voxel_weights
    real(8), dimension(ny,nz), intent(in)         :: col_weights
    real(8), dimension(nz), intent(in)         :: level_weights
    integer, dimension(nx,ny,nz), optional, intent(out)   :: option1

    ! Local variables
    integer :: i, numberOfAtmsPhotons, startPhoton,  ii, ij, ik
    real    :: RN!, test
!    real, dimension(0:nx*ny*nz)                      :: temp1, temp2
    
    if(present(option1))option1=0

    ! Checks
    if(numberOfPhotons <= 0) &
      call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")
       
    if(.not. stateIsFailure(status)) then
!PRINT *, "number of photons > 0"
      allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons),  photons%solarAzimuth(numberOfPhotons))

              
     ! divide the photons into sfc vs. atms sources 
      numberOfAtmsPhotons=MIN(numberOfPhotons,atms_photons)
      if ( numberOfAtmsPhotons .gt. 0)then
!PRINT *, "numberOfAtmsPhotons .gt. 0"
          startPhoton = 1
     
          do i = startPhoton,  numberOfAtmsPhotons ! Loop over photons from atms source 
!           do while (i .le. numberOfPhotoins) 
            ! Random initial positions
            RN = getRandomReal(randomNumbers)
!if(i .eq. 1)PRINT *, "i=", i, " RN=", RNi          
            DO ik=1,nz
              if(RN .le. level_weights(ik))then
                DO ij=1,ny
                  if(RN .le. col_weights(ij,ik))then
                    DO ii=1,nx
                      if(RN .le. voxel_weights(ii,ij,ik))then
                        exit
                      endif     
                    ENDDO
                    exit
                  endif
                ENDDO
                exit
              endif
            ENDDO

!PRINT *, RN, ii, ij, ik
!write(16,"(I5 ,2X, E30.20)") i, voxel_weights(ii,ij,ik)-RN

 !if(present(option1))then
 !   option1(ii,ij,ik)=option1(ii,ij,ik)+1
 !end if 

            photons%zPosition(i) = ((ik-1)*(1.0_8)/nz) + dble(getRandomReal(randomNumbers)/nz) ! The first term represents the fractional position of the layer bottom in the column, such that ik=1 corresponds to a position of 0. The second term respresents the position within the layer.
            if(ik .eq. 1 .and. photons%zPosition(i) .eq. 0.0_8) photons%zPosition(i)=0.0_8+spacing(1.0_8)
            if(ik .eq. nz .and. photons%zPosition(i) .gt. 1.0_8-2.0_8*spacing(1.0_8)) photons%zPosition(i)=photons%zPosition(i) - (2.0_8*spacing(1.0_8))
            photons%xPosition(i) = ((ii -1)*1.0_8/nx) + dble(getRandomReal(randomNumbers)*(1.0/nx)) 
            photons%yPosition(i) = ((ij -1)*1.0_8/ny) + dble(getRandomReal(randomNumbers)*(1.0/ny)) 
!if(i .eq. 1) PRINT *, 'ind= ', ind, 'i= ', ik, 'j= ', ij, 'k= ', ik, 'xPos= ', photons%xPosition(i)
            ! Random initial directions
!            test=getRandomReal(randomNumbers)
!PRINT *, 'i=', i, ' atms_rand=', test 
!            photons%solarMu(i) = 1-(2.*test)
          DO
            photons%solarMu(i) = 1-(2.*getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative transfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu no longer has to be negative. The name "solarMu" should really just be "sourceMu" 
            if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the layer it's initialized in
          END DO
            photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
            if (numberOfPhotons < numberOfAtmsPhotons)exit
!           end do
          end do  !loop over i

     end if

     if (numberOfPhotons > numberOfAtmsPhotons)then
      do i = numberOfAtmsPhotons+1, numberOfPhotons ! Loop over sfc source
        ! Random initial positions
        photons%xPosition(   i) = getRandomReal(randomNumbers) ! Assuming the surface temperature and emissivities are the same everywhere the x,y values don't need to be weighted
        photons%yPosition(   i) = getRandomReal(randomNumbers)
        ! Random initial directions
!        test = getRandomReal(randomNumbers)
!PRINT *, 'i=', i, ' sfc_rand=', test
!         photons%solarMu(     i) = sqrt(test)
       DO
        photons%solarMu(     i) = sqrt(getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative trasnfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu has to be positive for a sfc source. The name "solarMu" should really just be "sourceMu"
        if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the lowest layer. this is especially bad when there is no atmospheric extinction because then the photon will never leave the domain or hit the surface or become extinct. Thus we enter an infinite loop
       END DO
        photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
      end do  !loop over i
        photons%zPosition(numberOfAtmsPhotons+1:numberOfPhotons) = 0.0_8 ! must be from the sfc. Make sure this value makes sense for how the position is interpretted
     end if
!PRINT *,  ' photons%solarMu=', photons%solarMu(1)
      photons%currentPhoton = 1
      call setStateToSuccess(status) 
    end if    
  end function newPhotonStream_LWemission
!-------------------------------------------------------------------------------------------
  function newPhotonStream_BBEmission(theseWeights, iLambda, numberOfPhotons, randomNumbers, status, option1) result(photons)
	implicit none

	type(Weights), intent(in)                               :: theseWeights
	integer, intent(in)                                     :: numberOfPhotons, iLambda
	type(randomNumberSequence), intent(inout)               :: randomNumbers
	type(ErrorMessage), intent(inout)                       :: status
	integer, allocatable, dimension(:,:,:), optional, intent(out)        :: option1
	type(photonStream)                                      :: photons

	real(8)                                                 :: fracAtmsPower
	real(8), allocatable, dimension(:)                      :: levelWeights
	real(8), allocatable, dimension(:,:)                    :: colWeights
	real(8), allocatable, dimension(:,:,:)                  :: voxelWeights
	real                                                    :: RN
	integer                                                 :: numX, numY, numZ, i, ii, ij, ik, sfcTally, atmsTally 

    sfcTally = 0; atmsTally = 0
    ! Checks
        if(numberOfPhotons <= 0) &
          call setStateToFailure(status, "setIllumination: must ask for non-negative number of photons.")

        if(.not. stateIsFailure(status)) then
	  call getInfo_Weights(theseWeights=theseWeights, iLambda=iLambda, numX=numX, numY=numY, numZ=numZ, status=status)
	  if(present(option1))then
            allocate(option1(1:numX,1:numY,1:numZ))
            option1=0
          end if

          allocate(photons%xPosition(numberOfPhotons),   photons%yPosition(numberOfPhotons), &
               photons%zPosition(numberOfPhotons),                                       &
               photons%solarMu(numberOfPhotons),  photons%solarAzimuth(numberOfPhotons))

	! get relevant theseWeights values    !!!!!!!BUT i REALLY ONLY HAVE TO DO THIS FOR THE ATMOSPHERE SO SHOULD IT GO IN A DIFFERENT SUBROUTINE, MAYBE IN THE EMISSIONbbwEIGHTS MODULE?????
	  allocate(levelWeights(1:numX), colWeights(1:numX,1:numY), voxelWeights(1:numX,1:numY,1:numZ))
	  call getInfo_weights(theseWeights=theseWeights, iLambda=iLambda, fracAtmsPower=fracAtmsPower, levelWeights=levelWeights, colWeights=colWeights, voxelWeights=voxelWeights, status=status)
	  DO i = 1, numberOfPhotons
	  ! determine if it's emitting from the surface or atmosphere
	    RN = getRandomReal(randomNumbers)
	    if (RN .gt. fracAtmsPower) then ! must be from surface
		photons%xPosition(   i) = getRandomReal(randomNumbers) ! Assuming the surface temperature and emissivities are the same everywhere the x,y values don't need to be weighted
        	photons%yPosition(   i) = getRandomReal(randomNumbers)
		DO
        	  photons%solarMu(     i) = sqrt(getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative trasnfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu has to be positive for a sfc source. The name "solarMu" should really just be "sourceMu"
        	  if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the lowest layer. this is especially bad when there is no atmospheric extinction because then the photon will never leave the domain or hit the surface or become extinct. Thus we enter an infinite loop
       		END DO
		sfcTally = sfcTally + 1
		photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
		photons%zPosition(i) = 0.0_8
	    else  ! must be from atmos
		RN = getRandomReal(randomNumbers)
		ik = findCDFIndex(RN, levelWeights)
		ij = findCDFIndex(RN, colWeights(:,ik))
		ii = findCDFIndex(RN, voxelWeights(:,ij,ik))

		photons%zPosition(i) = ((ik-1)*(1.0_8)/numZ) + dble(getRandomReal(randomNumbers)/numZ) ! The first term represents the fractional position of the layer bottom in the column, such that ik=1 corresponds to a position of 0. The second term respresents the position within the layer.
            	if(ik .eq. 1 .and. photons%zPosition(i) .eq. 0.0_8) photons%zPosition(i)=0.0_8+spacing(1.0_8)
            	if(ik .eq. numZ .and. photons%zPosition(i) .gt. 1.0_8-2.0_8*spacing(1.0_8)) photons%zPosition(i)=photons%zPosition(i) - (2.0_8*spacing(1.0_8))
            	photons%xPosition(i) = ((ii -1)*1.0_8/numX) + dble(getRandomReal(randomNumbers)*(1.0/numX)) 
            	photons%yPosition(i) = ((ij -1)*1.0_8/numY) + dble(getRandomReal(randomNumbers)*(1.0/numY))
		DO
            	  photons%solarMu(i) = 1-(2.*getRandomReal(randomNumbers))    ! This formula is from section 10.5 of '3D radiative transfer in cloudy atmospheres'. These angles will stay the same...truly random numbers, since LW emission is isotropic. But the Mu no longer has to be negative. The name "solarMu" should really just be "sourceMu" 
            	  if(abs(photons%solarMu(i)) > 2 * tiny (photons%solarMu(i)))exit  ! This ensures that there is some vertical component to the photon trajectory, so the photon doesn't get permanently stuck in the layer it's initialized in
          	END DO
		photons%solarAzimuth(i) = getRandomReal(randomNumbers) * 2. * acos(-1.)
		atmsTally = atmsTally + 1
		if(present(option1))then
	   	  option1(ii,ij,ik)=option1(ii,ij,ik)+1
		end if
	    end if
	  END DO
PRINT *, atmsTally, sfcTally
	  photons%currentPhoton = 1
	  call setStateToSuccess(status)
 	end if
  end function newPhotonStream_BBEmission
!------------------------------------------------------------------------------------------
!  ACTUALLY I DON'T THINK I NEED A NEW FUNCTION FOR BB SOLAR. THE MONOCHROMATIC ONE SHOULD WORK FINE. I WONDER IF THAT MEANS i COULD REPLACE THE MONOCHROMATIC EMISSION ONE WITH ONE THAT WORKS FOR BOTH BB AND MONO....
!  function newPhotonStream_BBSolar(theseWeights, solarMu, solarAzimuth, numberOfPhotons, randomNumbers, status) result(photons)
!    implicit none

!    type(Weights), intent(in)             :: theseWeights
!    real, intent(in)                      :: solarMu, solarAzimuth
!    integer, intent(in)                   :: numberOfPhotons
!    type(randomNumberSequence), intent(inout) :: randomNumbers
!    type(ErrorMessage), intent(inout)     :: status



!  end function newPhotonStream_BBSolar 
  !------------------------------------------------------------------------------------------
  ! Are there more photons? Get the next photon in the sequence
  !------------------------------------------------------------------------------------------
  function morePhotonsExist(photons)
    type(photonStream), intent(inout) :: photons
    logical                           :: morePhotonsExist
    
    morePhotonsExist = photons%currentPhoton > 0 .and. &
                       photons%currentPhoton <= size(photons%xPosition)
  end function morePhotonsExist
  !------------------------------------------------------------------------------------------
  subroutine getNextPhoton(photons, xPosition, yPosition, zPosition, solarMu, solarAzimuth, status)
    type(photonStream), intent(inout) :: photons
    real(8),                 intent(  out) :: xPosition, yPosition, zPosition
    real,                 intent(  out) ::  solarMu, solarAzimuth
    type(ErrorMessage),   intent(inout) :: status
    
    ! Checks 
    ! Are there more photons?
    if(photons%currentPhoton < 1) &
      call setStateToFailure(status, "getNextPhoton: photons have not been initialized.")  
    if(.not. stateIsFailure(status)) then
      if(photons%currentPhoton > size(photons%xPosition)) &
        call setStateToFailure(status, "getNextPhoton: Ran out of photons")
    end if
      
    if(.not. stateIsFailure(status)) then
      xPosition    = photons%xPosition(photons%currentPhoton) 
      yPosition    = photons%yPosition(photons%currentPhoton)
      zPosition    = photons%zPosition(photons%currentPhoton)
      solarMu      = photons%solarMu(photons%currentPhoton)
!PRINT *, 'solarMu=', solarMu
      solarAzimuth = photons%solarAzimuth(photons%currentPhoton)
      photons%currentPhoton = photons%currentPhoton + 1
    end if
     
  end subroutine getNextPhoton
  !------------------------------------------------------------------------------------------
  ! Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_PhotonStream(photons)
    type(photonStream), intent(inout) :: photons
    ! Free memory and nullify pointers. This leaves the variable in 
    !   a pristine state
    if(associated(photons%xPosition))    deallocate(photons%xPosition)
    if(associated(photons%yPosition))    deallocate(photons%yPosition)
    if(associated(photons%zPosition))    deallocate(photons%zPosition)
    if(associated(photons%solarMu))      deallocate(photons%solarMu)
    if(associated(photons%solarAzimuth)) deallocate(photons%solarAzimuth)
    
    photons%currentPhoton = 0
  end subroutine finalize_PhotonStream

end module monteCarloIllumination
