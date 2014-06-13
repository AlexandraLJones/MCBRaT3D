
module emissionAndBBWeights

  use OpticalProperties
  use ErrorMessages
  use RandomNumbers
  use NumericUtilities
  use netcdf 

  implicit none
  private 

  type weights
    private
    real(8)                                  :: spectrIntgrFlux = 0.0_8
    real(8), dimension(:), pointer            :: fracAtmsPower  => null()
    real(8), dimension(:), pointer            :: totalPowerCDF  => null()

    real(8), dimension(:,:), pointer          :: levelWeights => null()
    real(8), dimension(:,:,:), pointer        :: colWeights => null()
    real(8), dimension(:,:,:,:), pointer      :: voxelWeights => null()
  end type weights


  
!  interface xxxxx

!  end interface xxxxx

  public :: weights

  public :: new_Weights, finalize_Weights, emission_weighting, getFrequencyDistr, getInfo_weights, read_SolarSource, solar_Weighting !, xxxxx

  contains

   function new_Weights(numX, numY, numZ, numLambda, status) result(theseWeights)
    integer, optional, intent(in)      :: numX, numY, numZ
    integer, intent(in)                :: numLambda
    type(ErrorMessage), intent(inout) :: status
    type(Weights)                      :: theseWeights

    allocate(theseWeights%fracAtmsPower(1:numLambda), theseWeights%totalPowerCDF(1:numLambda))
    theseWeights%fracAtmsPower = 0.0_8
    theseWeights%totalPowerCDF  = 0.0_8

    if(present(numX) .or. present(numY) .or. present(numZ))then
      if(present(numX) .and. present(numY) .and. present(numZ))then
	allocate(theseWeights%voxelWeights(1:numX, 1:numY, 1:numZ, 1:numLambda))
	theseWeights%voxelWeights(:,:,:,:) = 0.0_8
	theseWeights%colWeights => theseWeights%voxelWeights(numX, :, :, :)
	theseWeights%levelWeights => theseWeights%colWeights(numY, :, :)
      else
	call setStateToFailure(status, "new_Weights: must supply all physical dimensions for emission weighting arrays")
      end if
    end if
   end function new_Weights

   subroutine finalize_Weights(theseWeights)
     type(Weights), intent(inout)  :: theseWeights

     deallocate(theseWeights%fracAtmsPower, theseWeights%totalPowerCDF)
     if (associated(theseWeights%voxelWeights))then
        NULLIFY(theseWeights%levelWeights)
        NULLIFY(theseWeights%colWeights)
        NULLIFY(theseWeights%voxelWeights)
     end if
   end subroutine finalize_Weights

   subroutine getInfo_Weights(theseWeights, iLambda, numX, numY, numZ, fracAtmsPower, levelWeights, colWeights, voxelWeights, status)
     implicit none

     type(Weights), intent(in)                           :: theseWeights
     integer, optional, intent(in)                       :: iLambda
     integer, optional, intent(out)                      :: numX, numY, numZ
     real(8), optional, intent(out)                      :: fracAtmsPower
     real(8), optional, dimension(:),intent(out)         :: levelWeights
     real(8), optional, dimension(:,:),intent(out)       :: colWeights
     real(8), optional, dimension(:,:,:),intent(out)     :: voxelWeights
     type(ErrorMessage), intent(inout)                   :: status

     integer                                             :: j, k

     if(present(numX)) numX = size(theseWeights%voxelWeights,1)
     if(present(numY)) numY = size(theseWeights%voxelWeights,2)
     if(present(numZ)) numZ = size(theseWeights%voxelWeights,3)

     if(present(fracAtmsPower))then
	if(present(iLambda))then
	  fracAtmsPower = theseWeights%fracAtmsPower(iLambda)
	else
	  call setStateToFailure(status, "getInfo_Weights: need to specify iLambda for fracAtmsPower.")
	end if
     end if

     if (present(levelWeights) .and. present(iLambda))then
!PRINT *, size(levelWeights), size(theseWeights%levelWeights,1)
	if(size(levelWeights) .ne. size(theseWeights%voxelWeights,3)) then
	  call setStateToFailure(status, "getInfo_Weights: array for levelWeights is wrong dimensions.")
	else
	  levelWeights = theseWeights%levelWeights(:,ilambda)
	end if
     else
	call setStateToFailure(status, "getInfo_Weights: must supply iLambda for levelWeights.")
     end if
     if (present(colWeights) .and. present(iLambda))then
!PRINT *, size(colWeights,1), size(colWeights,2), size(theseWeights%colWeights,1), size(theseWeights%colWeights,2)
        if(size(colWeights,1) .ne. size(theseWeights%voxelWeights,2) .or. size(colWeights,2) .ne. size(theseWeights%voxelWeights,3)) then
          call setStateToFailure(status, "getInfo_Weights: array for colWeights is wrong dimensions.")
        else
          colWeights = theseWeights%colWeights(:,:,iLambda)
        end if
     else
	call setStateToFailure(status, "getInfo_Weights: must supply iLambda for colWeights.")
     end if
     if (present(voxelWeights) .and. present(iLambda))then
!PRINT *, size(voxelWeights,1), size(voxelWeights,2), size(voxelWeights,3), size(theseWeights%voxelWeights,1), size(theseWeights%voxelWeights,2), size(theseWeights%voxelWeights,3)
        if(size(voxelWeights,1) .ne. size(theseWeights%voxelWeights,1) .or. size(voxelWeights,2) .ne. size(theseWeights%voxelWeights,2) .or. size(voxelWeights,3) .ne. size(theseWeights%voxelWeights,3)) then
          call setStateToFailure(status, "getInfo_Weights: array for voxelWeights is wrong dimensions.")
        else
          voxelWeights = theseWeights%voxelWeights(:,:,:,iLambda)
if (ilambda .eq. 10)then
    write(32,"(36F12.8)") levelWeights(:)
    DO k = 1, size(theseWeights%voxelWeights,3)
      write(33,"(100F12.8)") colWeights(:,k)
      DO j = 1, size(theseWeights%voxelWeights,2)
        write(31,"(100F12.8)") voxelWeights(:,j,k)
      end do
    end do
    close(31)
    close(32)
    close(33)
end if
        end if
     else
	call setStateToFailure(status, "getInfo_Weights: must supply iLambda for voxelWeights.")
     end if

   end subroutine getInfo_Weights
!---------------------------------------------------------------------------------------------------
   subroutine solar_Weighting(theseWeights, nLambda, radianceFunction, lambdas, solarMu, totalFlux, status)
     type(Weights), intent(inout)                         :: theseWeights
     integer, intent(in)                                  :: nLambda
     real(8), dimension(1:nLambda), intent(in)            :: lambdas, radianceFunction
     real, intent(in)                                     :: solarMu
     real(8), intent(out)                                 :: totalFlux
     type(ErrorMessage), intent(inout)                    :: status

     real(8), parameter                                   :: Des = 1.496E8 ! [km] mean earth-sun distance
     real(8), parameter                                   :: Rs = 6.96E5   ! [km] mean radius of sun's photosphere
     real(8), parameter                                   :: Pi=4*DATAN(1.0_8)
     real(8), parameter                                   :: solarSolidAngle=2.0_8*Pi*(1.0_8 - DCOS(DASIN(Rs/(Rs+Des))))
     integer                                              :: i
     real(8)                                              :: dLambda

     dLambda = lambdas(2)-lambdas(1)
     theseWeights%totalPowerCDF(1) = dLambda*solarSolidAngle*solarMu*radianceFunction(1)

     DO i=2, nLambda
        if (i .gt. 1 .and. i .lt. nlambda) then
          dLambda = (lambdas(i+1)-lambdas(i-1))/2.0_8 ! half points between ilambda and the adjacent values
        elseif (i .eq. nlambda) then
          dLambda = lambdas(i)-lambdas(i-1)
        else ! should never end up here
          PRINT *, 'solar_weighting: ended up in error condition of the loop. i = ', i
        end if 
        theseWeights%totalPowerCDF(i) = theseWeights%totalPowerCDF(i-1)+ (dLambda*solarSolidAngle*solarMu*radianceFunction(i))
!PRINT *, 'cumFlux= ', theseWeights%totalPowerCDF(i), 'dLambda= ',dLambda, 'solid angle= ', solarSolidAngle, 'mu= ', solarMu, 'radiance= ', radianceFunction(i)
     END DO
     theseWeights%spectrIntgrFlux = theseWeights%totalPowerCDF(nLambda)
     totalFlux = theseWeights%spectrIntgrFlux
     theseWeights%totalPowerCDF = theseWeights%totalPowerCDF/theseWeights%totalPowerCDF(nLambda)
PRINT *, "solar_weighting: totalFlux and CDF ", totalFlux, theseWeights%totalPowerCDF
   end subroutine solar_Weighting
!-------------------------------------------------------------------------------------------
   subroutine emission_weighting(theseDomains, nLambda, theseWeights, sfcTemp, totalPhotons, atmsPhotons, totalFlux, status)
!Computes Planck Radiance for each surface and atmosphere pixel to determine the wieghting for the distribution of photons.
!Written by ALexandra Jones, University of Illinois, Urbana-Champaign, Fall 2011
! Updated Fall 2012 to remove predetermination of number of photons emitted per column
! Updated Spring 2013 to work with Broadband code
     implicit none

     type(Domain), dimension(nLambda),intent(in)        :: theseDomains
     type(Weights), intent(inout)                         :: theseWeights
     real(8), intent(in)                                  :: sfcTemp
     integer,  intent(in)                                 :: nLambda, totalPhotons
     integer,  intent(out)                             :: atmsPhotons
!     real(8), allocatable, dimension(:,:,:,:), intent(out)          :: voxel_weights
!     real(8), allocatable, dimension(:,:,:), intent(out)          :: col_weights
!     real(8), allocatable, dimension(:,:), intent(out)          :: level_weights
     real(8),                                intent(out)  :: totalFlux
     type(ErrorMessage), intent(inout)                         :: status
     !Local variables
     integer                                           :: ix, iy, iz, ilambda, nx, ny, nz, nComps !, last
     real(8)                                              ::  sfcPlanckRad, sfcPower,  atmsPower, totalPower, totalAbsCoef, b, lambda, albedo, emiss, dlambda
     real(8), dimension(1:3)                          :: lambda_u
     real(8)                                          :: previous, corr_contrib,corr,temp_sum, prev_exact, tempPower
     real(8), allocatable, dimension(:)                             :: dx, dy, dz
     real(8), allocatable, dimension(:)                             :: xPosition, yPosition, zPosition
     real(8), allocatable, dimension(:,:,:)                         :: cumExt, atmsTemp
     real(8), allocatable, dimension(:,:,:,:)                       :: ssas, ext
     real(8)                                               :: atmsPlanckRad
     
!     real, dimension(1:nx, 1:ny)                       :: atmsColumn_power

     real(8), parameter                                   :: h=6.62606957e-34 !planck's constant [Js]
     real(8), parameter                                   :: c=2.99792458e+8 !speed of light [ms^-1]
     real(8), parameter                                   :: k=1.3806488e-23 !boltzman constant [J/K molecule]
     real(8), parameter                                   :: a=2.0_8*h*c**2.0_8
     real(8), parameter                                   :: Pi=4*DATAN(1.0_8)


     call getInfo_Domain(theseDomains(1), numX=nx, numY=ny, numZ=nz, namelistNumLambda=nlambda, &
                         numberOfComponents=nComps, status=status)
!     allocate(voxel_weights(1:nx,1:ny,1:nz,1:nlambda), col_weights(1:ny,1:nz,1:nlambda), &
!              level_weights(1:nz,1:nlambda))
     allocate(xPosition(1:nx+1), yPosition(1:ny+1), zPosition(1:nz+1), dz(1:nz),         &
              dy(1:ny), dx(1:nx), atmsTemp(1:nx,1:ny,1:nz), ssas(1:nx,1:ny,1:nz,1:nComps),&
              ext(1:nx,1:ny,1:nz,1:nComps), cumExt(1:nx,1:ny,1:nz))
     call getInfo_Domain(theseDomains(nlambda), xPosition=xPosition, yPosition=yPosition,           &
                         temps=atmsTemp, zPosition=zPosition, status=status)
!PRINT *, 'emissionWeighting: min temp: ', MINVAL(atmsTemp)
!PRINT *, h, c, k, lambda, Pi, a, b
!calculate arrays of depths from the position arrays in km
     dz(1:nz)=zPosition(2:nz+1)-zPosition(1:nz)
     dy(1:ny)=yPosition(2:ny+1)-yPosition(1:ny)
     dx(1:nx)=xPosition(2:nx+1)-xPosition(1:nx)
     

     tempPower = 0.0_8
     lambda_u = 0.0_8
     
      !!$OMP DO ORDERED PRIVATE(ilambda, lambda, dlambda, albedo, ssas, ext, cumExt, &
      !!$OMP& emiss, b, sfcPower, sfcPlanckRad, atmsPower, previous, corr_contrib, &
      !!$OMP& temp_sum, corr, prev_exact, atmsPlanckRad, totalAbsCoef, ix, iy, iz) &
      !!$OMP& FIRSTPRIVATE(tempPower, lambda_u, status) LASTPRIVATE(status)
   !!THIS IS A COMMENT!!$OMP& SHARED(theseWeights, dz, dy, dx, nlambda, theseDomains, h, c, k, sfcTemp, a, Pi, xPosition, yPosition, atmsTemp, nx, ny, nz) 
     DO ilambda = 1, nlambda
         !!$OMP ORDERED
	if (ilambda .gt. 1 .and. ilambda .lt. nlambda) then ! this if block has to be executed in order
	  call getInfo_Domain(theseDomains(ilambda+1), lambda=lambda_u(3), status=status)
	  dlambda = (lambda_u(3)-lambda_u(1))/2.0_8 ! half points between ilambda and the adjacent values
	  lambda = lambda_u(2)
	  lambda_u(1) = lambda_u(2)
          lambda_u(2) = lambda_u(3)
	elseif (ilambda .eq. 1) then
	  call getInfo_Domain(theseDomains(ilambda), lambda=lambda_u(1), status=status)
	  call getInfo_Domain(theseDomains(ilambda+1), lambda=lambda_u(2), status=status)
	  dlambda = lambda_u(2)-lambda_u(1)
	  lambda = lambda_u(1)
	elseif (ilambda .eq. nlambda) then
!	  PRINT *, ilambda
	  dlambda = lambda_u(2)-lambda_u(1)
	  lambda = lambda_u(2)
	else ! should never end up here
	  PRINT *, 'emission_weighting: ended up in error condition of the loop. ilambda = ', ilambda
        end if
         !!$OMP END ORDERED

!PRINT *, lambda, dlambda

        call getInfo_Domain(theseDomains(ilambda), albedo=albedo, ssa=ssas, ext=ext, totalExt=cumExt, status=status)
          
        emiss = (1.0_8 - albedo)
        lambda=lambda/(10.0_8**6.0_8) ! convert lambda from micrometers to meters
	dlambda = dlambda/(10.0_8**6.0_8)
        b=h*c/(k*lambda)        
!dz(1:nz)= 0.04                                 ! be sure to remove this line after debugging FOR DIAGNOSTIC PURPOSES ONLY!

!     last=nx*ny*nz ! the value of the index of the last element of the voxel_weights array

!first compute atms planck radiances then combine algorithms from mcarWld_fMC_srcDist and mcarWld_fMC_srcProf to determine the  weights of each voxel taking into consideration the ones that would be emitted from the surface instead.
     if (emiss .eq. 0.0_8 .or. sfcTemp .eq. 0.0_8)then
        sfcPower=0.0_8
     else
        sfcPlanckRad=dlambda*(a/((lambda**5.0_8)*(exp(b/sfcTemp)-1.0_8))) ! for broadband quantity we have to factor in discretization width, dlambda
        sfcPower = Pi*emiss*sfcPlanckRad*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)     ! [W] factor of 1000^2 needed to convert area from km to m
     end if
!PRINT *, 'emission_weighting: ilambda= ', ilambda, ' sfcPower= ', sfcPower
     atmsPower = 0.0_8
     previous=0.0_8
     corr_contrib=0.0_8
     temp_sum=0.0_8
     corr=0.0_8
     prev_exact=0.0_8

!PRINT *, 'emission_weighting: how many atmsTemp < 0?', COUNT(atmsTemp .le. 0.0_8)

    if(COUNT(atmsTemp .le. 0.0_8) .eq. 0)then
     do iz = 1, nz
       do iy = 1, ny
         do ix = 1, nx
           atmsPlanckRad= dlambda*(a/((lambda**5.0_8)*(exp(b/atmsTemp(ix,iy,iz))-1.0_8))) ! for broadband quantity we have to factor in discretization width, dlambda
           totalAbsCoef=cumExt(ix,iy,iz)-sum(ssas(ix,iy,iz,:) * ext(ix,iy,iz,:))
!PRINT *, cumExt(ix,iy,iz),ssas(ix,iy,iz,:), ext(ix,iy,iz,:), sum(ssas(ix,iy,iz,:) * ext(ix,iy,iz,:)), totalAbsCoef
           corr_contrib = (4.0_8*Pi* atmsPlanckRad * totalAbsCoef*dz(iz))-corr     ! [Wm^-2]
           temp_sum = previous + corr_contrib
           corr = (temp_sum - previous)-corr_contrib
           previous = temp_sum
           theseWeights%voxelWeights(ix,iy,iz,ilambda) = previous
           prev_exact=prev_exact + dble(1.0_8/(nx*ny*nz))
!           write(11, "(6E30.20)") atmsTemp(ix,iy,iz), atmsPlanckRad, totalAbsCoef, 4.0*Pi* atmsPlanckRad * totalAbsCoef*dz(iz), dz(iz), voxel_weights(ix,iy,iz)
!            write(11, "(9E30.20)") atmsTemp(ix,iy,iz), atmsPlanckRad, totalAbsCoef, 4.0_8*Pi* atmsPlanckRad * totalAbsCoef*dz(iz), dz(iz), voxel_weights(ix,iy,iz), dble( ((iz-1)*nx*ny)+((iy-1)*nx)+ix  )/dble(nx*ny*nz), prev_exact,corr
         end do ! i loop
!         col_weights(iy,iz)= previous !!!!!!!!!!!!!!!!I THINK I DON'T NEED THIS LINE ANYMORE BECAUSE OF POINTING THESEWEIGHTS%COLWEIGHTS AT THE PROPER INDICIES OF VOXELWEIGHTS
!          write(10, "(3I5, A, E30.20, A, E30.20)" ) ix, iy, iz, 'voxel_weights= ', voxel_weights(ix-1,iy,iz), 'col_weights= ', col_weights(iy,iz)
       end do   ! j loop
!       level_weights(iz)= previous  !!!!!!!!!!!!!!!!!!!I THINK I DON'T NEED THIS LINE ANYMORE BECAUSE OF POINTING THESEWEIGHTS%LEVELWEIGHTS AT THE PROPER INDICIES OF COLWEIGHTS
!       write(10, "(3I5, A, E30.20, A, E30.20, A, E30.20)" ) ix, iy, iz, 'voxel_weights= ', voxel_weights(ix-1,iy-1,iz), 'col_weights= ', col_weights(iy-1,iz), 'level_weights= ', level_weights(iz)
     end do     ! k loop
    end if
!PRINT *, 'emission_weighting: temp, cumExt, SUM(ssas*ext), totalAbsoCoeff ', atmsTemp(nx,ny,nz), cumExt(nx,ny,nz), sum(ssas(nx,ny,nz,:) * ext(nx,ny,nz,:)), cumExt(nx,ny,nz)-sum(ssas(nx,ny,nz,:) * ext(nx,ny,nz,:))
!PRINT *, 'emission_weighting: final Radiance, SSA, voxelWeights ',atmsPlanckRad, totalAbsCoef, theseWeights%voxelWeights(nx,ny,nz,ilambda)
          if (theseWeights%voxelWeights(nx,ny,nz,ilambda) .gt. 0.0_8) then
               atmsPower = theseWeights%voxelWeights(nx,ny,nz,ilambda)*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)/dble(nx*ny)  ! [W] total power emitted by atmosphere. Factor of 1000^2 is to convert dx and dy from km to m
               theseWeights%voxelWeights(:,:,:,ilambda)=theseWeights%voxelWeights(:,:,:,ilambda)/theseWeights%voxelWeights(nx,ny,nz,ilambda)     ! normalized
!               do iz = 1, nz
!                  do iy = 1, ny
!                     PRINT *,  theseWeights%voxelWeights(:,iy,iz,ilambda)
!		     PRINT *,  theseWeights%colWeights(iy,iz,ilambda)
!                  end do
!		  PRINT *, theseWeights%levelWeights(iz,ilambda)
!               end do

!               col_weights(:,:)=col_weights(:,:)/col_weights(ny,nz)  I THINK I DON'T NEED THESE 
!               level_weights(:)=level_weights(:)/level_weights(nz)   LINES B/C OF POINTERS

               theseWeights%voxelWeights(nx,ny,nz,ilambda)=1.0_8     ! need this to be 1 for algorithm used to select emitting voxel
!               col_weights(ny,nz)=1.0_8                              I THINK I DON'T NEED THESE
!               level_weights(nz)=1.0_8                               LINES B/C OF POINTERS
		theseWeights%fracAtmsPower(ilambda) = atmsPower/(atmsPower + sfcPower)
          end if
!PRINT *, 'emission weighting: ilambda= ', ilambda,  ' atmosPower= ', atmsPower, ' fractAtmosPower= ', theseWeights%fracAtmsPower(ilambda)
	  theseWeights%totalPowerCDF(ilambda) = tempPower + atmsPower + sfcPower
	  tempPower = theseWeights%totalPowerCDF(ilambda)
!PRINT *, 'level_weights= ', level_weights, 'col_weights= ', col_weights
     END DO
     deallocate( zPosition, dz, dy, dx, atmsTemp, ssas, ext, cumExt)    

     if (theseWeights%totalPowerCDF(nlambda) .eq. 0.0_8)then
        CALL setStateToFailure(status, 'emission_weighting: Neither surface nor atmosphere will emitt photons since total power is 0. Not a valid solution')
     else
       totalFlux=theseWeights%totalPowerCDF(nlambda)/((xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8))  ! We want the units to be [Wm^-2] but the x and y positions are in km
       theseWeights%spectrIntgrFlux = totalFlux     

!PRINT *, 'atmsPower= ',atmsPower, 'sfcPower= ', sfcPower, ' totalFlux=', totalFlux, ' totalArea=', (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1)), &
!         ' average column area=', (SUM(dx)/dble(nx))*(SUM(dy)/dble(ny)), (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))/dble(nx*ny), ' expected radiance=', atmsPlanckRad*(1.0_8-exp(-1.0_8*totalAbsCoef*(zPosition(nz+1)-zPosition(1))))

       atmsPhotons=ceiling(SUM(totalPhotons * theseWeights%fracAtmsPower(:)))
       theseWeights%totalPowerCDF = theseWeights%totalPowerCDF/theseWeights%totalPowerCDF(nlambda)
!PRINT *, "emission_weighting: fraction of atmos power ", theseWeights%fracAtmsPower
!PRINT *, "emission_weighting: total power CDF ", theseWeights%totalPowerCDF
     end if
     deallocate(xPosition, yPosition)
   end subroutine emission_weighting

   subroutine getFrequencyDistr(theseWeights, totalPhotons, randomNumbers, distribution)
     implicit none

     type(Weights), intent(in)                     :: theseWeights
     integer, intent(in)                           :: totalPhotons
     type(randomNumberSequence), intent(inout)     :: randomNumbers
     integer, dimension(:), intent(out)            :: distribution

     integer                                       :: numLambda, i, n
     real                                          :: RN


     numLambda = size(theseWeights%totalPowerCDF)
!     allocate(distribution(1:numLambda))
     distribution = 0.0_8

     DO n = 1, totalPhotons
	RN = getRandomReal(randomNumbers)
	i = findCDFIndex(RN, theseWeights%totalPowerCDF)
	distribution(i) = distribution(i)+1
     END DO
   end subroutine getFrequencyDistr
!---------------------------------------------------------------------------------
   subroutine read_SolarSource(fileName, nLambda, sourceFunc, lambdas, status)
    character(len = *), intent(in   ) :: fileName
    integer, intent(in)               :: nLambda
    real(8), dimension(1:nLambda), intent(out)   :: sourceFunc, lambdas    
    type(ErrorMessage), intent(inout) :: status

    integer                           :: dims, ncFileID, ncDimID, ncVarID
    integer, dimension(16)            :: ncStatus
    
    ncStatus(:) = nf90_NoErr
    if(nf90_open(trim(fileName), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      call setStateToFailure(status, "read_SolarSource: Can't open file " // trim(fileName))
    end if

    if(.not. stateIsFailure(status)) then
      ncStatus( 1) = nf90_inq_dimid(ncFileId, "Lambdas", ncDimId)
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = dims) 
      if(dims .ne. nLambda) &
        call setStateToFailure(status, "read_SolarSource: " // trim(fileName) // &
                               " dimension of solar source function does not match numLambdas from namelist.")
      ncStatus( 3) = nf90_inq_varid(ncFileId, "Lambdas", ncVarId)
      ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, lambdas)
      ncStatus( 5) = nf90_inq_varid(ncFileId, "SourceFunction", ncVarId)
      ncStatus( 6) = nf90_get_var(ncFileId, ncVarId, sourceFunc)
PRINT *, 'Radmax=', MAXVAL(sourceFunc), 'Radmin=', MINVAL(sourceFunc)
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "read_SolarSource: " // trim(fileName) // &
                               " doesn't look a solar source function file.")
    end if

! IT WOULD BE NICE TO BE ABLE TO CHECK THE LAMBDA VALUES HERE AGAINST THE DOMAIN FILE VALUES, BUT THAT MEANS IT WOULD HAVE TO BE CALLED LATER AFTER THE DOMAINS ARE ALL READ IN

    if(.not. stateIsFailure(status)) call setStateToSuccess(status) 
      

   end subroutine read_SolarSource

end module emissionAndBBWeights
