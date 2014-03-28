
module emissionAndBBWeights

  use OpticalProperties
  use ErrorMessages
  use RandomNumbers
  use NumericUtilities

  implicit none
  private 

  type weights
    private
    real(8)                                  :: spectrIntgrFlux = 0.0_8
    real(8)                                  :: spectrIntgrFractAtmsPower = 0.0_8
    real(8), dimension(:), pointer            :: atmsPowerCDF  => null()
    real(8), dimension(:), pointer            :: sfcPowerCDF  => null()

    real(8), dimension(:,:), pointer          :: levelWeights => null()
    real(8), dimension(:,:,:), pointer        :: colWeights => null()
    real(8), dimension(:,:,:,:), pointer      :: voxelWeights => null()
  end type weights


  
!  interface xxxxx

!  end interface xxxxx

  public :: weights

  public :: new_Weights, emission_weighting!, getInfo_weights, xxxxx

  contains

   function new_Weights(numX, numY, numZ, numLambda, status) result(theseWeights)
    integer, optional, intent(in)      :: numX, numY, numZ
    integer, intent(in)                :: numLambda
    type(ErrorMessage), intent(inout) :: status
    type(Weights)                      :: theseWeights

    allocate(theseWeights%atmsPowerCDF(1:numLambda), theseWeights%sfcPowerCDF(1:numLambda))
    theseWeights%atmsPowerCDF = 0.0_8
    theseWeights%sfcPowerCDF  = 0.0_8

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

   subroutine emission_weighting(thisDomain, theseWeights, sfcTemp, totalPhotons, atmsPhotons, voxel_weights, col_weights, level_weights, totalFlux, status)
!Computes Planck Radiance for each surface and atmosphere pixel to determine the wieghting for the distribution of photons.
!Written by ALexandra Jones, University of Illinois, Urbana-Champaign, Fall 2011
! Updated Fall 2012 to remove predetermination of number of photons emitted per column
! Updated Spring 2013 to work with Broadband code
     implicit none

     type(Domain), intent(in)                             :: thisDomain
     type(Weights), intent(inout)                         :: theseWeights
     real(8), intent(in)                                  :: sfcTemp
     integer,  intent(in)                                 :: totalPhotons
     integer,  intent(out)                             :: atmsPhotons
     real(8), allocatable, dimension(:,:,:,:), intent(out)          :: voxel_weights
     real(8), allocatable, dimension(:,:,:), intent(out)          :: col_weights
     real(8), allocatable, dimension(:,:), intent(out)          :: level_weights
     real(8),                                intent(out)  :: totalFlux
     type(ErrorMessage), intent(inout)                         :: status
     !Local variables
     integer                                           :: ix, iy, iz, ilambda, nx, ny, nz, nlambda, nComps !, last
     real(8)                                              ::  sfcPlanckRad, sfcPower,  atmsPower, totalPower, totalAbsCoef, b, lambda, lambda_u, albedo, emiss
     real(8)                                          :: previous, corr_contrib,corr,temp_sum, prev_exact
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


     call getInfo_Domain(thisDomain, numX=nx, numY=ny, numZ=nz, domainNumLambda=nlambda, &
                         numberOfComponents=nComps, status=status)
     allocate(voxel_weights(1:nx,1:ny,1:nz,1:nlambda), col_weights(1:ny,1:nz,1:nlambda), &
              level_weights(1:nz,1:nlambda))
     allocate(xPosition(1:nx+1), yPosition(1:ny+1), zPosition(1:nz+1), dz(1:nz),         &
              dy(1:ny), dx(1:nx), atmsTemp(1:nx,1:ny,1:nz), ssas(1:nx,1:ny,1:nz,1:nComps),&
              ext(1:nx,1:ny,1:nz,1:nComps), cumExt(1:nx,1:ny,1:nz))
     call getInfo_Domain(thisDomain, xPosition=xPosition, yPosition=yPosition,           &
                         temps=atmsTemp, zPosition=zPosition, status=status)
!PRINT *, h, c, k, lambda, Pi, a, b
!calculate arrays of depths from the position arrays in km
     dz(1:nz)=zPosition(2:nz+1)-zPosition(1:nz)
     dy(1:ny)=yPosition(2:ny+1)-yPosition(1:ny)
     dx(1:nx)=xPosition(2:nx+1)-xPosition(1:nx)

     DO ilambda = 1, nlambda
        call getInfo_Domain(thisDomain, lambda=lambda_u, lambdaIndex=ilambda, &
                            albedo=albedo, ssa=ssas, ext=ext, totalExt=cumExt, status=status)
        emiss = (1.0_8 - albedo)
        lambda=lambda_u/(10.0_8**6.0_8) ! convert lambda from micrometers to meters
        b=h*c/(k*lambda)        
!dz(1:nz)= 0.04                                 ! be sure to remove this line after debugging FOR DIAGNOSTIC PURPOSES ONLY!

!     last=nx*ny*nz ! the value of the index of the last element of the voxel_weights array

!first compute atms planck radiances then combine algorithms from mcarWld_fMC_srcDist and mcarWld_fMC_srcProf to determine the  weights of each voxel taking into consideration the ones that would be emitted from the surface instead.
     if (emiss .eq. 0.0_8 .or. sfcTemp .eq. 0.0_8)then
        theseWeights%sfcPowerCDF(ilambda)=0.0_8
     else
        sfcPlanckRad=(a/((lambda**5.0_8)*(exp(b/sfcTemp)-1.0_8)))/(10.0_8**6.0_8)
        theseWeights%sfcPowerCDF(ilambda)= Pi*emiss*sfcPlanckRad*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)     ! [W] factor of 1000^2 needed to convert area from km to m
     end if

     theseWeights%atmsPowerCDF(ilambda)=0.0_8
     previous=0.0_8
     corr_contrib=0.0_8
     temp_sum=0.0_8
     corr=0.0_8
     prev_exact=0.0_8

    if(COUNT(atmsTemp .le. 0.0_8) .eq. 0)then
     do iz = 1, nz
       do iy = 1, ny
         do ix = 1, nx
           atmsPlanckRad= (a/((lambda**5.0_8)*(exp(b/atmsTemp(ix,iy,iz))-1.0_8)))/(10.0_8**6.0_8) ! the 10^-6 factor converts it from Wsr^-1m^-3 to Wm^-2sr^-1micron^-1
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
          if (theseWeights%voxelWeights(nx,ny,nz,ilambda) .gt. 0.0_8) then
               theseWeights%atmsPowerCDF(ilambda) = theseWeights%voxelWeights(nx,ny,nz,ilambda)*(xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8)/dble(nx*ny)  ! [W] total power emitted by atmosphere. Factor of 1000^2 is to convert dx and dy from km to m
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
          end if

!PRINT *, 'level_weights= ', level_weights, 'col_weights= ', col_weights
     END DO

     voxel_weights = theseWeights%voxelWeights
     col_weights = theseWeights%colWeights
     level_weights = theseWeights%levelWeights

     atmsPower=SUM(theseWeights%atmsPowerCDF)
     sfcPower=SUM(theseWeights%sfcPowerCDF)
     totalPower=atmsPower + sfcPower
     if (totalPower .eq. 0.0_8)then
        CALL setStateToFailure(status, 'emission_weighting: Neither surface nor atmosphere will emitt photons since total power is 0. Not a valid solution')
     end if
     totalFlux=totalPower/((xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))*(1000.0_8**2.0_8))  ! We want the units to be [Wm^-2] but the x and y positions are in km
     theseWeights%spectrIntgrFlux = totalFlux     

!PRINT *, 'atmsPower= ',atmsPower, 'sfcPower= ', sfcPower, ' totalFlux=', totalFlux, ' totalArea=', (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1)), &
!         ' average column area=', (SUM(dx)/dble(nx))*(SUM(dy)/dble(ny)), (xPosition(nx+1)-xPosition(1))*(yPosition(ny+1)-yPosition(1))/dble(nx*ny), ' expected radiance=', atmsPlanckRad*(1.0_8-exp(-1.0_8*totalAbsCoef*(zPosition(nz+1)-zPosition(1))))
     if (atmsPower .eq. 0.0_8)then
         atmsPhotons=0
     else
        theseWeights%spectrIntgrFractAtmsPower = atmsPower / totalPower
!	PRINT *, theseWeights%spectrIntgrFractAtmsPower, atmsPower, totalPower
!	PRINT *, theseWeights%spectrIntgrFlux
!	PRINT *, theseWeights%atmsPowerCDF
!	PRINT *, theseWeights%sfcPowerCDF
        atmsPhotons=ceiling(totalPhotons * theseWeights%spectrIntgrFractAtmsPower)
     end if
   end subroutine emission_weighting

end module emissionAndBBWeights
