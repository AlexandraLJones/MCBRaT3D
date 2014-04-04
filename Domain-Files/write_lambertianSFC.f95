
program write_lambertianSFC

use ErrorMessages
use userInterface
use surfaceProperties

implicit none

real(8) :: domainSize=10.0_8 
integer :: nx=20, ny=20, nParam=1
real(8) :: a1=0.0_8, a2 =0.4_8 
character*100 :: ncpath="lambertian20x20.sfc"
type(ErrorMessage) :: status
real(8) :: deltaX, deltaY
integer :: i, j
real(8), allocatable, dimension(:) :: xPosition
real(8), allocatable, dimension(:) :: yPosition
real(8), allocatable, dimension(:,:,:) :: albedo
type(surfaceDescription) :: surface


allocate(xPosition(nx+1))
allocate(yPosition(ny+1))
allocate(albedo(1,nx,ny)) 


deltaX=domainSize/dble(nx)
deltaY=domainSize/dble(ny)

xPosition = deltaX * (/0.0_8, (dble(i), i=1,nx) /)
yPosition = deltaY * (/0.0_8, (dble(i), i=1,ny) /)

forall (i=1:nx)
 forall (j=1:ny)
  albedo (1,i,j) = dble((i*j))/dble((nx*ny))
 end forall
end forall

surface = new_SurfaceDescription(albedo, xPosition, yPosition, status)
call printStatus(status)
call write_SurfaceDescription(surface, trim(ncpath),status)
call printStatus(status)
print*, "Successfully wrote surface file"


end program write_lambertianSFC
