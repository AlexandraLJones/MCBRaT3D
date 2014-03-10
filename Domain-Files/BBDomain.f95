! Broadband Domain creation code
! Written by Alexandra Jones Feb. 2014

program BBDomain
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  use netcdf  
  implicit none
  
  character(len=256) :: LGfile,  domfile, component_string
  integer ::  ncid, varid, status, nx, ny, nz, i, pEnt,maxz=20, nLG, nlambda, ilambda
  integer :: nLegendreCoefficients = 180
  real(8)    :: re, dx, dy, dz, isoT, concen, isoSSA, beta_e, s_lambda, dlambda
  real, allocatable	 ::  LG(:)
  real(8), allocatable    :: Ext(:,:,:,:), SSA(:,:,:,:),  Temp(:,:,:), lambdas(:)
  integer, allocatable	::	pInd(:,:,:,:), LGind(:)
  integer, dimension(3) ::  dims
  real, dimension(1000) :: thetas
  real, dimension(1000,1) :: values
  real, dimension(2) 			:: LegendreCoefs
  logical								::	HGphase=.false.
    
  type(ErrorMessage)              :: domstatus
  type(phaseFunction), allocatable :: PhaseFunc(:,:)
  type(phaseFunctionTable)        :: table, OnePhaseFuncTable
  type(domain)                    :: cloud

  re = 10.0

! read data from input
  read(*,'(1A)') LGfile
  read(*,'(1I)') nLG
  read(*,'(1I)') nx
  read(*,'(1I)') ny
  read(*,'(1I)') nz
  read(*,'(1I)') nlambda
  read(*,'(1F)') isoT
  read(*,'(1F)') dx
  read(*,'(1F)') dy
  read(*,'(1F)') dz
  read(*,'(1F)') dlambda
  read(*,'(1F)') s_lambda
  read(*,'(1E12.0)') isoSSA
  read(*,'(1E12.0)') isobeta_e
  read(*,'(1I)') isoIndex
  read(*,'(1A)') domfile
  read(*,'(1A)') component_string


! read data from files
  allocate(LG(1:nLG-1))
  allocate(LGind(1:nLG-1))
  open(unit=2,file=LGfile,status='OLD')
    read (2, '(1I12, 1G24.19)')
    DO i=1, nLG-1
       read(2,'(1I12, 1G24.19)') LGind(i), LG(i)
!       PRINT *, LGind(i), LG(i)
    END DO        
  close(2)

  !--Create 3D array of temps
   allocate(Temp(1:nx,1:ny,1:nz))
   Temp(:,:,:)=isoT

   allocate(lambdas(1:nlambda))
   lambdas = s_lambda + (dlambda * (/(dble(i), i = 0, nlambda-1)/))

   allocate(PhaseFunc(size(re),1:nlambda))

  !--Create Domain
  cloud =                                                             &
    new_Domain(xPosition = dx * (/ 0., (real(i), i = 1, nx) /), &
               yPosition = dy * (/ 0., (real(i), i = 1, ny) /), &
               zPosition = dz * (/ 0., (real(i), i = 1, nz) /), &
               lambdas = lambdas, &
               temps = Temp, status = domstatus)
   call printStatus(domstatus)

  !--create arrays of cloud optical properties
  allocate(Ext(1:nx,1:ny,1:nz,1:nlambda))
  allocate(SSA(1:nx,1:ny,1:nz,1:nlambda))
  allocate(pInd(1:nx,1:ny,1:nz,1:nlambda))
  Ext=isobeta_e
  SSA=isoSSA
  pInd=1

  ! --Create Phase Function Table
  DO ilambda= 1, nlambda
  PhaseFunc(ilambda) = new_PhaseFunction(LG, status=domstatus)
  call printStatus(domstatus)
  END DO
  OnePhaseFuncTable = new_PhaseFunctionTable (PhaseFunc, re_key=(/re/), lambda_key=lambdas, status=domstatus)
  call printStatus(domstatus)

  !--Add the cloud
	print*,'Adding cloud to domain.'
  call addOpticalComponent(cloud, component_string,  &
                           Ext, SSA, pInd, OnePhaseFuncTable, status = domstatus)
  call printStatus(domstatus)
  
    
  !--Write it to the file

  call write_Domain(cloud, trim(domfile), status = domstatus)
  print *, "Writing domain: ",trim(domfile)
  call printStatus(domstatus)
  if(stateIsSuccess(domstatus)) call setStateToCompleteSuccess(domstatus)
  
  
end program BBDomain

