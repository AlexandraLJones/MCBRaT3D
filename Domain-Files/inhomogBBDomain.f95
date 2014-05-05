! Broadband Domain creation code
! Written by Alexandra Jones Feb. 2014

program homogBBDomain
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  use netcdf  
  implicit none
  
  character(len=256) :: LGfile,  domfile, domfile_str, component_string, file_list
  integer ::  ncid, varid, status, nx, ny, nz, i, pEnt,maxz=20, nLG, nlambda, ilambda, isoIndex
  integer :: nLegendreCoefficients = 180
  integer :: nxc,nyc,nzc,xs,ys,zs, xe, ye, ze, i
  real(8)    ::  dx, dy, dz, Tb, Tt, concen, SSAs, SSAe, beta_eE, beta_eS, s_lambda, dlambda, albedo
  real     :: re, g
  real, allocatable	 ::  LG(:)
  real(8), allocatable    :: Ext(:,:,:), SSA(:,:,:),  Temp(:,:,:), lambdas(:)
  integer, allocatable	::	pInd(:,:,:), LGind(:)
  integer, dimension(3) ::  dims
  real, dimension(1000) :: thetas
  real, dimension(1000,1) :: values
  real, dimension(2) 			:: LegendreCoefs
  logical								::	HGphase=.false.
    
  type(ErrorMessage)              :: domstatus
  type(phaseFunction)		  :: PhaseFunc
  type(phaseFunctionTable)        :: table, OnePhaseFuncTable
  type(domain)                    :: cloud
  
  nLG = 12
  re = 10.0
  file_list = './namelist/inhomogDomfiles_list.in'

! read data from input
  
  read(*,'(1F)') albedo
  read(*,'(1F)') g
  read(*,'(1I)') nx
  read(*,'(1I)') ny
  read(*,'(1I)') nz
  read(*,'(1I)') nxc
  read(*,'(1I)') nyc
  read(*,'(1I)') nzc
  read(*,'(1I)') xs
  read(*,'(1I)') ys
  read(*,'(1I)') zs
  read(*,'(1I)') nlambda
  read(*,'(1F)') Tb
  read(*,'(1F)') Tt
  read(*,'(1F)') dx
  read(*,'(1F)') dy
  read(*,'(1F)') dz
  read(*,'(1F)') dlambda
  read(*,'(1F)') s_lambda
  read(*,'(1E12.0)') SSAs
  read(*,'(1E12.0)') SSAe
  read(*,'(1E12.0)') beta_eS
  read(*,'(1E12.0)') beta_eE
  read(*,'(1I)') isoIndex
  read(*,'(1A)') domfile_str
  read(*,'(1A)') component_string


  allocate(LG(1:nLG))
  allocate(LGind(1:nLG))
  LG(:)=(g**(/(i, i=1, nLG)/)) !! HG phase function
  LGind=isoIndex

!--Create 3D array of temps
   allocate(Temp(1:nx,1:ny,1:nz))
   DO i = 0, nz-1
     Temp(:,:,i+1)=(Tb *(dble(nz-1-i)/(nz-1))) + (Tt *(i/dble(nz-1)))
   END DO

   allocate(lambdas(1:nlambda))
   lambdas = s_lambda + (dlambda * (/(dble(i), i = 0, nlambda-1)/))


  !--create arrays of cloud optical properties
  allocate(Ext(1:nx,1:ny,1:nz))
  allocate(SSA(1:nx,1:ny,1:nz))
  allocate(pInd(1:nx,1:ny,1:nz))
  Ext=0.0_8
  SSA=0.0_8
  pInd = isoIndex
  xe=xs+nxc-1
  ye=ys+nyc-1
  ze=zs+nzc-1
  DO i=0, nxc-1
     Ext(xs+i,ys:ye,zs:ze) = (beta_eS *(dble(nxc-1-i)/(nxc-1))) + (beta_eE *(i/dble(nxc-1)))
     SSA(xs+i,ys:ye,zs:ze) = (SSAs *(dble(nxc-1-i)/(nxc-1))) + (SSAe *(i/dble(nxc-1)))
  END DO

  open(unit=22, file=TRIM(file_list), status='UNKNOWN')
  ! --Create Phase Function Table
  DO ilambda= 1, nlambda
    write(domfile, '(A,1ES11.5,A)') TRIM(domfile_str), lambdas(ilambda) , '.dom'
    write(22,"(A)") TRIM(domfile)
    cloud =                                                             &
    new_Domain(xPosition = dx * (/ 0., (real(i), i = 1, nx) /), &
               yPosition = dy * (/ 0., (real(i), i = 1, ny) /), &
               zPosition = dz * (/ 0., (real(i), i = 1, nz) /), &
               lambda = lambdas(ilambda), lambdaI = ilambda, nlambda=nlambda, &
               albedo = albedo, temps = Temp, status = domstatus)
   call printStatus(domstatus)

    PhaseFunc = new_PhaseFunction(LG, status=domstatus)
    call printStatus(domstatus)
 
    OnePhaseFuncTable = new_PhaseFunctionTable ((/PhaseFunc/), key=(/re/), status=domstatus)
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
  END DO
close(22)  
end program homogBBDomain

