! SitCom to Domain Generator
! Daeven Jackson
! Sept. 7, 2012
! Modified by Alexandra Jones Nov. 2013

program cloudsdale
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  use netcdf  
  implicit none
  
  character(len=256) :: LGfile,  domfile, component_string
  integer ::  ncid, varid, status, nx, ny, nz, i, pEnt,maxz=20, nxc,nyc,nzc,xs,ys,zs,xe,ye,ze, nLG
  integer :: nLegendreCoefficients = 180
  real(8)    :: dx, dy, dz, Tb, Tt, concen, isoSSA, beta_e, beta_e86
  real 	 ::  tauScale, re, g, tau86 ,sigma_e, sigma_s, tau
  real, allocatable	 ::  Ql(:,:,:), p(:), T(:), LWC(:,:,:), LG(:)
  real(8), allocatable    :: Ext(:,:,:), SSA(:,:,:), atmSSA(:), RaylExt(:), Temp(:,:,:)
  real	::RaylWavelen,raylcoef
  integer, allocatable	::	pInd(:,:,:), atmpInd(:), LGind(:)
  integer, dimension(3) ::  dims
  real, dimension(1000) :: thetas
  real, dimension(1000,1) :: values
  real, dimension(2) 			:: LegendreCoefs
  logical								::	HGphase=.false.
    
  type(ErrorMessage)              :: domstatus
  type(phaseFunction)             :: phase, PhaseFunc
  type(phaseFunctionTable)        :: table, OnePhaseFuncTable
  type(domain)                    :: cloud

! read data from input
  read(*,'(1A)') LGfile
  read(*,'(1I)') nLG
  read(*,'(1I)') nx
  read(*,'(1I)') ny
  read(*,'(1I)') nz
  read(*,'(1I)') nxc
  read(*,'(1I)') nyc
  read(*,'(1I)') nzc
  read(*,'(1I)') xs
  read(*,'(1I)') ys
  read(*,'(1I)') zs
  read(*,'(1F)') Tb
  read(*,'(1F)') Tt
  read(*,'(1F)') tau86
  read(*,'(1E12.0)') beta_e86
  read(*,'(1F)') dx
  read(*,'(1F)') dy
  read(*,'(1F)') dz
  read(*,'(1E12.0)') isoSSA
  read(*,'(1E12.0)') beta_e
  read(*,'(1E12.0)') re
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
   if(Tb .eq. Tt)then
        Temp(:,:,:)=Tb
   else
        Temp(:,:,1)=Tb
        Temp(:,:,nz)=Tt
        do i=2,nz-1
          Temp(:,:,i)=Tb + ((i-1) * ((Tt-Tb)/(nz-1))) 
        end do
   end if

  
  !--Create Domain
  cloud =                                                             &
    new_Domain(xPosition = dx * (/ 0., (real(i), i = 1, nx) /), &
               yPosition = dy * (/ 0., (real(i), i = 1, ny) /), &
               zPosition = dz * (/ 0., (real(i), i = 1, nz) /), &
               temps = Temp, status = domstatus)
   call printStatus(domstatus)

  !--create arrays of cloud optical properties
  allocate(Ext(1:nx,1:ny,1:nz))
  allocate(SSA(1:nx,1:ny,1:nz))
  allocate(pInd(1:nx,1:ny,1:nz))
  Ext=0.0_8
  SSA=0.0_8
  xe=xs+nxc-1
  ye=ys+nyc-1
  ze=zs+nzc-1
  if(beta_e86 .gt. 0.0_8)then
    tau=tau86*beta_e/beta_e86
  else
    tau = 0.0_8
  end if
 PRINT *, tau86, beta_e86, beta_e, tau
  Ext(xs:xe,ys:ye,zs:ze)=tau/(nzc*dz)
 PRINT *, Ext(xs,ys,zs), tau86/(nzc*dz), nzc, dz
  SSA(xs:xe,ys:ye,zs:ze)=isoSSA
  pInd=1

  ! --Create Phase Function Table
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
  
  
end program cloudsdale

