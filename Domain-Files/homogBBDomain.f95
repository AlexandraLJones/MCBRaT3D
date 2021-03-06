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
  
  character(len=256) :: LGfile,  domfile, domfile_str, component_string, file_list, sourcefile_str
  integer ::  ncid, varid, status, nx, ny, nz, i, pEnt,maxz=20, nLG, nlambda, ilambda, isoIndex, source, ncFileId, DimId, ncVarId
  integer :: nLegendreCoefficients = 180
  real(8)    ::  dx, dy, dz, isoT, concen, isoSSA, isobeta_e, s_lambda, dlambda, albedo, Rad
  real     :: re, g
  real, allocatable	 ::  LG(:)
  real(8), allocatable    :: Ext(:,:,:), SSA(:,:,:),  lambdas(:), sourceFunc(:)!, Temp(:,:,:)
  integer, allocatable	::	pInd(:,:,:), LGind(:)
  integer, dimension(3) ::  dims
  integer, dimension(16) :: ncStatus
  real, dimension(1000) :: thetas
  real, dimension(1000,1) :: values
  real, dimension(2) 			:: LegendreCoefs
  logical								::	HGphase=.false.
    
  type(ErrorMessage)              :: domstatus, srcstatus
  type(phaseFunction)		  :: PhaseFunc
  type(phaseFunctionTable)        :: table, OnePhaseFuncTable
  type(commonDomain)              :: commonD
  type(domain)                    :: cloud

  re = 10.0
  file_list = './namelist/homogDomfiles_list.in'

! read data from input
  
  read(*,'(1I)') source
  read(*,'(1F)') Rad
  read(*,'(1A)') sourcefile_str
  read(*,'(1F)') albedo
  read(*,'(1F)') g
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
  read(*,'(1A)') domfile_str
  read(*,'(1A)') component_string


! read data from files
  nLG = 12
  allocate(LG(1:nLG))
  allocate(LGind(1:nLG))
  LG(:)=(g**(/(i, i=1, nLG)/)) !! HG phase function
  LGind=isoIndex

!--Create 3D array of temps
   allocate(commonD%temps(1:nx,1:ny,1:nz))
   commonD%temps(:,:,:)=isoT

   allocate(lambdas(1:nlambda))
   lambdas = s_lambda + (dlambda * (/(dble(i), i = 0, nlambda-1)/))


  !--create arrays of cloud optical properties
  allocate(Ext(1:nx,1:ny,1:nz))
  allocate(SSA(1:nx,1:ny,1:nz))
  allocate(pInd(1:nx,1:ny,1:nz))
  Ext=isobeta_e
  SSA=isoSSA
  pInd=1

  open(unit=22, file=TRIM(file_list), status='UNKNOWN')
  ! --Create Phase Function Table
  allocate(commonD%xPosition(1:nx), commonD%yPosition(1:ny), commonD%zPosition(1:nz))
  DO ilambda= 1, nlambda
    write(domfile, '(A,1ES11.5,A)') TRIM(domfile_str), lambdas(ilambda) , '.dom'
    write(22,"(A)") TRIM(domfile)
    
    commonD%xPosition=dx * (/ 0., (real(i), i = 1, nx) /)
    commonD%yPosition=dy * (/ 0., (real(i), i = 1, ny) /)
    commonD%zPosition=dz * (/ 0., (real(i), i = 1, nz) /)

    cloud = new_Domain(commonD, lambda = lambdas(ilambda), lambdaI = ilambda, nlambda=nlambda, &
                        albedo = albedo, status = domstatus)
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

  if(source .lt. 0)then  ! solar source; create solar source file
     allocate(sourceFunc(1:nLambda))
     sourceFunc = Rad
     ncStatus(:) = nf90_NoErr
     ncStatus( 1) = nf90_create(trim(sourcefile_str), nf90_Clobber, ncFileId)
     ncStatus( 2) = nf90_def_dim(ncFileId, "Lambdas", nLambda, DimId)
     ncStatus( 3) = nf90_def_var(ncFileId, "Lambdas", nf90_double, DimId, ncVarId)
     ncStatus( 4) = nf90_def_var(ncFileId, "SourceFunction", nf90_double, DimId, ncVarId)

     if(any(ncStatus(:) /= nf90_NoErr)) &
         call setStateToFailure(srcstatus, "homogBBDomain: error writing source function file")
     if(.not. stateIsFailure(srcstatus))then
         ncStatus( 1) = nf90_EndDef(ncFileId)
         ncStatus( 2) = nf90_inq_varid(ncFileId, "Lambdas", ncVarID)
         ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, lambdas)
         ncStatus( 4) = nf90_inq_varid(ncFileId, "SourceFunction", ncVarID)
         ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, sourceFunc)
         ncStatus( 6) = nf90_close(ncFileId)
     end if

  end if
end program homogBBDomain

