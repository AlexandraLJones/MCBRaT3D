! Broadband Domain creation code
! Written by Alexandra Jones Feb. 2014

program inhomogBBDomain
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  use netcdf  
  implicit none
  
  character(len=256) :: LGfile,  domfile, domfile_str, component_string, file_list, sourcefile_str, int_str, SRFfile_str
  integer ::  ncid, varid, status, nx, ny, nz, i, pEnt,maxz=20, nLG, nlambda, ilambda, isoIndex, source, shapeCode, SRFcode
  integer :: nLegendreCoefficients = 180
  integer :: nxc,nyc,nzc,xs,ys,zs, xe, ye, ze, i, ncFileId, DimId, ncVarId, nLines, center
  real(8)    ::  dx, dy, dz, Tb, Tt, concen, SSAs, SSAe, beta_eE, beta_eS 
  real(8) :: s_lambda, dlambda, albedo, Rad1, Rad2, W0, alpha, delta, SRho, y, x, Srange, v0
  real(8), parameter  :: Pi=4*DATAN(1.0_8)
  real     :: re, g
  real, allocatable	 ::  LG(:)
  real(8), allocatable    :: Ext(:,:,:), SSA(:,:,:), lambdas(:), sourceFunc(:), SRF(:) !, Temp(:,:,:)
  integer, allocatable	::	pInd(:,:,:), LGind(:)
  integer, dimension(3) ::  dims
  integer, dimension(16) :: ncStatus
  real, dimension(1000) :: thetas
  real, dimension(1000,1) :: values
  real, dimension(2) 			:: LegendreCoefs
  logical								::	HGphase=.false.
    
  type(commonDomain)              :: commonD
  type(ErrorMessage)              :: domstatus, srcstatus
  type(phaseFunction)		  :: PhaseFunc
  type(phaseFunctionTable)        :: table, OnePhaseFuncTable
  type(domain)                    :: cloud
  
  nLG = 12
  re = 10.0
  file_list = './namelist/inhomogDomfiles_list.in'

! read data from input

  read(*,'(1I)') source  
  read(*,'(1F)') alpha
  read(*,'(1F)') W0
  read(*,'(1I)') nLines
  read(*,'(1F)') SRho
  read(*,'(1F)') Rad1
  read(*,'(1F)') Rad2
  read(*,'(1I)') shapeCode
  read(*,'(1A)') sourcefile_str
  read(*,'(1I)') SRFcode
  read(*,'(1A)') SRFfile_str
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

PRINT *, SRFfile_str

  allocate(LG(1:nLG))
!PRINT *, 'allocated LG'
  allocate(LGind(1:nLG))
!PRINT *, 'allocated LGind'
  LG(:)=(g**(/(i, i=1, nLG)/)) !! HG phase function
  LGind=isoIndex

!--Create 3D array of temps
   allocate(commonD%temps(1:nx,1:ny,1:nz))
!PRINT *, 'allocated commonD%temps'
   if(nz .eq. 1)then
     commonD%temps(:,:,:) = Tb
PRINT *, 'nz=1', commonD%temps(1,1,i+1)
   else
     DO i = 0, nz-1
        commonD%temps(:,:,i+1)=(Tb *(dble(nz-1-i)/(nz-1))) + (Tt *(i/dble(nz-1)))
PRINT *, 'z=', i+1, commonD%temps(1,1,i+1)
     END DO
   end if

   allocate(lambdas(1:nlambda))
!PRINT *, 'allocated lambdas'
   lambdas = s_lambda + (dlambda * (/(dble(i), i = 0, nlambda-1)/))
   Srange = SUM((lambdas(3:nlambda)-lambdas(1:nlambda-2))/2.0_8) + (lambdas(2)-lambdas(1)) + (lambdas(nlambda)-lambdas(nlambda-1))  ! this includes the buffers at each end because that is how the I3RC handles it. Essentially eash lambda is considered to be the center of the spectral bin
   v0 = ((3.0_8*lambdas(1))-lambdas(2)+Srange)/2.0_8  ! center of spectral interval

  !--create arrays of cloud optical properties
  allocate(Ext(1:nx,1:ny,1:nz))
!PRINT *, 'allocated Ext'
  allocate(SSA(1:nx,1:ny,1:nz))
!PRINT *, 'allocated SSA'
  allocate(pInd(1:nx,1:ny,1:nz))
!PRINT *, 'allocated pInd'
  Ext=0.0_8
  SSA=0.0_8
  pInd = isoIndex
  xe=xs+nxc-1
  ye=ys+nyc-1
  ze=zs+nzc-1

  if (nxc .eq. 1)then
    Ext(xs,ys,zs) = beta_eS
    SSA(xs,ys,zs) = SSAs
  else
    DO i=0, nxc-1
       Ext(xs+i,ys:ye,zs:ze) = (beta_eS *(dble(nxc-1-i)/(nxc-1))) + (beta_eE *(i/dble(nxc-1)))
       SSA(xs+i,ys:ye,zs:ze) = (SSAs *(dble(nxc-1-i)/(nxc-1))) + (SSAe *(i/dble(nxc-1)))
    END DO
  end if

  open(unit=22, file=TRIM(file_list), status='UNKNOWN')
  ! --Create Phase Function Table
  allocate(commonD%xPosition(1:nx), commonD%yPosition(1:ny), commonD%zPosition(1:nz))
!PRINT *, 'allocated commonD% xPos, yPos, zPos'
  DO ilambda= 1, nlambda
!    write(domfile, '(A,1ES11.5,A)') TRIM(domfile_str), lambdas(ilambda) , '.dom'
    write(int_str, '(I0.3)') ilambda
    write(domfile, '(3A)') TRIM(domfile_str), TRIM(int_str) , '.dom'
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
    SELECT CASE (shapeCode)
      CASE (1,2) ! homogeneous spectral properties or linear inc/dec solar spectra with homogeneous extinction 
	! nothing needs to be done

      CASE (3) ! isolated square line shape
	SSA = 0.0_8
	Ext = 0.0_8
!	W0 = dlambda * nlambda /10.0_8  ! Play with this denominator until you reach a desired width
	PRINT *, 'W0=', W0, 'spectral range=', Srange 
!	PRINT *, 'central lambda = ', ((lambdas(nlambda)+lambdas(1))/2)
!	PRINT *, 'min, lambda, max', (((lambdas(nlambda)+lambdas(1))/2)-(W0/2)), lambdas(ilambda), (((lambdas(nlambda)+lambdas(1))/2)+(W0/2))
	if (lambdas(ilambda) .gt. (v0-(W0/2)) .and. lambdas(ilambda) .lt. (v0+(W0/2))) Ext = SRho / W0 

      CASE (4) ! isolated Lorentz line shape
	SSA = 0.0_8
	Ext = SRho*alpha/(Pi*((lambdas(ilambda)-v0)**2 + alpha**2))
	PRINT *, 'Ext=', Ext(1,1,1), 'spectral range=', Srange 

      CASE (5) ! Elsasser band model (regualry spaced, equal strength lorentz lines)
	SSA = 0.0_8
	delta = Srange/nLines
	y = alpha/delta
	x = lambdas(ilambda)/delta
	Ext = SRho*DSINH(2.0_8*pi*y)/(delta*(DCOSH(2*Pi*y)-DCOS(2*Pi*x)))
        PRINT *, 'y=', y, 'delta=', delta, 'spectral range=', Srange

      CASE (6,7) ! absorption increases linearly with wavelength
	SSA = 0.0_8
	Ext = beta_eS * lambdas(ilambda)
	PRINT *, 'Ext=', Ext(1,1,1), 'spectral range=', Srange

      CASE DEFAULT
          PRINT *, 'no matching case for shapeCode = ', shapeCode
    END SELECT
	PRINT*,'Adding cloud to domain.'
    call addOpticalComponent(cloud, component_string,  &
                           Ext, SSA, pInd, OnePhaseFuncTable, status = domstatus)
    call printStatus(domstatus)
  
    
  !--Write it to the file

    call write_Domain(cloud, trim(domfile), status = domstatus)
    PRINT *, "Writing domain: ",trim(domfile)
    call printStatus(domstatus)
    if(stateIsSuccess(domstatus)) call setStateToCompleteSuccess(domstatus)
  END DO
  close(22)  

  if(source .lt. 0)then  ! solar source; create solar source file
     allocate(sourceFunc(1:nLambda))
!PRINT *, 'allocated sourceFunc'
     SELECT CASE (shapeCode)
	CASE (1,3,4,5,6) !white Spectra 
	  sourceFunc = Rad1
	CASE (2,7) ! linearly increasing or decreasing spectral properties
	  DO i = 0, nLambda-1
     	     sourceFunc(i+1)=(Rad1 *(dble(nLambda-1-i)/(nLambda-1))) + (Rad2 *(i/dble(nLambda-1)))
   	  END DO  
!	CASE (3) ! isolated square line; uniform source function, SSA=0, special beta_e
	  
!	CASE (4) ! isolated lorentz line; uniform source function, SSA=0, special beta_e

!	CASE (5) ! Elsasser band model (regularly spaced, equal strength lines); uniform source function, SSA=0, special beta_e
	  
	CASE DEFAULT
	  PRINT *, 'no matching case for shapeCode = ', shapeCode
     END SELECT

     ncStatus(:) = nf90_NoErr
     ncStatus( 1) = nf90_create(trim(sourcefile_str), nf90_Clobber, ncFileId)
     ncStatus( 2) = nf90_def_dim(ncFileId, "Lambdas", nLambda, DimId)
     ncStatus( 3) = nf90_def_var(ncFileId, "Lambdas", nf90_double, DimId, ncVarId)
     ncStatus( 4) = nf90_def_var(ncFileId, "SourceFunction", nf90_double, DimId, ncVarId)

     if(any(ncStatus(:) /= nf90_NoErr)) &
         call setStateToFailure(srcstatus, "inhomogBBDomain: error writing source function file")
     if(.not. stateIsFailure(srcstatus))then
         ncStatus( 1) = nf90_EndDef(ncFileId)
         ncStatus( 2) = nf90_inq_varid(ncFileId, "Lambdas", ncVarID)
         ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, lambdas)
         ncStatus( 4) = nf90_inq_varid(ncFileId, "SourceFunction", ncVarID)
         ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, sourceFunc)
         ncStatus( 6) = nf90_close(ncFileId)
     end if
  end if

  if (SRFcode .gt. 0)then
PRINT *, "writing to file: ", SRFfile_str
     allocate(SRF(1:nLambda))
     SRF=0.0_8
     SELECT CASE (SRFcode)
       CASE(1) ! homogeneous
	SRF = 1.0_8
       CASE(2) ! rectangular
	DO ilambda =1,nLambda
	   if(DBLE(ilambda) .gt. nLambda/4.0_8 .and. DBLE(ilambda) .le. 3*nLambda/4.0_8)SRF(ilambda)=1.0_8
	END DO
       CASE(3) ! increasing
	DO ilambda =1,nLambda
	   SRF(ilambda) = DBLE(iLambda)/DBLE(nLambda)
	END DO
       CASE(4) ! decreasing
	DO ilambda =1,nLambda
	   SRF(ilambda) = DBLE(nLambda-ilambda+1)/DBLE(nLambda)
	END DO
       CASE(5) ! triangle
	center = nLambda/2
	DO ilambda =1,nLambda
	   if(ilambda .lt. center)then
	     SRF(ilambda) = DBLE(ilambda)/DBLE(center)
	   else if(ilambda .gt. center)then
	     SRF(ilambda) = DBLE(nlambda-ilambda +1)/DBLE(center)
	   else ! ilambda = center
	     SRF(ilambda) = 1.0_8
	   end if
        END DO
       CASE DEFAULT
	PRINT *, 'no matching case for SRFcode=', SRFcode
     END SELECT
     ncStatus(:) = nf90_NoErr
     ncStatus( 1) = nf90_create(trim(SRFfile_str), nf90_Clobber, ncFileId)
     ncStatus( 2) = nf90_def_dim(ncFileId, "Lambdas", nLambda, DimId)
     ncStatus( 3) = nf90_def_var(ncFileId, "Lambdas", nf90_double, DimId, ncVarId)
     ncStatus( 4) = nf90_def_var(ncFileId, "SRF", nf90_double, DimId, ncVarId)

     if(any(ncStatus(:) /= nf90_NoErr)) &
         call setStateToFailure(srcstatus, "inhomogBBDomain: error writing SRF file")
     if(.not. stateIsFailure(srcstatus))then
         ncStatus( 1) = nf90_EndDef(ncFileId)
         ncStatus( 2) = nf90_inq_varid(ncFileId, "Lambdas", ncVarID)
         ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, lambdas)
         ncStatus( 4) = nf90_inq_varid(ncFileId, "SRF", ncVarID)
         ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, SRF)
         ncStatus( 6) = nf90_close(ncFileId)
     end if
  end if

end program inhomogBBDomain

