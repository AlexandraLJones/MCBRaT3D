! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program monteCarloDriver
  ! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
  ! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Example-Drivers/monteCarloDriver.f95 $
  !Last edited 2011 by Maxwell Smith, Dept of Atmospheric Sciences, University of Illinois at Urbana-Champaign

  ! Monte Carlo radiative transfer program that computes top and bottom
  ! of domain pixel level outgoing fluxes, domain average absorption
  ! profile, 3D absorption fields, and domain top radiances from an 
  ! input "domain" file describing the optical properties for the medium.  
  ! Does monochromatic solar radiative transfer only.
  ! 
  ! The input parameters are specified with a namelist file. 
  !
  ! Outputs four types of results: 
  !    Pixel level and domain average fluxes
  !    Pixel level radiances for the specified directions
  !    Domain average absorption profiles (flux/km)
  !    3D absorption field (flux/km)
  ! Each output has the mean value and standard error of the mean
  ! over numBatches of photon batches (to estimate the Monte Carlo noise).
  !
  ! Assumptions: 
  !   Monochromatic solar radiative transfer.
  !   Lambertian surface reflection.
  !   Periodic horizontal boundary conditions.
  !   Uniform optical properties in each grid cell.

  !    Frank Evans    University of Colorado     December 2005
  !      Modifications by Robert Pincus, Climate Diagnostics Center, January 2006
  !      Modifications by Robert Pincus and Frank Evans, June-July 2006
  !      Adapted for MPI by Robert Pincus, August 2007
  !      Merged with single processor version by Robert Pincus, January 2009

  ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages
  use MultipleProcesses
  use RandomNumbers
  use scatteringPhaseFunctions
  use opticalProperties
  use monteCarloIllumination
  use monteCarloRadiativeTransfer
  use UserInterface
  use surfaceProperties
  use emissionAndBBWeights

  implicit none
 
  real, parameter :: Pi = acos(-1.)

  ! Input parameters
  !   Radiative transfer 
  real(8)                 :: solarFlux = 1., surfaceAlbedo = 0.
  real                 :: solarMu = 1., solarAzimuth = 0.,  LW_flag = -1. 
  real(8)                 :: lambda = 6.0, surfaceTemp = 300.0       ! added by Alexandra Jones Fall 2011. lambda in microns and surface temp in K
  integer, parameter   :: maxNumRad = 648 !18mus*36phis , oldVal:72
  real                 :: intensityMus(maxNumRad)  = 0., &
                          intensityPhis(maxNumRad) = 0.
  logical              :: angleFill = .false.
  real, dimension(3)   :: thetaFill = -1., phiFill = -1.
  real, allocatable, dimension(:) :: mus, phis
  integer              :: nMu, nPhi, numLambda, lambdaI
  !   Monte Carlo
  integer              :: numPhotonsPerBatch = 0, numBatches = 100, &
                          iseed = 10, nPhaseIntervals = 10001
                          
  !   Monte Carlo algorithmic choices 
  logical              :: useRayTracing = .true., useRussianRoulette = .true. 
  logical              :: useHybridPhaseFunsForIntenCalcs = .false. 
  real                 :: hybridPhaseFunWidth = 7. 
  integer              :: numOrdersOrigPhaseFunIntenCalcs = 0
  logical              :: useRussianRouletteForIntensity = .true. 
  real                 :: zetaMin = 0.3 
  logical              :: limitIntensityContributions = .false.
  real                 :: maxIntensityContribution    = 77.
  
  ! Control over output
  logical              :: reportVolumeAbsorption = .false., &
                          reportAbsorptionProfile = .false. 
  
  ! File names
  character(len=256)   :: domainFileName = "", domainFileList = "", solarSourceFile = ""
  character(len=256)   :: outputFluxFile = "", outputRadFile = "",  &
                          outputAbsProfFile = "", outputAbsVolumeFile = "", &
                          outputNetcdfFile = ""

  !Aux hist-- auxiliary output
  !auxhist01-fluxes organized by the order of scattering, and intensities where
  !   appropriate
  logical              :: recScatOrd = .false.
  integer              :: numRecScatOrd = 0
  character(len=256)   :: auxhist01_radFile=""
  character(len=256)   :: auxhist01_fluxFile=""

  namelist /radiativeTransfer/ solarFlux, solarMu, solarAzimuth, surfaceTemp, &
                               intensityMus, intensityPhis, angleFill, thetaFill, phiFill, LW_flag, numLambda
  
  namelist /monteCarlo/        numPhotonsPerBatch, numBatches, iseed, nPhaseIntervals
  
  namelist /algorithms/        useRayTracing, useRussianRoulette,                    &
                               useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                               numOrdersOrigPhaseFunIntenCalcs,                      &
                               useRussianRouletteForIntensity, zetaMin,              &
                               limitIntensityContributions, maxIntensityContribution
                        
  namelist /output/            reportVolumeAbsorption, reportAbsorptionProfile, &
                               recScatOrd, numRecScatOrd, &
                               auxhist01_fluxFile, auxhist01_radFile
  
  namelist /fileNames/         domainFileList, solarSourceFile, &
                               outputRadFile, outputFluxFile, &
                               outputAbsProfFile, outputAbsVolumeFile, outputNetcdfFile
  

   ! Local variables
  character(len=256)   :: namelistFileName, voxel_file, voxel_file2, horiz_file, level_file, col_file, row_file, diff_file, photon_file
  integer              :: nX, nY, nZ, phtn
  integer              :: i, j, k, batch, ix, iy, iz
  integer              :: numRadDir
  integer              :: numberOfComponents
  logical              :: computeIntensity
  real                 :: cpuTime0, cpuTime1, cpuTime2, cpuTimeTotal, cpuTimeSetup
  real                 :: meanFluxUp, meanFluxDown, meanFluxAbsorbed
  real(8)                 :: emittedFlux
  real                 :: meanFluxUpStats(2), meanFluxDownStats(2), meanFluxAbsorbedStats(2)
  real(8), allocatable    :: xPosition(:), yPosition(:), zPosition(:)
  real(8), allocatable    :: solarSourceFunction(:), centralLambdas(:)
  real, allocatable    :: fluxUp(:, :), fluxDown(:, :), fluxAbsorbed(:, :)
  real, allocatable    :: fluxUpStats(:, :, :), fluxDownStats(:, :, :), fluxAbsorbedStats(:, :, :)
  real, allocatable    :: absorbedProfile(:), absorbedProfileStats(:, :)
  real, allocatable    :: absorbedVolume(:, :, :), absorbedVolumeStats(:, :, :, :)
  real, allocatable    :: Radiance(:, :, :), RadianceStats(:, :, :, :)
  real, allocatable    :: meanFluxUpByScatOrd(:), meanFluxDownByScatOrd(:)
  real, allocatable    :: meanFluxUpByScatOrdStats(:,:), meanFluxDownByScatOrdStats(:,:)
  real, allocatable    :: fluxUpByScatOrd(:,:,:), fluxDownByScatOrd(:,:,:)
  real, allocatable    :: fluxUpByScatOrdStats(:,:,:,:), fluxDownByScatOrdStats(:,:,:,:)
  real, allocatable    :: intensityByScatOrd(:,:,:,:), intensityByScatOrdStats(:,:,:,:,:)
  integer              ::  N, atms_photons
  real(8), allocatable    :: voxel_weights(:,:,:,:), col_weights(:,:,:),level_weights(:,:)
  integer, allocatable   :: voxel_tallys1(:,:,:), voxel_tallys2(:,:,:), voxel_tallys1_sum(:,:,:), voxel_tallys2_sum(:,:,:), voxel_tallys1_total(:,:,:), voxel_tallys2_total(:,:,:), freqDistr(:)

   ! I3RC Monte Carlo code derived type variables
  type(domain)               :: thisDomain
  type(domain), allocatable, dimension(:) :: BBDomain
  type(ErrorMessage)         :: status
  type(randomNumberSequence) :: randoms
  type(photonStream)         :: incomingPhotons
  type(photonStream), allocatable, dimension(:) :: incomingBBPhotons
  type(integrator)           :: mcIntegrator
  type(Weights)              :: theseWeights

  ! Variables related to splitting up the job across processors
  integer            :: thisProc            ! ID OF CURRENT PROCESSOR; default 0
  integer            :: numProcs            ! TOTAL NUMBER OF PROCESSORS USED; default 1
  integer            :: batchesPerProcessor

  ! -----------------------------------------
  ! Start communications among multiple processes.
  !
!PRINT*, 'about to call initializeProcesses' 
  call initializeProcesses(numProcs, thisProc)
!PRINT*, 'returned from initializeProcesses'

! temporary output file for the purposes of debugging
  write(voxel_file, '(A,I0.4)') "voxel.out.",thisProc
  write(photon_file, '(A,I0.4)') "photon.out.",thisProc
  write(col_file, '(A,I0.4)') "col.out.",thisProc
  write(level_file, '(A,I0.4)') "level.out.",thisProc
!  write(horiz_file, '(A,I0.4)') "horiz.out.",thisProc
!  write(diff_file, '(A,I0.4)') "diff.out.",thisProc
!  write(voxel_file2, '(A,I0.4)') "voxel2.out.",thisProc

  open(unit=31, file=trim(voxel_file) , status='UNKNOWN')
  open(unit=32, file=trim(level_file) , status='UNKNOWN')
  open(unit=33, file=trim(col_file) , status='UNKNOWN')
!  open(unit=14, file=trim(row_file) , status='UNKNOWN')
!  open(unit=15, file=trim(horiz_file) , status='UNKNOWN')
!  open(unit=16, file=trim(diff_file) , status='UNKNOWN')
!  open(unit=17, file=trim(voxel_file2) , status='UNKNOWN')


  ! -----------------------------------------
  ! Get the input variables from the namelist file
  !
  call cpu_time(cpuTime0)
  namelistFileName = getOneArgument()
  open (unit = 1, file = trim(namelistFileName), status='OLD')
  read (1, nml = radiativeTransfer); rewind(1)
  read (1, nml = monteCarlo);        rewind(1)
  read (1, nml = algorithms);        rewind(1)
  read (1, nml = output);            rewind(1)
  read (1, nml = fileNames);         rewind(1)
  close (1)

! if Solar simulation we want to read in the arrays to store solar spectrum from the netCDF file
  if (LW_flag .lt. 0) then
     allocate(solarSourceFunction(1:numLambda))
     allocate(centralLambdas(1:numLambda))
     call read_SolarSource(solarSourceFile, numLambda, solarSourceFunction, centralLambdas, status=status)
PRINT *, 'Driver: Radmax= ', MAXVAL(solarSourceFunction), ' Radmin= ', MINVAL(solarSourceFunction)
     call printStatus(status)
  end if
!automatic assignment of observation angles. If angleFill
!mu and phi must be specified
  if(angleFill) then
  if((phiFill(3) >= 0.) .and. (thetaFill(3) >= 0. )  .and. &
     (thetaFill(2) .ge. thetaFill(1)) .and. (phiFill(2) .ge. phiFill(1)) )   then 
    nMu = int( (thetaFill(2) - thetafill(1))/thetafill(3)) +1
    nPhi = int((phiFill(2) - phiFill(1))/phiFill(3)) +1 
    
!    if(thetaFill(2).eq.thetaFill(1)) nMu=1
!    if(phiFill(2) .eq. phiFill(1)) nPhi = 1

    if (nMu >= 1 .and. nPhi >= 1) then
    allocate(mus(nMu), phis(nPhi))
    mus = thetaFill(1)+((/(REAL(N), N=0, nMu-1)/))*thetaFill(3) 
!     FORALL (i=0:nMu-1)
!        mus(i)=thetaFill(1) + (i*thetaFill(3))
!     END FORALL
      mus = cos(Pi/180. * mus)
    phis= phiFill(1)  +((/(REAL(N),N=0,nPhi-1)/)*phiFill(3)) 
!     FORALL (i=0:nPhi-1)
!        phis(i)=phiFill(1) + (i*phiFill(3))
!     END FORALL
    end if
 
    do i = 1,nMu
      do j = 1,nPhi
         k = (i-1)*nPhi+j
         intensityMus(k) = mus(i)
         intensityPhis(k) = phis(j)
      end do
    end do
  end if
  end if 
  numRadDir = count(abs(intensityMus(:)) > 0.) 
  computeIntensity = numRadDir > 0 .and. &
                     (len_trim(outputRadFile) > 0 .or. len_trim(outputNetcdfFile) > 0)
  if(.not. computeIntensity) outputRadFile = ""
  
  !set recScatOrd to false if the number requested is negative
  !no need to do the opposite because the integrator 
  !will listen only to numRecScatOrd if it is present
  if(numRecScatOrd < 0) recScatOrd=.false.

  ! -----------------------------------------
  !  Read the domain file
  !
  allocate(BBDomain(numLambda))
PRINT *, 'namelist numLambda= ', numLambda
!!!!!!!!!! Open list of domain files, then loop over reading them into BBDomain
  open(unit=22, file=TRIM(domainFileList), status='UNKNOWN')
! first we need to get the dimensions needed to allocate position arrays
  read(22,'(1A)') domainFileName
  call read_Domain(domainFileName, thisDomain, status)
  call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, namelistNumLambda=numLambda, domainNumLambda=numLambda, status = status)
  REWIND(22)
  allocate(xPosition(nx+1), yPosition(ny+1), zPosition(nz+1))

  DO i = 1, numLambda  !!!!!!!!!!! CAN I PARALLELIZE THIS LOOP? !!!!!!!!!!!!!!!!!!!!!!!!
     read(22,'(1A)') domainFileName
     call read_Domain(domainFileName, thisDomain, status)
PRINT *, 'Driver: Read Domain'  
  !call printStatus(status)
  !call write_Domain(thisDomain, 'test_write.dom', status)
     call printStatus(status)
   
     call getInfo_Domain(thisDomain, lambda = lambda, lambdaIndex = lambdaI,  &
                      xPosition = xPosition, yPosition = yPosition, &
                      zPosition = zPosition, numberOfComponents=numberOfComponents, status = status) 
!PRINT *, lambda, lambdaI
     BBDomain(lambdaI)=thisDomain
print *, 'Driver: got info for Domain: ', i
      call printStatus(status)
  END DO
  call getInfo_Domain(BBDomain(numLambda), lambda = lambda, lambdaIndex = lambdaI, status=status)
PRINT *, 'retrieved Domain info for ', lambda, lambdaI, numLambda
  call printStatus(status)

  close(22)
!
!  I MAY WANT TO IMPLEMENT AN EXPLICIT CHECK HERE TO MAKE SURE CERTAIN COMPONENTS OF ALL THE DOMAINS ARE CONSISTENT WITH EACHOTHER
!

  ! Set up the integrator object. It only grabs information from the domain that is consistent between all of them so I only need to pass any one domain
  mcIntegrator = new_Integrator(thisDomain, status = status)     
PRINT *, 'initialized Integrator'
  call printStatus(status)					 ! CAN I MOVE THIS SECTION INITIALIZING THE INTEGRATOR TO BEFORE THE LOOP RADING IN ALL THE DOMAINS BUT AFTER READING IN THE FIRST DOMAIN SO THAT I CAN JUST PARALLELIZE THE REST? 
!  call finalize_Domain(thisDomain)                               !

   ! Set the surface albedo, table sizes, and maybe the radiance directions
  call specifyParameters (mcIntegrator,                          &
                          minInverseTableSize = nPhaseIntervals, &
                          LW_flag = LW_flag,                     &
                          status = status)
PRINT *, 'Specified photon source'
  call printStatus(status) 

  if (computeIntensity) then
    call specifyParameters (mcIntegrator, &
                            minForwardTableSize=nPhaseIntervals, &
                            intensityMus=intensityMus(1:numRadDir), &
                            intensityPhis=intensityPhis(1:numRadDir), &
                            computeIntensity=computeIntensity,numComps=numberOfComponents, status=status)
PRINT *, 'Specified Intensity Calcs'
    call printStatus(status) 
  endif

  if(recScatOrd) then
    call specifyParameters (mcIntegrator, &
                            recScatOrd=recScatOrd, &
                            numRecScatOrd=numRecScatOrd, & 
                            status = status)
PRINT *, 'SPecified recScatOrder'
    call printStatus(status)
  end if 

  !
  ! Make the algorithmic choices
  !
  call specifyParameters(mcIntegrator,                              &
                         useRayTracing      = useRayTracing,        &
                         useRussianRoulette = useRussianRoulette,   &
                         status = status)
PRINT *, 'Specfied ray tracing'
  call printStatus(status) 
  
  !
  ! Algorithmic choices for intensity calculations
  !
  if (computeIntensity) then
    call specifyParameters(mcIntegrator,                            &
                         useHybridPhaseFunsForIntenCalcs =          &
                                 useHybridPhaseFunsForIntenCalcs,   &
                         hybridPhaseFunWidth = hybridPhaseFunWidth, &
                         numOrdersOrigPhaseFunIntenCalcs =          &
                                 numOrdersOrigPhaseFunIntenCalcs,   &
                         useRussianRouletteForIntensity =           &
                                 useRussianRouletteForIntensity,    &
                         zetaMin = zetaMin,                         &
                         limitIntensityContributions =              &
                                   limitIntensityContributions,     &
			 numComps = numberOfComponents,		    &
                         maxIntensityContribution =                 &
                                   maxIntensityContribution,        &
                         status = status)
PRINT *, 'Specified variance reduction'
    call printStatus(status) 
  end if

!PRINT *, 'Driver: Specified Parameters'
   ! Allocate and zero the arrays for radiative quantities and moments 
  allocate (voxel_tallys1(nX, nY, nZ), voxel_tallys1_sum(nX, nY, nZ), voxel_tallys1_total(nX, nY, nZ))
  allocate (voxel_tallys2(nX, nY, nZ), voxel_tallys2_sum(nX, nY, nZ), voxel_tallys2_total(nX, nY, nZ))
  allocate (fluxUp      (nX, nY), fluxUpStats      (nX, nY, 2))
  allocate (fluxDown    (nX, nY), fluxDownStats    (nX, nY, 2))
  allocate (fluxAbsorbed(nX, nY), fluxAbsorbedStats(nX, nY, 2))
  allocate (absorbedProfile(nZ), absorbedProfilestats(nZ, 2))
  allocate (absorbedVolume(nX, nY, nZ), absorbedVolumeStats(nX, nY, nZ, 2))
  voxel_tallys1(:,:,:)=0 ; voxel_tallys1_sum(:,:,:) = 0 ; voxel_tallys1_total(:,:,:) = 0
!  voxel_tallys2(:,:,:)=0 ; voxel_tallys2_sum(:,:,:) = 0 ; voxel_tallys2_total(:,:,:) = 0
  meanFluxUpStats(:) = 0.0  ; meanFluxDownStats(:) = 0.0  ; meanFluxAbsorbedStats(:) = 0.0
  fluxUpStats(:, :, :) = 0.0  ; fluxDownStats(:, :, :) = 0.0  ; fluxAbsorbedStats(:, :, :) = 0.0
  absorbedProfilestats(:, :) = 0.0 ;  absorbedVolumeStats(:, :, :, :) = 0.0
  if (computeIntensity) then
    allocate (Radiance(nX, nY, numRadDir), RadianceStats(nX, nY, numRadDir, 2))
    RadianceStats(:, :, :, :) = 0.0
  endif  
  
  if(recScatOrd .and. numRecScatOrd >=0) then 
    allocate(meanFluxUpByScatOrd(0:numRecScatOrd), meanFluxDownByScatOrd(0:numRecScatOrd),&
             fluxUpByScatOrd(nX,nY,0:numRecScatOrd),fluxDownByScatOrd(nX,nY,0:numRecScatOrd), &
             meanFluxUpByScatOrdStats(0:numRecScatOrd, 2), meanFluxDownByScatOrdStats(0:numRecScatOrd,2), &
             fluxUpByScatOrdStats(nX,nY, 0:numRecScatOrd, 2), fluxDownByScatOrdStats(nX,nY,0:numRecScatOrd,2))
    meanFluxUpByScatOrdStats=0.0;meanFluxDownByScatOrdStats=0.0;
    fluxUpByScatOrdStats=0.0;fluxDownByScatOrdStats=0.0;
    
    if(computeIntensity) then
      allocate(intensityByScatOrd(nX,nY,numRadDir,0:numRecScatOrd), &
               intensityByScatOrdStats(nX,nY,numRadDir, 0:numRecScatOrd, 2))
      intensityByScatOrdStats=0.0;
    end if
  end if

  ! --------------------------------------------------------------------------
  ! Compute radiative transfer with a trivial number of photons. 
  !   This checks to see if the integrator is properly set up before 
  !   running all the batches, and also allows the integrator to 
  !   do any internal pre-computation. 

  ! Seed the random number generator.
  randoms = new_RandomNumberSequence(seed = (/ iseed, 0 /) )
  allocate (freqDistr(1:numLambda), incomingBBPhotons(1:numLambda))
   ! The initial direction and position of the photons are precomputed and 
   !   stored in an "illumination" object.
!PRINT *, 'solarFlux=', solarFlux
  if(LW_flag >= 0.0)then
     theseWeights=new_Weights(numX=nX, numY=nY, numZ=nZ, numlambda=numLambda, status=status)
PRINT *, 'initialized thermal weights'
call printStatus(status)
     allocate (voxel_weights(nX,nY,nZ,numLambda),col_weights(nY,nZ,numLambda), level_weights(nZ,numLambda))
!     call getInfo_Domain(thisDomain, temps=temps, ssa=ssa, totalExt=cumExt, ext=ext, status=status)
!PRINT *, 'retrieved fields for emission weighting'
!call printStatus(status)

     call emission_weighting(BBDomain, numLambda, theseWeights, surfaceTemp, numPhotonsPerBatch, atms_photons,voxel_weights, col_weights=col_weights, level_weights=level_weights, totalFlux=emittedFlux, status=status) 
PRINT *, 'returned from emission_weighting'

!    write(32,"(36F12.8)") level_weights(:,1)
!    DO k = 1, nZ
!      write(33,"(100F12.8)") col_weights(:,k,1)
!      DO j = 1, nY
!        write(31,"(100F12.8)") voxel_weights(:,j,k,1)
!      end do
!    end do
!    close(31)
!    close(32)
!    close(33)

!     solarFlux=31.25138117141156822262   !!!MAKE SURE TO COMMENT OUT THIS LINE. DIAGNOSTICE PURPOSES ONLY!!!
!PRINT *, 'total atms photons=', atms_photons)
PRINT *, 'emittedFlux=', emittedFlux, ' solarFlux=', solarFlux
call printStatus(status)
!PRINT *, 'Driver: calculated emission weighting'
     call getFrequencyDistr(theseWeights, numPhotonsPerBatch*numBatches,randoms, freqDistr) ! TECHNICALLY NUMPHOTONSPERBATCH is NOT BE CORRECT but this is just for the set up step so i think its OK
PRINT*, freqDistr
     DO i = 1, numLambda
        incomingBBPhotons(i) = new_PhotonStream (theseWeights=theseWeights, iLambda=i, numberOfPhotons=freqDistr(i),randomNumbers=randoms, status=status, option1=voxel_tallys1)
	voxel_tallys1_sum = voxel_tallys1_sum + voxel_tallys1 
     END DO
PRINT *, 'initialized thermal BB photon stream'
call printStatus(status)
    open(unit=12, file=trim(photon_file) , status='UNKNOWN')
    DO k = 1, nZ
      DO j = 1, nY
        write(12,"(100I12)") voxel_tallys1_sum(:,j,k)
      end do
    end do
    close(12)
     incomingPhotons = new_PhotonStream (numberOfPhotons=1, atms_photons=atms_photons, voxel_weights=voxel_weights(:,:,:,1), col_weights=col_weights(:,:,1), level_weights=level_weights(:,1), nX=nX, nY=nY, nZ=nZ, randomNumbers=randoms, status=status)  
!PRINT *, 'LW', ' incomingPhotons%SolarMu=', incomingPhotons%solarMu(1)
PRINT *, 'initialized single photon'
call printStatus(status)
!PRINT *, 'Driver: initialized single photon'
  else
     theseWeights=new_Weights(numLambda=numLambda, status=status)
PRINT *, 'initialized solar weights'
call printStatus(status)
     call solar_Weighting(theseWeights, numLambda, solarSourceFunction, centralLambdas, solarMu, emittedFlux, status=status)   ! convert the solar source function to CDF and total Flux
PRINT *, 'filled solar weights'
call printStatus(status)
     deallocate(solarSourceFunction)
     deallocate(centralLambdas)
     call getFrequencyDistr(theseWeights, numPhotonsPerBatch*numBatches,randoms, freqDistr)
PRINT *, freqDistr
     DO i = 1, numLambda
       incomingBBPhotons(i) = new_PhotonStream (solarMu, solarAzimuth, &
                                      numberOfPhotons = freqDistr(i),   &
                                      randomNumbers = randoms, status=status)
     END DO
!PRINT *, 'not LW', 'incomingPhotons%SolarMu=', incomingPhotons%solarMu(1)
PRINT *, 'initialized SW photon stream'
call printStatus(status)
     incomingPhotons = new_PhotonStream (solarMu, solarAzimuth, &
                                      numberOfPhotons = 1,   &
                                      randomNumbers = randoms, status=status)
call printStatus(status)
PRINT *, 'initialized SW test photon'     
  end if
  solarFlux = emittedFlux
!PRINT *, 'incomingPhotons%solarMu=', incomingPhotons%solarMu(1)
!  call finalize_Domain(thisDomain)
STOP

  ! Now we compute the radiative transfer for a single photon 
  if(.not. isReady_Integrator (mcIntegrator)) stop 'Integrator is not ready.'
  call computeRadiativeTransfer (mcIntegrator,thisDomain, randoms, incomingPhotons, status, voxel_tallys2)
PRINT *, 'Computed RT for a single photon'
  call printStatus(status) 
  call finalize_PhotonStream (incomingPhotons)
!PRINT *, 'Driver: succesfully tested photon initialization'

  call cpu_time(cpuTime1)
!PRINT *, "called cpu_time"
  call synchronizeProcesses
!PRINT *, "called synchronizeProcesses"
  cpuTimeSetup = sumAcrossProcesses(cpuTime1 - cpuTime0) 
PRINT *, "cpuTimeSetup=", cpuTimeSetup
STOP
  if (MasterProc) &
    print *, "Setup CPU time (secs, approx): ", int(cpuTimeSetup)
!PRINT *, "setup completed"

  ! --------------------------------------------------------------------------

  ! The  loop over batches is for estimating the uncertainty in the flux and
  !   radiance from the variance between numBatches independent calculations. 
  numBatches = max(numBatches,2)
  batchesPerProcessor = numBatches/numProcs
  ! If the number of batches doesn't divide among the processors evenly increase the 
  !   number until it does. 
  if(mod(numBatches, numProcs) /= 0) then 
    batchesPerProcessor = batchesPerProcessor + 1
    numBatches = batchesPerProcessor * numProcs
  end if 
  if (MasterProc) &
    print *, "Doing ", batchesPerProcessor, " batches on each of ", numProcs, " processors." 
!PRINT *, "entering batch loop"
  batches: do batch = thisProc*batchesPerProcessor + 1, thisProc*batchesPerProcessor + batchesPerProcessor
    ! Seed the random number generator.
    !   Variable randoms holds the state of the random number generator. 
    randoms = new_RandomNumberSequence(seed = (/ iseed, batch /) )

    ! The initial direction and position of the photons are precomputed and 
    !   stored in an "illumination" object. 
    if(LW_flag >= 0.0)then
       incomingPhotons = new_PhotonStream (numberOfPhotons=numPhotonsPerBatch, atms_photons=atms_photons, voxel_weights=voxel_weights(:,:,:,1), col_weights=col_weights(:,:,1),&
level_weights=level_weights(:,1), nX=nX, nY=nY, nZ=nZ, randomNumbers=randoms, status=status, option1=voxel_tallys1)  
!DO phtn=1,numPhotonsPerBatch
!   PRINT *, incomingPhotons%xPosition(phtn), incomingPhotons%yPosition(phtn), incomingPhotons%zPosition(phtn)
!ENDDO      
    else
       incomingPhotons = new_PhotonStream (solarMu, solarAzimuth,                &
                                        numberOfPhotons = numPhotonsPerBatch, &
                                        randomNumbers = randoms, status = status)
    end if
    call printStatus(status)
!PRINT *, 'Driver: Initialized photons for batch', batch
    ! Now we compute the radiative transfer for this batch of photons. 
    call computeRadiativeTransfer (mcIntegrator, thisDomain,randoms, incomingPhotons, status, voxel_tallys2)
!PRINT *, 'Driver: COmputed RT for batch', batch
     ! Get the radiative quantities:
     !   This particular integrator provides fluxes at the top and bottom 
     !   of the domain for both the domain mean and pixel level fluxes,
     !   the absorbed flux profile, 3D field of absorbed flux, and
     !   the pixel level radiances at top and/or bottom of domain.
    call reportResults (mcIntegrator, &
           meanFluxUp=meanFluxUp, meanFluxDown=meanFluxDown, meanFluxAbsorbed=meanFluxAbsorbed,       &
           fluxUp=fluxUp(:, :), fluxDown=fluxDown(:, :), fluxAbsorbed=fluxAbsorbed(:, :), &
           absorbedProfile=absorbedProfile(:), volumeAbsorption=absorbedVolume(:, :, :), status = status)

     ! Accumulate the first and second moments of each quantity over the batches 
!    voxel_tallys1_sum = voxel_tallys1_sum + voxel_tallys1
!    voxel_tallys2_sum = voxel_tallys2_sum + voxel_tallys2
    meanFluxUpStats(1)       = meanFluxUpStats(1)       + solarFlux*meanFluxUp
    meanFluxUpStats(2)       = meanFluxUpStats(2)       + (solarFlux*meanFluxUp)**2
    meanFluxDownStats(1)     = meanFluxDownStats(1)     + solarFlux*meanFluxDown
    meanFluxDownStats(2)     = meanFluxDownStats(2)     + (solarFlux*meanFluxDown)**2
    meanFluxAbsorbedStats(1) = meanFluxAbsorbedStats(1) + solarFlux*meanFluxAbsorbed
    meanFluxAbsorbedStats(2) = meanFluxAbsorbedStats(2) + (solarFlux*meanFluxAbsorbed)**2
          fluxUpStats(:, :, 1) =       fluxUpStats(:, :, 1) + solarFlux*fluxUp(:, :)
          fluxUpStats(:, :, 2) =       fluxUpStats(:, :, 2) + (solarFlux*fluxUp(:, :))**2
        fluxDownStats(:, :, 1) =     fluxDownStats(:, :, 1) + solarFlux*fluxDown(:, :)
        fluxDownStats(:, :, 2) =     fluxDownStats(:, :, 2) + (solarFlux*fluxDown(:, :))**2
    fluxAbsorbedStats(:, :, 1) = fluxAbsorbedStats(:, :, 1) + solarFlux*fluxAbsorbed(:, :)
    fluxAbsorbedStats(:, :, 2) = fluxAbsorbedStats(:, :, 2) + (solarFlux*fluxAbsorbed(:, :))**2
    absorbedProfileStats(:, 1) = absorbedProfileStats(:, 1) + solarFlux*absorbedProfile(:)
    absorbedProfileStats(:, 2) = absorbedProfileStats(:, 2) + (solarFlux*absorbedProfile(:))**2
    absorbedVolumeStats(:, :, :, 1) = absorbedVolumeStats(:, :, :, 1) + solarFlux*absorbedVolume(:, :, :)
    absorbedVolumeStats(:, :, :, 2) = absorbedVolumeStats(:, :, :, 2) + (solarFlux*absorbedVolume(:, :, :))**2

    if (computeIntensity) then
      call reportResults(mcIntegrator, intensity = Radiance(:, :, :), status = status)
      RadianceStats(:, :, :,1) = RadianceStats(:, :, :,1) + solarFlux*Radiance(:, :, :)
!PRINT *, "mean RadianceStats =", sum (RadianceStats(:, :, 1,1))/real (nX * nY)
      RadianceStats(:, :, :,2) = RadianceStats(:, :, :,2) + (solarFlux*Radiance(:, :, :))**2
    endif

   if(recScatOrd) then 
      call reportResults(mcIntegrator, &
             meanFluxUpByScatOrd=meanFluxUpByScatOrd, meanFluxDownByScatOrd=meanFluxDownByScatOrd, &
             fluxUpByScatOrd=fluxUpByScatOrd, fluxDownByScatOrd=fluxDownByScatOrd, status=status)
      !accumulate quantities as above
      meanFluxUpByScatOrdStats(:,1) = meanFluxUpByScatOrdStats(:,1) + solarFlux*meanFluxUpByScatOrd(:)
      meanFluxUpByScatOrdStats(:,2) = meanFluxUpByScatOrdStats(:,2) + (solarFlux*meanFluxUpByScatOrd(:))**2
      meanFluxDownByScatOrdStats(:,1) = meanFluxDownByScatOrdStats(:,1) +  solarFlux*meanFluxDownByScatOrd(:)
      meanFluxDownByScatOrdStats(:,2) = meanFluxDownByScatOrdStats(:,2) + (solarFlux*meanFluxUpByScatOrd(:))**2
      meanFluxDownByScatOrdStats(:,1) = meanFluxDownByScatOrdStats(:,1) +  (solarFlux*meanFluxDownByScatOrd(:))**2
            fluxUpByScatOrdStats(:,:,:,1) = fluxUpByScatOrdStats(:,:,:,1) +  solarFlux*fluxUpByScatOrd(:,:,:)
            fluxUpByScatOrdStats(:,:,:,2) = fluxUpByScatOrdStats(:,:,:,2) + ( solarFlux*fluxUpByScatOrd(:,:,:))**2
          fluxDownByScatOrdStats(:,:,:,1) = fluxDownByScatOrdStats(:,:,:,1) +  solarFlux*fluxDownByScatOrd(:,:,:)
          fluxDownByScatOrdStats(:,:,:,2) = fluxDownByScatOrdStats(:,:,:,2) + ( solarFlux*fluxDownByScatOrd(:,:,:))**2
      if(computeIntensity) then
        call reportResults(mcIntegrator, intensityByScatOrd=intensityByScatOrd, status=status)
        intensityByScatOrdStats(:,:,:,:,1) = intensityByScatOrdStats(:,:,:,:,1) +  solarFlux*intensityByScatOrd(:,:,:,:)
        intensityByScatOrdStats(:,:,:,:,2) = intensityByScatOrdStats(:,:,:,:,2) + ( solarFlux*intensityByScatOrd(:,:,:,:))**2
      end if
   end if
!PRINT *, 'Driver: Reported results for batch', batch
     ! Release the photon "illumination" object memory
    call finalize_PhotonStream (incomingPhotons)
    call printStatus(status)
  end do batches

  if (allocated(voxel_weights)) deallocate (voxel_weights)
  if (allocated(col_weights)) deallocate (col_weights)
  if (allocated(level_weights)) deallocate (level_weights)  
  !
  ! Accumulate statistics from across all the processors
  !
! voxel_tallys1_total = sumAcrossProcesses(voxel_tallys1_sum)
! voxel_tallys2_total = sumAcrossProcesses(voxel_tallys2_sum)

  ! Domain-mean fluxes
  
  meanFluxUpStats(:)       = sumAcrossProcesses(meanFluxUpStats)
  meanFluxDownStats(:)     = sumAcrossProcesses(meanFluxDownStats)
  meanFluxAbsorbedStats(:) = sumAcrossProcesses(meanFluxAbsorbedStats)

  ! Pixel-by-pixel fluxes
  fluxUpStats(:, :, :)       = sumAcrossProcesses(fluxUpStats)
  fluxDownStats(:, :, :)     = sumAcrossProcesses(fluxDownStats)
  fluxAbsorbedStats(:, :, :) = sumAcrossProcesses(fluxAbsorbedStats)
  
  ! Absorption (mean profile and cell-by-cell)
  absorbedProfileStats(:, :)      = sumAcrossProcesses(absorbedProfileStats)
  absorbedVolumeStats(:, :, :, :) = sumAcrossProcesses(absorbedVolumeStats)
  
  ! Radiance
  if (computeIntensity) &
    RadianceStats(:, :, :, :) = sumAcrossProcesses(RadianceStats)
!PRINT *, "mean total RadianceStats =", sum(RadianceStats(:, :, 1,1))/real (nX * nY)  

   if(recScatOrd) then
     meanFluxUpByScatOrdStats(:,:) = sumAcrossProcesses(meanFluxUpByScatOrdStats)
     meanFluxDownByScatOrdStats(:,:) = sumAcrossProcesses(meanFluxDownByScatOrdStats)
     fluxUpByScatOrdStats(:,:,:,:) = sumAcrossProcesses(fluxUpByScatOrdStats)
     fluxDownByScatOrdStats(:,:,:,:) = sumAcrossProcesses(fluxDownByScatOrdStats)
    
     if(computeIntensity) then
     intensityByScatOrdStats(:,:,:,:,:) = sumAcrossProcesses(intensityByScatOrdStats)
     end if
   end if

!PRINT *, 'Driver: accumulated results'
!  close(11)
!  close(12)
!  close(13)
!  close(14)
!  close(15)
!  close(16)
!  close(17)

  call synchronizeProcesses
  call cpu_time(cpuTime2)
  cpuTimeTotal = sumAcrossProcesses(cpuTime2 - cpuTime0)
  call finalizeProcesses

  if (MasterProc) print *, "Total CPU time (secs, approx): ", int(cpuTimeTotal)

   ! Calculate the mean and standard error of the radiative quantities from the two moments
  meanFluxUpStats(:)       = meanFluxUpStats(:)/numBatches
  meanFluxUpStats(2)       = sqrt( max(0.0, meanFluxUpStats(2) - meanFluxUpStats(1)**2) /(numBatches-1))
  meanFluxDownStats(:)     = meanFluxDownStats(:)/numBatches
  meanFluxDownStats(2)     = sqrt( max(0.0, meanFluxDownStats(2) - meanFluxDownStats(1)**2) /(numBatches-1))
  meanFluxAbsorbedStats(:) = meanFluxAbsorbedStats(:)/numBatches
  meanFluxAbsorbedStats(2) = sqrt( max(0.0, meanFluxAbsorbedStats(2) - meanFluxAbsorbedStats(1)**2) /(numBatches-1))
  fluxUpStats(:, :, :)       =fluxUpStats(:, :, :)/numBatches
  fluxUpStats(:, :, 2)       = sqrt( max(0.0, fluxUpStats(:, :,2) - fluxUpStats(:, :,1)**2) /(numBatches-1))
  fluxDownStats(:, :, :)     =fluxDownStats(:, :, :)/numBatches
  fluxDownStats(:, :, 2)     = sqrt( max(0.0, fluxDownStats(:, :,2) - fluxDownStats(:, :,1)**2) /(numBatches-1))
  fluxAbsorbedStats(:, :, :) =fluxAbsorbedStats(:, :, :)/numBatches
  fluxAbsorbedStats(:, :, 2) = sqrt( max(0.0, fluxAbsorbedStats(:, :,2) - fluxAbsorbedStats(:, :,1)**2) /(numBatches-1))
  absorbedProfileStats(:, :) =absorbedProfileStats(:, :)/numBatches
  absorbedProfileStats(:, 2) = sqrt( max(0.0, absorbedProfileStats(:,2) - absorbedProfileStats(:,1)**2) /(numBatches-1))
  absorbedVolumeStats(:, :, :, :) = absorbedVolumeStats(:, :, :, :)/numBatches
  absorbedVolumeStats(:, :, :, 2) = sqrt( max(0.0, absorbedVolumeStats(:, :, :, 2) - absorbedVolumeStats(:, :, :,1)**2) / &
                                    (numBatches-1))
  if (computeIntensity) then
    RadianceStats(:, :, :, :) = RadianceStats(:, :, :, :)/numBatches
!PRINT *, "mean Radiance stats including solarflux =", sum (RadianceStats(:, :, 1,1))/real (nX * nY)
    RadianceStats(:, :, :,2) = sqrt( max(0.0, RadianceStats(:, :, :,2) - RadianceStats(:, :, :,1)**2) /(numBatches-1))
  endif

  if(recScatOrd) then 
    meanFluxUpByScatOrdStats(:,:) =meanFluxUpByScatOrdStats(:,:)/numBatches
    meanFluxUpByScatOrdStats(:,2) = sqrt(max(0.0, meanFluxUpByScatOrdStats(:,2) - meanFluxUpByScatOrdStats(:,1)**2) / (numBatches-1))
    meanFluxDownByScatOrdStats(:,:) = meanFluxDownByScatOrdStats(:,:)/numBatches
    meanFluxDownByScatOrdStats(:,2) = sqrt(max(0.0, meanFluxDownByScatOrdStats(:,2) - meanFluxDownByScatOrdStats(:,1)**2) / (numBatches-1))
    fluxUpByScatOrdStats(:,:,:,:) = fluxUpByScatOrdStats(:,:,:,:)/numBatches
    fluxUpByScatOrdStats(:,:,:,2) = sqrt(max(0.0, fluxUpByScatOrdStats(:,:,:,2) - fluxUpByScatOrdStats(:,:,:,1)**2) / (numBatches-1))
    fluxDownByScatOrdStats(:,:,:,:) =fluxDownByScatOrdStats(:,:,:,:)/numBatches
    fluxDownByScatOrdStats(:,:,:,2) = sqrt(max(0.0, fluxDownByScatOrdStats(:,:,:,2) - fluxDownByScatOrdStats(:,:,:,1)**2) / (numBatches-1))
    
    if(computeIntensity) then
      intensityByScatOrdStats(:,:,:,:,:) = intensityByScatOrdStats(:,:,:,:,:)/numBatches
      intensityByScatOrdStats(:,:,:,:,2) = sqrt(max(0.0, intensityByScatOrdStats(:,:,:,:,2)-intensityByScatOrdStats(:,:,:,:,1)**2) / &
                                           (numBatches-1))
    end if
  end if

!PRINT *, 'Driver: calculated radiative quantities'
  if(MasterProc) then ! Write a single output file. 
!    open(unit=12, file=trim(photon_file) , status='UNKNOWN')

!    DO k = 1, nZ
!      DO j = 1, nY
!        write(12,"(100I12)") voxel_tallys1_total(:,j,k)
!      end do
!    end do
!    close(12)

    if(any( (/ len_trim(outputFluxFile),      len_trim(outputAbsProfFile), &
               len_trim(outputAbsVolumeFile), len_trim(outputRadFile)      /) > 0)) then 
      call writeResults_ASCII(domainFileName,  numPhotonsPerBatch, numBatches,      &
                              useRayTracing, useRussianRoulette,                    &
                              useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                              solarFlux, solarMu, solarAzimuth, surfaceAlbedo,      &
                              xPosition, yPosition, zPosition,                      &
                              outputFluxFile, meanFluxUpStats, meanFluxDownStats,   &
                              meanFluxAbsorbedStats,                                &
                              fluxUpStats, fluxDownStats, fluxAbsorbedStats,        &
                              outputAbsProfFile, absorbedProfileStats,              &
                              outputAbsVolumeFile, absorbedVolumeStats,             &
                              outputRadFile, intensityMus, intensityPhis, RadianceStats)
      print *, "Wrote ASCII results"  
    end if 
  
    if(len_trim(outputNetcdfFile) > 0) then
      if(computeIntensity .and. .not. recScatOrd) then 
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats)      
      elseif(computeIntensity .and. recScatOrd) then
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats,      &
                                 numRecScatOrd,                                   &
                                 fluxUpByScatOrdStats, fluxDownByScatOrdStats,    &
                                 intensityByScatOrdStats                          )
      elseif(.not. computeIntensity .and. recScatOrd) then
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 numRecScatOrd=numRecScatOrd,                     &
                                 fluxUpByScatOrdStats=fluxUpByScatOrdStats,       &
                                 fluxDownByScatOrdStats=fluxDownByScatOrdStats      )
      else
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats)
  
      end if 
      print *, "Wrote netcdf results"                         
    end if
  end if

  !
  ! Release all the memory. We should be able to finalize the derived types before we write 
  !   the output but this fails intermittently, so we want to be sure to get our results 
  !   before we take a chance on blowing up. 
  ! 
  deallocate (fluxUp, fluxUpStats, fluxDown, fluxDownStats, fluxAbsorbed, fluxAbsorbedStats)
  deallocate (absorbedProfile, absorbedProfileStats, absorbedVolume, absorbedVolumeStats)
  if (computeIntensity) deallocate (Radiance, RadianceStats)

  call finalize_RandomNumberSequence(randoms)
  call finalize_Integrator (mcIntegrator)

contains

! -------------------------------------------------------------------------------
  subroutine writeResults_ASCII(domainFileName,  numPhotonsPerBatch, numBatches,&
                          useRayTracing, useRussianRoulette,                    &
                          useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                          solarFlux, solarMu, solarAzimuth, surfaceAlbedo,      &
                          xPosition, yPosition, zPosition,                      &
                          outputFluxFile, meanFluxUpStats, meanFluxDownStats,   &
                          meanFluxAbsorbedStats,                                &
                          fluxUpStats, fluxDownStats, fluxAbsorbedStats,        &
                          outputAbsProfFile, absorbedProfileStats,              &
                          outputAbsVolumeFile, absorbedVolumeStats,             &
                          outputRadFile, intensityMus, intensityPhis, RadianceStats)
    !
    ! Writes Monte Carlo results to ASCII files. 
    !   Fluxes, absorption profiles, 3D absorption, and intensities are written to separate files.
    !
    ! Variables describing the problem 
    character(len=*),   intent(in) :: domainFileName
    integer,            intent(in) :: numPhotonsPerBatch, numBatches
    logical,            intent(in) :: useRayTracing, useRussianRoulette, useHybridPhaseFunsForIntenCalcs
    real,               intent(in) :: hybridPhaseFunWidth
    real,               intent(in) ::  solarMu, solarAzimuth
    real(8),    intent(in) :: solarFlux,surfaceAlbedo
    real(8), dimension(:), intent(in) :: xPosition, yPosition, zPosition
    ! Flux variables
    character(len = *),       intent(in) :: outputFluxFile
    real, dimension(:),       intent(in) :: meanFluxUpStats, meanFluxDownStats, meanFluxAbsorbedStats
    real, dimension(:, :, :), intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    character(len = *),       intent(in) :: outputAbsProfFile
    real, dimension(:, :),    intent(in) :: absorbedProfileStats
    character(len = *),       intent(in) :: outputAbsVolumeFile
    real, dimension(:, :, :, :), &
                              intent(in) :: absorbedVolumeStats
    ! Intensity variable
    character(len = *),       intent(in) :: outputRadFile
    real, dimension(:),       intent(in) :: intensityMus, intensityPhis
    real, dimension(:, :, :, :), &
                              intent(in) :: RadianceStats

    !
    ! Local variables
    !
    integer  :: i, j, k, nx, ny, nz, numRadDir

    nx = size(xPosition) - 1; ny = size(yPosition) - 1; nz = size(zPosition) - 1

    ! -----------------------------------------
    !  Output the flux results to an ASCII file
    !
    if(len_trim(outputFluxFile) > 0) then 
      open (unit=2, file=outputFluxFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Flux'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', INT(numPhotonsPerBatch * numBatches, KIND=8)
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Pixel Flux'
      write (2,'(A,F7.3,A,F7.3)') '!  Upwelling_Level=', zPosition(nZ+1), &
                                  '   Downwelling_level=', zPosition(1)
      write (2,'(A)') '!   X      Y           Flux_Up             Flux_Down            Flux_Absorbed '
      write (2,'(A)') '!                  Mean     StdErr       Mean     StdErr       Mean     StdErr'
      write (2,'(A14,3(1X,2(1X,F9.4)))') '!  Average:   ', meanFluxUpStats(1:2), &
                             meanFluxDownStats(1:2), meanFluxAbsorbedStats(1:2)
      do j = 1, nY
        do i = 1, nX  
          write (2,'(2(F7.3),3(1X,2(1X,F9.4)))') &
               sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., &
               fluxUpStats(i,j,1:2), fluxDownStats(i,j,1:2), fluxAbsorbedStats(i,j,1:2)
        enddo
      enddo  
      close (2)
      end if 


    ! -----------------------------------------
    !  Output the absorption profile results to a file
    !
    if(len_trim(outputAbsProfFile) > 0) then 
      open (unit=2, file=outputAbsProfFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Absorption Profile'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', INT(numPhotonsPerBatch * numBatches, KIND=8)
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Absorption Profile'
      write (2,'(A)') '!   Z    Absorbed_Flux (flux/km) '
      write (2,'(A)') '!          Mean     StdErr '
      do k = 1, nZ
        write (2,'(F7.3,1X,2(1X,F9.4))') 0.5*(zPosition(k)+zPosition(k+1)), absorbedProfileStats(k,1:2) 
      enddo  
      close (2)
    end if


    ! -----------------------------------------
    !  Output the volume absorption results to a file
    !
    if(len_trim(outputAbsVolumeFile) > 0) then 
      open (unit=2, file=outputAbsVolumeFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: 3D Absorption Field'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', INT(numPhotonsPerBatch * numBatches, KIND=8)
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Volume Absorption '
      write (2,'(A)') '!    X       Y        Z       Absorbed_Flux (flux/km)'
      write (2,'(A)') '!                               Mean     StdErr '
      do i = 1, nX  
        do j = 1, nY
          do k = 1, nZ
            write (2,'(3(F7.3,1X),2(1X,F9.4))') &
              sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., sum(zPosition(k:k+1))/2., &
              absorbedVolumeStats(i,j,k,1:2) 
          enddo  
        enddo  
      enddo  
      close (2)
    end if

    ! -----------------------------------------
    !  Output the radiance results to a file
    !
    if (len_trim(outputRadFile) > 0) then
      numRadDir = count(abs(intensityMus(:)) > 0.)
      open (unit=2, file=outputRadFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Radiance'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', INT(numPhotonsPerBatch * numBatches, KIND=8)
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,L1,A,F5.2)') '!  Intensity_uses_Russian_Roulette=', useRussianRouletteForIntensity, &
                                '   Intensity_Russian_Roulette_zeta_min=',zetaMin
      write (2,'(A,L1,A,F5.2)') '!  limited_intensity_contributions=', limitIntensityContributions, &
                                '   max_intensity_contribution=',maxIntensityContribution
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                        '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Pixel Radiance'
      write (2,'(A,F7.3,3(A,I4))') '!  RADIANCE AT Z=',  zPosition(nZ+1), &
              '   NXO=',nX, '   NYO=',nY, '   NDIR=',numRadDir
      write (2,'(A)') '!   X      Y         Radiance (Mean, StdErr)'
      do k = 1, numRadDir
        WRITE (2,'(A,1X,F8.5,1X,F6.2,2X,A)') &
            '! ', intensityMus(k), intensityPhis(k), '<- (mu,phi)'
        do j = 1, nY
          do i = 1, nX  
            write (2,'(2(F7.3),2(1X,F9.4))') sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., &
               RadianceStats(i,j,k,1:2)
          enddo
        enddo  
      enddo  
      close (2)
    endif
  end subroutine writeResults_ASCII


! -------------------------------------------------------------------------------
  subroutine writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputFileName,                                  &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats,      &
                                 numRecScatOrd,                                   &
                                 fluxUpByScatOrdStats, fluxDownByScatOrdStats,    &
                                 intensityByScatOrdStats                          )
    use netcdf
    !
    ! Writes Monte Carlo results to a netCDF file. 
    !   Writes top-of-dmain fluxes, absorption profiles, and intensities 
    !   Monte Carlo noise for each quantitiy is estimated as the std. dev. across the nBatches. 
    !
    ! Variables describing the problem 
    character(len = *), intent(in) :: domainFileName
    character(len = *), intent(in) :: outputFileName
    integer,            intent(in) :: numPhotonsPerBatch, numBatches
    real,               intent(in) ::  solarMu, solarAzimuth
    real(8),               intent(in) :: solarFlux, surfaceAlbedo
    ! The position variables mark the edges of the cells - we'll report the 
    !   results at the midpoints. 
    real(8), dimension(:), intent(in) :: xPosition, yPosition, zPosition

    ! Flux variables
    real, dimension(:, :, :),     intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    real, dimension(:, :),        intent(in) :: absorbedProfileStats
    real, dimension(:, :, :, :),  intent(in) :: absorbedVolumeStats
    ! Intensity variable
    real, dimension(:), optional, intent(in) :: intensityMus, intensityPhis
    real, dimension(:, :, :, :), &
                        optional, intent(in) :: RadianceStats
    integer,            optional, intent(in) :: numRecScatOrd
    real, dimension(:,:,:,:), &
                        optional, intent(in) :: fluxUpByScatOrdStats, fluxDownByScatOrdStats
    real, dimension(:,:,:,:,:), &
                        optional, intent(in) :: intensityByScatOrdStats

    !
    ! Local variables
    !
    logical                :: computeIntensity
    integer, dimension(32) :: ncStatus = nf90_NoErr
    integer                :: ncFileId, xDimId, yDimId, zDimId, dirDimId, ncVarId, &
                              numX, numY, numZ, numIntensityDirs, scatDimId
    logical                :: recScatOrd
    integer                :: i, N

!    real, dimension(numRecScatOrd+1)  :: scatOrderHolder = (/(REAL(i),i=0,numRecScatOrd)/)
    ! ---------------------------------------------------------------------------
    numX = size(xPosition) - 1; numY = size(yPosition) - 1; numZ = size(zPosition) - 1
    computeIntensity = present(intensityMus) .and. present(intensityPhis) .and. &
                       present(RadianceStats)
    recScatOrd       = present(numRecScatOrd) .and. present(fluxUpByScatOrdStats) .and. &
                       present(fluxDownByScatOrdStats)

    ncStatus( 1) = nf90_create(trim(outputFileName), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncFileId)

    ! Store the simulation parameters as global attributes
    !
    ncStatus( 5) = nf90_put_att(ncFileId, NF90_Global, "description", &
                                                       "Output from I3RC Community Monte Carlo Model")
    ncStatus( 6) = nf90_put_att(ncFileId, NF90_Global, "Domain_filename", trim(domainFileName))
    ncStatus( 7) = nf90_put_att(ncFileId, NF90_Global, "Surface_albedo", surfaceAlbedo)
    ncStatus( 8) = nf90_put_att(ncFileId, NF90_Global, "Total_number_of_photons", &
                                                       INT(numPhotonsPerBatch * numBatches, KIND=8))
    ncStatus( 9) = nf90_put_att(ncFileId, NF90_Global, "Number_of_batches", numBatches)
    ncStatus(10) = nf90_put_att(ncFileId, NF90_Global, "Solar_flux", SolarFlux)
    ncStatus(11) = nf90_put_att(ncFileId, NF90_Global, "Solar_mu", SolarMu)
    ncStatus(12) = nf90_put_att(ncFileId, NF90_Global, "Solar_phi", SolarAzimuth)

    ncStatus(13) = nf90_put_att(ncFileId, NF90_Global, "Random_number_seed", iseed)
    ncStatus(14) = nf90_put_att(ncFileId, NF90_Global, "Phase_function_table_sizes", &
                                                       nPhaseIntervals)
    if(useRayTracing) then
      ncStatus(15) = nf90_put_att(ncFileId, NF90_Global, "Algorithm", "Ray_tracing")
    else
      ncStatus(15) = nf90_put_att(ncFileId, NF90_Global, "Algorithm", "Max_cross_section")
    end if

    if(useHybridPhaseFunsForIntenCalcs) then 
      ncStatus(16) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_hyrbid_phase_functions", 1)
      ncStatus(17) = nf90_put_att(ncFileId, NF90_Global, "Hybrid_phase_function_width", &
                                                          hybridPhaseFunWidth)
    else
      ncStatus(16) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_hyrbid_phase_functions", 0)
      ncStatus(17) = nf90_put_att(ncFileId, NF90_Global, "Hybrid_phase_function_width", 0.)
    end if 

    if(useRussianRouletteForIntensity) then 
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_Russian_roulette", 1)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "Intensity_Russian_roulette_zeta_min",  zetaMin)
    else
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_Russian_roulette", 0)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "Intensity_Russian_roulette_zeta_min", 0.)
    end if 
    
    if(limitIntensityContributions) then 
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "limited_intensity_contributions", 1)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "max_intensity_contribution",  maxIntensityContribution)
    else
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "limited_intensity_contributions", 0)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "max_intensity_contribution", 0.)
    end if 
    ncStatus(20) = nf90_put_att(ncFileId, NF90_Global, "Cpu_time_total", CpuTimeTotal)
    ncStatus(21) = nf90_put_att(ncFileId, NF90_Global, "Cpu_time_setup", CpuTimeSetup)
    ncStatus(22) = nf90_put_att(ncFileId, NF90_Global, "Number_of_processors_used", numProcs)
    
    if(recScatOrd) then
    ncStatus(23) = nf90_put_att(ncFileId, NF90_GLOBAL, "Highest_recorded_scattering_order", numRecScatOrd)
    end if

    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Attributes ", ncStatus(:19)

    !
    ! Dimensions
    !
    ncStatus( 1) = nf90_def_dim(ncFileId, "x", numX, xDimId)
    ncStatus( 2) = nf90_def_dim(ncFileId, "y", numY, yDimId)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) & 
      ncStatus( 3) = nf90_def_dim(ncFileId, "z", numZ, zDimId)
    !
    ! Dimension variables
    !
    ncStatus( 4) = nf90_def_var(ncFileId, "x", NF90_DOUBLE, xDimId, ncVarId)
    ncStatus( 5) = nf90_def_var(ncFileId, "y", NF90_DOUBLE, yDimId, ncVarId)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) & 
      ncStatus( 6) = nf90_def_var(ncFileId, "z", NF90_DOUBLE, zDimId, ncVarId)
    !
    ! Flux variables
    !
    ncStatus( 7) = nf90_def_var(ncFileId, "fluxUp", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus( 8) = nf90_def_var(ncFileId, "fluxDown", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus( 9) = nf90_def_var(ncFileId, "fluxAbsorbed", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(10) = nf90_def_var(ncFileId, "fluxUp_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(11) = nf90_def_var(ncFileId, "fluxDown_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(12) = nf90_def_var(ncFileId, "fluxAbsorbed_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    !
    ! Absorption profile
    !
    if(reportAbsorptionProfile) then 
      ncStatus(13) = nf90_def_var(ncFileId, "absorptionProfile", &
                                  nf90_float, zDimId, ncVarId)
      ncStatus(14) = nf90_def_var(ncFileId, "absorptionProfile_StdErr", &
                                  nf90_float, zDimId, ncVarId)
    end if 
    
    !
    ! Volume absorption
    !
    if(reportVolumeAbsorption) then 
      ncStatus(15) = nf90_def_var(ncFileId, "absorbedVolume", &
                                  nf90_float, (/ xDimId, yDimId, zDimID /), ncVarId)
      ncStatus(16) = nf90_def_var(ncFileId, "absorbedVolume_StdErr", &
                                  nf90_float, (/ xDimId, yDimId, zDimId /), ncVarId)
    end if 
    
    !
    ! Intensity
    !
    if(computeIntensity) then
      numIntensityDirs = size(RadianceStats, 3)
      ncStatus(17) = nf90_def_dim(ncFileId, "direction",     numIntensityDirs,      dirDimId)
      ncStatus(19) = nf90_def_var(ncFileId, "intensityMus",  nf90_float, dirDimId, ncVarId)
      ncStatus(19) = nf90_def_var(ncFileId, "intensityPhis", nf90_float, dirDimId, ncVarId)
      ncStatus(20) = nf90_def_var(ncFileId, "intensity", &
                                  nf90_float, (/ xDimId, yDimId, dirDimId /), ncVarId)
      ncStatus(21) = nf90_def_var(ncFileId, "intensity_StdErr", &
                                  nf90_float, (/ xDimId, yDimId, dirDimId /), ncVarId)
    end if
 
    !recScatOrd
    if(recScatOrd) then 
      ncStatus(23) = nf90_def_dim(ncFileId, "numRecScatOrd", numRecScatOrd+1, scatDimId)
      ncstatus(24) = nf90_def_var(ncFileId, "Scattering_Order", nf90_float, scatDimId, ncVarId)
      ncStatus(25) = nf90_def_var(ncFileId, "fluxUpByScatOrd", nf90_float, &
                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
      ncStatus(26) = nf90_def_var(ncFileId, "fluxDownByScatOrd", nf90_float, &
                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
      ncStatus(27) = nf90_def_var(ncFileId, "fluxUpByScatOrd_StdErr", nf90_float, &
                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
      ncStatus(28) = nf90_def_var(ncFileId, "fluxDownByScatOrd_StdErr", nf90_float, &
                                   (/ xDimId, yDimId, scatDimId /), ncVarId)

      if(computeIntensity .and. present(intensityByScatOrdStats)) then
      ncStatus(29) = nf90_def_var(ncFileId, "intensityByScatOrd", nf90_float, &
                                   (/ xDimId, yDimId, dirDimId, scatDimId /), ncVarId) 
      ncStatus(30) = nf90_def_var(ncFileId, "intensityByScatOrd_StdErr", nf90_float, &
                                   (/ xDimId, yDimId, dirDimId, scatDimId /), ncVarId)
      end if
    end if

    ncStatus(31) = nf90_EndDef(ncFileId)
    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:31)
    ! ----------------------
    ! File is set up - now write each variable
    ! 
    ! Dimension variables
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "x", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, (xPosition(:numX) + xPosition(2:))/2)
    ncStatus( 3) = nf90_inq_varid(ncFileId, "y", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, (yPosition(:numY) + yPosition(2:))/2)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) then 
      ncStatus( 5) = nf90_inq_varid(ncFileId, "z", ncVarId)
      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, (zPosition(:numZ) + ZPosition(2:))/2)
    end if 
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Position", ncStatus(:6)
    ncstatus(:) = nf90_NoErr
    !
    ! Upward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxUp", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxUp_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux up", ncStatus(:4)

    !
    ! Downward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxDown", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxDown_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux down", ncStatus(:4)

    !
    ! Absorbed flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxAbsorbed", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxAbsorbed_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux absorbed", ncStatus(:4)

    !
    ! Absorption profile
    !
    if(reportAbsorptionProfile) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorptionProfile", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorptionProfile_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Absorption profile", ncStatus(:4)
    end if 
    
    !
    ! Volume absorption
    !
    if(reportVolumeAbsorption) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorbedVolume", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorbedVolume_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux convergence", ncStatus(:4)
    end if 
    
    !
    ! Intensity
    !
    if(computeIntensity) then
      ncStatus( 1) = nf90_inq_varid(ncFileId, "intensityMus", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, intensityMus(:numIntensityDirs))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "intensityPhis", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, intensityPhis(:numIntensityDirs))
      ncStatus( 5) = nf90_inq_varid(ncFileId, "intensity", ncVarId)
      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, RadianceStats(:, :, :,1))
      ncStatus( 7) = nf90_inq_varid(ncFileId, "intensity_StdErr", ncVarId)
      ncStatus( 8) = nf90_put_var(ncFileId, ncVarId, RadianceStats(:, :, :,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Intensity", ncStatus(:8)
    end if

    if(recScatOrd) then
      ncStatus( 1) = nf90_inq_varId(ncFileId, "Scattering_Order", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId,(/(REAL(N), N=0,numRecScatOrd)/) )
!       ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, scatOrderHolder)
      ncStatus( 3) = nf90_inq_varId(ncFileId, "fluxUpByScatOrd", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxUpByScatOrdStats(:,:,:,1))
      ncStatus( 5) = nf90_inq_varId(ncFileId, "fluxUpByScatOrd_StdErr", ncVarId)
      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, fluxUpByScatOrdStats(:,:,:,2))
      ncStatus( 7) = nf90_inq_varId(ncFileId, "fluxDownByScatOrd", ncVarId)
      ncStatus( 8) = nf90_put_var(ncFileId, ncVarId, fluxDownByScatOrdStats(:,:,:,1))
      ncStatus( 9) = nf90_inq_varId(ncFileId, "fluxDownByScatOrd_StdErr", ncVarId)
      ncStatus(10) = nf90_put_var(ncFileId, ncVarId, fluxDownByScatOrdStats(:,:,:,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print*, "Flux*ByScatOrd", ncStatus(:10)
      
     if(computeIntensity .and. present(intensityByScatOrdStats)) then
       ncStatus( 1) = nf90_inq_varId(ncFileId, "intensityByScatOrd", ncVarId)
       ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, intensityByScatOrdStats(:,:,:,:,1))
       ncStatus( 3) = nf90_inq_varId(ncFileId, "intensityByScatOrd_StdErr", ncVarId)
       ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, intensityByScatOrdStats(:,:,:,:,2))
       if(any(ncStatus(:) /= nf90_NoErr)) print*, "IntensityByScatOrd", ncStatus(:4)
     end if
    end if

    ncStatus( 1) = nf90_close(ncFileId)
    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:1)

  end subroutine writeResults_netcdf
! -------------------------------------------------------------------------------
end program monteCarloDriver
