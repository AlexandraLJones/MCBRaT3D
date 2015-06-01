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
  use netcdf
  use numericUtilities
  implicit none
 
  real, parameter :: Pi = dacos(-1.0_8)

  ! Input parameters
  !   Radiative transfer 
  real(8)                 :: surfaceAlbedo = 0.
  real                 :: solarMu = 1., solarAzimuth = 0.,  LW_flag = -1. 
  real(8)                 :: surfaceTemp = 300.0       ! added by Alexandra Jones Fall 2011. lambda in microns and surface temp in K
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
  character(len=256)   :: domainFileName = "", solarSourceFile = "", &
			  instrResponseFile="", SSPfilename="", physDomainFile=""
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

  namelist /radiativeTransfer/ solarMu, solarAzimuth, surfaceTemp, &
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
  
  namelist /fileNames/         solarSourceFile, instrResponseFile, &
			       SSPfilename, physDomainFile, &
                               outputRadFile, outputFluxFile, &
                               outputAbsProfFile, outputAbsVolumeFile, outputNetcdfFile
  

   ! Local variables
  character(len=256)   :: namelistFileName!, photon_file, col_file
!  character(len=256)   :: voxel_file, voxel_file2, horiz_file, level_file, col_file, row_file, diff_file, photon_file
  character(len=256)   :: batch_file, checkpointFile
  integer              :: nX, nY, nZ, phtn, checkFreq = 100
  integer              :: i, j, k, batch, ix, iy, iz
  integer              :: numRadDir
  integer              :: numberOfComponents
  logical              :: computeIntensity, ompParallel
  real                 :: cpuTime0, cpuTime1, cpuTime2, cpuTimeTotal, cpuTimeSetup
  real                 :: prevTotal, total, writeFreq=3600.0
  real                 :: meanFluxUp, meanFluxDown, meanFluxAbsorbed
  real(8)                 :: solarFlux, emittedFlux, lambda, lambdaAbove, lambdaBelow, dLambda, corr, tempSum, corrContr
  real(8)                 :: meanFluxUpStats(2), meanFluxDownStats(2), meanFluxAbsorbedStats(2)
  real(8), allocatable    :: xPosition(:), yPosition(:), zPosition(:), tempExt(:,:,:)
  real(8), allocatable    :: solarSourceFunction(:), centralLambdas(:), fluxCDF(:)
  real, allocatable    :: fluxUp(:, :), fluxDown(:, :), fluxAbsorbed(:, :), absorbedProfile(:), absorbedVolume(:, :, :), Radiance(:, :, :)
  real(8), allocatable    :: fluxUpStats(:, :, :), fluxDownStats(:, :, :), fluxAbsorbedStats(:, :, :), &
				absorbedProfileStats(:, :), absorbedVolumeStats(:, :, :, :), RadianceStats(:, :, :, :),&
				dummy2D(:,:), dummy3D(:,:,:), dummy4D1(:,:,:,:), dummy4D2(:,:,:,:)
!  real, allocatable    :: meanFluxUpByScatOrd(:), meanFluxDownByScatOrd(:)
!  real, allocatable    :: meanFluxUpByScatOrdStats(:,:), meanFluxDownByScatOrdStats(:,:)
!  real, allocatable    :: fluxUpByScatOrd(:,:,:), fluxDownByScatOrd(:,:,:)
!  real, allocatable    :: fluxUpByScatOrdStats(:,:,:,:), fluxDownByScatOrdStats(:,:,:,:)
!  real, allocatable    :: intensityByScatOrd(:,:,:,:), intensityByScatOrdStats(:,:,:,:,:)
  integer              ::  N, atms_photons, maxThreads, availProcs, numPhotonsProcessed, n, ierr
  integer              :: eindex, nonempty, counts, totalPhotonsProcessed, batchsum, tempIndex
  integer, allocatable :: indexes(:), tempFreqs(:), left(:)
  integer              :: batchesAssigned, batchesCompleted
  integer(8)           ::  totalNumPhotons = 0, counter
  integer              :: currentFreq(1)
  integer              :: SSPfileID, err
!  real(8), allocatable    :: voxel_weights(:,:,:,:), col_weights(:,:,:),level_weights(:,:)
!  integer, allocatable   :: voxel_tallys1(:,:,:), voxel_tallys2(:,:,:), voxel_tallys1_sum(:,:,:), voxel_tallys2_sum(:,:,:), voxel_tallys1_total(:,:,:), voxel_tallys2_total(:,:,:)
  integer, allocatable   ::  freqDistr(:), photonCDF(:)
  integer, dimension(MPI_STATUS_SIZE)  :: mpiStatus
  integer, parameter     :: EXIT_TAG = 99, NEW_FREQ = 1, SAME_FREQ = 2, FIRST_FREQ=0, PHOTONS = 4

   ! I3RC Monte Carlo code derived type variables
  type(commonDomain)         :: commonPhysical!, commonConcen
  type(domain)               :: thisDomain, tempDomain
  type(ErrorMessage)         :: status
  type(randomNumberSequence), allocatable, dimension(:) :: randoms
!  type(photonStream)         :: incomingPhotons
!  type(photonStream), allocatable, dimension(:) :: incomingBBPhotons
  type(photonStream)         :: incomingPhotons
  type(integrator)           :: mcIntegrator
  type(Weights)              :: theseWeights

  ! Variables related to splitting up the job across processors
  integer            :: thisProc, thisThread            ! ID OF CURRENT PROCESSOR; default 0
  integer            :: numProcs, numThreads            ! TOTAL NUMBER OF PROCESSORS USED; default 1
  integer            :: batchesPerProcessor, lambdaPerProc


  ! -----------------------------------------
  ! Start communications among multiple processes.
  !
!PRINT*, 'about to call initializeProcesses' 
  call initializeProcesses(numProcs, thisProc)

 !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thisThread)
    
 availProcs = OMP_GET_NUM_PROCS()
! maxThreads = OMP_GET_THREAD_LIMIT()
 ompParallel =  OMP_IN_PARALLEL()
!PRINT *, 'in parallel region? ', ompParallel, 'available procs:', availProcs

  numThreads = OMP_GET_NUM_THREADS()
!PRINT *, 'numThreads=', numThreads
  !!$OMP SINGLE
  allocate(randoms(0:numThreads-1))
   
!PRINT*, 'returned from initializeProcesses'

! temporary output file for the purposes of debugging
!  write(voxel_file, '(A,I0.4)') "voxel.out.",thisProc
!  write(photon_file, '(A,I0.4)') "photon.out.",thisProc
!  write(col_file, '(A,I0.4)') "col.out.",thisProc
!  write(level_file, '(A,I0.4)') "level.out.",thisProc
!  write(horiz_file, '(A,I0.4)') "horiz.out.",thisProc
!  write(diff_file, '(A,I0.4)') "diff.out.",thisProc
!  write(voxel_file2, '(A,I0.4)') "voxel2.out.",thisProc
!  write(batch_file, '(A,I0.4)') "batch.out.", thisProc

!   open(unit=34, file=trim(photon_file) , status='UNKNOWN')
!  open(unit=33, file=trim(voxel_file) , status='UNKNOWN')
!  open(unit=32, file=trim(level_file) , status='UNKNOWN')
!  open(unit=31, file=trim(col_file) , status='UNKNOWN')
!  open(unit=14, file=trim(row_file) , status='UNKNOWN')
!  open(unit=15, file=trim(horiz_file) , status='UNKNOWN')
!  open(unit=16, file=trim(diff_file) , status='UNKNOWN')
!  open(unit=17, file=trim(voxel_file2) , status='UNKNOWN')
!  open(unit=51, file=trim(batch_file), status='UNKNOWN')

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
  !!$OMP END SINGLE

 !!$OMP DO PRIVATE(i, j, k) 
    do i = 1,nMu
      do j = 1,nPhi
         k = (i-1)*nPhi+j
         intensityMus(k) = mus(i)
         intensityPhis(k) = phis(j)
      end do
    end do
   end if
  end if
  !!$OMP SINGLE
 if (allocated(mus)) deallocate(mus)
 if (allocated(phis)) deallocate(phis) 

  numRadDir = count(abs(intensityMus(:)) > 0.) 
  computeIntensity = numRadDir > 0 .and. &
                     (len_trim(outputRadFile) > 0 .or. len_trim(outputNetcdfFile) > 0)
  if(.not. computeIntensity) outputRadFile = ""
  
  !set recScatOrd to false if the number requested is negative
  !no need to do the opposite because the integrator 
  !will listen only to numRecScatOrd if it is present
!  if(numRecScatOrd < 0) recScatOrd=.false.

  ! -----------------------------------------
  !  Read the domain file
  !
!  allocate(BBDomain(numLambda))
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'namelist numLambda= ', numLambda
  call read_Common(physDomainFile, commonPhysical, status)
  call printStatus(status)
  allocate(freqDistr(1:numLambda), photonCDF(1:numLambda))
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Driver: Read common domain'

  if (LW_flag .ge. 0.0)then
     allocate(fluxCDF(1:numLambda))
     fluxCDF=0
     lambdaPerProc = CEILING(DBLE(numLambda/numProcs))
     err = nf90_open(trim(SSPfileName), nf90_NoWrite, SSPFileID)
     if(err /= nf90_NoErr) then
	PRINT *, "Driver: error opening file ", err, trim(nf90_strerror(err))
	STOP
     end if
     DO i = thisProc*lambdaPerProc+1, MIN(numLambda, thisProc*lambdaPerProc+lambdaPerProc), 1  
PRINT *, "in LW loop at iteration ", i
        call read_SSPTable(SSPfileID, i, commonPhysical, thisDomain, status) ! domain is initialized within this routine
        call printStatus(status)
	call getInfo_Domain(thisDomain, lambda=lambda, status=status)
        call printStatus(status)

        if (i .eq. thisProc*lambdaPerProc+1)then
           call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, status=status)
           call printStatus(status)

	   call read_SSPTable(SSPfileID, i+1, commonPhysical, tempDomain, status) ! domain is initialized within this routine   
	   call getInfo_Domain(tempDomain, lambda=lambdaAbove, status=status)
           call printStatus(status)
	   if (i .gt. 1) then ! we have a more accurate way to calculate dLambda
	      call read_SSPTable(SSPfileID, i-1, commonPhysical, tempDomain, status) ! domain is initialized within this routine
              call getInfo_Domain(tempDomain, lambda=lambdaBelow, status=status)
              call printStatus(status)
	      dLambda=ABS((lambdaAbove-lambdaBelow)/2.0_8)
	   else
	       dLambda = ABS(lambdaAbove-lambda) ! best we can do for i = 1
	   end if
	   lambdaBelow=lambda
	   lambda=lambdaAbove
	elseif(i .gt. thisProc*lambdaPerProc+1 .and. i .lt. numLambda)then
	   call finalize_Domain(tempDomain)
	   call read_SSPTable(SSPfileID, i+1, commonPhysical, tempDomain, status) ! domain is initialized within this routine
           call getInfo_Domain(tempDomain, lambda=lambdaAbove, status=status)
           call printStatus(status)
	   dLambda=ABS((lambdaAbove-lambdaBelow)/2.0_8)
	   lambdaBelow=lambda
           lambda=lambdaAbove
	else !i must equal numLambda
	   dLambda=ABS(lambda-lambdaBelow)
	   call finalize_Domain(tempDomain)
        end if

        theseWeights = new_Weights(numX=nX, numY=nY, numZ=nZ, numlambda=1, status=status)
        call printStatus(status)
        call emission_weighting(thisDomain, numLambda, i, theseWeights, surfaceTemp, instrResponseFile, dLambda=dLambda, totalFlux=emittedFlux, status=status)
        call printStatus(status)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'returned from emission_weighting'

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

!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'emittedFlux=', emittedFlux

 !!$OMP DO SCHEDULE(static, 1) PRIVATE(thisThread)
!     DO i = 1, numLambda
!PRINT *, 'numThreads=', numThreads, 'thisThread=', thisThread
!        incomingBBPhotons(i) = new_PhotonStream (theseWeights=theseWeights, iLambda=i, numberOfPhotons=freqDistr(i),randomNumbers=randoms(thisThread), status=status)
!       voxel_tallys1_sum = voxel_tallys1_sum + voxel_tallys1
!     END DO
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'initialized thermal BB photon stream'
!call printStatus(status)
!    open(unit=12, file=trim(photon_file) , status='UNKNOWN')
!    DO k = 1, nZ
 !     DO j = 1, nY
  !      write(12,"(100I12)") voxel_tallys1_sum(:,j,k)
 !     end do
!    end do
!    close(12)
	fluxCDF(i)=emittedFlux
	call finalize_Weights(theseWeights)
     END DO
     err=nf90_close(SSPfileId)
!     fluxCDF(:) = sumAcrossProcesses(fluxCDF)
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,fluxCDF,numLambda,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!     DO i = numLambda, 1, -1
!	fluxCDF(i)=SUM(fluxCDF(1:i))
!     END DO

     corr = 0.0_8
     DO i = 2,numLambda
	corrContr = fluxCDF(i)-corr
	tempSum = fluxCDF(i-1)+corrContr
	corr = (tempSum-fluxCDF(i-1))-corrContr
	fluxCDF(i)=tempSum
     END DO

     solarFlux = fluxCDF(numLambda)
     fluxCDF=fluxCDF/solarFlux
     fluxCDF(numLambda)=1.0_8

     thisThread = OMP_GET_THREAD_NUM()
!PRINT *, 'thisThread=', thisThread
  !randoms(thisThread) = new_RandomNumberSequence(seed = (/ iseed, thisProc, thisThread /) )
     randoms(thisThread) = new_RandomNumberSequence(seed = (/ iseed,thisProc, thisThread /) ) ! I think for the set up step every processor should have an identical set of random numbers so that when work is being divided up the photonstreams and photon frequency distributions are consistent with one another.
     if (MasterProc)then
	call getFrequencyDistr(numLambda, fluxCDF, &
	     MOD(numPhotonsPerBatch*numBatches, numProcs-1),&
	     randoms(thisThread), freqDistr)
     else
	call getFrequencyDistr(numLambda, fluxCDF, numPhotonsPerBatch*numBatches/(numProcs-1),randoms(thisThread), freqDistr)
     end if

     deallocate(fluxCDF)
     call synchronizeProcesses
     freqDistr(:) = sumAcrossProcesses(freqDistr)
  else
if(MasterProc .and. thisThread .eq. 0)PRINT *, "in SW part about to read ", SSPfilename
     err = nf90_open(trim(SSPfileName), nf90_NoWrite, SSPFileID)
     if(err /= nf90_NoErr) then
        PRINT *, "Driver: error opening file ", err, trim(nf90_strerror(err))
        STOP
     end if
     call read_SSPTable(SSPfileID, 1, commonPhysical, thisDomain, status) ! domain is initialized within this routine
     call printStatus(status)
     err=nf90_close(SSPfileId)
     call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, status=status)
     call printStatus(status)
     theseWeights=new_Weights(numLambda=numLambda, status=status)
     call printStatus(status)
     allocate(solarSourceFunction(1:numLambda))
     allocate(centralLambdas(1:numLambda))
     call read_SolarSource(solarSourceFile, numLambda, solarSourceFunction, centralLambdas, status=status)
     call printStatus(status) 
     call solar_Weighting(theseWeights, numLambda, solarSourceFunction, centralLambdas, solarMu, instrResponseFile, emittedFlux, status=status)   ! convert the solar source function to CDF and total Flux
!     if(MasterProc .and. thisThread .eq. 0)PRINT *, 'filled solar weights'
     call printStatus(status)
     solarFlux = emittedFlux
     deallocate(centralLambdas)
     deallocate(solarSourceFunction)

     thisThread = OMP_GET_THREAD_NUM()
!PRINT *, 'thisThread=', thisThread
  !randoms(thisThread) = new_RandomNumberSequence(seed = (/ iseed, thisProc, thisThread /) )
     randoms(thisThread) = new_RandomNumberSequence(seed = (/ iseed, thisProc, thisThread /) ) ! I think for the set up step every processor should have an identical set of random numbers so that when work is being divided up the photonstreams and photon frequency distributions are consistent with one another.
     if (MasterProc)then
        call getFrequencyDistr(theseWeights, &
             MOD(numPhotonsPerBatch*numBatches, numProcs-1),&
             randoms(thisThread), freqDistr)
     else
        call getFrequencyDistr(theseWeights, numPhotonsPerBatch*numBatches/(numProcs-1),randoms(thisThread), freqDistr)
     end if
     call finalize_Weights(theseWeights)
     call synchronizeProcesses
     freqDistr(:) = sumAcrossProcesses(freqDistr)
  end if

  if(MasterProc .and. thisThread .eq. 0)PRINT *, 'solarFlux=', solarFlux
!  if(MasterProc .and. thisThread .eq. 0)PRINT*, 'frequency distribution:', freqDistr
  if(MasterProc .and. thisThread .eq. 0)PRINT*, 'number of non-zero frequencies:', COUNT(freqDistr .gt. 0)
  
  photonCDF(1) = freqDistr(1)
!if(MasterProc .and. thisThread .eq. 0) write(34,"(2E30.20, 2I12)") centralLambdas(1), solarSourceFunction(1), photonCDF(1), freqDistr(1)
  DO i = 2, numLambda
     photonCDF(i)= photonCDF(i-1)+freqDistr(i)
!if(MasterProc .and. thisThread .eq. 0) write(34,"(2E30.20, 2I12)") centralLambdas(i), solarSourceFunction(i), photonCDF(i), freqDistr(i) 
  END DO
!if(MasterProc .and. thisThread .eq. 0) close(34)

  allocate(xPosition(nx+1), yPosition(ny+1), zPosition(nz+1))
  call getInfo_Domain(thisDomain, xPosition = xPosition, yPosition = yPosition, &
                      zPosition = zPosition, numberOfComponents=numberOfComponents, status = status)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'retrieved Domain info for ', lambda, lambdaI, numLambda
  call printStatus(status)

  ! Set up the integrator object. It only grabs information from the domain that is consistent between all of them so I only need to pass any one domain
  mcIntegrator = new_Integrator(thisDomain, status = status)     
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'initialized Integrator'
  call finalize_Domain(thisDomain)
  call printStatus(status)					 ! CAN I MOVE THIS SECTION INITIALIZING THE INTEGRATOR TO BEFORE THE LOOP RADING IN ALL THE DOMAINS BUT AFTER READING IN THE FIRST DOMAIN SO THAT I CAN JUST PARALLELIZE THE REST? 
!  call finalize_Domain(thisDomain)                               !

   ! Set the surface albedo, table sizes, and maybe the radiance directions
  call specifyParameters (mcIntegrator,                          &
                          minInverseTableSize = nPhaseIntervals, &
                          LW_flag = LW_flag,                     &
                          status = status)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Specified photon source'
  call printStatus(status) 

  if (computeIntensity) then
    call specifyParameters (mcIntegrator, &
                            minForwardTableSize=nPhaseIntervals, &
                            intensityMus=intensityMus(1:numRadDir), &
                            intensityPhis=intensityPhis(1:numRadDir), &
                            computeIntensity=computeIntensity,numComps=numberOfComponents, status=status)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Specified Intensity Calcs'
    call printStatus(status) 
  endif

!  if(recScatOrd) then
!    call specifyParameters (mcIntegrator, &
!                            recScatOrd=recScatOrd, &
!                            numRecScatOrd=numRecScatOrd, & 
!                            status = status)
!PRINT *, 'SPecified recScatOrder'
!    call printStatus(status)
!  end if 

  !
  ! Make the algorithmic choices
  !
  call specifyParameters(mcIntegrator,                              &
                         useRayTracing      = useRayTracing,        &
                         useRussianRoulette = useRussianRoulette,   &
                         status = status)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Specfied ray tracing'
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
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Specified variance reduction'
    call printStatus(status) 
  end if

!PRINT *, 'Driver: Specified Parameters'
   ! Allocate and zero the arrays for radiative quantities and moments 
!  allocate (voxel_tallys1(nX, nY, nZ), voxel_tallys1_sum(nX, nY, nZ), voxel_tallys1_total(nX, nY, nZ))
!  allocate (voxel_tallys2(nX, nY, nZ), voxel_tallys2_sum(nX, nY, nZ), voxel_tallys2_total(nX, nY, nZ))
  allocate (fluxUp      (nX, nY), fluxUpStats      (nX, nY, 2))
  allocate (fluxDown    (nX, nY), fluxDownStats    (nX, nY, 2))
  allocate (fluxAbsorbed(nX, nY), fluxAbsorbedStats(nX, nY, 2))
  allocate (absorbedProfile(nZ), absorbedProfilestats(nZ, 2))
  allocate (absorbedVolume(nX, nY, nZ), absorbedVolumeStats(nX, nY, nZ, 2))
!  voxel_tallys1(:,:,:)=0 ; voxel_tallys1_sum(:,:,:) = 0 ; voxel_tallys1_total(:,:,:) = 0
!  voxel_tallys2(:,:,:)=0 ; voxel_tallys2_sum(:,:,:) = 0 ; voxel_tallys2_total(:,:,:) = 0
  meanFluxUpStats(:) = 0.0  ; meanFluxDownStats(:) = 0.0  ; meanFluxAbsorbedStats(:) = 0.0
  fluxUpStats(:, :, :) = 0.0  ; fluxDownStats(:, :, :) = 0.0  ; fluxAbsorbedStats(:, :, :) = 0.0
  absorbedProfilestats(:, :) = 0.0 ;  absorbedVolumeStats(:, :, :, :) = 0.0
  if (computeIntensity) then
    allocate (Radiance(nX, nY, numRadDir), RadianceStats(nX, nY, numRadDir, 2))
    RadianceStats(:, :, :, :) = 0.0
  endif  
  !!$OMP END SINGLE
  
!  if(recScatOrd .and. numRecScatOrd >=0) then 
!    allocate(meanFluxUpByScatOrd(0:numRecScatOrd), meanFluxDownByScatOrd(0:numRecScatOrd),&
!             fluxUpByScatOrd(nX,nY,0:numRecScatOrd),fluxDownByScatOrd(nX,nY,0:numRecScatOrd), &
!             meanFluxUpByScatOrdStats(0:numRecScatOrd, 2), meanFluxDownByScatOrdStats(0:numRecScatOrd,2), &
!             fluxUpByScatOrdStats(nX,nY, 0:numRecScatOrd, 2), fluxDownByScatOrdStats(nX,nY,0:numRecScatOrd,2))
!    meanFluxUpByScatOrdStats=0.0;meanFluxDownByScatOrdStats=0.0;
!    fluxUpByScatOrdStats=0.0;fluxDownByScatOrdStats=0.0;
!    
!    if(computeIntensity) then
!      allocate(intensityByScatOrd(nX,nY,numRadDir,0:numRecScatOrd), &
!               intensityByScatOrdStats(nX,nY,numRadDir, 0:numRecScatOrd, 2))
!      intensityByScatOrdStats=0.0;
!    end if
!  end if

  call cpu_time(cpuTime1)
!PRINT *, "called cpu_time"
  call synchronizeProcesses
!PRINT *, "called synchronizeProcesses"
  cpuTimeSetup = sumAcrossProcesses(cpuTime1 - cpuTime0) 
!if(MasterProc .and. thisThread .eq. 0)PRINT *, "cpuTimeSetup=", cpuTimeSetup
!STOP
  if (MasterProc) &
    print *, "Setup CPU time (secs, approx): ", int(cpuTimeSetup)
!PRINT *, "setup completed"

  ! --------------------------------------------------------------------------

  ! The  loop over batches is for estimating the uncertainty in the flux and
  !   radiance from the variance between numBatches independent calculations. 
!  numBatches = max(numBatches,2)
!  batchesPerProcessor = numBatches/numProcs
!  ! If the number of batches doesn't divide among the processors evenly increase the 
!  !   number until it does. 
!  if(mod(numBatches, numProcs) /= 0) then 
!    batchesPerProcessor = batchesPerProcessor + 1
!    numBatches = batchesPerProcessor * numProcs
!  end if 
 
!  if (MasterProc) &
!    print *, "Doing ", batchesPerProcessor, " batches on each of ", numProcs, " processors." 
!PRINT *, "entering batch loop"
!  batches: do batch = thisProc*batchesPerProcessor + 1, thisProc*batchesPerProcessor + batchesPerProcessor
  if(MasterProc)then  ! must be master process
    currentFreq = 1								! This is a new variable. Should only range from 1 to nlambda
!    n = 1
!    allocate(startingPhoton(1:numlambda))               ! This is a new variable. Needs to be initialized to 1
!    startingPhoton = 0
!    WHERE(freqDistr .ne. 0)startingPhoton = 1 
!PRINT *, startingPhoton
    batchesAssigned = 0
    batchesCompleted = 0
    ! hand out initial work to processes and update photons and currentFreq. !!numProcs must be .le. numLambda!!
    DO n = 1, numProcs-1
!PRINT *, "beginning of loop", currentFreq(1), freqDistr(currentFreq(1))
	DO WHILE(freqDistr(currentFreq(1)) .lt. 1) ! this do while loop let's us skip over empty frequencies
	  currentFreq = currentFreq + 1
	  if (currentFreq(1) .gt. numLambda) EXIT
	END DO	
!PRINT *, "after skipping", currentFreq(1), freqDistr(currentFreq(1))
	if (freqDistr(currentFreq(1)).ge. numPhotonsPerBatch)then
	  CALL MPI_SEND(numPhotonsPerBatch, 1, MPI_INTEGER, n, FIRST_FREQ, MPI_COMM_WORLD, ierr) ! first send starting index of photons to work on(this always should be the first type of message received along with tag indicating what type of messages will follow)
	  CALL MPI_SEND(currentFreq, 1, MPI_INTEGER, n, FIRST_FREQ, MPI_COMM_WORLD, ierr) ! then send frequency to work on
	  freqDistr(currentFreq(1))=freqDistr(currentFreq(1))-numPhotonsPerBatch 
 	else
	  ! update the variables for the next proc
!	  freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch
!	  startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!	  batchesAssigned = batchesAssigned + 1
!	  PRINT *, 'MasterProc=', thisProc, 'sending photons of frequency index ', currentFreq, 'which has', freqDistr(currentFreq), 'photons remaining, to rank ', n
!	else
!          EXIT
	  eindex = findCDFIndex(numPhotonsPerBatch, photonCDF)
	  nonempty = COUNT(freqDistr(currentFreq(1):eindex).gt.0) ! how many indeses in that range actually have photons?
	  if(nonempty .eq. 0)CYCLE 
	  allocate(indexes(1:nonempty),tempFreqs(1:nonempty))
	  batchSum = 0
	  tempIndex = 0
!PRINT *, "rank: currentFreq, ending index, nonempty indexes: ", n, currentFreq(1), eindex, nonempty
	  DO i = currentFreq(1), eindex
	     if(freqDistr(i).gt.0)then
		tempIndex = tempIndex + 1
		batchsum = batchsum + freqDistr(i)
		indexes(tempIndex) = i
		tempFreqs(tempIndex) = freqDistr(i)
		freqDistr(i)=freqDistr(i)-tempFreqs(tempIndex)
	     end if  
	  END DO
	  ! correction for overestimate in last bin
	  tempFreqs(nonempty)=tempFreqs(nonempty)-(batchsum-numPhotonsPerBatch)
	  freqDistr(indexes(nonempty))=freqDistr(indexes(nonempty))+(batchsum-numPhotonsPerBatch)

	  CALL MPI_SEND(tempFreqs, nonempty, MPI_INTEGER, n, FIRST_FREQ, MPI_COMM_WORLD, ierr)
	  CALL MPI_SEND(indexes, nonempty, MPI_INTEGER, n, FIRST_FREQ, MPI_COMM_WORLD, ierr)
	  deallocate(indexes,tempFreqs)
	  currentFreq=eindex
        end if
	batchesAssigned = batchesAssigned + 1
	photonCDF=photonCDF-numPhotonsPerBatch
!	n = n + 1
!	currentFreq = currentFreq + 1			! so that next proc will be assigned a dif freq
    END DO

!    DO WHILE (n .le. numProcs-1) ! if there are any procs left that don't have tasks
!	currentFreq = MAXLOC(freqDistr(:))
!	CALL MPI_SEND(freqDistr(currentFreq), 1, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr) ! first send starting index of photons to work on(this always should be the first type of message received along with tag indicating what type of messages will follow)
!        CALL MPI_SEND(currentFreq, 1, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr) ! then send frequency to work on
        ! update the variables for the next proc
!        freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch
!        startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!        batchesAssigned = batchesAssigned + 1
!        PRINT *, 'MasterProc=', thisProc, 'sending photons of frequency index ', currentFreq, 'which has', freqDistr(currentFreq), 'photons remaining, to rank ', n
!	n = n + 1
!    END DO

		! after the initial round of batches to processors there probably be work left to do. Allow workers to come to you for  more tasks. Listen for a query message
    counter = 1
    prevTotal = cpuTime0 
    DO WHILE(batchesCompleted .lt. batchesAssigned)
	CALL MPI_RECV(numPhotonsProcessed, 1, MPI_INTEGER, MPI_ANY_SOURCE, PHOTONS, MPI_COMM_WORLD, mpiStatus, ierr)
	n =  mpiStatus(MPI_SOURCE)
	totalNumPhotons = totalNumPhotons + numPhotonsProcessed
	batchesCompleted = batchesCompleted + 1
	PRINT *, 'have completed', batchesCompleted, 'of the', batchesAssigned, 'batches assigned'
!	PRINT *, 'MasterProc recieved completion message ', mpiStatus(MPI_TAG), 'from rank', mpiStatus(MPI_SOURCE), 'asking if there are more photons in frequency index ', currentFreq
!PRINT *, MOD(batchesCompleted, checkFreq)
!	if(MOD(batchesCompleted, checkFreq) .eq. 0)then
!	     call cpu_time(total) 
!PRINT *, "time to check write Frequency", total-prevTotal
!	   if(total-prevTotal .gt. writeFreq)then
!		prevTotal = total
!          ! initiate intermediate sum across processes and write
!!PRINT *, LEN(outputNetcdfFile), LEN(TRIM(outputNetcdfFile))
!PRINT *, "time to write checkpoint file"
!		if(MOD(counter, 2) .eq. 0)then
!               	   write(checkpointFile, '(A,I1)') TRIM(outputNetcdfFile), 2
!             	else
!             	   write(checkpointFile, '(A,I1)') TRIM(outputNetcdfFile), 1
!             	end if
!          	counter = counter + 1
!PRINT *, "summing across processes"          
!          	! Pixel-by-pixel fluxes
!          	fluxUpStats(:, :, :)       = sumAcrossProcesses(fluxUpStats)
!          	fluxDownStats(:, :, :)     = sumAcrossProcesses(fluxDownStats)
!          	fluxAbsorbedStats(:, :, :) = sumAcrossProcesses(fluxAbsorbedStats)
!
!          	! Absorption (mean profile and cell-by-cell)
!          	absorbedProfileStats(:, :)      = sumAcrossProcesses(absorbedProfileStats)
!          	absorbedVolumeStats(:, :, :, :) = sumAcrossProcesses(absorbedVolumeStats)
!
!          	! Radiance
!          	if (computeIntensity)then
!                    RadianceStats(:, :, :, :) = sumAcrossProcesses(RadianceStats)
!PRINT *, "writing checkpoint file"
!               	    call writeResults_netcdf(domainFileName,  totalNumPhotons, batchesCompleted, &
!                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
!                                 xPosition, yPosition, zPosition,                 &
!                                 checkpointFile,                                &
!                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
!                                 absorbedProfileStats, absorbedVolumeStats,       &
!                                 intensityMus, intensityPhis, RadianceStats)
!		    RadianceStats(:, :, :, :) = 0.0
!          	else
!PRINT *, "writing checkpoint file"
!                    call writeResults_netcdf(domainFileName,  totalNumPhotons, batchesCompleted, &
!                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
!                                 xPosition, yPosition, zPosition,                 &
!                                 checkpointFile,                                &
!                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
!                                 absorbedProfileStats, absorbedVolumeStats)
!          	end if
!         	print *, "Wrote to checkpoint file"
!		fluxUpStats(:, :, :) = 0.0  ; fluxDownStats(:, :, :) = 0.0  ; fluxAbsorbedStats(:, :, :) = 0.0
!		absorbedProfilestats(:, :) = 0.0 ;  absorbedVolumeStats(:, :, :, :) = 0.0
!     	   end if
!	end if



	if (photonCDF(numLambda).gt.0)then ! if there is still work left to do...
	   if (freqDistr(currentFreq(1)).ge. numPhotonsPerBatch)then
        	CALL MPI_SEND(numPhotonsPerBatch, 1, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr) ! first send starting index of photons to work on(this always should be the first type of message received along with tag indicating what type of messages will follow)
          	CALL MPI_SEND(currentFreq, 1, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr) ! then send frequency to work on
		freqDistr(currentFreq(1))=freqDistr(currentFreq(1))-numPhotonsPerBatch
           else
          ! update the variables for the next proc
!         freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch
!         startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!         batchesAssigned = batchesAssigned + 1
!         PRINT *, 'MasterProc=', thisProc, 'sending photons of frequency index ', currentFreq, 'which has', freqDistr(currentFreq), 'photons remaining, to rank ', n
!       else
!          EXIT
          	eindex = findCDFIndex(numPhotonsPerBatch, photonCDF)
          	nonempty = COUNT(freqDistr(currentFreq(1):eindex).gt.0) ! how many indeses in that range actually have photons?
		if(nonempty .eq. 0)CYCLE
          	allocate(indexes(1:nonempty),tempFreqs(1:nonempty))
		batchsum = 0
		tempIndex = 0
!PRINT *, "rank, currentFreq, ending index, nonempty indexes, photons: ", n, currentFreq(1), eindex, nonempty, freqDistr(currentFreq(1):eindex) 
          	DO i = currentFreq(1), eindex
!PRINT *, "inside do loop: rank, current index, photons: ", n, i, freqDistr(i)
             	   if(freqDistr(i).gt.0)then
			tempIndex = tempIndex +1
!PRINT *, "inside if statement: rank, current index, relative index, photons: ", n, i,tempIndex, freqDistr(i)
                   	indexes(tempIndex) = i
                   	tempFreqs(tempIndex) = freqDistr(i)
			batchsum = batchsum + freqDistr(i)
			freqDistr(i)=freqDistr(i)-tempFreqs(tempIndex)
             	   end if
          	END DO
		! correction for overestimate in last bin
          	tempFreqs(nonempty)=tempFreqs(nonempty)-(batchsum-numPhotonsPerBatch)
          	freqDistr(indexes(nonempty))=freqDistr(indexes(nonempty))+(batchsum-numPhotonsPerBatch)
          	CALL MPI_SEND(tempFreqs, nonempty, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr)
          	CALL MPI_SEND(indexes, nonempty, MPI_INTEGER, n, NEW_FREQ, MPI_COMM_WORLD, ierr)
          	currentFreq=eindex
		deallocate(indexes, tempFreqs)
           end if
           batchesAssigned = batchesAssigned + 1
           photonCDF=photonCDF-numPhotonsPerBatch   
!	     if(freqDistr(currentFreq(1)) .gt. 0)then ! use flag SAME_FREQ. ! are more photons available at the frequency it's currently working on?
!PRINT *, freqDistr(currentFreq(1)), 'Photons left  on frequency index', currentFreq, 'Telling rank ', mpiStatus(MPI_SOURCE), 'to keep going starting with ', startingPhoton(currentFreq)
!		CALL MPI_SEND(freqDistr(currentFreq), 1, MPI_INTEGER, mpiStatus(MPI_SOURCE), SAME_FREQ, MPI_COMM_WORLD, ierr) 
!		 freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch ! update book keeping variables
!		startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!		batchesAssigned = batchesAssigned + 1
!		PRINT *, 'frequency index ',currentFreq, 'has ', freqDistr(currentFreq), 'remaining photons'
!	     elseif (ANY(startingPhoton(:) .eq. 1))THEN ! There are unassigned frequencies left. use flag NEW_FREQ
!		currentFreq = MINLOC(startingPhoton(:), MASK=startingPhoton .gt. 0) ! Returns first location of unassigned frequency with actual photons to run
!PRINT *, currentFreq, "frequency was previously unassigned with", freqDistr(currentFreq), "photons left. Starting with photon ", startingPhoton(currentFreq) 
!		CALL MPI_SEND(freqDistr(currentFreq), 1, MPI_INTEGER, mpiStatus(MPI_SOURCE), NEW_FREQ, MPI_COMM_WORLD, ierr) ! first send tag so it knows what kind of data to expect
!		CALL MPI_SEND(currentFreq, 1, MPI_INTEGER, mpiStatus(MPI_SOURCE), NEW_FREQ, MPI_COMM_WORLD, ierr)   ! send the messages containing the needed info
!		freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch ! update book keeping variables
!		startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!		batchesAssigned = batchesAssigned + 1
!	    else ! assign it to the frequency with the most work left to do. use flag NEW_FREQ
!PRINT *, "minimum starting photon greater than 0 is ", MINVAL(startingPhoton,  MASK=startingPhoton .gt. 0), "and MINVAL(startingPhoton(:)) .eq. 1 evaluates as ", MINVAL(startingPhoton(:)) .eq. 1
!		currentFreq = MAXLOC(freqDistr(:))
!PRINT *, 'assigned rank ', mpiStatus(MPI_SOURCE), 'to work on frequency with most remaining photons', currentFreq, freqDistr(currentFreq), 'starting with photon', startingPhoton(currentFreq)
!		CALL MPI_SEND(freqDistr(currentFreq), 1, MPI_INTEGER, mpiStatus(MPI_SOURCE), NEW_FREQ, MPI_COMM_WORLD, ierr) ! first send tag so it knows what kind of data to expect
!		CALL MPI_SEND(currentFreq, 1, MPI_INTEGER, mpiStatus(MPI_SOURCE), NEW_FREQ, MPI_COMM_WORLD, ierr) ! send the messages containing the needed info
!		freqDistr(currentFreq) = freqDistr(currentFreq) - numPhotonsPerBatch ! update the book keeping variables
!		startingPhoton(currentFreq) = startingPhoton(currentFreq) + numPhotonsPerBatch
!		batchesAssigned = batchesAssigned + 1
!	    END IF
!	    PRINT *, 'assigned rank ', mpiStatus(MPI_SOURCE), 'to work on frequency', currentFreq, 'starting with photon', startingPhoton(currentFreq)-numPhotonsPerBatch, 'which now has', freqDistr(currentFreq), 'photons remaining'
	END IF
    END DO

    DO n = 1,numProcs-1	! send exit message to procs				
	! pass to worker the needed info
	CALL MPI_SEND(freqDistr(1), 1, MPI_INTEGER, n, EXIT_TAG, MPI_COMM_WORLD, ierr)
    END DO
PRINT *, "Done with batches"
!
!
!
!
!
!
!
!
  else ! must be a worker process
    err = nf90_open(trim(SSPfileName), nf90_NoWrite, SSPFileID)
    if(err /= nf90_NoErr) then
       PRINT *, "Driver:initial batches: error opening file ", err, trim(nf90_strerror(err))
       STOP
    end if
     DO i=0,numThreads-1
       randoms(thisThread)=new_RandomNumberSequence(seed = (/ iseed, thisProc, thisThread /) )
    END DO
    DO WHILE (mpiStatus(MPI_TAG) .ne. EXIT_TAG)
	CALL MPI_PROBE(0,MPI_ANY_TAG, MPI_COMM_WORLD, mpiStatus, ierr)
	CALL MPI_GET_COUNT(mpiStatus,MPI_INTEGER,counts,ierr)
	if(mpiStatus(MPI_TAG) .eq. EXIT_TAG)EXIT
       allocate(left(1:counts), indexes(1:counts))
       CALL MPI_RECV(left, counts, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpiStatus, ierr) ! recieve initial work assignment
       CALL MPI_RECV(indexes, counts, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpiStatus, ierr)
!    PRINT *, 'Rank', thisProc, 'Received initial message ', mpiStatus(MPI_TAG), 'frequency index ', currentFreq
       totalPhotonsProcessed = 0
       DO i= 1, counts
	  if(mpiStatus(MPI_TAG) .eq. NEW_FREQ)then
	     CALL finalize_Domain(thisDomain)
	    if(LW_flag >= 0.0)CALL finalize_Weights(theseWeights)
	  end if
	  CALL read_SSPTable(SSPFileID, indexes(i), commonPhysical, thisDomain, status)
          if(LW_flag >= 0.0)then   ! need to reconstruct a domain and weighting array
             theseWeights=new_Weights(numX=nX, numY=nY, numZ=nZ, numLambda=1, status=status)
	     CALL printStatus(status)
	     CALL emission_weighting(thisDomain, numLambda, indexes(i), theseWeights, surfaceTemp, instrResponseFile, status=status)
	     CALL printStatus(status)
	  end if

	  if(LW_flag >= 0.0)then 
	     incomingPhotons=new_PhotonStream(theseWeights=theseWeights,iLambda=1, numberOfPhotons=MIN(numPhotonsPerBatch, left(i)), randomNumbers=randoms(thisThread), status=status) ! monochromatic thermal source call
!PRINT *, 'For its first work, Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons. Tag is: ', mpiStatus(MPI_TAG)
	     CALL printStatus(status)
          else
	     incomingPhotons=new_PhotonStream(solarMu, solarAzimuth, numberOfPhotons=MIN(numPhotonsPerBatch, left(i)), randomNumbers=randoms(thisThread), status=status) ! monochromatic solar source call
!	PRINT *, 'For its first work, Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons. Tag is: ', mpiStatus(MPI_TAG)
             CALL printStatus(status)  
          end if


!    if(LW_flag >= 0.0)then   ! need to reconstruct a domain and weighting array
!        theseWeights=new_Weights(numX=nX, numY=nY, numZ=nZ, numLambda=1, status=status)
!        CALL printStatus(status)
!        CALL emission_weighting(thisDomain, numLambda, currentFreq(1), theseWeights, surfaceTemp, instrResponseFile, status=status)
!        CALL printStatus(status)
!	incomingPhotons=new_PhotonStream(theseWeights=theseWeights,iLambda=1, numberOfPhotons=MIN(numPhotonsPerBatch, left), randomNumbers=randoms(thisThread), status=status) ! monochromatic thermal source call
!PRINT *, 'For its first work, Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons. Tag is: ', mpiStatus(MPI_TAG)
!        CALL printStatus(status)
!    else
!	incomingPhotons=new_PhotonStream(solarMu, solarAzimuth, numberOfPhotons=MIN(numPhotonsPerBatch, left), randomNumbers=randoms(thisThread), status=status) ! monochromatic solar source call
!	PRINT *, 'For its first work, Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons. Tag is: ', mpiStatus(MPI_TAG)
!        CALL printStatus(status)  
!    end if
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thisThread)
         prevTotal = cpuTime0
!    DO WHILE (mpiStatus(MPI_TAG) .ne. EXIT_TAG)
!	PRINT *, 'Rank ', thisProc, 'recieved this tag:', mpiStatus(MPI_TAG), 'and will process photons of frequency', currentFreq, 'which has', left, 'photons left'
	! set the current photon of the appropriate frequency to the starting value
!	call setCurrentPhoton(incomingBBPhotons(currentFreq(1)), start, status)
!	call printStatus(status)
	! Do the work on the sent data
!	DO i=0,numThreads-1
!	    ! seed the random number generator
!	    randoms(i) = new_RandomNumberSequence(seed = (/ iseed, thisProc, thisThread /) )
!	END DO
	! Now we compute the radiative transfer for this batch of photons.
        call computeRadiativeTransfer (mcIntegrator, thisDomain, randoms(thisThread), incomingPhotons, numPhotonsPerBatch, numPhotonsProcessed, status)
        call printStatus(status)
        totalPhotonsProcessed=totalPhotonsProcessed+numPhotonsProcessed
	! get contribution arrays from integrator
	!   This particular integrator provides fluxes at the top and bottom
     	!   of the domain for both the domain mean and pixel level fluxes,
     	!   the absorbed flux profile, 3D field of absorbed flux, and
     	!   the pixel level radiances at top and/or bottom of domain.
        call reportResults (mcIntegrator, &
	 meanFluxUp=meanFluxUp, meanFluxDown=meanFluxDown, meanFluxAbsorbed=meanFluxAbsorbed,   &
	 fluxUp=fluxUp(:, :), fluxDown=fluxDown(:, :), fluxAbsorbed=fluxAbsorbed(:, :), &
	 absorbedProfile=absorbedProfile(:), volumeAbsorption=absorbedVolume(:, :, :), status = status)
        call printStatus(status)
	!
	! Accumulate the first and second moments of each quantity over the batches
	!
!	totalNumPhotons = totalNumPhotons + numPhotonsProcessed

	meanFluxUpStats(1) = meanFluxUpStats(1) + meanFluxUp*numPhotonsProcessed
	meanFluxUpStats(2) = meanFluxUpStats(2) + numPhotonsProcessed*(meanFluxUp**2.0_8)

	meanFluxDownStats(1) = meanFluxDownStats(1) + meanFluxDown*numPhotonsProcessed
        meanFluxDownStats(2) = meanFluxDownStats(2) + numPhotonsProcessed*(meanFluxDown**2.0_8)

	meanFluxAbsorbedStats(1) = meanFluxAbsorbedStats(1) + meanFluxAbsorbed*numPhotonsProcessed
        meanFluxAbsorbedStats(2) = meanFluxAbsorbedStats(2) + numPhotonsProcessed*(meanFluxAbsorbed**2.0_8)

	FluxUpStats(:,:,1) = FluxUpStats(:,:,1) + FluxUp*numPhotonsProcessed
        FluxUpStats(:,:,2) = FluxUpStats(:,:,2) + numPhotonsProcessed*(FluxUp**2.0_8)

        FluxDownStats(:,:,1) = FluxDownStats(:,:,1) + FluxDown*numPhotonsProcessed
        FluxDownStats(:,:,2) = FluxDownStats(:,:,2) + numPhotonsProcessed*(FluxDown**2.0_8)

        FluxAbsorbedStats(:,:,1) = FluxAbsorbedStats(:,:,1) + FluxAbsorbed*numPhotonsProcessed
        FluxAbsorbedStats(:,:,2) = FluxAbsorbedStats(:,:,2) + numPhotonsProcessed*(FluxAbsorbed**2.0_8)

	AbsorbedProfileStats(:,1) = AbsorbedProfileStats(:,1) + AbsorbedProfile*numPhotonsProcessed
        AbsorbedProfileStats(:,2) = AbsorbedProfileStats(:,2) + numPhotonsProcessed*(AbsorbedProfile**2.0_8)

	AbsorbedVolumeStats(:,:,:,1) = AbsorbedVolumeStats(:,:,:,1) + AbsorbedVolume*numPhotonsProcessed
        AbsorbedVolumeStats(:,:,:,2) = AbsorbedVolumeStats(:,:,:,2) + numPhotonsProcessed*(AbsorbedVolume**2.0_8)

	if (computeIntensity) then
	    call reportResults(mcIntegrator, intensity = Radiance(:, :, :), status = status)
	    RadianceStats(:, :, :,1) = RadianceStats(:, :, :,1) + Radiance(:, :, :)*numPhotonsProcessed
	    RadianceStats(:, :, :,2) = RadianceStats(:, :, :,2) + numPhotonsProcessed*(Radiance(:, :, :)**2.0_8)
!PRINT *, numPhotonsProcessed, 'numPhotons', solarFlux*FluxUp, '=FluxUp', solarFlux*FluxDown, '=FluxDown', solarFlux*AbsorbedVolume, '=AbsVolume', solarFlux*Radiance(:, :, :), '=Intensity'
	end if
     END DO
!WRITE(51, '(E26.16, I10, I4, 4E26.16 )') solarFlux, numPhotonsProcessed, currentFreq, Radiance(1,1,1), RadianceStats(1,1,1,1:2)

     CALL MPI_SEND(totalPhotonsProcessed, 1, MPI_INTEGER, 0, PHOTONS, MPI_COMM_WORLD, ierr)
     if(allocated(left))deallocate(left)
     if(allocated(indexes))deallocate(indexes)
!	CALL cpu_time(total)
!	if(total-prevTotal .gt. writeFreq)then
!PRINT *, 'Rank', thisProc, 'entered the sumAcrossprocesses block', total-prevTotal
!	   allocate(dummy4D1(nx,ny,nz,2), dummy2D(nz,2), dummy3D(nx,ny,2))
!	   ! Pixel-by-pixel fluxes
!           dummy3D(:, :, :)       = sumAcrossProcesses(fluxUpStats)
!           dummy3D(:, :, :)     = sumAcrossProcesses(fluxDownStats)
!           dummy3D(:, :, :) = sumAcrossProcesses(fluxAbsorbedStats)
!
!          ! Absorption (mean profile and cell-by-cell)
!           dummy2D(:, :)      = sumAcrossProcesses(absorbedProfileStats)
!           dummy4D1(:, :, :, :) = sumAcrossProcesses(absorbedVolumeStats)
!	  deallocate(dummy4D1, dummy2D, dummy3D)
!          ! Radiance
!          if (computeIntensity) then
!	     allocate(dummy4D2(nx,ny,numRadDir, 2))
!             dummy4D2(:, :, :, :) = sumAcrossProcesses(RadianceStats) 
!	     deallocate(dummy4D2)
!	  end if
!	  prevTotal = total
!	end if
  END DO
  err = nf90_close(SSPFileID)
  if(err /= nf90_NoErr) then
     PRINT *, "Driver:initial batches: error closing file ", err, trim(nf90_strerror(err))
     STOP
  end if
!	CALL MPI_RECV(left, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpiStatus, ierr) ! recieve initial work assignment
!	if(mpiStatus(MPI_TAG) .eq. EXIT_TAG)EXIT
!	if(mpiStatus(MPI_TAG) .eq. NEW_FREQ) then
!	    CALL MPI_RECV(currentFreq, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, mpiStatus, ierr)
!	    CALL finalize_Domain(thisDomain)
!	    err = nf90_open(trim(SSPfileName), nf90_NoWrite, SSPFileID)
 !    	    if(err /= nf90_NoErr) then
  !      	PRINT *, "Driver:batches: error opening file ", err, trim(nf90_strerror(err))
   !     	STOP
    ! 	    end if
!	    CALL read_SSPTable(SSPfileID, currentFreq(1), commonPhysical, thisDomain, status) !!!! TODO !!!!
!	    err = nf90_close(SSPfileID)
!	    if(LW_flag >= 0.0)then   ! need to reconstruct a domain and weighting array
!		CALL finalize_Weights(theseWeights)
!		theseWeights=new_Weights(numX=nX, numY=nY, numZ=nZ, numLambda=1, status=status)
!	        CALL printStatus(status)
 !       	CALL emission_weighting(thisDomain, numLambda, currentFreq(1), theseWeights, surfaceTemp, instrResponseFile, status=status)
!	        CALL printStatus(status)
!	    end if
!	end if

!	if(LW_flag >= 0.0)then
!	   incomingPhotons=new_PhotonStream(theseWeights=theseWeights, iLambda=1, numberOfPhotons=MIN(numPhotonsPerBatch, left), randomNumbers=randoms(thisThread), status=status) ! monochromatic thermal source call
!	   PRINT *, 'Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons and recieved tag: ', mpiStatus(MPI_TAG)
!	   CALL printStatus(status)
!	else
!	   incomingPhotons=new_PhotonStream(solarMu, solarAzimuth, numberOfPhotons=MIN(numPhotonsPerBatch, left), randomNumbers=randoms(thisThread), status=status) ! monochromatic solar source call
!	   PRINT *, 'Rank ', thisProc, 'initialized ', MIN(numPhotonsPerBatch, left), 'photons and recieved tag: ', mpiStatus(MPI_TAG)
!	   CALL printStatus(status)
!	end if    
!	PRINT *, 'rank ', thisProc, 'recieved message ', mpiStatus(MPI_TAG)
 !   END DO
!PRINT *, 'thisProc=', thisProc,  'RadianceStatsSums=', RadianceStats(1,1,1,1:3)
  end if !  end do batches
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Driver: finished tracing photons'
  if(allocated(freqDistr))deallocate(freqDistr)
!  if (allocated(startingPhoton)) deallocate(startingPhoton)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Driver: about to finalize photon stream'
!close(51)
  call finalize_PhotonStream (incomingPhotons)
  !
  ! Accumulate statistics from across all the processors
  !
! voxel_tallys1_total = sumAcrossProcesses(voxel_tallys1_sum)
! voxel_tallys2_total = sumAcrossProcesses(voxel_tallys2_sum)

!  totalNumPhotons = sumAcrossProcesses(totalNumPhotons)

  ! Domain-mean fluxes
!if(MasterProc .and. thisThread .eq. 0)PRINT *, "about to do final sumAcrossProcesses" 
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

if(MasterProc .and. thisThread .eq. 0)then
!  PRINT *, 'Driver: accumulated results.' 
!  PRINT *, 'RadianceStats=', RadianceStats(1,1,1,1:3)
end if
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
  meanFluxUpStats(:)       = solarFlux*meanFluxUpStats(:)/totalNumPhotons
  meanFluxUpStats(2)       = solarFlux*meanFluxUpStats(2)
  meanFluxUpStats(2)       = sqrt(max(0.0, meanFluxUpStats(2)-(meanFluxUpStats(1)**2.0_8))/(batchesCompleted-1))

  meanFluxDownStats(:)     = solarFlux*meanFluxDownStats(:)/totalNumPhotons
  meanFluxDownStats(2)     = solarFlux*meanFluxDownStats(2)
  meanFluxDownStats(2)     = sqrt(max(0.0, meanFluxDownStats(2)-(meanFluxDownStats(1)**2.0_8))/(batchesCompleted-1))

  meanFluxAbsorbedStats(:) = solarFlux*meanFluxAbsorbedStats(:)/totalNumPhotons
  meanFluxAbsorbedStats(2) = solarFlux*meanFluxAbsorbedStats(2)
  meanFluxAbsorbedStats(2) = sqrt(max(0.0, meanFluxAbsorbedStats(2)-(meanFluxAbsorbedStats(1)**2.0_8))/(batchesCompleted-1))

  fluxUpStats(:, :, :)       =solarFlux*fluxUpStats(:, :, :)/totalNumPhotons
  fluxUpStats(:, :, 2)       = solarFlux*fluxUpStats(:, :, 2)
  fluxUpStats(:, :, 2)       = sqrt(max(0.0,fluxUpStats(:, :, 2)-(fluxUpStats(:, :,1)**2.0_8))/(batchesCompleted-1))

  fluxDownStats(:, :, :)     =solarFlux*fluxDownStats(:, :, :)/totalNumPhotons
  fluxDownStats(:, :, 2)     = solarFlux*fluxDownStats(:, :, 2)
  fluxDownStats(:, :, 2)     = sqrt(max(0.0, fluxDownStats(:, :, 2)-(fluxDownStats(:, :,1)**2.0_8))/(batchesCompleted-1))

  fluxAbsorbedStats(:, :, :) =solarFlux*fluxAbsorbedStats(:, :, :)/totalNumPhotons
  fluxAbsorbedStats(:, :, 2) = solarFlux*fluxAbsorbedStats(:, :, 2)
  fluxAbsorbedStats(:, :, 2) = sqrt(max(0.0, fluxAbsorbedStats(:, :, 2)-(fluxAbsorbedStats(:, :,1)**2.0_8))/(batchesCompleted-1))

  absorbedProfileStats(:, :) =solarFlux*absorbedProfileStats(:, :)/totalNumPhotons
  absorbedProfileStats(:, 2) = solarFlux*absorbedProfileStats(:, 2)
  absorbedProfileStats(:, 2) = sqrt(max(0.0, absorbedProfileStats(:, 2)-(absorbedProfileStats(:,1)**2.0_8))/(batchesCompleted-1))

  absorbedVolumeStats(:, :, :, :) = solarFlux*absorbedVolumeStats(:, :, :, :)/totalNumPhotons
  absorbedVolumeStats(:, :, :, 2) = solarFlux*absorbedVolumeStats(:, :, :, 2)
  absorbedVolumeStats(:, :, :, 2) = sqrt(max(0.0, absorbedVolumeStats(:, :, :, 2)-(absorbedVolumeStats(:, :, :,1)**2.0_8))/(batchesCompleted-1))

  if (computeIntensity) then
    RadianceStats(:, :, :, :) = solarFlux*RadianceStats(:, :, :, :)/totalNumPhotons
    RadianceStats(:, :, :,2) = solarFlux*RadianceStats(:, :, :, 2)
!PRINT *, "mean Radiance stats including solarflux =", sum (RadianceStats(:, :, 1,1))/real (nX * nY)
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'RadianceStats=', RadianceStats(:, :, :,1:3)
    RadianceStats(:, :, :,2) = sqrt(max(0.0, RadianceStats(:, :, :,2)-(RadianceStats(:, :, :,1)**2.0_8))/(batchesCompleted-1))
!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'RadianceStats=', RadianceStats(:, :, :,1:3)
  endif

!if(MasterProc .and. thisThread .eq. 0)PRINT *, 'Driver: calculated radiative quantities. Mean and error of FluxUp'
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
      call writeResults_ASCII(domainFileName,  totalNumPhotons, batchesCompleted,      &
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
      if(computeIntensity) then 
        call writeResults_netcdf(domainFileName,  totalNumPhotons, batchesCompleted, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats)      
!      elseif(computeIntensity .and. recScatOrd) then
!        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
!                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
!                                 xPosition, yPosition, zPosition,                 &
!                                 outputNetcdfFile,                                &
!                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
!                                 absorbedProfileStats, absorbedVolumeStats,       &
!                                 intensityMus, intensityPhis, RadianceStats,      &
!                                 numRecScatOrd,                                   &
!                                 fluxUpByScatOrdStats, fluxDownByScatOrdStats,    &
!                                 intensityByScatOrdStats                          )
!      elseif(.not. computeIntensity .and. recScatOrd) then
!        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
!                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
!                                 xPosition, yPosition, zPosition,                 &
!                                 outputNetcdfFile,                                &
!                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
!                                 absorbedProfileStats, absorbedVolumeStats,       &
!                                 numRecScatOrd=numRecScatOrd,                     &
!                                 fluxUpByScatOrdStats=fluxUpByScatOrdStats,       &
!                                 fluxDownByScatOrdStats=fluxDownByScatOrdStats      )
      else
        call writeResults_netcdf(domainFileName,  totalNumPhotons, batchesCompleted, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats)
  
      end if 
      print *, "Wrote netcdf results"                         
    end if
  end if
  deallocate(xPosition,yPosition,zPosition)
  !
  ! Release all the memory. We should be able to finalize the derived types before we write 
  !   the output but this fails intermittently, so we want to be sure to get our results 
  !   before we take a chance on blowing up. 
  ! 
  deallocate (fluxUp, fluxUpStats, fluxDown, fluxDownStats, fluxAbsorbed, fluxAbsorbedStats)
  deallocate (absorbedProfile, absorbedProfileStats, absorbedVolume, absorbedVolumeStats)
  if (computeIntensity) deallocate (Radiance, RadianceStats)

  call finalize_RandomNumberSequence(randoms(thisThread))
  call finalize_Integrator (mcIntegrator)

contains

! -------------------------------------------------------------------------------
  subroutine writeResults_ASCII(domainFileName,  totalNumPhotons, numBatches,&
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
    integer,            intent(in) :: numBatches
    integer(8),          intent(in) :: totalNumPhotons
    logical,            intent(in) :: useRayTracing, useRussianRoulette, useHybridPhaseFunsForIntenCalcs
    real,               intent(in) :: hybridPhaseFunWidth
    real,               intent(in) ::  solarMu, solarAzimuth
    real(8),    intent(in) :: solarFlux,surfaceAlbedo
    real(8), dimension(:), intent(in) :: xPosition, yPosition, zPosition
    ! Flux variables
    character(len = *),       intent(in) :: outputFluxFile
    real(8), dimension(:),       intent(in) :: meanFluxUpStats, meanFluxDownStats, meanFluxAbsorbedStats
    real(8), dimension(:, :, :), intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    character(len = *),       intent(in) :: outputAbsProfFile
    real(8), dimension(:, :),    intent(in) :: absorbedProfileStats
    character(len = *),       intent(in) :: outputAbsVolumeFile
    real(8), dimension(:, :, :, :), &
                              intent(in) :: absorbedVolumeStats
    ! Intensity variable
    character(len = *),       intent(in) :: outputRadFile
    real, dimension(:),       intent(in) :: intensityMus, intensityPhis
    real(8), dimension(:, :, :, :), &
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
      write (2,'(A,I10)')  '!  Num_Photons=', totalNumPhotons
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
      write (2,'(A,I10)')  '!  Num_Photons=', totalNumPhotons
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
      write (2,'(A,I10)')  '!  Num_Photons=', totalNumPhotons
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
      write (2,'(A,I10)')  '!  Num_Photons=', totalNumPhotons
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
  subroutine writeResults_netcdf(domainFileName,  totalNumPhotons, numBatches, &
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
    integer,            intent(in) :: numBatches
    integer(8),          intent(in) :: totalNumPhotons
    real,               intent(in) ::  solarMu, solarAzimuth
    real(8),               intent(in) :: solarFlux, surfaceAlbedo
    ! The position variables mark the edges of the cells - we'll report the 
    !   results at the midpoints. 
    real(8), dimension(:), intent(in) :: xPosition, yPosition, zPosition

    ! Flux variables
    real(8), dimension(:, :, :),     intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    real(8), dimension(:, :),        intent(in) :: absorbedProfileStats
    real(8), dimension(:, :, :, :),  intent(in) :: absorbedVolumeStats
    ! Intensity variable
    real, dimension(:), optional, intent(in) :: intensityMus, intensityPhis
    real(8), dimension(:, :, :, :), &
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
    ncStatus( 8) = nf90_put_att(ncFileId, NF90_Global, "Total_number_of_photons", totalNumPhotons)
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
    
!    if(recScatOrd) then
!    ncStatus(23) = nf90_put_att(ncFileId, NF90_GLOBAL, "Highest_recorded_scattering_order", numRecScatOrd)
!    end if

!    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Attributes ", ncStatus(:19)

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
!    if(recScatOrd) then 
!      ncStatus(23) = nf90_def_dim(ncFileId, "numRecScatOrd", numRecScatOrd+1, scatDimId)
!      ncstatus(24) = nf90_def_var(ncFileId, "Scattering_Order", nf90_float, scatDimId, ncVarId)
!      ncStatus(25) = nf90_def_var(ncFileId, "fluxUpByScatOrd", nf90_float, &
!                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
!      ncStatus(26) = nf90_def_var(ncFileId, "fluxDownByScatOrd", nf90_float, &
!                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
!      ncStatus(27) = nf90_def_var(ncFileId, "fluxUpByScatOrd_StdErr", nf90_float, &
!                                   (/ xDimId, yDimId, scatDimId /), ncVarId)
!      ncStatus(28) = nf90_def_var(ncFileId, "fluxDownByScatOrd_StdErr", nf90_float, &
!                                   (/ xDimId, yDimId, scatDimId /), ncVarId)

!      if(computeIntensity .and. present(intensityByScatOrdStats)) then
!      ncStatus(29) = nf90_def_var(ncFileId, "intensityByScatOrd", nf90_float, &
!                                   (/ xDimId, yDimId, dirDimId, scatDimId /), ncVarId) 
!      ncStatus(30) = nf90_def_var(ncFileId, "intensityByScatOrd_StdErr", nf90_float, &
!                                   (/ xDimId, yDimId, dirDimId, scatDimId /), ncVarId)
!      end if
!    end if

    ncStatus(31) = nf90_EndDef(ncFileId)
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:31)
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
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Position", ncStatus(:6)
    ncstatus(:) = nf90_NoErr
    !
    ! Upward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxUp", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxUp_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,2))
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux up", ncStatus(:4)

    !
    ! Downward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxDown", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxDown_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,2))
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux down", ncStatus(:4)

    !
    ! Absorbed flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxAbsorbed", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxAbsorbed_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,2))
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux absorbed", ncStatus(:4)

    !
    ! Absorption profile
    !
    if(reportAbsorptionProfile) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorptionProfile", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorptionProfile_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,2))
!      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Absorption profile", ncStatus(:4)
    end if 
    
    !
    ! Volume absorption
    !
    if(reportVolumeAbsorption) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorbedVolume", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorbedVolume_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,2))
!      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux convergence", ncStatus(:4)
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
!      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Intensity", ncStatus(:8)
    end if

!    if(recScatOrd) then
!      ncStatus( 1) = nf90_inq_varId(ncFileId, "Scattering_Order", ncVarId)
!      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId,(/(REAL(N), N=0,numRecScatOrd)/) )
!!       ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, scatOrderHolder)
!      ncStatus( 3) = nf90_inq_varId(ncFileId, "fluxUpByScatOrd", ncVarId)
!      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxUpByScatOrdStats(:,:,:,1))
!      ncStatus( 5) = nf90_inq_varId(ncFileId, "fluxUpByScatOrd_StdErr", ncVarId)
!      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, fluxUpByScatOrdStats(:,:,:,2))
!      ncStatus( 7) = nf90_inq_varId(ncFileId, "fluxDownByScatOrd", ncVarId)
!      ncStatus( 8) = nf90_put_var(ncFileId, ncVarId, fluxDownByScatOrdStats(:,:,:,1))
!      ncStatus( 9) = nf90_inq_varId(ncFileId, "fluxDownByScatOrd_StdErr", ncVarId)
!      ncStatus(10) = nf90_put_var(ncFileId, ncVarId, fluxDownByScatOrdStats(:,:,:,2))
!      if(any(ncStatus(:) /= nf90_NoErr)) print*, "Flux*ByScatOrd", ncStatus(:10)
      
!     if(computeIntensity .and. present(intensityByScatOrdStats)) then
!       ncStatus( 1) = nf90_inq_varId(ncFileId, "intensityByScatOrd", ncVarId)
!       ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, intensityByScatOrdStats(:,:,:,:,1))
!       ncStatus( 3) = nf90_inq_varId(ncFileId, "intensityByScatOrd_StdErr", ncVarId)
!       ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, intensityByScatOrdStats(:,:,:,:,2))
!       if(any(ncStatus(:) /= nf90_NoErr)) print*, "IntensityByScatOrd", ncStatus(:4)
!     end if
!    end if

    ncStatus( 1) = nf90_close(ncFileId)
!    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:1)

  end subroutine writeResults_netcdf
! -------------------------------------------------------------------------------
end program monteCarloDriver
