! performs the steps neccessary to produce a PDF and CDF of photons frequency
program BBsetup
! module usage statements
 use UserInterface
 use RandomNumbers
 use MultipleProcesses
 use ErrorMessages
 use opticalProperties
 use emissionAndBBWeights
 use netcdf
 use numericUtilities
 use scatteringPhaseFunctions
 
implicit none

! input parameters from namelist
  !   Radiative transfer
  real(8)                 :: surfaceAlbedo = 0.
  real                 :: solarMu = 1., solarAzimuth = 0.,  LW_flag = -1.
  real(8)                 :: surfaceTemp = 300.0       ! added by Alexandra Jones Fall 2011. lambda in microns and surface temp in K
  integer, parameter   :: maxNumRad = 648 !18mus*36phis , oldVal:72
  real                 :: intensityMus(maxNumRad)  = 0., &
                          intensityPhis(maxNumRad) = 0.
  logical              :: angleFill = .false., calcRayl = .True.
  real, dimension(3)   :: thetaFill = -1., phiFill = -1.
  real, allocatable, dimension(:) :: mus, phis
  integer              :: nMu, nPhi, numLambda, lambdaIi, nx, ny, nz
  !   Monte Carlo
  integer(8)           :: numPhotonsPerBatch = 0
  integer              :: numBatches = 100, &
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
                          instrResponseFile="", physDomainFile=""
  character(len=256), dimension(4)  :: SSPfilename=""
  character(len=256)   :: outputFluxFile = "", outputRadFile = "",  &
                          outputAbsProfFile = "", outputAbsVolumeFile = "", &
                          outputNetcdfFile = "", distrFileName = ""

  !Aux hist-- auxiliary output
  !auxhist01-fluxes organized by the order of scattering, and intensities where
  !   appropriate
  logical              :: recScatOrd = .false.
  integer              :: numRecScatOrd = 0
  character(len=256)   :: auxhist01_radFile=""
  character(len=256)   :: auxhist01_fluxFile=""

  namelist /radiativeTransfer/ solarMu, solarAzimuth, surfaceTemp, intensityMus, intensityPhis, &
                               angleFill, thetaFill, phiFill, LW_flag, numLambda, calcRayl

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
                               SSPfilename, physDomainFile, distrFileName, &
                               outputRadFile, outputFluxFile, &
                               outputAbsProfFile, outputAbsVolumeFile, outputNetcdfFile


! Local variables
character(len=256)   :: namelistFileName
real(8), allocatable    :: solarSourceFunction(:), centralLambdas(:), fluxCDF(:)
integer                 :: n, i, err, ierr
integer, dimension(4):: SSPFileID=-9999   ! must have number of dimensions equal to SSPfilename
real(8)                 :: solarFlux, emittedFlux, lambda, lambdaAbove, lambdaBelow, dLambda, corr, tempSum, corrContr
integer(8), allocatable   ::  freqDistr(:), photonCDF(:)

!derived type variables
type(randomNumberSequence) :: randoms
type(Domain)               :: thisDomain,tempDomain
type(ErrorMessage)         :: status
type(Weights)              :: theseWeights
type(commonDomain), TARGET         :: commonPhysical

! Variables related to splitting up the job across processors
  integer            :: thisProc, numProcs, lambdaPerProc    
  real		     :: cpuTime0, cpuTime1, cpuTimeSetup


  ! -----------------------------------------
  ! Start communications among multiple processes.
  !
  call initializeProcesses(numProcs, thisProc)

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


! Read the physical domain file
  call read_Common(physDomainFile, commonPhysical, status)
  call printStatus(status)

  if (LW_flag .ge. 0.0)then
     allocate(fluxCDF(1:numLambda))
     fluxCDF=0
     lambdaPerProc = CEILING(DBLE(numLambda/numProcs))
     n = 1
     DO
        err = nf90_open(trim(SSPfilename(n)), nf90_NoWrite, SSPFileID(n))
        if(err /= nf90_NoErr) then
           PRINT *, "Driver: error opening file ", SSPfilename(n), trim(nf90_strerror(err))
           STOP
        end if
        n=n+1
        if(len_trim(SSPfilename(n)).le.0)EXIT
     END DO
     DO i = thisProc*lambdaPerProc+1, MIN(numLambda, thisProc*lambdaPerProc+lambdaPerProc), 1
	call read_SSPTable(SSPFileID, i, commonPhysical, thisDomain,.True.,.False., status) ! domain is initialized within this routine
        call getInfo_Domain(thisDomain, lambda=lambda, status=status)
        call printStatus(status)
        if (i .eq. thisProc*lambdaPerProc+1)then
           call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, status=status)
           call printStatus(status)

           call read_SSPTable(SSPFileID, i+1, commonPhysical, tempDomain,.True.,.False., status) ! domain is initialized within this routine
           call printStatus(status)
           if (i .gt. 1) then ! we have a more accurate way to calculate dLambda
              call finalize_Domain(tempDomain)
              call read_SSPTable(SSPFileID, i-1, commonPhysical, tempDomain, .True.,.False., status) ! domain is initialized within this routine
              call printStatus(status)
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
           call read_SSPTable(SSPFileID, i+1, commonPhysical, tempDomain, .True.,.False., status) ! domain is initialized within this routine
           call printStatus(status)
           call getInfo_Domain(tempDomain, lambda=lambdaAbove, status=status)
           call printStatus(status)
           dLambda=ABS((lambdaAbove-lambdaBelow)/2.0_8)
           lambdaBelow=lambda
           lambda=lambdaAbove
        else !i must equal numLambda
           dLambda=ABS(lambda-lambdaBelow)
        end if
        call finalize_Domain(tempDomain)
        theseWeights = new_Weights(numX=nX, numY=nY, numZ=nZ, numlambda=1, status=status)
        call printStatus(status)
        call emission_weighting(thisDomain, numLambda, i, theseWeights, surfaceTemp, &
                instrResponseFile, dLambda=dLambda, totalFlux=emittedFlux, status=status)
        call printStatus(status)
        if (i .lt. MIN(numLambda, thisProc*lambdaPerProc+lambdaPerProc)) call finalize_Domain(thisDomain)
        fluxCDF(i)=emittedFlux
        call finalize_Weights(theseWeights)
     END DO
     n = 1
     DO
        err = nf90_close(SSPFileID(n))
        if(err /= nf90_NoErr) then
           PRINT *, "Driver: error closing file ", SSPfilename(n), trim(nf90_strerror(err))
           STOP
        end if
        n=n+1
        if(len_trim(SSPfilename(n)).le.0)EXIT
     END DO
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,fluxCDF,numLambda,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

! do the cumulative sum
     corr = 0.0_8
     DO i = 2,numLambda
        corrContr = fluxCDF(i)-corr
        tempSum = fluxCDF(i-1)+corrContr
        corr = (tempSum-fluxCDF(i-1))-corrContr
        fluxCDF(i)=tempSum
     END DO
! normalize
     solarFlux = fluxCDF(numLambda)
     fluxCDF=fluxCDF/solarFlux
     fluxCDF(numLambda)=1.0_8
!get random number stream
     randoms = new_RandomNumberSequence(seed = (/ iseed,thisProc/))
! assign frequencies to photons
     allocate(freqDistr(1:numLambda), photonCDF(1:numLambda))
     if (MasterProc)then
        call getFrequencyDistr(numLambda, fluxCDF, &
             MOD(numPhotonsPerBatch*numBatches, numProcs-1),&
             randoms, freqDistr)
     else
        call getFrequencyDistr(numLambda, fluxCDF, numPhotonsPerBatch*numBatches/(numProcs-1),randoms, freqDistr)
     end if
     deallocate(fluxCDF)
     call synchronizeProcesses
     freqDistr(:) = sumAcrossProcesses(freqDistr)
  else ! SW simulations
     n = 1
     DO
        err = nf90_open(trim(SSPfilename(n)), nf90_NoWrite, SSPFileID(n))
        if(err /= nf90_NoErr) then
           PRINT *, "Driver: error opening file ", SSPfilename(n), trim(nf90_strerror(err))
           STOP
        end if
        n=n+1
        if(len_trim(SSPfilename(n)).le.0)EXIT
     END DO
     call read_SSPTable(SSPFileID, 1, commonPhysical, thisDomain, .True.,.False., status) ! domain is initialized within this routine
     call printStatus(status)
     n = 1
     DO
        err = nf90_close(SSPFileID(n))
        if(err /= nf90_NoErr) then
           PRINT *, "Driver: error closing file ", SSPfilename(n), trim(nf90_strerror(err))
           STOP
        end if
        n=n+1
        if(len_trim(SSPfilename(n)).le.0)EXIT
     END DO
     call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, status=status)
     call printStatus(status)
     theseWeights=new_Weights(numLambda=numLambda, status=status)
     call printStatus(status)
     allocate(solarSourceFunction(1:numLambda))
     allocate(centralLambdas(1:numLambda))
     call read_SolarSource(solarSourceFile, numLambda, solarSourceFunction, centralLambdas, status=status)
     call printStatus(status)
     call solar_Weighting(theseWeights, numLambda, solarSourceFunction, centralLambdas, solarMu, &
                instrResponseFile, emittedFlux, status=status)   ! convert the solar source function to CDF and total Flux
     call printStatus(status)
     solarFlux = emittedFlux
     deallocate(centralLambdas)
     deallocate(solarSourceFunction)
     allocate(freqDistr(1:numLambda), photonCDF(1:numLambda))
     randoms = new_RandomNumberSequence(seed = (/ iseed, thisProc/)) 
     if (MasterProc)then
        call getFrequencyDistr(theseWeights, &
             MOD(numPhotonsPerBatch*numBatches, numProcs-1),&
             randoms, freqDistr)
     else
        call getFrequencyDistr(theseWeights, numPhotonsPerBatch*numBatches/(numProcs-1),randoms, freqDistr)
     end if
     call finalize_Weights(theseWeights)
     call synchronizeProcesses
     freqDistr(:) = sumAcrossProcesses(freqDistr)
  end if


  if(MasterProc)PRINT *, 'solarFlux=', solarFlux
  if(MasterProc)PRINT*, 'number of non-zero frequencies:', COUNT(freqDistr .gt. 0)

! create photon CDF from freq distribution
  if(MasterProc)then
      photonCDF(1) = freqDistr(1)
      DO i = 2, numLambda
	photonCDF(i)= photonCDF(i-1)+freqDistr(i)
      END DO
  end if

  call cpu_time(cpuTime1)
  call synchronizeProcesses
  cpuTimeSetup = sumAcrossProcesses(cpuTime1 - cpuTime0)
  if (MasterProc) &
    print *, "Setup CPU time (secs, approx): ", int(cpuTimeSetup)

! write photon CDF and freq Distr to file along with RN seed
  call writePhotonDistributions(numLambda,solarFlux,freqDistr, iseed, photonCDF, distrFileName)

contains

  subroutine writePhotonDistributions(nsize, BBFlux, PDF, seed, CDF, outputFileName)
  use netcdf
  implicit none
  integer, intent(in)				:: nsize, seed
  real(8), intent(in)				:: BBFlux
  integer(8), dimension(1:nsize), intent(in)	:: PDF, CDF
  character(len=256), intent(in)			:: outputFileName
! Local variables
  integer                                       :: ncFileId, DimId, cdfVarId, pdfVarId
  integer, dimension(10) 			:: ncStatus = nf90_NoErr
   
  ncStatus( 1) = nf90_create(trim(outputFileName), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncFileId) 
ncStatus(2) = nf90_put_att(ncFileId, NF90_Global, "Random_number_seed", seed)
ncStatus(10) = nf90_put_att(ncFileId, NF90_Global, "totalFlux", BBFlux)
ncStatus( 3) = nf90_def_dim(ncFileId, "numSpectralPoints", nsize, DimId)
ncStatus( 4) = nf90_def_var(ncFileId, "photonFrequencyDistr", NF90_INT64, DimId, pdfVarId)
ncStatus( 5) = nf90_def_var(ncFileId, "photonCDF", NF90_INT64, DimId, cdfVarId)
ncStatus(6) = nf90_EndDef(ncFileId)
ncStatus( 7) = nf90_put_var(ncFileId, pdfVarId,PDF)
ncStatus( 8) = nf90_put_var(ncFileId, cdfVarId,CDF)
ncStatus( 9) = nf90_close(ncFileId)

  end subroutine writePhotonDistributions

end program BBsetup
