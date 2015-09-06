! new outline for parallelization of MieSSPTableCreate
Program MakeMieSSPTable

  use ErrorMessages
  use scatteringPhaseFunctions
  use UserInterface
  use netcdf
  use multipleProcesses

  IMPLICIT NONE
  !
  ! Input parameters
  !
  INTEGER :: NRETAB = 0
  REAL    :: WAVELEN1 = 0.0, WAVELEN2 = 0.0, DELTAWAVE = 0., PARDENS = 1.
  REAL    :: SRETAB = 0.0, ERETAB = 0.0, ALPHA = 0.0, MAXRADIUS = 0.0
!  REAL(8) :: albedo = 0.164
  COMPLEX :: RINDEX
  CHARACTER(LEN=1)   :: PARTYPE = "W", AVGFLAG = "C", DISTFLAG = "G"
  CHARACTER(LEN=256) :: phaseFunctionTableFile = "phaseFunctionTable.pft", frequencyFile
  NAMELIST /mie_table_input/ DISTFLAG, ALPHA, NRETAB, SRETAB, ERETAB, phaseFunctionTableFile, frequencyFile

  ! Local variable
  !
  INTEGER :: NSIZE, MAXLEG, I, J, L, NL, n, m, work, ierr, src, count
  INTEGER :: ncDimId, nfreq, ncVarId, ncFileId, entryDimId, fDimId, coeffDimId, thisProc, numProcs
  INTEGER, dimension(23) :: ncStatus
  INTEGER, dimension(8) :: VarId
  LOGICAL :: LOGSPACEDREFF
  REAL    :: PI, WAVELENCEN, XMAX, SCATTER
  CHARACTER(LEN=256) :: namelistFileName, thisPrefix="Component1_"
  character(len=256) :: tableDescription
  integer, parameter :: ind=1, ext=2, ssa=3, len=4, st=5, LG=6

  INTEGER, ALLOCATABLE :: NLEG1(:), NLEG(:), start(:,:), length(:,:), tempStart(:)
  REAL,    ALLOCATABLE :: RADII(:), ND(:)
  REAL(8),    ALLOCATABLE :: EXTINCT1(:), SCATTER1(:),EXTINCT(:), SSALB(:), lambdas(:)
   REAL(8),    ALLOCATABLE :: freqs(:), TOTALEXT(:,:), TOTALSSA(:,:)
  REAL,    ALLOCATABLE :: REFF(:),LEGEN1(:,:), LEGCOEF(:,:), TOTALLG(:,:), buffr(:)
  integer, dimension(MPI_STATUS_SIZE)  :: mpiStatus
   
  !
  ! Temperature at which to evaluate index of refraction for water and ice
  !
  real, parameter :: waterTemperature = 283., iceTemperature = 243.

  type(ErrorMessage)               :: status

! initiate procs
  call initializeProcesses(numProcs, thisProc)

! all procs read in freqs and calculate lambdas
  PI = ACOS(-1.0)
  RINDEX = CMPLX(0., 0.)


 ! Read namelist containing input filename, paramters describing effective radius, shape parameter, distribution type, and output file name
  namelistFileName = getOneArgument()
  OPEN (UNIT=1, FILE=trim(namelistFileName), STATUS='OLD')
  READ (1,NML=mie_table_input)
  CLOSE (1)
  if(thisProc .eq. 0)PRINT *,  frequencyFile
  !
  ! Check for valid input values and consistency among parameters
  !
  if(nReTab == 0) stop "MakeMieTable: Must specify at least one effective radius"
  if(alpha <= 0.) stop "MakeMieTable: must specify parameter alpha for size distribution"
  if(sretab <= 0) stop "MakeMieTable: must specify a starting effective radius (sretab)"
  if(eretab <= 0.) then
    eretab = sretab
    if(nretab > 1) print *, "MakeMieTable: specified a range of effective radii but requested only 1. " // &
                            "Single scattering parameters will be computed at only one value."
    nretab = 1
  end if
  LOGSPACEDREFF = NRETAB < 0
  NRETAB = ABS(NRETAB)

 ! Max radius should be set to twice the largest Re
  MAXRADIUS = 2*ERETAB

!Find proper file in the following directory: /mnt/b/projects/sciteam/jq0/LookupTables/final/500

  ncStatus(:) = nf90_NoErr
  if(nf90_open(trim(frequencyFile), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      call setStateToFailure(status, "Can't open file " // trim(frequencyFile))
  end if
    ! read in freq_grid(f_grid_nelem)
  if(.not. stateIsFailure(status)) then
        ncstatus( 1)=nf90_inq_dimid(ncFileId, "f_grid_nelem", ncDimId)
        ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nfreq)
        ALLOCATE(freqs(1:nfreq), lambdas(1:nfreq))
        ncStatus( 3) = nf90_inq_varid(ncFileId, "freq_grid", ncVarId)
        ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, freqs)
        ! convert to wavelength
        lambdas = 2.99792458E14/freqs  ! in microns
        if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "MieSSPTable: " // trim(frequencyFile) // &
                               " problem reading frequencies")
  end if			

! allocate and fill REFF
  ALLOCATE(REFF(NRETAB))
  ! At which values of effective radius do we build the table?
  IF (NRETAB == 1) THEN
        REFF(:nReTab) = SRETAB
  ELSE
        IF (LOGSPACEDREFF) THEN
            REFF(:) = SRETAB*(ERETAB/SRETAB)**(FLOAT( (/ (i, i = 0, nReTab - 1)/))/(NRETAB-1))
        ELSE
            REFF(:) = SRETAB + (ERETAB - SRETAB) * FLOAT((/ (i, i = 0, nReTab - 1) /) )/(NRETAB-1)
        ENDIF
  ENDIF


  if (thisProc .gt. 0) deallocate(freqs)

! loop over frequencies
  DO n = nfreq-thisProc, 1, -1*(numProcs-1)
	PRINT *, "now computing frequency index: ", n
	! only proceed if first iteration or thisProc .gt.0
     if(n .eq. nfreq .or. thisProc .gt. 0)then
		! determine MAXLEG
	 WAVELEN1 = lambdas(n)
       	WAVELEN2 = WAVELEN1

       	CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)

        ! Maximum number of Legendre coefficients to compute (for largest
        !    size parameter)
        ! Reference Wiscombe, W.J., 1980: Improved Mie scattering algorithms
        !   Appl. Opt. Vol 19, pg 1505-1509.
        !
        XMAX = 2*PI*MAXRADIUS/WAVELENCEN
        MAXLEG = NINT(2*(XMAX + 4.0*XMAX**0.3334 + 2))
	if(thisProc .eq. 0)PRINT *, "maximum LG coefficients needed: ", MAXLEG
		 
		! get RINDEX
	RINDEX =  GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2)
		
		! determine number of radii needed, nsize
	nsize = GET_NSIZE(SRETAB, MAXRADIUS, WAVELENCEN)
        !PRINT *, "total number of individual radii needed: ", NSIZE
		
		! allocate vectors of length nsize
		
        ! Drop radius and number concentration (for all drops)
        ALLOCATE (RADII(NSIZE), ND(NSIZE))
        ! Extinction, single scattering albedo, number of Legendre coefficents, and
        !   the value of those coefficients
        ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), NLEG1(NSIZE))
		
		! allocate vectors of length NRETAB (except REFF which was already done and including start)
	 ALLOCATE (EXTINCT(NRETAB), SSALB(NRETAB), NLEG(NRETAB), tempStart(NRETAB))
		 
		! allocate arrays of size (MAXLEG+1, NSIZE)
	ALLOCATE(LEGCOEF(0:MAXLEG,NRETAB), LEGEN1(0:MAXLEG,NSIZE))
		
		! get discreet particle radii
	radii(:NSIZE) = GET_SIZES (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE)
		
		! do the mie comps for each radius
	CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXLEG, &
                              EXTINCT1, SCATTER1, NLEG1, LEGEN1)
							  
		! loop over effective radii
	DO I = 1, NRETAB
		! Set tabulated effective radius

		! Calculate the discrete size number concentrations (ND), which vary
		!   according to a truncated gamma or lognormal distribution,
		!   that gives the desired effective radius (REFF) and LWC (1 g/m^3).
		CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, RADII, REFF(I), ALPHA, ND(:))

		! Sum the scattering properties over the discrete size distribution
		extinct(i) = dot_product(nd(:nsize), extinct1(:nsize))
		scatter    = dot_product(nd(:nsize), scatter1(:nsize))
		LEGCOEF(:,I) = 0.0
		NL = 1
		do J = 1, NSIZE
		   NL = max(NL, nLeg1(J))
		   LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I) + ND(J) * LEGEN1(0:NL,J)
		END do

		LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I)/SCATTER
		tempStart(1)=1
		DO L = 0, NL
		   IF (LEGCOEF(L,I) .GT. 0.5E-5)then
 			NLEG(I) = L
			if(I .lt. NRETAB) tempStart(I+1)=tempStart(I)+nLeg(I)
		   END IF	
		ENDDO

    !
			! Sanity check - the first Legendre coefficient should be identically 1
    !
		IF (ABS(LEGCOEF(0,I)-1.0) > 0.0001) THEN
		   PRINT *,'Phase function not normalized for Reff=',REFF(I),LEGCOEF(0,I)
		   STOP
		ENDIF
		IF (EXTINCT(I) > 0.0) THEN
			SSALB(I) = SCATTER/EXTINCT(I)
		ENDIF
		LegCoef(0:Nleg(i), i) = LegCoef(0:Nleg(i), i) / (/ (2*l+1, l=0, Nleg(i)) /)
	ENDDO  ! end of effective radius loop

	if(thisProc .eq. 0)EXIT
		
		! put LegCoef into buffer to send
	ALLOCATE(buffr(1:SUM(nLeg)))
	DO m = 1, NRETAB
		buffr(tempStart(m):tempStart(m)+nLeg(m)-1)=LegCoef(1:nLeg(m),m)
	ENDDO
	DEALLOCATE(LegCoef)
		
		! send arrays to root (EXTINCT, SSALB, nLeg, start, LegCoeff)
	CALL MPI_SEND(n, 1, MPI_INT, 0, ind, MPI_COMM_WORLD, ierr)
	CALL MPI_SEND(EXTINCT, NRETAB, MPI_REAL8, 0, ext, MPI_COMM_WORLD, ierr)
	CALL MPI_SEND(SSALB, NRETAB, MPI_REAL8, 0, ssa, MPI_COMM_WORLD, ierr)
	CALL MPI_SEND(nLeg, NRETAB, MPI_INT, 0, len, MPI_COMM_WORLD, ierr)
	CALL MPI_SEND(tempStart, NRETAB, MPI_INT, 0, st, MPI_COMM_WORLD, ierr)
	CALL MPI_SEND(buffr, SUM(nLeg), MPI_REAL, 0, LG, MPI_COMM_WORLD, ierr)
		
		
		! deallocate variables allocated within freq loop
	DEALLOCATE (RADII, ND)
        DEALLOCATE (EXTINCT1, SCATTER1, NLEG1)
	DEALLOCATE (EXTINCT, SSALB, NLEG, tempStart)
	DEALLOCATE (LEGEN1, buffr)
     end if
  ENDDO ! end of frequency loop

	
! root process
  if(thisProc .eq. 0)then
	! deallocate lambdas
	DEALLOCATE(lambdas)
	
	! sum nLeg
	MAXLEG=SUM(nLeg)+5
	
	! allocate total arrays (TOTALEXT, TOTALSSA, TOTALLG, start, length)
	ALLOCATE(TOTALEXT(NRETAB,nfreq), TOTALSSA(NRETAB,nfreq),start(NRETAB,nfreq), length(NRETAB,nfreq), TOTALLG(MAXLEG,nfreq))
	
	! fill with results from first iteration
	TOTALEXT(:,nfreq)=EXTINCT
	TOTALSSA(:,nfreq)=SSALB
	length(:,nfreq)=nLeg
	start(:,nfreq)=tempStart
	DO m = 1, NRETAB
		TOTALLG(tempStart(m):tempStart(m)+nLeg(m)-1,nfreq)=LegCoef(1:nLeg(m),m)
	ENDDO
	
	! DEALLOCATE single frequency vars that we no longer need or need to reallocate later
	DEALLOCATE (RADII, ND, EXTINCT, SSALB, nLeg, tempStart)
    	DEALLOCATE (EXTINCT1, SCATTER1, NLEG1)
	DEALLOCATE(LEGCOEF, LEGEN1)
	
	! initialize output file
	tableDescription = "Mie phase function table for spheres made of water at a concentration of 1 g/m^3. Key is in microns. "
    	select case (trim(distflag))
      	case('g', 'G')
       	tableDescription = trim(tableDescription) // " Gamma size distribution. "
      	case('l', 'L')
       	tableDescription = trim(tableDescription) // " Lognormal size distribution. "
    	end select
	
	ncStatus(:) = nf90_NoErr
    	ncStatus(1)=  nf90_create(trim(phaseFunctionTableFile), NF90_NETCDF4, ncFileID) ! need to use netcdf4 when file size will be over 2GB
    	if(ANY(ncStatus == nf90_NoErr)) then
        ! define dimensions
        	ncStatus(2) = nf90_def_dim(ncFileId, "f_grid_nelem", nfreq, fDimId)
        	ncStatus(3) = nf90_def_dim(ncFileId, trim(thisPrefix) // "phaseFunctionNumber", NRETAB, entryDimId)
        	ncStatus(4) = nf90_def_dim(ncFileId, trim(thisPrefix) // "maxCoefficients", size(TOTALLG(:,1)), coeffDimId)
        ! define variables
        	ncStatus( 5) = nf90_def_var(ncFileId, "f_grid", nf90_double, fDimId, VarId(1))
!       ncStatus( 6) = nf90_def_var(ncFileId, "surfaceAlbedo", nf90_double, fDimId, VarId(2))
        	ncStatus( 7) = nf90_def_var(ncFileId, trim(thisPrefix) // "ExtinctionT", nf90_double, (/entryDimId,fDimId/), VarId(3))
        	ncStatus( 8) = nf90_def_var(ncFileId, trim(thisPrefix) // "SingleScatterAlbedoT", nf90_double, (/entryDimId,fDimId/), VarId(4))
        	ncStatus( 9) = nf90_def_var(ncFileId, trim(thisPrefix) // "phaseFunctionKeyT", nf90_float, entryDimId, VarId(5))
        	ncStatus(10) = nf90_def_var(ncFileId, trim(thisPrefix) // "start", nf90_int, (/entryDimId,fDimId/), VarId(6))
        	ncStatus(11) = nf90_def_var(ncFileId, trim(thisPrefix) // "length", nf90_int, (/entryDimId,fDimId/), VarId(7))
        	ncStatus(12) = nf90_def_var(ncFileId, trim(thisPrefix) // "legendreCoefficients", nf90_float, (/coeffDimId,fDimId/), VarId(8))
        ! put attributes
		ncStatus(13) = nf90_put_att(ncFileID,nf90_Global, "numberOfComponents", 1)
        	ncStatus(14) = nf90_put_att(ncFileID,nf90_Global, "title", "SSP Table based on Mie Properties from I3RC's Mie Tool")
        	ncStatus(15) = nf90_put_att(ncFileID,nf90_Global, trim(thisPrefix) // "Name", "Water Droplets")
        	ncStatus(16) = nf90_put_att(ncFileID,nf90_Global, trim(thisPrefix) // "description", tableDescription)
        	ncStatus(17) = nf90_put_att(ncFileID,nf90_Global, trim(thisPrefix) // "zLevelBase", 1)
        	ncStatus(18) = nf90_put_att(ncFileID,nf90_Global, trim(thisPrefix) // "extType", "volExt")
        	ncStatus(19) = nf90_put_att(ncFileID,nf90_Global, trim(thisPrefix) // "phaseFunctionStorageType", "LegendreCoefficients")
        	ncStatus(20) = nf90_put_att(ncFileID,nf90_Global, "freqUnits", "Hz")
        	ncStatus(21) = nf90_put_att(ncFileID,nf90_Global, "starting_freq", freqs(1))
        	ncStatus(22) = nf90_put_att(ncFileID,nf90_Global, "ending_freq", freqs(nfreq))
    	end if
    	ncStatus(23) = nf90_EndDef(ncFileID)
    		
	work = nfreq-1
	DO WHILE (work .gt. 0)
		! listen for and recv arrays from other procs then put them into their proper place
		CALL MPI_RECV(n, 1, MPI_INT, MPI_ANY_SOURCE, ind, MPI_COMM_WORLD, mpiStatus, ierr)
		src = mpiStatus(MPI_SOURCE)
		CALL MPI_RECV(TOTALEXT(:,n), NRETAB, MPI_REAL8, src, ext, MPI_COMM_WORLD, mpiStatus, ierr)
		CALL MPI_RECV(TOTALSSA(:,n), NRETAB, MPI_REAL8, src, ssa, MPI_COMM_WORLD, mpiStatus, ierr)
		CALL MPI_RECV(length(:,n), NRETAB, MPI_INT, src, len, MPI_COMM_WORLD, mpiStatus, ierr)
		CALL MPI_RECV(start(:,n), NRETAB, MPI_INT, src, st, MPI_COMM_WORLD, mpiStatus, ierr)
		CALL MPI_PROBE(src, LG, MPI_COMM_WORLD, mpiStatus, ierr)
		CALL MPI_GET_COUNT(mpiStatus, MPI_REAL, count)
		CALL MPI_RECV(TOTALLG(1:count, n), count, MPI_REAL, src, LG, MPI_COMM_WORLD, mpiStatus, ierr)
		
		! decrement remaining work counter
		work = work - 1
	END DO	
	! do any final transformations on variables
	! convert to units of extinction in gm^-3
	TOTALEXT=TOTALEXT * 0.001
	! fill output file and deallocate variables as you go
	if(any(ncStatus(:) /= nf90_noErr)) call setStateToFailure(status, "MIESSPTableCreate: problem opening or defining SSP file.")
    	if(.not. stateIsFailure(status)) then
        	ncStatus = nf90_noErr
		ncStatus(1)=nf90_put_var(ncFileId, varId(1), freqs)
		DEALLOCATE(freqs)
!		ncStatus(2)=nf90_put_var(ncFileId, varId(2), albedo)
!		DEALLOCATE(albedo)
		ncStatus(3)=nf90_put_var(ncFileId, varId(3), TOTALEXT)
		DEALLOCATE(TOTALEXT)
        	ncStatus(4)=nf90_put_var(ncFileId, varId(4), TOTALSSA)
		DEALLOCATE(TOTALSSA)
        	ncStatus(5)=nf90_put_var(ncFileId, varId(5), REFF)
		DEALLOCATE(REFF)
        	ncStatus(6)=nf90_put_var(ncFileId, varId(6), start)
		DEALLOCATE(start)
        	ncStatus(7)=nf90_put_var(ncFileId, varId(7), length)
		DEALLOCATE(length)
        	ncStatus(8)=nf90_put_var(ncFileId, varId(8), TOTALLG)
		DEALLOCATE(TOTALLG)
	end if
	ncStatus(23) = nf90_close(ncFileId)
  end if
		
! DEALLOCATE remaining allocated variables if allocated (freqs, lambdas)
 IF(ALLOCATED(freqs))DEALLOCATE(freqs)
 IF(ALLOCATED(lambdas))DEALLOCATE(lambdas)
 IF(ALLOCATED(REFF))DEALLOCATE(REFF)
 
! finalize processes
  call finalizeProcesses


	CONTAINS
	
! ------------------------------------------------------------------------------

  elemental function planckRadiation(wavelength, temperature)
    implicit none
    real, intent(in) :: wavelength, temperature
    real             :: planckRadiation
    real             :: h, c, k

    h = 6.62606957e-34   ! [J*s]
    c = 2.99792458e+8    ! [m/s]
    k = 1.3806488e-23    ! [J/K]
  !
  ! Computes the Planck blackbody radiation (W/m2-str) as a function of
  !   wavelength (microns) and black body temperature (K)
  !

    planckRadiation = ((h*c**2)/wavelength**5)/(EXP((h*c)/(wavelength * temperature * k))-1.)

  end function planckRadiation

! ------------------------------------------------------------------------------

   function effectiveBlackBodyTemp(wavelength1, wavelength2)
     implicit none
     real, intent(in) :: wavelength1, wavelength2 ! in microns
     real             :: effectiveBlackBodyTemp   ! K
     !
     ! Computes the black body temperature at the center
     !   of a wavelength interval. In the
     !   shortwave (wavelengths < 3 microns) we use the sun's
     !   temperature (5800K); in the longwave (> 5 microns) we use a temperaure
     !   of 270K, and in between we set the BB temperature = -1.
     !
     real, parameter :: solarBBTemp = 5800., terrestrialBBTemp = 270., &
                        maxSolarBBWavelen = 3.0, minTerrestBBWavelen = 5.0
     real :: WAVECEN

     WAVECEN = 0.5*(wavelength1+wavelength2)
     IF (WAVECEN < maxSolarBBWavelen) THEN
       effectiveBlackBodyTemp = solarBBTemp
     ELSE IF (WAVECEN > minTerrestBBWavelen) THEN
       effectiveBlackBodyTemp = terrestrialBBTemp
     ELSE
       effectiveBlackBodyTemp = -1.0
     ENDIF

   end function effectiveBlackBodyTemp

! -----------------------------------------------------------------------------

  function computeNumPlanckWeightingWaves(wavelength1, wavelength2)
    implicit none
    real, intent(in) :: wavelength1, wavelength2
    integer          :: computeNumPlanckWeightingWaves
    !
    ! Compute how many indivdual wavelengths are required in order
    !   to weight a quantity by the Planck function across a wavelength
    !   interval
    !
    real :: wavecen, delwave

    IF (wavelength1 == wavelength2) THEN
      computeNumPlanckWeightingWaves = 1
    ELSE
      WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
      DELWAVE = MIN(WAVECEN/100., 0.1*ABS(WAVELEN2-WAVELEN1))
      DELWAVE = MAX(DELWAVE, WAVECEN*1.0E-5)
      computeNumPlanckWeightingWaves = int(ABS(WAVELEN2-WAVELEN1)/delwave)
    END IF

  end function computeNumPlanckWeightingWaves

! ------------------------------------------------------------------------------

  subroutine planckWeightingWavelengths(wavelength1, wavelength2, wavelengths)
    real, intent(in)  :: wavelength1, wavelength2
    real, dimension(0:), &
          intent(out) :: wavelengths
    !
    ! Returns the wavelengths requiredto weight a quantity by the Planck
    !  function across a wavelength interval
    !

    integer :: nWavelengths, i

    nWavelengths = size(wavelengths) - 1
    if(computeNumPlanckWeightingWaves(wavelength1, wavelength2) /= nWavelengths) &
      stop 'planckWeightingWavelengths: wrong number of wavelengths requested'

    IF (wavelength1 == wavelength2) THEN
      wavelengths(0:) = wavelength1
    ELSE
      wavelengths(0:) = wavelength1 + (wavelength2 - wavelength1) * &
                                      float((/ (i, i = 0, nWavelengths) /))/ nWavelengths
    END IF
  end subroutine planckWeightingWavelengths

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------


  SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  !  Returns the Planck weighted center wavelength averaged over the
  ! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
  ! blackbody temperature of 5800 K is used for the Planck weighting
  ! if the average wavelength is less than 3 microns, no Planck weighting
  ! is done for an average wavelength between 3 and 5 microns, and a
  ! blackbody temperature of 270 K is done for an average wavelength
  ! greater than 5 microns.
    IMPLICIT NONE
    REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
    REAL, INTENT(OUT) :: WAVELENCEN

    ! Local variables
    REAL              :: BBTEMP
    integer           :: nSteps
    real, allocatable :: wavelengths(:), planckValues(:)

    IF (WAVELEN1 == WAVELEN2) THEN
      WAVELENCEN = WAVELEN1
    ELSE
      BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)

      nSteps = computeNumPlanckWeightingWaves(WAVELEN1, WAVELEN2)
      allocate(wavelengths(0:nSteps), planckValues(0:nSteps))
      call planckWeightingWavelengths(WAVELEN1, WAVELEN2, wavelengths)

      if(bbtemp > 0) then
        planckValues(:) = planckRadiation(wavelengths(:), bbtemp)
      else
        planckValues(:) = 1.
      end if
      wavelencen = .001 * int(1000 * dot_product(planckValues(:), wavelengths(:)) / &
                                     sum(planckValues(:)) )
      deallocate(wavelengths, planckValues)
    ENDIF
  END SUBROUTINE GET_CENTER_WAVELEN

! ------------------------------------------------------------------------------
  function GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2) result(RINDEX)
   ! Returns the index of refraction for water or ice averaged over
   ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
   ! averaging is done at 0.05 micron intervals and is weighted by
   ! a Planck function at 5800 K temperature for central wavelengths
   ! less than 3 microns, a flat function between 3 and 5 microns, and
   ! 270 K Planck function beyond 5 microns.
   ! The index of refraction is using parameters iceTemperature and
   !   waterTemperature in the main program
   ! (the temperature dependence is important in the microwave).
    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
    REAL,             INTENT(IN) :: WAVELEN1, WAVELEN2
    COMPLEX                      :: RINDEX

    REAL    :: BBTEMP
    REAL    :: MRE, MIM, A
    integer :: i, nsteps
    real, allocatable :: M_real(:), M_Complex(:), wavelengths(:), planckValues(:)

    BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)

    nSteps = computeNumPlanckWeightingWaves(WAVELEN1, WAVELEN2)
    allocate(wavelengths(0:nSteps), planckValues(0:nSteps), &
             M_real(0:nSteps), M_Complex(0:nSteps))
    call planckWeightingWavelengths(WAVELEN1, WAVELEN2, wavelengths)

    if(BBTEMP > 0) then
      planckValues = planckRadiation(wavelengths, BBTEMP)
    else
      planckValues = 1.
    end if

    IF (PARTYPE == 'I') THEN
      do i = 0, nSteps
        CALL REFICE (0, wavelengths(i),   iceTemperature, M_real(I), M_Complex(i), A, A)
      end do
    ELSE
      do i = 0, nSteps
        CALL REFWAT (0, wavelengths(i), waterTemperature, M_real(I), M_Complex(i), A, A)
      end do
    end if

    MRE = dot_product(planckValues(:), M_real(:)   ) / sum(planckValues(:))
    MIM = dot_product(planckValues(:), M_Complex(:)) / sum(planckValues(:))
    RINDEX = CMPLX(MRE, -MIM)
  END function GET_REFRACT_INDEX

! ------------------------------------------------------------------------------
  function GET_NSIZE (SRETAB, MAXRADIUS, WAVELEN) result(nsize)
   ! Calculates the number of radii for which the Mie computation will be run.
   ! The formula and spacing in size parameter can be changed to trade
   ! off size distribution integration accuracy vs. computer time.
    IMPLICIT NONE
    REAL,    INTENT(IN)  :: SRETAB, MAXRADIUS, WAVELEN
    INTEGER              :: NSIZE

    ! Local variables
    REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD

    TWOPI = 2.0*ACOS(-1.0)
    RADMIN = 0.02*SRETAB
    RAD = RADMIN
    NSIZE = 1
    DO WHILE (RAD < MAXRADIUS)
      X = TWOPI*RAD/WAVELEN
      DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
  !    DELX = 0.1                     ! One alternative method
      DELRAD = DELX*WAVELEN/TWOPI
      RAD = RAD + DELRAD
      NSIZE = NSIZE + 1
    ENDDO
  END function GET_NSIZE

! ------------------------------------------------------------------------------
  function GET_SIZES (SRETAB, MAXRADIUS, WAVELEN, NSIZE) result(radii)
   ! Calculates the radii for which the Mie computation will be run and
   ! from which all the size distributions will be computed.
   ! The formula and spacing in size parameter can be changed to trade
   ! off size distribution integration accuracy vs. computer time.
    IMPLICIT NONE
    INTEGER, INTENT(IN ) :: NSIZE
    REAL,    INTENT(IN ) :: SRETAB, MAXRADIUS, WAVELEN
    REAL                 :: RADII(NSIZE)

    ! Local variables
    INTEGER :: N
    REAL    :: TWOPI, RAD, X, DELX, DELRAD

    TWOPI = 2.0*ACOS(-1.0)
    RAD = 0.02*SRETAB
    RADII(1) = RAD
    DO N = 2, NSIZE
      X = TWOPI*RAD/WAVELEN
      DELX = MAX(0.01,0.03*sqrt(X))    ! coarser spacing at large size parameters
  !    DELX = 0.1                     ! One alternative method
      DELRAD = DELX*WAVELEN/TWOPI
      RAD = RAD + DELRAD
      RADII(N) = RAD
    ENDDO
  END function GET_SIZES

! ------------------------------------------------------------------------------
  SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                    PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                    MAXLEG, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
   ! Does a Mie computation for each particle radius in RADII and returns the
   ! optical properties in arrays EXTINCT1, SCATTER1, NLEG1, and LEGEN1.
   ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
   ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
   ! integration of the Mie properties over wavelength is performed for
   ! each radius.  For each wavelength, with spacing DELTAWAVE, the water
   ! or ice (depending on PARTYPE) index of refraction is obtained and
   ! used in the Mie computation for that wavelength, and the Mie optical
   ! properties are averaged with Planck function weighting (blackbody
   ! temperature depends on wavelength).  The Legendre coefficients are
   ! returned with the product of the phase function times the scattering
   ! coefficient.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NSIZE, MAXLEG
    REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
    REAL,    INTENT(IN) :: RADII(NSIZE)
    COMPLEX, INTENT(IN) :: RINDEX
    CHARACTER(LEN=1), &
             INTENT(IN) :: AVGFLAG, PARTYPE
    INTEGER, INTENT(OUT) :: NLEG1(NSIZE)
    REAL(8),    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
    REAL,    INTENT(OUT) :: LEGEN1(0:MAXLEG,NSIZE)

    ! Local variables
    INTEGER :: I, NL
    REAL    :: WAVE, BBTEMP, PLANCK, SUMP, A
    REAL    :: MRE, MIM, EXT, SCAT, LEG(0:MAXLEG)
    COMPLEX :: REFIND

    IF (AVGFLAG == 'C') THEN
       ! For using one central wavelength: just call Mie routine for each radius
      DO I = 1, NSIZE
        CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXLEG, &
                      EXTINCT1(I), SCATTER1(I), NLEG1(I), LEGEN1(0,I) )
      ENDDO

    ELSE
      ! For averaging over wavelength range:
      BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)

      EXTINCT1(:) = 0.0_8
      SCATTER1(:) = 0.0_8
      NLEG1(:) = 1
      LEGEN1(:,:) = 0.0
      SUMP = 0.0
      WAVE = WAVELEN1
      DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
        IF (BBTEMP > 0) PLANCK = planckRadiation(WAVE, BBTEMP)
        SUMP = SUMP + PLANCK
        IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
          CALL REFICE (0, WAVE, iceTemperature,   MRE, MIM, A, A)
        ELSE
          CALL REFWAT (0, WAVE, waterTemperature, MRE, MIM, A, A)
        ENDIF
        REFIND = CMPLX(MRE,-MIM)
        DO I = 1, NSIZE
          CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXLEG, EXT, SCAT, NL, LEG)
          EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
          SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
          NLEG1(I) = MAX(NLEG1(I),NL)
          LEGEN1(0:NL,I) = LEGEN1(0:NL,I) + PLANCK*LEG(0:NL)
        ENDDO
        WAVE = WAVE + DELTAWAVE
      ENDDO
      EXTINCT1(:) = EXTINCT1(:)/SUMP
      SCATTER1(:) = SCATTER1(:)/SUMP
      LEGEN1(:,:) = LEGEN1(:,:)/SUMP
    ENDIF

  END SUBROUTINE COMPUTE_MIE_ALL_SIZES

! ------------------------------------------------------------------------------
  SUBROUTINE MAKE_SIZE_DIST (DISTFLAG, PARDENS, RADII, REFF, ALPHA, ND)
   ! Calculates the number concentrations (ND in cm^-3) at
   ! discrete particle RADII (micron) of a gamma or lognormal size distribution
   ! with an effective radius of REFF (micron), gamma shape parameter or
   ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN ) :: DISTFLAG
    REAL,             INTENT(IN ) :: RADII(:), REFF, ALPHA, PARDENS
    REAL,             INTENT(OUT) :: ND(:)

    REAL,    PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
    integer, parameter :: maxIterations = 8
    INTEGER :: I, NSIZE
    REAL    :: TRUERE, F, REHI, RELO, REMID

    NSIZE = size(radii)
    if(size(nd) /= NSIZE) &
      stop 'MAKE_SIZE_DIST: vectors RADII and ND must be the same length'

     ! See if the true effective radius is already close enough
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REFF, RADII, ND, TRUERE)
    IF (ABS(TRUERE-REFF) < TOL*REFF) RETURN
    F = REFF/TRUERE

    IF (TRUERE < REFF) THEN
      ! Find Reff that gives true Reff above desired value
      RELO = REFF
      REHI = F*REFF
      I = 0
      TRUERE = REFF/F
      DO WHILE (TRUERE <= REFF .AND. I < maxIterations)
        REHI = F*REHI
        I = I + 1
        CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, REHI, RADII, ND, TRUERE)
      ENDDO
      IF (TRUERE <= REFF) THEN
        PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
        STOP
      ENDIF
    ELSE
      ! Find Reff that gives true Reff below desired value
      REHI = REFF
      RELO = F*REFF
      I = 0
      TRUERE = REFF/F
      DO WHILE (TRUERE >= REFF .AND. I < maxIterations)
        RELO = F*RELO
        I = I + 1
        CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, RELO, RADII, ND, TRUERE)
      ENDDO
      IF (TRUERE >= REFF) THEN
        PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
        STOP
      ENDIF
    ENDIF
    ! Do bisection to get correct effective radius
    DO WHILE (ABS(TRUERE-REFF) > TOL*REFF)
      REMID = 0.5*(RELO+REHI)
      CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REMID, RADII, ND, TRUERE)
      IF (TRUERE < REFF) THEN
        RELO = REMID
      ELSE
        REHI = REMID
      ENDIF
    ENDDO
  END SUBROUTINE MAKE_SIZE_DIST

! ------------------------------------------------------------------------------
  SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, RE, RADII, ND, TRUERE)
   ! For the input effective radius (RE) [um], returns the number concentrations
   ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a
   ! gamma or lognormal size distribution with mass content of 1 g/m^3.
    IMPLICIT NONE
    REAL,    INTENT(IN) :: PARDENS
    CHARACTER(LEN=*), &
             INTENT(IN)  :: DISTFLAG
    REAL,    INTENT(IN)  :: ALPHA, RE, RADII(:)
    REAL,    INTENT(OUT) :: ND(:), TRUERE

    INTEGER :: nSize
    REAL    :: PI, A, B, LWC, SUM2, SUM3
    real, dimension(:) &
            :: deltaR(size(RADII))
    ! External function
    real    :: GAMMLN

    nSize = size(radii)
    if(size(nd) /= nSize) &
      stop 'DO_SIZE_DIST: vectors RADII and ND must be the same length'
    PI = ACOS(-1.0)

    deltaR(2:nSize-1) = sqrt(radii(2:nSize-1) * radii(3:nSize)) - &
                        sqrt(radii(2:nSize-1) * radii(1:nSize-2))
    deltaR(1        ) = sqrt(radii(2    )     * radii(3   )) - &
                        radii(1)
    deltaR(nSize    ) = radii(nsize) - &
                        sqrt(radii(nsize) * radii(nSize-1))

    IF (DISTFLAG == 'G') THEN       ! Gamma distribution
      B = (ALPHA+3)/RE
      A = 1.E6/( (4*PI/3.)*PARDENS *B**(-ALPHA-4) *EXP(GAMMLN(ALPHA+4.)) )
      ND(:) = A * RADII(:)**ALPHA *EXP(-B*RADII(:)) * deltaR(:)
    ELSE                            ! Lognormal distribution
      B = RE*EXP(-2.5*ALPHA**2)
      A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
      ND(:) = (A/RADII(:))*EXP(-0.5*(LOG(RADII(:)/B))**2/ALPHA**2) * deltaR(:)
    ENDIF

    SUM2 = dot_product(nd(:), RADII(:)**2)
    SUM3 = dot_product(nd(:), RADII(:)**3)
    TRUERE = SUM3/SUM2
    LWC = 1.0E-6 * PARDENS * (4.*PI/3.) * SUM3
    ND(:) = (1.0/LWC) * ND(:)

  END SUBROUTINE DO_SIZE_DIST
! ------------------------------------------------------------------------------		

	End Program MakeMieSSPTable
