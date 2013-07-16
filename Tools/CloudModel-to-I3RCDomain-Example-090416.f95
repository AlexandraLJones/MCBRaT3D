!
! This program provides an example of creating and writing a domain file for the 
!   I3RC Community Monte Carlo model using existing files and Fortran calls. It's meant 
!   as a demonstration, to help you write your own adapter; you should use it as a 
!   flexible guide and not a rigid formula. 
!
! The cloud field is taken from a large-eddy simulation 
!   snapshot in a netcdf format. Effective radius is determined from the liquid water
!   content and (specified) drop number concentration, and the appropriate optical properties 
!   are taken from a phase function table (pre-computed using Tools/MakeMieTable). 
! The domain includes horizontally-uniform aerosols whose optical properties come from a 
!   separate pre-computed phase function table. 
! The domain also includes Rayleigh scattering using the same parameterization used in 
!   Tools/PhysicalPropertiesToDomain. 
!
! The code should be linked with all the object files in the Code/ directory. 
!
! The resulting domain file can be used as input to Example-Drivers/MonteCarloDriver to calculate
!   fluxes, flux divergence, and intensity
!
program RicoDomain
  use netcdf
  use ErrorMessages, only: ErrorMessage
  use UserInterface, only: printStatus
  use numericUtilities, only : findIndex
  use scatteringPhaseFunctions
  use opticalProperties
  implicit none
  
  real, parameter :: Pi = 3.141593, R_air = 287.05, cp = 1004., k = R_air/cp
  real, parameter :: wavelength = 0.67, DropNumConc = 70., &
                     totalAerosolTau = .233, aerosolLayerDepth = 3. ! in km
  !
  ! Aerosol optical depth is at 0.675; needs to be scaled where it's computed otherwise
  !
                     
  integer :: iz, ix, iy
  integer :: nx = -1, ny = -1, nz
  real    :: dx, dy, dz
  real    :: re, maxRe, minRe, fraction
  real,    dimension(:),       allocatable :: xPos, yPos, zPos, & 
                                              lwc, pressure, temperature, pressureMean, temperatureMean
  real,    dimension(:, :, :), allocatable :: extinction, ssa
  integer, dimension(:, :, :), allocatable :: reIndex
  real,    dimension(:),       allocatable :: extinction1D, ssa1D
  integer, dimension(:),       allocatable :: reIndex1D
  
  !
  ! File names
  !
  character(len = 256) :: cloudPhaseFunction = "cloud_w0.67_mie.phasetab", &
                          aerosolPhaseFunction = "cu_aero_w0.67_mie.phasetab", &
                          RicoFileName = "rico.384x384x200.nc", &
                          domainFileName = "cu_w0.67.dom"

  !
  ! I3RC variables 
  !
  type(ErrorMessage)       :: status
  type(phaseFunction), dimension(1) &
                           :: phaseFuncs
  type(phaseFunctionTable) :: ppTable
  type(domain)             :: Rico
  real, dimension(:), allocatable &
                           :: tabulatedRe, tabulatedExtinction, tabulatedSSA
  integer                  :: nEntries
  
  ! Netcdf vars
  integer :: ncFileId, ncDimId, ncVarId, pressureVarId, waterVarId, thetaVarId
  integer, dimension(32) :: ncStatus
  integer, dimension(2)  :: start = (/ 1, 1 /) 
  
  
  ! ---------------------------------------------------------------------------
  !
  ! Read phase function table; find values at which it's been tabulated
  !   and extinction, ssa values at each size
  !
  call read_PhaseFunctionTable(cloudPhaseFunction, table = ppTable, status = status)
  call printStatus(status)
  call getInfo_PhaseFunctionTable(ppTable, nEntries = nEntries,  status = status)
  call printStatus(status)
  allocate(tabulatedExtinction(nEntries), tabulatedSSA(nEntries), tabulatedRe(nEntries)) 
  call getInfo_PhaseFunctionTable(ppTable, key = tabulatedRe,       &
                                  extinction = tabulatedExtinction, &
                                  singleScatteringAlbedo = tabulatedSSA, status = status)
  call printStatus(status)
  minRe = minval(tabulatedRe); maxRe = maxval(tabulatedRe)
  print *, "Range of effective radius in table", minRe, maxRe
  
  ! -----------------------------------------------------------------
  !
  ! Domain size in each dimension
  !   Defaults to the whole domain but can be changed by specifying nx, ny, and start. 
  !
  nx = 96; ny = nx; start = (/ 257, 236 /)
  
  ncStatus( :) = nf90_NoErr
  ncStatus( 1) = nf90_open(trim(RIcoFileName), nf90_NoWrite, ncFileID)
  if(nx < 0) then 
    ncstatus( 2) = nf90_inq_dimid(ncFileId, "xt", ncDimId)
    ncstatus( 3) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nx)
  end if 
  if (ny < 0) then 
    ncstatus( 4) = nf90_inq_dimid(ncFileId, "yt", ncDimId)
    ncstatus( 5) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = ny)
  end if
  ncstatus( 6) = nf90_inq_dimid(ncFileId, "zt", ncDimId)
  ncstatus( 7) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nz)
  nz = nz - 1 ! Bottom layer mirrors across lower boundary
  !
  ! Here's where we restrict ourselves to the bottom 3 km, knowing that the vertical resolution is 20 m
  !
  nz = 150
  print *, "Domain size (x, y, z)", nx, ny, nz
  allocate(xPos(nx+1), yPos(ny+1), zPos(nz+1)) 

  !
  ! Read x, y, z positions (m), then liquid water mixing ratio (kg/kg)
  !
  ncstatus( 8) = nf90_inq_varid(ncFileId, "xt", ncVarId) 
  ncstatus( 9) = nf90_get_var(ncFileId, ncVarId, xPos(:nx), start = start(1:1))
  ncstatus(10) = nf90_inq_varid(ncFileId, "yt", ncVarId) 
  ncstatus(11) = nf90_get_var(ncFileId, ncVarId, yPos(:ny), start = start(2:2))
  ncstatus(12) = nf90_inq_varid(ncFileId, "zt", ncVarId) 
  ncstatus(13) = nf90_get_var(ncFileId, ncVarId, zPos(:nz), start = (/ 2 /))
  if(any(ncstatus(:13) /= nf90_NoErr)) stop "Can't read the positions!" 
  !
  ! Extinction in the phase function tables is in inverse km. 
  !
  xPos(:) = xPos(:)/1000.
  yPos(:) = yPos(:)/1000.
  zPos(:) = zPos(:)/1000.
  
  
  !
  ! LES fields are ordered (z, x, y); I3RC wants (x, y, z) 
  !   We'll read one column at a time from the LES results
  !
  allocate(extinction(nx, ny, nz), ssa(nx, ny, nz), reIndex(nx, ny, nz))
  allocate(lwc(nz), pressure(nz), temperature(nz), pressureMean(nz), temperatureMean(nz)) 
  allocate(extinction1D(nz), ssa1D(nz), reIndex1D(nz))                                           
  
  ! ---------------------------------------------------------------------------
  !
  ! Compute phase function index; look up extinction and single scattering albedo
  !
  pressureMean(:) = 0.; temperatureMean(:) = 0. 
  ncstatus(1) = nf90_inq_varid(ncFileId, "rl", waterVarId) 
  ncstatus(2) = nf90_inq_varid(ncFileId, "pressure", pressureVarId) 
  ncstatus(3) = nf90_inq_varid(ncFileId, "theta", thetaVarId) 
  if(any(ncstatus(:3) /= nf90_NoErr)) stop "Can't find the variables in LES file " // trim(RicoFileName)
  do iy = 1, ny
    do ix = 1, nx
      !
      ! Read liquid water content, pressure, and potential temperature one column at a time. 
      !
      ncstatus(1) = nf90_get_var(ncFileId, waterVarId,    lwc,      start = (/ 2, ix + start(1) - 1, iy + start(2) - 1 /))
      ncstatus(2) = nf90_get_var(ncFileId, pressureVarId, pressure, start = (/ 2, ix + start(1) - 1, iy + start(2) - 1 /))
      ncstatus(3) = nf90_get_var(ncFileId, thetaVarId, temperature, start = (/ 2, ix + start(1) - 1, iy + start(2) - 1 /))
      if(any(ncstatus(:3) /= nf90_NoErr)) stop "Can't read the variables" 
      pressureMean(:)    = pressureMean(:) + pressure(:)
      temperatureMean(:) = temperatureMean(:) + temperature(:) * (Pressure(:) / 100000.0)**k
      
      
      !
      ! Convert mixing ratio (kg/kg) to liquid water content (g/m3); 
      !   compute effective radius (microns)
      !   Temperature has to be computed from potential temperature 
      !   (we'll ignore vapor effects on density)
      !   Pressure is in Pa
      !
      lwc(:) = 1000. * lwc(:) * pressure(:) / (R_air * temperature(:) * (Pressure(:) / 1000.0)**k)
      do iz = 1, nz
        if(lwc(iz) > 0.) then 
          !
          ! The factor of 1.3889 comes from PhysicalPropertiesToDomain. Quoth Frank Evans
          !   "(alpha+3)^2/((alpha+1)*(alpha+2)) for gamma distribution parameter
          !    alpha=7.  Isn't this obvious? ;-) " 
          ! 
          re = 100. * (lwc(iz) * 0.75 * 1.3889/(Pi * DropNumConc))**(1./3)
          if(re < minRe .or. re > maxRe) then 
            print *,  "Re out of bounds at ix, iy, iz=", ix, iy, iz, ": lwc ", lwc(iz), "g/m3, re ", re, "microns"
            reIndex(ix, iy, iz) = 0 
            extinction(ix, iy, iz) = 0.  
            ssa(ix, iy, iz) = 0. 
          else  
            !
            ! Find re index at or below current value; interpolate extinction and single scattering 
            !   albedo linearly to this value of re
            !
            reIndex(ix, iy, iz) = findIndex(re, tabulatedRe)
            fraction = (re                                 - tabulatedRe(reIndex(ix, iy, iz))) / &
                       (tabulatedRe(reIndex(ix, iy, iz)+1) - tabulatedRe(reIndex(ix, iy, iz)))
            extinction(ix, iy, iz) = lwc(iz) * & 
                                     ((1. - fraction) * tabulatedExtinction(reIndex(ix, iy, iz)) +  & 
                                      fraction        * tabulatedExtinction(reIndex(ix, iy, iz)+1))
            ssa       (ix, iy, iz) = ((1. - fraction) * tabulatedSSA(reIndex(ix, iy, iz)) + &
                                      fraction        * tabulatedSSA(reIndex(ix, iy, iz)+1))
            !
            ! Use phase function from closest re. 
            !
            if(abs( re - tabulatedRe(reIndex(ix, iy, iz)) ) >     &
               abs( re - tabulatedRe(reIndex(ix, iy, iz) + 1) ))  &
              reIndex(ix, iy, iz) = reIndex(ix, iy, iz) + 1
          end if 
        else
          reIndex(ix, iy, iz) = 0 
          extinction(ix, iy, iz) = 0.  
          ssa(ix, iy, iz) = 0. 
        end if
      end do
    end do
  end do
  ncStatus(1) = nf90_close(ncFileId) 
  print *, "Cloudy volume fraction is ", count(reIndex(:, :, :) > 0)/real(size(reIndex))
  print *, "Projected cloud fraction is ", count(any(reIndex(:, :, :) > 0, dim = 3)) / real(nx * ny) 
  
  !
  ! Mean temperature and pressure profiles for computing Rayleigh extinction
  !
  pressureMean(:)    = pressureMean(:)/real(nx * ny) 
  temperatureMean(:) = temperatureMean(:)/real(nx * ny) 
  
  deallocate(lwc, pressure, temperature, tabulatedRe, tabulatedSSA, tabulatedExtinction)
  ! ---------------------------------------------------------------------------
  !
  ! Create the domain and add clouds
  !   Convert positions from grid centers to edges (assume constant spacing)
  !
  dx = xPos(2) - xPos(1); xPos(:) = (/ xPos(1) - dx/2,  xPos(:) + dx/2 /)
  dy = yPos(2) - yPos(1); yPos(:) = (/ yPos(1) - dy/2,  yPos(:) + dy/2 /)
  dz = zPos(2) - zPos(1); zPos(:) = (/ zPos(1) - dz/2,  zPos(:) + dz/2 /)
  Rico = new_Domain (xPos, yPos, zPos, status)
  call printStatus(status)
  print *, "Created the domain" 

  call addOpticalComponent (Rico, "Clouds", &
                            extinction, ssa, reIndex, ppTable,  &
                            status = status)
  call printStatus(status)
  print *, "Added the clouds" 
  deallocate(extinction, ssa, reindex) 
  call finalize_PhaseFunctionTable(ppTable)
  ! ---------------------------------------------------------------------------
  !
  ! Aerosols, which are all the same size
  !
  allocate(tabulatedExtinction(1), tabulatedSSA(1)) 
  call read_PhaseFunctionTable(aerosolPhaseFunction, table = ppTable, status = status)
  call printStatus(status)
  call getInfo_PhaseFunctionTable(ppTable, extinction = tabulatedExtinction, &
                                  singleScatteringAlbedo = tabulatedSSA, status = status)
  call printStatus(status)
  
  extinction1D(:) = 0; ssa1D(:) = 1. 
  reIndex1D(:) = 1
  where(zPos <= aerosolLayerDepth) 
    ! 
    ! We know that grid is spaced evenly in the vertical 
    !
    extinction1D(:) = totalAerosolTau/aerosolLayerDepth
    ssa1D(:) = tabulatedSSA(1) 
  end where 
  call addOpticalComponent (Rico, "aerosols", &
                            extinction1D, ssa1D, reIndex1D, ppTable,  &
                            status = status)
  print *, "Added the aerosol" 
  print *, "  Total optical depth is ", dz * sum(extinction1D) 
  call finalize_PhaseFunctionTable(ppTable)
  
  ! ---------------------------------------------------------------------------
  !
  ! Scattering by molecules. Parameterization from PhysicalPropertiesToDomain
  !
  ssa1D(:)     = 1.0
  reIndex1D(:) = 1
  ! 
  ! Rayleigh scattering phase function - the I3RC convention is *not* to include the 2l + 1 term in 
  !   the phase function expansion coefficients (DISORT convention, not SHDOM convention)
  !
  PhaseFuncs(1) = new_PhaseFunction ((/ 0.0, 0.5 /) / (/ 2.*1. + 1., 2.*2. + 1. /) , status=status)
  call printStatus(status)
  ppTable = new_PhaseFunctionTable (PhaseFuncs(1:1), key=(/ 0.0 /),&
                                    tableDescription = "Rayleigh scattering", &
                                    status=status)
  extinction1D(:) = (2.97E-4 * wavelength**(-4.15+0.2*wavelength)) * &
                     pressureMean(:) / 100. * 1./temperatureMean(:) 
 
  call addOpticalComponent (Rico, "Rayleigh", &
                            extinction1D, ssa1D, reIndex1D, ppTable,  &
                            status = status)
  call finalize_PhaseFunctionTable(ppTable)
  print *, "Added Rayleigh extinction" 
  print *, "  Total optical depth is ", dz * sum(extinction1D) 
  print *, "  Scales to ", dz * sum(extinction1D) * (maxval(PressureMean) - 0.) / (maxval(pressureMean) - minval(pressureMean))
  ! ---------------------------------------------------------------------------
  !
  ! Write out the domain 
  !
  
  call write_Domain(Rico, domainFileName, status)
  call printStatus(status)
  print *, "Wrote domain file ", trim(domainFileName) 
  call finalize_Domain(Rico) 

end program RicoDomain

