! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/I3RC-Examples/i3rcStepCloud.f95 $
program planeParallel
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  implicit none
  !
  ! Write the domains corresponding to a plane parallel cloud
  !   The domain is .5 km wide and 32 columns across. Optical
  !   depth is defined to be 1 everywhere. 
  !   We use a Henyey-Greenstein phase function with g=0.85. 
  !   We write two domains, one with single scattering albedo = 0.99 
  !   and the other entirely non-absorbing. 
  ! Two parameters we need that aren't specified are the number of 
  !    Legendre moments and the physical thickness of the cloud. 
  !
  
  ! I3RC parameters
  real,    parameter :: tau = 0.5
  real,    parameter :: domainSize = 500, g = 0.85
  integer, parameter :: nColumns   = 1, nSSAs = 2
  real, dimension(nSSAs), &
           parameter :: SSAs = (/ 1., 0.99 /) 
  

  ! Other parameters
  integer, parameter :: nLegendreCoefficients = 64, nLayers = 32
  real,    parameter :: physicalThickness = 250
  real               :: deltaX, deltaZ 
  character(len = 40), dimension(2), &
           parameter :: fileNames = (/ "Data/planeParallel_NonAbsorbing.dom", &
                                       "Data/planeParallel_Absorbing.dom   " /)
  
  ! Variables
  integer                           :: i
    
  real,    dimension(nColumns, 1, nLayers) :: extinction, singleScatteringAlbedo, temps
  integer, dimension(nColumns, 1, nLayers) :: phaseFunctionIndex
  
  type(ErrorMessage)              :: status
  type(phaseFunction)             :: phase
  type(phaseFunctionTable)        :: table
  ! type(InversePhaseFunctionTable) :: inverseTable
  type(domain)                    :: stepCloud

  ! ------------------------------------------
  phase = &
    new_PhaseFunction(g**(/ (real(i), i = 1, nLegendreCoefficients )/), status = status)
  table = &
    new_PhaseFunctionTable((/ phase /), key = (/ 1. /), status = status)
  if(.not. stateIsFailure(status)) call setStateToCompleteSuccess(status)
  
  deltaX = domainSize/real(nColumns)
  deltaZ = physicalThickness/real(nLayers)
  extinction(:, 1, :) = spread((/ (tau, i = 1, nColumns) /), &
                               dim = 2, nCopies = nLayers) / physicalThickness
  phaseFunctionIndex(:, :, :) = 1
  
  singleScatteringAlbedo(:, :, :) = SSAs(1)
  
  temps(:,:,:) =  0.0

  !
  ! Define the domain
  !
  stepCloud =                                                            &
    new_Domain(xPosition = deltaX * (/ 0., (real(i), i = 1, nColumns) /), &
               yPosition = (/ 0., 500.0 /),                               &
               zPosition = deltaZ * (/ 0., (real(i), i = 1, nLayers) /) , &
               temps=temps, status = status)
  !
  ! Add the clouds
  !
  call addOpticalComponent(stepCloud, "cloud: non-absorbing",  &
                           extinction, singleScatteringAlbedo, &
                           phaseFunctionIndex, table, status = status)
  call printStatus(status)
  !
  ! Write it to the file
  !
  call write_Domain(stepCloud, trim(fileNames(1)), status = status)
  print *, "Writing non-absorbing domain"; call printStatus(status)
  if(stateIsSuccess(status)) call setStateToCompleteSuccess(status)
  
  !
  ! Now write out the same domain but with SSA = =.99
  !
  singleScatteringAlbedo(:, :, :) = SSAs(2)
  call replaceOpticalComponent(stepCloud, 1, "cloud: absorbing", &
                               extinction, singleScatteringAlbedo,   &
                               phaseFunctionIndex, table, status = status)
  call write_Domain(stepCloud, trim(fileNames(2)), status = status)
  print *, "Writing absorbing domain"; call printStatus(status)
  if(stateIsSuccess(status)) call setStateToCompleteSuccess(status)
end program planeParallel
