!! Maxwell smith, 28 Jun 2011
!! Department of Atmospheric Sciences, University of Illinois at Urbana-Champaign
!Edited from i3rcStepCloud.f95, from the I3RC Monte Carlo Community Radiative Transfer Model
! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/I3RC-Examples/i3rcStepCloud.f95 $

program cubic
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  implicit none
  !
  ! Write the domains corresponding to the I3RC step cloud. 
  !   The domain is .5 km wide and 32 columns across. The first 16 
  !   columns have optical depth 2, the second 16 are optical depth 18. 
  !   We use a Henyey-Greenstein phase function with g=0.85. 
  !   We write two domains, one with single scattering albedo = 0.99 
  !   and the other entirely non-absorbing. 
  ! Two parameters weed need that aren't specified are the number of 
  !    Legendre moments and the physical thickness of the cloud. 
  !
  
  ! I3RC parameters
  real,    parameter :: domainSize = 2250, g = 0.85
  integer, parameter :: nColumns   = 225,nRows = 225, nSSAs = 2
  real, dimension(nSSAs), &
           parameter :: SSAs = (/ 1., 0.99 /) 
  

  ! Other parameters
  integer, parameter :: nLegendreCoefficients = 180, nLayers = 50
  real,    parameter :: physicalThickness = 500
  real               :: deltaX, deltaY, deltaZ 
  integer            :: D = nLayers*2 
  real               :: tau = 100
  character(len = 60), dimension(2), &
           parameter :: fileNames = (/ "Data/cubic_NonAbsorbing.dom", &
                                       "Data/cubic_Absorbing.dom   " /)
  
  ! Variables
  integer                           :: i,j,k
    
  real,    dimension(nColumns, nRows, nLayers) :: extinction=0., singleScatteringAlbedo=0., temps=0.
  integer, dimension(nColumns, nRows, nLayers) :: phaseFunctionIndex=0
  
  type(ErrorMessage)              :: status
  type(phaseFunction)             :: phase
  type(phaseFunctionTable)        :: table
  type(domain)                    :: stepCloud

  ! ------------------------------------------
  phase = &
    new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
  table = &
    new_PhaseFunctionTable((/ phase /), key = (/ 1. /), status = status)
  if(.not. stateIsFailure(status)) call setStateToCompleteSuccess(status)
  
  deltaX = domainSize/real(nColumns)
  deltaY = domainSize/real(nRows)
  deltaZ = physicalThickness/real(nLayers)
  !extinction(:, 1, :) = spread((/ (2, i = 1, nColumns/2), (18, i = 1, nColumns/2) /), &
  !                             dim = 2, nCopies = nLayers) / physicalThickness
  
  forall (i=D:(D+nLayers/2), j=D:(D+nLayers/2), k=(nLayers/2):nLayers)
     extinction(i,j,k) = tau*2/physicalThickness
  end forall
  
  where(extinction > 0.)
    phaseFunctionIndex(:, :, :) = 1
    singleScatteringAlbedo(:, :, :) = SSAs(1)
  end where
  
  !
  ! Define the domain
  !
  stepCloud =                                                             &
    new_Domain(xPosition = deltaX * (/ 0., (real(i), i = 1, nColumns) /), &
               yPosition = deltaY * (/ 0., (real(i), i = 1, nRows)    /), &
               zPosition = deltaZ * (/ 0., (real(i), i = 1, nLayers)  /), &
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
  where(extinction >0.)
    singleScatteringAlbedo(:, :, :) = SSAs(2)
  end where
  call replaceOpticalComponent(stepCloud, 1, "cloud: absorbing", &
                               extinction, singleScatteringAlbedo,   &
                               phaseFunctionIndex, table, status = status)
  call write_Domain(stepCloud, trim(fileNames(2)), status = status)
  print *, "Writing absorbing domain"; call printStatus(status)
  if(stateIsSuccess(status)) call setStateToCompleteSuccess(status)
end program cubic
