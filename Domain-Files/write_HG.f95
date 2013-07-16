
!Maxwell Smith 05 Jul 2011
!Plane Parallel step cloud with 2 optical depth portions
!semi-infinite from the center of the domain
!Adapted from the I3RC Step cloud example code, credits below:



! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/I3RC-Examples/i3rcStepCloud.f95 $

program stepcloud_2tau
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
  real,    parameter :: g = 0.85
  

  ! Other parameters
  integer, parameter :: nLegendreCoefficients = 120
     
  ! Variables
  integer                           :: i
      
  type(ErrorMessage)              :: status
  type(phaseFunction)             :: phase
  type(phaseFunctionTable)        :: table

  integer, parameter :: nAngles = 1801
  real, dimension(nAngles) :: angles=0.
  real, dimension(nAngles) :: values=0.

  ! ------------------------------------------
  phase = &
    new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
  table = &
    new_PhaseFunctionTable((/ phase /), key = (/ 1. /), status = status)
  call write_PhaseFunctionTable(table, 'HG_i3rc_table.nc', status)

  angles=acos(-1.)/real(nAngles-1) * (/ (real(i),i = 0, nAngles-1) /)
  call getPhaseFunctionValues(phase, angles, values, status)
  if(.not. stateIsFailure(status)) call setStateToCompleteSuccess(status)

  open(66, file="HG_values.txt", form="formatted")
  do i = 1,nAngles
     write(66, 67) angles(i)*180./acos(-1.), values(i)
  end do
67 format(F5.1, F35.30)

  

end program stepcloud_2tau
