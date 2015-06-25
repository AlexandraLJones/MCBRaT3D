
program RayleighTrans

  use opticalProperties
  use netcdf
  use ErrorMessages
  use scatteringPhaseFunctions
  use characterUtils
  use inversePhaseFunctions
  use numericUtilities
  
  implicit none

!  Variable declarations
  real(8), parameter :: Na = 6.02214129E23 ![mol^-1]
  real(8), parameter :: Rstar = 8.3144621 ![J K^-1 mol^-1]
  integer, parameter :: nlevels=60
  type(commonDomain)         :: commonPhysical
  type(ErrorMessage)         :: status
  type(phaseFunctionTable)   :: table
  integer                    :: I, ncFileID, ncStatus, ncVarID
  real(8)                    :: lambda(1000),absx(nlevels-1),zPos(1:nlevels), depths(nlevels-1)
  real(8)                    :: trans
  real(8)                    :: temps(1,1,nlevels-1),prssr(1,1,nlevels-1), N(1,1,1,nlevels-1)
  real(8)                    :: ext(1,1,nlevels-1),ssa(1,1,nlevels-1)
  integer                    :: phaseInd(1,1,nlevels-1)
  
 
  lambda = (/(I/1000.0_8, I=201, 1200, 1)/)
  absx=0.0_8
!  call read_Common("/u/sciteam/aljones4/ARTS/Cases/CIRCcase7/input/CIRCcase7PhysDomain.dom", commonPhysical, status)
  ncStatus=nf90_open("/u/sciteam/aljones4/ARTS/Cases/CIRCcase7/input/CIRCcase7PhysDomain.dom", nf90_NoWrite, ncFileID)
  ncStatus=nf90_inq_varid(ncFileId, "z-edges", ncVarId)
  ncStatus=nf90_get_var(ncFileId, ncVarId, zPos)
  ncStatus = nf90_inq_varid(ncFileId, "Temperatures", ncVarId)
  ncStatus = nf90_get_var(ncFileId, ncVarId, temps)
  ncStatus = nf90_inq_varid(ncFileId,"Pressures", ncVarID)
  ncStatus = nf90_get_var(ncFileId, ncVarId, prssr)
  N(1,1,1,:) = (prssr(1,1,:)*100.0*Na)/(Rstar*temps(1,1,:))
  depths = zPos(2:nlevels)-zPos(1:nlevels-1)
!  PRINT *, depths
  open(unit=34, file='RayleighTrans.txt', status='UNKNOWN')

  DO i = 1, size(lambda)
	CALL calc_RayleighScattering(absx,lambda(i),N(1,1,1,:),ext(1,1,:),ssa(1,1,:),phaseInd(1,1,:),table,status)
	trans = exp(-1.0_8*SUM(depths*ext(1,1,:)))
!	PRINT *, lambda(i), trans
	write(34, "(5F20.16)") lambda(i), absx(1), ext(1,1,1), ssa(1,1,1), trans
  END DO

!  PRINT *, ext(1,1,:)
  
  close(34)

end program RayleighTrans
