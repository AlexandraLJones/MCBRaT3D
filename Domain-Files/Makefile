# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/I3RC-Examples/Makefile $
# ------------------------------------------------------------------
# Programs for creating domain files for the I3RC test cases
# Run the radiative transfer using the driver in Example-Drivers

include ../Makefile 

Modules += $(ModuleFlag)$(NetcdfHome)/include $(ModuleFlag)$(CodeDir) $(ModuleFlag)$(IntegDir)


rtObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
          $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
          $(CodeDir)/$(multipleProcCode) $(CodeDir)/scatteringPhaseFunctions.o \
          $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o \
          $(CodeDir)/RandomNumbersForMC.o $(CodeDir)/monteCarloIllumination.o \
          $(CodeDir)/surfaceProperties.o $(IntegDir)/monteCarloRadiativeTransfer.o
 

Targets = i3rcStepCloud i3rcRadarCloud i3rcLandsatCloud planeParallel cubic stepcloud_2tau write_HG

all: $(Targets)

i3rcStepCloud: i3rcStepCloud.o  $(rtObjs)
i3rcStepCloud.o: i3rcStepCloud.f95 $(rtObjs)

i3rcRadarCloud: i3rcRadarCloud.o  $(rtObjs)
i3rcRadarCloud.o: i3rcRadarCloud.f95 $(rtObjs)


i3rcLandsatCloud: i3rcLandsatCloud.o  $(rtObjs)
i3rcLandsatCloud.o: i3rcLandsatCloud.f95 $(rtObjs)

planeParallel: planeParallel.o $(rtObjs)
planeParallel.o: planeParallel.f95 $(rtObjs)

cubic: cubic.o $(rtObjs)
cubic.o: cubic.f95 $(rtObjs)

stepcloud_2tau: stepcloud_2tau.o $(rtObjs)
stepcloud_2tau.o: stepcloud_2tau.f95 $(rtObjs)

write_HG: write_HG.o $(rtObjs)
write_HG.o: write_HG.f95 $(rtObjs)

test: all
	./i3rcStepCloud
	./i3rcRadarCloud
	./i3rcLandsatCloud

clean::
	/bin/rm *.lst *.o $(Targets)
        
