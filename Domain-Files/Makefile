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
          $(CodeDir)/RandomNumbersForMC.o $(CodeDir)/emissionAndBroadBandWeights.o \
          $(CodeDir)/monteCarloIllumination.o $(CodeDir)/surfaceProperties.o \
          $(IntegDir)/monteCarloRadiativeTransfer.o
 

Targets =  #homogBBDomain inhomogBBDomain

all: $(Targets)




#homogBBDomain: homogBBDomain.o $(rtObjs)
#homogBBDomain.o: homogBBDomain.f95 $(rtObjs)

#inhomogBBDomain: inhomogBBDomain.o $(rtObjs)
#inhomogBBDomain.o: inhomogBBDomain.f95 $(rtObjs)

test: all
#	./homogBBDomain
#	./inhomogBBDomain

clean::
	/bin/rm *.lst *.o *.mod *.pdb *.inst.f95 $(Targets)
        
