#! $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
#! $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Code/Makefile $

include ../Makefile

all_obs = characterUtils.o ErrorMessages.o RandomNumbersForMC.o numericUtilities.o \
      scatteringPhaseFunctions.o inversePhaseFunctions.o \
      surfaceProperties.o opticalProperties.o userInterface_Unix.o \
      emissionAndBroadBandWeights.o monteCarloIllumination.o ${multipleProcCode}

all: ${all_obs} 

ErrorMessages.o: ErrorMessages.f95

characterUtils.o: characterUtils.f95

numericUtilities.o: numericUtilities.f95

RandomNumbersForMC.o: RandomNumbersForMC.f95

multipleProcesses_nompi.o: multipleProcesses_nompi.f95

multipleProcesses_mpi.o: multipleProcesses_mpi.f95

userInterface_Unix.o: userInterface_Unix.f95 ErrorMessages.o ${multipleProcCode}

kDistribution.o: kDistribution.f95 ErrorMessages.o 

surfaceProperties.o: surfaceProperties.f95 ErrorMessages.o numericUtilities.o

scatteringPhaseFunctions.o: scatteringPhaseFunctions.f95 ErrorMessages.o numericUtilities.o characterUtils.o

opticalProperties.o: opticalProperties.f95 scatteringPhaseFunctions.o \
                     characterUtils.o ErrorMessages.o 
                     
inversePhaseFunctions.o: inversePhaseFunctions.f95 scatteringPhaseFunctions.o \
                         numericUtilities.o ErrorMessages.o

emissionAndBroadBandWeights.o: emissionAndBroadBandWeights.f95 ErrorMessages.o \
                               RandomNumbersForMC.o numericUtilities.o \
                               opticalProperties.o

monteCarloIllumination.o: monteCarloIllumination.f95 ErrorMessages.o RandomNumbersForMC.o \
                          numericUtilities.o emissionAndBroadBandWeights.o characterUtils.o

clean::
	-/bin/rm *.lst *.mod *.o *.pdb *.inst.f95

