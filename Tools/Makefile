# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Tools/Makefile $
include ../Makefile
Modules += $(ModuleFlag)$(CodeDir)

mieObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
           $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
           $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/$(multipleProcCode)
           
optObjs = $(mieObjs) $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o

all: MakeMieTable 

# Object files
#
MakeMieTable.o: MakeMieTable.f95  $(mieObjs)







# Executables
#
MakeMieTable: MakeMieTable.o mieindsub.o RefractiveIndex-IceAndWater.o $(mieObjs)







clean::
	-/bin/rm *.lst *.mod *.o MakeMieTable 
