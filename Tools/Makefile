# $Revision: 6 $, $Date: 2009-03-10 20:13:07 +0000 (Tue, 10 Mar 2009) $
# $URL: http://i3rc-monte-carlo-model.googlecode.com/svn/trunk/Tools/Makefile $
include ../Makefile
Modules += $(ModuleFlag)$(CodeDir) $(ModuleFlag)$(NetcdfHome)/include

mieObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
           $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
           $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/$(multipleProcCode)
           
optObjs = $(mieObjs) $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o

rayObjs = $(optObjs) $(CodeDir)/characterUtils.o $(CodeDir)/ErrorMessages.o \
	  $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/numericUtilities.o

Targets = MakeMieTable MieSSPTableCreate PhysicalPropertiesToDomainNew molecScattTrans

all: $(Targets)

# Object files
#
MakeMieTable.o: MakeMieTable.f95  $(mieObjs)

MieSSPTableCreate.o: MieSSPTableCreate.f95 $(mieObjs)

PhysicalPropertiesToDomainNew.o: PhysicalPropertiesToDomainNew.f95 $(optObjs)

molecScattTrans.o: molecScattTrans.f95 $(rayObjs)



# Executables
#
MakeMieTable: MakeMieTable.o mieindsub.o RefractiveIndex-IceAndWater.o $(mieObjs)

MieSSPTableCreate: MieSSPTableCreate.o mieindsub.o RefractiveIndex-IceAndWater.o $(mieObjs)

PhysicalPropertiesToDomainNew: PhysicalPropertiesToDomainNew.o $(optObjs)

molecScattTrans: molecScattTrans.o $(rayObjs)



test: all
	./MakeMieTable
	./MieSSPTableCreate
	./PhysicalPropertiesToDomainNew
	./molecScattTrans
	

clean::
	-/bin/rm *.lst *.mod *.o $(Targets) 
