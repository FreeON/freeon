#================================================================================
#
#     Matt Challacombe
#     Los Alamos National Laboratory
#     Copyright 2000, The University of California
#
#================================================================================
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions \
                -I$(MONDO_HOME)/MondoMods/MMA/Functionals
#
MModObjs =    GlobalScalars.o                                         \
              GlobalCharacters.o                                      \
              DerivedTypes.o                                          \
              GlobalObjects.o                                         \
              BasisSetParameters.o                                    \
              MPIInclude.o                                            \
              ProcessControl.o Trap.o                                 \
              MemMan.o                                                \
              MondoMPI.o                                              \
              Clock.o Seconds.o                                       \
              ParsingConstants.o                                      \
              Parse.o                                                 \
              Indexing.o                                              \
              SetEq.o SetEq77.o                                       \
              InOut.o                                                 \
              Order.o Order77.o Random.o                              \
              PrettyPrint.o                                           \
              McMurch.o                                               \
              LinAlg.o LinAlg77.o MatFunk.o SMWrap.o                  \
              SpecFun.o                                               \
              Thresholding.o                                          \
              Functionals.o                                           \
              PBC.o AtomPairs.o                                       \
              BoundingBox.o BraKetBloks.o                             \
              Macros.o                                                \
              HDF5Wrap.o
#
all:	Mods DGEMM SM 
#
small:	Mods
#
clean:	CMods CDGEMM CSM
#
purge:	clean
	rm -f $(REMOVEALL)
	rm -f $(MONDO_HOME)/Libs/libMondoMods.a
	ls -l $(MONDO_HOME)/Libs
#
Mods:	$(MModObjs)
	$(AR) $(ARFLAGS) $(MONDO_HOME)/Libs/libMondoMods.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_HOME)/Libs/libMondoMods.a	
#
CMods:
	rm -f $(REMOVEMISC)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
#
DGEMM:	
	$(MAKE) -C $(MONDO_DGEMMS)
#
CDGEMM:	
	$(MAKE) -i -C $(MONDO_DGEMMS) clean
	rm -f $(MONDO_HOME)/MondoMods/DGEMMs
#
SMTestObjs=  SMTest.o 
#
SM:	$(SMTestObjs)
	$(LD) $(LDFLAGS) SMTest $(SMTestObjs) $(LibOpts) $(F90_LIBRARIES) 
	./SMTest>SMPerf.log
	rm -f $(SMTestObjs) ./SMTest
	cat SMPerf.log
#
CSM:	
	rm -f $(SMTestObjs)
	rm -f ./SMTest
#---------------------------------------------------------------------------------------------
#       DEPENDENCIES
#
AtomPairs.o:		Thresholding.o PBC.o BoundingBox.o \
			GlobalObjects.o GlobalCharacters.o \
			GlobalScalars.o DerivedTypes.o 
#
BasisSetParameters.o: 	GlobalCharacters.o 
#
BoundingBox.o: 		DerivedTypes.o 
#
BraKetBloks.o: 		McMurch.o PrettyPrint.o GlobalScalars.o \
			DerivedTypes.o 
#
Clock.o: 		MondoMPI.o DerivedTypes.o 
#
DerivedTypes.o: 	GlobalCharacters.o GlobalScalars.o 
#
Functionals.o: 		Parse.o ProcessControl.o DerivedTypes.o \
                        MMA/Functionals/B88x.Inc  \
                        MMA/Functionals/LYPc.Inc  \
                        MMA/Functionals/PBEx.Inc  \
                        MMA/Functionals/PW91x.Inc \
                        MMA/Functionals/VWN3c.Inc \
                        MMA/Functionals/XAx.Inc   \
                        MMA/Functionals/PBEc.Inc  \
                        MMA/Functionals/PW91c.Inc \
                        MMA/Functionals/SDx.Inc   \
                        MMA/Functionals/VWN5c.Inc                                        
#
GlobalCharacters.o: 	GlobalScalars.o 
#
GlobalObjects.o: 	DerivedTypes.o 
#
Indexing.o: 		DerivedTypes.o 
#
InOut.o:   		SetEq.o Parse.o Indexing.o MemMan.o \
                        ProcessControl.o GlobalCharacters.o \
                        GlobalScalars.o DerivedTypes.o 
#
LinAlg.o: 		Thresholding.o MondoMPI.o Clock.o Order.o \
			InOut.o SetEq.o MemMan.o PrettyPrint.o \
			ProcessControl.o GlobalObjects.o GlobalScalars.o \
			DerivedTypes.o 
#
Macros.o: 		AtomPairs.o Functionals.o InOut.o Clock.o \
			Parse.o PrettyPrint.o ParsingConstants.o \
			GlobalCharacters.o GlobalObjects.o GlobalScalars.o \
			DerivedTypes.o 
#
MatFunk.o: 		PrettyPrint.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o 
#
McMurch.o:	 	MemMan.o GlobalScalars.o DerivedTypes.o 
#
MemMan.o: 		ProcessControl.o GlobalObjects.o GlobalCharacters.o\
			GlobalScalars.o DerivedTypes.o 
#
MondoMPI.o: 		MemMan.o MPIInclude.o ProcessControl.o \
			GlobalObjects.o GlobalScalars.o DerivedTypes.o 
#
Order.o: 		Order77.o MemMan.o ProcessControl.o GlobalObjects.o \
			GlobalScalars.o GlobalCharacters.o DerivedTypes.o 
#
Parse.o: 		MemMan.o ProcessControl.o ParsingConstants.o \
			GlobalCharacters.o GlobalScalars.o DerivedTypes.o 
#
ParsingConstants.o: 	GlobalScalars.o 
#
PBC.o: 			Thresholding.o PrettyPrint.o GlobalObjects.o \
			GlobalCharacters.o GlobalScalars.o DerivedTypes.o 
#
PrettyPrint.o: 		Order.o InOut.o SetEq.o Parse.o Clock.o \
			BasisSetParameters.o GlobalObjects.o GlobalScalars.o\
			GlobalCharacters.o DerivedTypes.o 
#
ProcessControl.o: 	MPIInclude.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o 
#
SetEq.o: 		MondoMPI.o MemMan.o GlobalObjects.o \
			GlobalScalars.o DerivedTypes.o 
#
SMTest.o: 		Order.o MemMan.o Parse.o PrettyPrint.o \
			InOut.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o 
#
SpecFun.o: 		MMA/Functions/Exp.Inc MMA/Functions/Erf.Inc \
			DerivedTypes.o 
#
Thresholding.o: 	GlobalObjects.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o InOut.o SpecFun.o
#------------------------------------------------------------------------------
