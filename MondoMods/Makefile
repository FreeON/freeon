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
              Thresholding.o                                          \
              LinAlg.o LinAlg77.o MatFunk.o SMWrap.o                  \
              SpecFun.o                                               \
              Functionals.o                                           \
              PBC.o AtomPairs.o                                       \
              BoundingBox.o BraKetBloks.o                             \
              Macros.o                                                \
              HDF5Wrap.o
#
#
#
all:	 Mods  DGEMM  SM
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
	make -C $(MONDO_DGEMMS)
CDGEMM:	
	make -i -C $(MONDO_DGEMMS) clean
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
#
# nag_depend95 *.f90 > make.tmp
#
AtomPairs.o : thresholding.mod cellsets.mod boundingbox.mod globalobjects.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
atompairs.mod : AtomPairs.o 
BasisSetParameters.o : globalcharacters.mod 
basissetparameters.mod : BasisSetParameters.o 
BoundingBox.o : derivedtypes.mod 
boundingbox.mod : BoundingBox.o 
BraKetBloks.o : mcmurchie.mod prettyprint.mod globalscalars.mod derivedtypes.mod 
braketbloks.mod : BraKetBloks.o 
Clock.o : mondompi.mod derivedtypes.mod 
clock.mod : Clock.o 
DerivedTypes.o : globalcharacters.mod globalscalars.mod 
derivedtypes.mod : DerivedTypes.o 
Functionals.o : parse.mod processcontrol.mod derivedtypes.mod 
functionals.mod : Functionals.o 
GlobalCharacters.o : globalscalars.mod 
globalcharacters.mod : GlobalCharacters.o 
GlobalObjects.o : derivedtypes.mod 
globalobjects.mod : GlobalObjects.o 
globalscalars.mod : GlobalScalars.o 
Indexing.o : derivedtypes.mod 
indexing.mod : Indexing.o 
InOut.o :   setxyz.mod parse.mod indexing.mod memman.mod processcontrol.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
inout.mod : InOut.o 
LinAlg.o : thresholding.mod mondompi.mod clock.mod order.mod inout.mod setxyz.mod memman.mod prettyprint.mod processcontrol.mod globalobjects.mod globalscalars.mod derivedtypes.mod 
linalg.mod common_debug.mod : LinAlg.o 
Macros.o : atompairs.mod functionals.mod inout.mod clock.mod parse.mod prettyprint.mod parsingconstants.mod globalcharacters.mod globalobjects.mod globalscalars.mod derivedtypes.mod 
macros.mod : Macros.o 
MatFunk.o : prettyprint.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
matfunk.mod : MatFunk.o 
McMurch.o : memman.mod globalscalars.mod derivedtypes.mod 
mcmurchie.mod : McMurch.o 
MemMan.o : processcontrol.mod globalobjects.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
memman.mod : MemMan.o 
MondoMPI.o : memman.mod mpiinclude.mod processcontrol.mod globalobjects.mod globalscalars.mod derivedtypes.mod 
mondompi.mod : MondoMPI.o 
Order.o : memman.mod processcontrol.mod globalobjects.mod globalscalars.mod globalcharacters.mod derivedtypes.mod 
order.mod : Order.o 
Parse.o : memman.mod processcontrol.mod parsingconstants.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
parse.mod : Parse.o 
ParsingConstants.o : globalscalars.mod 
parsingconstants.mod : ParsingConstants.o 
PBC.o : thresholding.mod prettyprint.mod globalobjects.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
cellsets.mod : PBC.o 
PrettyPrint.o : order.mod inout.mod setxyz.mod parse.mod clock.mod basissetparameters.mod globalobjects.mod globalscalars.mod globalcharacters.mod derivedtypes.mod 
prettyprint.mod : PrettyPrint.o 
ProcessControl.o : mpiinclude.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
processcontrol.mod : ProcessControl.o 
SetEq.o : mondompi.mod memman.mod globalobjects.mod globalscalars.mod derivedtypes.mod 
setxyz.mod : SetEq.o 
SMTest.o : order.mod memman.mod parse.mod prettyprint.mod inout.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
SpecFun.o : MMA/Functions/Exp.Inc MMA/Functions/Erf.Inc derivedtypes.mod 
specfun.mod : SpecFun.o 
Thresholding.o : globalobjects.mod globalcharacters.mod globalscalars.mod derivedtypes.mod 
thresholding.mod : Thresholding.o 
