#------------------------------------------------------------------------------
#--  This code is part of the MondoSCF suite of programs for linear scaling 
#    electronic structure theory and ab initio molecular dynamics.
#
#--  Copyright (c) 2001, the Regents of the University of California.  
#    This SOFTWARE has been authored by an employee or employees of the 
#    University of California, operator of the Los Alamos National Laboratory 
#    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
#    The U.S. Government has rights to use, reproduce, and distribute this 
#    SOFTWARE.  The public may copy, distribute, prepare derivative works 
#    and publicly display this SOFTWARE without charge, provided that this 
#    Notice and any statement of authorship are reproduced on all copies.  
#    Neither the Government nor the University makes any warranty, express 
#    or implied, or assumes any liability or responsibility for the use of 
#    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
#    such modified SOFTWARE should be clearly marked, so as not to confuse 
#    it with the version available from LANL.  The return of derivative works
#    to the primary author for integration and general release is encouraged. 
#    The first publication realized with the use of MondoSCF shall be
#    considered a joint work.  Publication of the results will appear
#    under the joint authorship of the researchers nominated by their
#    respective institutions. In future publications of work performed
#    with MondoSCF, the use of the software shall be properly acknowledged,
#    e.g. in the form "These calculations have been performed using MondoSCF, 
#    a suite of programs for linear scaling electronic structure theory and
#    ab initio molecular dynamics", and given appropriate citation.  
#------------------------------------------------------------------------------
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

CPPMISC = -DOLD_CHR_VECT #-Ddebug_interface -Ddebug_all
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/LookUpTables \
                -I$(MONDO_HOME)/MondoMods/MMA/Functionals  \
                -I$(MONDO_HOME)/MondoMods/MMA/PAC
#
MModObjs = GlobalScalars.o GlobalCharacters.o DerivedTypes.o    \
           GlobalObjects.o BasisSetParameters.o  MPIInclude.o   \
           ProcessControl.o Trap.o MemMan.o Clock.o  Seconds.o TimerWrapper.o \
           ParsingConstants.o Parse.o MondoMPI.o                \
           Indexing.o SetEq.o SetEq77.o InOut.o HDF5Wrap.o      \
           Anneal.o Order.o Order77.o Random.o PrettyPrint.o McMurch.o   \
           LinAlg.o LinAlg77.o MatFunk.o AInv.o SMWrap.o SpecFun.o     \
           D1Mach.o DCSEvl.o DErf.o DErfc.o DMach.o InitDS.o    \
           Mechanics.o Thresholding.o Functionals.o PBC.o AtomPairs.o       \
           BoundingBox.o BraBloks.o Macros.o FastMat.o \
           CholFact.o Slatec.o
#
MModRel  = GlobalScalars.x GlobalCharacters.x DerivedTypes.x    \
           GlobalObjects.x BasisSetParameters.x MPIInclude.x    \
           ProcessControl.x Trap.x MemMan.x Clock.x  Seconds.x  TimerWrapper.o\
           ParsingConstants.x Parse.x  MondoMPI.x               \
           Indexing.x SetEq.x SetEq77.x InOut.x HDF5Wrap.x      \
           Anneal.x Order.x Order77.x Random.x PrettyPrint.x McMurch.x   \
           LinAlg.x LinAlg77.x MatFunk.x AInv.x SMWrap.x SpecFun.x     \
           D1Mach.x DCSEvl.x DErf.x DErfc.x DMach.x InitDS.x    \
           Mechanics.x Thresholding.x Functionals.x PBC.x AtomPairs.x       \
           BoundingBox.x BraBloks.x Macros.x FastMat.x EllToInclude.x     \
           FTest.x SMTest.x CholFact.x Slatec.x
#
small:	Mods
#
all:	Mods ToEll DGEMM 
#
test:	SM F
#
clean:	CMods CDGEMM CSM CF
#
purge:	clean
	rm -f $(MONDO_LIB)/libMondoMods.a
	rm -f $(REMOVESRCE)
#
release:$(MModRel)
	rm -rf MMA/Functions
#
Mods:	$(MModObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libMondoMods.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libMondoMods.a	
#
CMods:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
#-------------------------------------------------------------
ToEll:	EllToInclude.o
	$(LD) $(LDFLAGS) EllInc EllToInclude.o $(LibOpts) $(F90_LIBRARIES) 
	./EllInc
	rm -f EllToInclude.o ./EllInc
#-------------------------------------------------------------
DGEMM:	
	$(MAKE) -C $(MONDO_DGEMMS)
CDGEMM:	
	$(MAKE) -i -C $(MONDO_DGEMMS) clean
	rm -f $(MONDO_HOME)/MondoMods/DGEMMs
#-------------------------------------------------------------
SMTestObjs=  SMTest.o 
SM:	$(SMTestObjs)
	$(LD) $(LDFLAGS) SMTest $(SMTestObjs)  $(F90_LIBRARIES) 
	./SMTest>Mods.log
	rm -f $(SMTestObjs) ./SMTest
	cat Mods.log
CSM:	
	rm -f $(SMTestObjs)
	rm -f ./SMTest
#-------------------------------------------------------------
FTestObjs=  FTest.o 
F:	$(FTestObjs)
	$(LD) $(LDFLAGS) FTest $(FTestObjs) $(F90_LIBRARIES) 
	./FTest>>Mods.log
	rm -f $(FTestObjs) ./FTest
	cat Mods.log
CF:	
	rm -f $(FTestObjs)
	rm -f ./FTest
#---------------------------------------------------------------------------------------------
#       DEPENDENCIES
#
AInv.o: 		PrettyPrint.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o MemMan.o LinAlg.o MatFunk.o 
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
			DerivedTypes.o Mechanics.o

FastMat.o:		GlobalScalars.o DerivedTypes.o GlobalObjects.o \
                        MemMan.o Order.o MondoMPI.o Thresholding.o
#
MatFunk.o: 		PrettyPrint.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o
#
McMurch.o:	 	MemMan.o GlobalScalars.o DerivedTypes.o 
#
MemMan.o: 		ProcessControl.o GlobalObjects.o GlobalCharacters.o\
			GlobalScalars.o DerivedTypes.o 
#
MondoMPI.o: 		MemMan.o MPIInclude.o ProcessControl.o Parse.o \
			GlobalObjects.o GlobalScalars.o DerivedTypes.o 
#
Anneal.o:               PrettyPrint.o MemMan.o
#
Order.o: 		PrettyPrint.o Order77.o MemMan.o ProcessControl.o GlobalObjects.o \
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
PrettyPrint.o: 		InOut.o SetEq.o Parse.o Clock.o \
			BasisSetParameters.o GlobalObjects.o GlobalScalars.o\
			GlobalCharacters.o DerivedTypes.o 
#
ProcessControl.o: 	MPIInclude.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o 
#
SetEq.o: 		MondoMPI.o MemMan.o GlobalObjects.o        \
			GlobalScalars.o DerivedTypes.o 
#
SMTest.o: 		Order.o MemMan.o Parse.o PrettyPrint.o     \
			InOut.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o 
#
SpecFun.o: 		MMA/LookUpTables/Exp.Inc \
			MMA/LookUpTables/Erf.Inc \
			MMA/LookUpTables/F0.Inc  \
			MMA/LookUpTables/F1.Inc  \
			MMA/LookUpTables/F2.Inc  \
			MMA/LookUpTables/F3.Inc  \
			MMA/LookUpTables/F4.Inc  \
			MMA/LookUpTables/F5.Inc  \
			MMA/LookUpTables/F7.Inc  \
			MMA/LookUpTables/F8.Inc  \
			MMA/LookUpTables/F9.Inc  \
			MMA/LookUpTables/F10.Inc \
			MMA/LookUpTables/F11.Inc \
			DerivedTypes.o 
#
Thresholding.o: 	GlobalObjects.o GlobalCharacters.o GlobalScalars.o \
			DerivedTypes.o InOut.o SpecFun.o Mechanics.o
#
CholFact.o:             DerivedTypes.o GlobalScalars.o GlobalObjects.o \
                        GlobalCharacters.o InOut.o MemMan.o  \
                        PrettyPrint.o ParsingConstants.o  
InCoords.o:             DerivedTypes.o GlobalScalars.o GlobalObjects.o InOut.o \
	                MemMan.o SetEq.o ProcessControl.o PrettyPrint.o ParsingConstants.o 
#------------------------------------------------------------------------------
