##############################################################################
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
##############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
EXTRA_INCLUDES=-I$(MONDO_HOME)/OneE/MMA/
#
ObjsM     = MBlok.o     MakeM.o  
ObjsU     = UBlok.o     MakeU.o  
ObjsS     = SBlok.o     MakeS.o 
ObjsT     = TBlok.o     MakeT.o 
ObjsR     = RhoTools.o  RhoBlok.o  MakeRho.o 
ObjsParaR = RhoTools.o  RhoBlok.o  ParaMakeRho.o
ObjSF = BlokTrWdS.o SForce.o 
ObjTF = BlokTrPdT.o TForce.o

ReleaseFiles = SBlok.x MakeS.x TBlok.x MakeT.x  UBlok.x MakeU.x \
               RhoTools.x   RhoBlok.x MakeRho.x ParaMakeRho.x \
               BlokTrWdS.x SForce.x BlokTrPdT.x \
               TForce.x MakeM.x MBlok.x
#
all: 	 MakeM MakeU MakeS MakeT MakeRho ParaMakeRho SForce TForce UForce 
#
MakeM:	 $(ObjsM)
	 $(LD) $(LDFLAGS) MakeM $(ObjsM) $(F90_LIBRARIES)
	 mv MakeM $(MONDO_EXEC)/MakeM

MakeU:	 $(ObjsU)
	 $(LD) $(LDFLAGS) MakeU $(ObjsU) $(F90_LIBRARIES)
	 mv MakeU $(MONDO_EXEC)/MakeU
#
MakeS:	 $(ObjsS)
	 $(LD) $(LDFLAGS) MakeS $(ObjsS) $(F90_LIBRARIES)
	 mv MakeS $(MONDO_EXEC)/MakeS
#
MakeT:	 $(ObjsT)
	 $(LD) $(LDFLAGS) MakeT $(ObjsT) $(F90_LIBRARIES)
	 mv MakeT $(MONDO_EXEC)/MakeT
#
ParaMakeRho: $(ObjsParaR)
	 $(LD) $(LDFLAGS) ParaMakeRho $(ObjsParaR) $(F90_LIBRARIES)
	 mv ParaMakeRho $(MONDO_EXEC)/ParaMakeRho
#
MakeRho: $(ObjsR)
	 $(LD) $(LDFLAGS) MakeRho $(ObjsR) $(F90_LIBRARIES)
	 mv MakeRho $(MONDO_EXEC)/MakeRho
#
SForce: $(ObjSF)
	$(LD) $(LDFLAGS) SForce $(ObjSF) $(F90_LIBRARIES)
	mv SForce $(MONDO_EXEC)/SForce
#
TForce:	$(ObjTF)
	$(LD) $(LDFLAGS) TForce $(ObjTF) $(F90_LIBRARIES)
	mv TForce $(MONDO_EXEC)/TForce
#
UForce:	UForce.o
	$(LD) $(LDFLAGS) UForce UForce.o UBlok.o $(F90_LIBRARIES)
	mv UForce $(MONDO_EXEC)/UForce
#
release:clean $(ReleaseFiles) 
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/MakeS
	rm -f $(MONDO_EXEC)/MakeT
	rm -f $(MONDO_EXEC)/MakeRho
	rm -f $(MONDO_EXEC)/SForce
	rm -f $(MONDO_EXEC)/TForce
	rm -f $(MONDO_EXEC)/ParaMakeRho
	ln -s /dev/null core
	ls -l
#
MakeU.o  : UBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeS.o  : SBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeT.o  : TBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeRho.o: RhoBlok.o  RhoTools.o  $(MONDO_LIB)/libMondoMods.a
SForce.o : BlokTrWdS.o            $(MONDO_LIB)/libMondoMods.a
TForce.o : BlokTrPdT.o            $(MONDO_LIB)/libMondoMods.a
