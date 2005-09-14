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
EXTRA_INCLUDES=-I$(MONDO_HOME)/OneE/MMA
#
Objs1 = Int1E.o     Make1X.o   Make1E.o
ObjsM     = MBlok.o     MakeM.o  
ObjsU     = UBlok.o     MakeU.o  
ObjsS     = SBlok.o     MakeS.o 
ObjsT     = TBlok.o     MakeT.o 
ObjsR     = RhoTools.o  RhoBlok.o  MakeRho.o 
ObjSF = BlokTrWdS.o SForce.o 
ObjTF = BlokTrPdT.o TForce.o
ObjMF = BlokTrPdM.o MForce.o
ObjTH = BlokTrPd2T.o THess.o
ObjUF = UBlok.o UForce.o

ReleaseFiles = SBlok.x MakeS.x TBlok.x MakeT.x  UBlok.x MakeU.x \
               RhoTools.x   RhoBlok.x MakeRho.x \
               BlokTrWdS.x SForce.x BlokTrPdT.x  TForce.x       \
               UForce.x MakeM.x MBlok.x BlokTrPd2T.x THess.x
#
all: 	 Make1E MakeM MakeU MakeS MakeT MakeRho SForce TForce UForce MForce THess
#
Make1E:	$(Objs1)
	$(LD) $(LDFLAGS) Make1E $(Objs1) $(F90_LIBRARIES) $(LAPACK)
	mv Make1E $(MONDO_EXEC)/Make1E

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
MakeRho: $(ObjsR)
	 $(LD) $(LDFLAGS) MakeRho $(ObjsR) $(F90_LIBRARIES) $(LAPACK)
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
MForce:	$(ObjMF)
	$(LD) $(LDFLAGS) MForce $(ObjMF) $(F90_LIBRARIES) $(LAPACK)
	mv MForce $(MONDO_EXEC)/MForce
#
UForce:	$(ObjUF)
	$(LD) $(LDFLAGS) UForce $(ObjUF) $(F90_LIBRARIES)
	mv UForce $(MONDO_EXEC)/UForce
#
THess:	$(ObjTH)
	$(LD) $(LDFLAGS) THess $(ObjTH) $(F90_LIBRARIES) $(LAPACK)
	mv THess $(MONDO_EXEC)/THess
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
	rm -f $(MONDO_EXEC)/Make1E
	rm -f $(MONDO_EXEC)/MakeS
	rm -f $(MONDO_EXEC)/MakeT
	rm -f $(MONDO_EXEC)/MakeRho
	rm -f $(MONDO_EXEC)/SForce
	rm -f $(MONDO_EXEC)/TForce
	rm -f $(MONDO_EXEC)/THess
	ln -s /dev/null core
	ls -l
#
MakeU.o  : UBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeS.o  : SBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeT.o  : TBlok.o                $(MONDO_LIB)/libMondoMods.a
MakeRho.o: RhoBlok.o  RhoTools.o  $(MONDO_LIB)/libMondoMods.a
SForce.o : BlokTrWdS.o            $(MONDO_LIB)/libMondoMods.a
TForce.o : BlokTrPdT.o            $(MONDO_LIB)/libMondoMods.a
MForce.o : BlokTrPdM.o            $(MONDO_LIB)/libMondoMods.a
THess.o  : BlokTrPd2T.o           $(MONDO_LIB)/libMondoMods.a
