###################################################################################
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
####################################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions 

all: ConvertTables 

LD=$(F90)

Objects1 = Convert.o GDump.o EDump.o HDump.o ADump.o GetGammaTable.o \
           ConvertG0.o  ConvertG1.o  ConvertG2.o  ConvertG3.o ConvertG4.o \
           ConvertG5.o  ConvertG6.o  ConvertG7.o  ConvertG8.o ConvertG9.o \
           ConvertG10.o ConvertG11.o ConvertG12.o \
           ConvertE.o 

Objects2 = MakeGD.o GCMaker.o

ReleaseFiles= \
           Convert.x GDump.x EDump.x HDump.x ADump.x GetGammaTable.x \
           ConvertG0.x  ConvertG1.x  ConvertG2.x  ConvertG3.x ConvertG4.x \
           ConvertG5.x  ConvertG6.x  ConvertG7.x  ConvertG8.x ConvertG9.x \
           ConvertG10.x ConvertG11.x ConvertG12.x \
           ConvertE.x MakeGD.x GCMaker.x

ConvertTables:	$(Objects1)
	$(LD) $(LDFLAGS) ConvertTables $(Objects1) $(LibOpts)
	./ConvertTables

MakeGDrivers:	$(Objects2)
	$(LD) $(LDFLAGS) MakeGDrivers $(Objects2) $(LibOpts)
	./MakeGDrivers

clean:
	rm -f $(REMOVEALL)
	rm -f ConvertTables MakeGDrivers
	rm -f $(MONDO_HOME)/MondoMods/MMA/Functions/*.bin
	ls

release:clean $(ReleaseFiles)
