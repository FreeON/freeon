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

EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE

all: ONX2

Objects = QuickSortDis.o PartDrv.o ONXParameters.o ONXInit.o ONXGet.o   \
	   ONXFillOut.o Stats.o ONXRng.o \
	   ONX2F77.o ONX2DataType.o ONXCtrSclg.o \
	   IntSSSS.o ONX2ComputK.o ONX2.o
ReleaseFiles= 
#
Libs =  $(F90_LIBRARIES) $(LAPACK) -lTwoE
#
ONX2:	$(Objects)
	$(LD) $(LDFLAGS) ONX2 $(Objects) $(Libs)
	mv ONX2 $(MONDO_EXEC)/ONX2
#
release:purge $(ReleaseFiles)
	rm -f HRR.tar
	rm -f Notes.txt
	rm -f README
#
purge:	clean
	rm -f $(REMOVESRCE)
#
clean:	
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/ONX
	rm -f $(Objects)
	ls
#
doc.txt:
	cat *.F90 | grep "^!H" | sed "s/^!H/  /" > doc.txt;




