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

EXTRA_INCLUDES = -I$(MONDO_HOME)/ONX/NONX

Objects = Int1111.o  Int2122.o  Int2262.o  Int6211.o  Int6622.o \
          Int1121.o  Int2161.o  Int2266.o  Int6221.o  Int6661.o \
          Int1122.o  Int2162.o  Int6111.o  Int6222.o  Int6662.o \
          Int1161.o  Int2166.o  Int6121.o  Int6261.o  Int6666.o \
          Int1162.o  Int2211.o  Int6122.o  Int6262.o \
          Int1166.o  Int2221.o  Int6161.o  Int6266.o \
          Int2111.o  Int2222.o  Int6162.o  Int6611.o \
          Int2121.o  Int2261.o  Int6166.o  Int6621.o \

ReleaseFiles= \
          Int1111.x  Int2122.x  Int2262.x  Int6211.x  Int6622.x \
          Int1121.x  Int2161.x  Int2266.x  Int6221.x  Int6661.x \
          Int1122.x  Int2162.x  Int6111.x  Int6222.x  Int6662.x \
          Int1161.x  Int2166.x  Int6121.x  Int6261.x  Int6666.x \
          Int1162.x  Int2211.x  Int6122.x  Int6262.x \
          Int1166.x  Int2221.x  Int6161.x  Int6266.x \
          Int2111.x  Int2222.x  Int6162.x  Int6611.x \
          Int2121.x  Int2261.x  Int6166.x  Int6621.x 


all:	$(Objects)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libVRR.a $(?:.o=.o)
	$(RANLIB) $(MONDO_LIB)/libVRR.a

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Libs/libVRR.a

release:clean $(ReleaseFiles)
	rm -f VRR-Maker.tar
	rm -f VRR.tar
