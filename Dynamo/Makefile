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
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/LookUpTables \
                -I$(MONDO_HOME)/MondoMods/MMA/Functionals  \
                -I$(MONDO_HOME)/MondoMods/MMA/PAC
#
DynamoObjs = dynamo.o definitions.o constants.o files.o io_units.o status.o \
        atoms.o sequence.o symmetry.o elements.o parsing.o \
        string.o pdb_io.o mm_system.o connectivity.o sort.o mm_file_data.o\
        mm_file_io.o mm_terms.o coordinate_io.o constraint.o  \
        linear_algebra.o 
#
DynamoRel = dynamo.x definitions.x constants.x files.x io_units.x status.x \
        atoms.x sequence.x symmetry.x  elements.x parsing.x\
        string.o pdb_io.x mm_system.x connectivity.x sort.x mm_file_data.x\
        mm_file_io.x mm_terms.x coordinate_io.x constraint.x  \
        linear_algebra.x 
#
small:	Dynamo
#
all:	Dynamo
#
clean:	CDynamo               
#
purge:	clean
	rm -f $(REMOVESRCE)
#
release:$(DynamoRel)
	rm -rf MMA/Functions
#
Dynamo:	$(DynamoObjs)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libMondoMods.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_LIB)/libMondoMods.a	
#
CDynamo:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	ln -s /dev/null core
	ls -l
#---------------------------------------------------------------------------------------------
#       DEPENDENCIES
#
dynamo.o: definitions.o constants.o files.o io_units.o status.o \
        atoms.o sequence.o symmetry.o elements.o parsing.o \
        string.o pdb_io.o mm_system.o connectivity.o sort.o mm_file_data.o\
        mm_file_io.o mm_terms.o coordinate_io.o constraint.o  \
        linear_algebra.o 
constants.o: definitions.o
elements.o: definitions.o
files.o: io_units.o definitions.o
io_units.o: definitions.o
status.o: io_units.o
sort.o: definitions.o
symmetry.o: definitions.o io_units.o parsing.o status.o
connectivity.o: definitions.o elements.o io_units.o sort.o
mm_file_data.o: definitions.o io_units.o atoms.o sequence.o
parsing.o: definitions.o io_units.o status.o string.o
atoms.o: constants.o definitions.o elements.o io_units.o
sequence.o: definitions.o files.o io_units.o status.o
pdb_io.o: constants.o definitions.o files.o io_units.o status.o atoms.o sequence.o symmetry.o
mm_system.o: definitions.o elements.o files.o io_units.o parsing.o status.o atoms.o connectivity.o mm_file_data.o mm_file_io.o mm_terms.o sequence.o symmetry.o
mm_file_io.o: constants.o definitions.o files.o io_units.o status.o mm_file_data.o sequence.o
mm_terms.o: definitions.o io_units.o atoms.o
coordinate_io.o: constants.o definitions.o  elements.o files.o io_units.o parsing.o status.o atoms.o sequence.o symmetry.o 
constraint.o: constants.o definitions.o files.o io_units.o linear_algebra.o status.o atoms.o
linear_algebra.o: definitions.o
#------------------------------------------------------------------------------
