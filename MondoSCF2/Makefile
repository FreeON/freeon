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
#    e.g. in the form <<These calculations have been performed using MondoSCF, 
#    a suite of programs for linear scaling electronic structure theory and
#    ab initio molecular dynamics>>, and given appropriate citation.  
##############################################################################
include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll
#
CPPMISC=#-DFULL_ON_FRONT_END_DEBUG
#
ObjectFiles = GetPID.o Spawn.o  ControlStructures.o NEB.o Numerics.o ParseCommands.o  \
              OptionKeys.o ConflictDetect.o ParseOptions.o DynamicsKeys.o             \
              ParseDynamics.o GeometryKeys.o ParseProperties.o InCoords.o             \
              ParseGeometries.o PeriodicKeys.o ParsePeriodic.o ParseBasis.o           \
              ParallelKeys.o ParseParallel.o Massage.o GeomOptKeys.o ParseGeomOpt.o   \
              ParseExtraCoords.o ParseInput.o PunchHDF.o Overlay.o PrintParsed.o      \
              SCFKeys.o SCFs.o MLP.o Response.o Hessian.o GDIIS.o Optimizer.o         \
	      MDynamics.o ZippyQuote.o MondoSCF.o
#
ReleaseFiles= GetPID.x Spawn.x  ControlStructures.x NEB.x Numerics.x ParseCommands.x  \
              OptionKeys.x ConflictDetect.x ParseOptions.x DynamicsKeys.x             \
              ParseDynamics.x GeometryKeys.x ParseProperties.x InCoords.x             \
              ParseGeometries.x PeriodicKeys.x ParsePeriodic.x ParseBasis.x           \
              ParallelKeys.x ParseParallel.x Massage.x GeomOptKeys.x ParseGeomOpt.x   \
              ParseExtraCoords.x ParseInput.x PunchHDF.x Overlay.x PrintParsed.x      \
              SCFKeys.x SCFs.x MLP.x Response.x Hessian.x GDIIS.x Optimizer.x         \
	      MDynamics.x ZippyQuote.x MondoSCF.x
#
all:	$(ObjectFiles) 
	$(LD) $(LDFLAGS) MondoSCF $(ObjectFiles) $(F90_LIBRARIES) $(LAPACK)
	mv MondoSCF $(MONDO_EXEC)/MondoSCF
#
release:$(ReleaseFiles) 
#
clean:
	rm -f $(REMOVEMISC) $(REMOVEMODS)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_EXEC)/MondoSCF
	ln -s /dev/null core
	ls -l

purge:	clean
	rm -f $(REMOVESRCE)
