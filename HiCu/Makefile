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
CPPMISC = -DRULE5 -DEXPLICIT_SOURCE
#
EXTRA_INCLUDES= -I$(MONDO_HOME)/MondoMods/MMA/Functions \
                -I$(MONDO_HOME)/MondoMods/MMA/Functionals
#
HObjs  = Grid.o RhoTree.o CubeTree.o KxcGen.o HiCu.o 
GObjs  = Grid.o RhoTree.o CubeTree.o GradXCBlok.o GradXC.o
NObjs  = Grid.o RhoTree.o CubeTree.o GradXCBlok_nonlinear.o GradXC_nonlinear.o
#
VObjs  = VRhoMD.o 
#
all:	VR HC GH NH
#
clean:	CHC CGH CNH
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	ls -lS
#
purge:	CHC CVR CGH CNH
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	ls -l
#
HC:	$(HObjs)
	$(LD) $(LDFLAGS) HiCu $(HObjs) $(LibOpts) $(F90_LIBRARIES)
	mv HiCu $(MONDO_HOME)/Exec
#
CHC:	
	rm -f $(MONDO_HOME)/Exec/HiCu
#
VR:	$(VObjs)
	$(LD) $(LDFLAGS) VRho $(VObjs) $(LibOpts) $(F90_LIBRARIES)
	./VRho
	rm -rf VRho VRho.o VRho.f90
#
GH:     $(GObjs)
	$(LD) $(LDFLAGS) GradXC $(GObjs) $(LibOpts) $(F90_LIBRARIES)
	mv GradXC $(MONDO_HOME)/Exec
#
CGH:
	rm -f $(MONDO_HOME)/Exec/GradXC
#
NH:     $(NObjs)
	$(LD) $(LDFLAGS) GradXCNL $(NObjs) $(LibOpts) $(F90_LIBRARIES)
	mv GradXCNL $(MONDO_HOME)/Exec
#
CNH:
	rm -f $(MONDO_HOME)/Exec/GradXCNL
#
CVR:	
	rm -rf Explicit*.Inc
#
CubeTree.o:	Grid.o                      \
		GeneralLeafPopulation.Inc   \
		GeneralLeafContribution.Inc   $(MONDO_HOME)/Libs/libMondoMods.a
Grid.o:		MMA/CubeRules/Rule5.Inc \
                MMA/CubeRules/Rule7.Inc        $(MONDO_HOME)/Libs/libMondoMods.a
HiCu.o:		KxcGen.o CubeTree.o RhoTree.o  $(MONDO_HOME)/Libs/libMondoMods.a
RhoTree.o:	                               $(MONDO_HOME)/Libs/libMondoMods.a
VRhoMD.o:	                               $(MONDO_HOME)/Libs/libMondoMods.a
