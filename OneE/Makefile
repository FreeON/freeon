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
Libs =  $(F90_LIBRARIES)
#
ObjsS = MakeS.o SBlok.o 
ObjsT = MakeT.o TBlok.o   
ObjsR = MakeRho.o RhoBlok.o 
ObjGS = GradSBlok.o GradS.o
ObjGT = GradTBlok.o GradT.o
#
all:	MakeS MakeT MakeRho GradS GradT
#
MakeS:	 $(ObjsS)
	 $(LD) $(LDFLAGS) MakeS $(ObjsS) $(LibOpts) $(Libs) 
	 mv MakeS $(MONDO_HOME)/Exec/MakeS

MakeT:	 $(ObjsT)
	 $(LD) $(LDFLAGS) MakeT $(ObjsT) $(LibOpts) $(Libs) 
	 mv MakeT $(MONDO_HOME)/Exec/MakeT

MakeRho: $(ObjsR)
	 $(LD) $(LDFLAGS) MakeRho $(ObjsR) $(LibOpts) $(Libs) 
	 mv MakeRho $(MONDO_HOME)/Exec/MakeRho

GradS:   $(ObjGS)
	$(LD) $(LDFLAGS) GradS $(ObjGS) $(LibOpts) $(Libs)
	mv GradS $(MONDO_HOME)/Exec/GradS

GradT:   $(ObjGT)
	$(LD) $(LDFLAGS) GradT $(ObjGT) $(LibOpts) $(Libs) 
	mv GradT $(MONDO_HOME)/Exec/GradT

clean:
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Exec/MakeS
	rm -f $(MONDO_HOME)/Exec/MakeT
	rm -f $(MONDO_HOME)/Exec/MakeRho
	rm -f $(MONDO_HOME)/Exec/GradS
	rm -f $(MONDO_HOME)/Exec/GradT
	ln -s /dev/null core
	ls -l

MakeS.o   : SBlok.o     $(MONDO_HOME)/Libs/libMondoMods.a
MakeT.o   : TBlok.o     $(MONDO_HOME)/Libs/libMondoMods.a
MakeRho.o : RhoBlok.o   $(MONDO_HOME)/Libs/libMondoMods.a

