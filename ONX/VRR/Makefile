include $(MONDO_COMPILER)
include $(MONDO_HOME)/Includes/Suffixes
include $(MONDO_HOME)/Includes/RemoveAll

EXTRA_INCLUDES = -I$(MONDO_HOME)/ONX/NONX

Objects = Int1111.o Int2111.o Int1121.o Int1122.o Int2211.o Int2121.o

all:	$(Objects)
	$(AR) $(ARFLAGS) $(MONDO_HOME)/Libs/libVRR.a $(?:.F90=.o)
	$(RANLIB) $(MONDO_HOME)/Libs/libVRR.a

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Libs/libVRR.a
