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

          

all:	$(Objects)
	$(AR) $(ARFLAGS) $(MONDO_LIB)/libVRR.a $(?:.o=.o)
	$(RANLIB) $(MONDO_LIB)/libVRR.a

clean:	
	rm -f $(REMOVEALL)
	rm -f \#*
	rm -f *~
	rm -f $(MONDO_HOME)/Libs/libVRR.a
