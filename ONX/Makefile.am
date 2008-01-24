include ../make.inc

bin_PROGRAMS = ONX GONX

ONX_SOURCES = \
  QSDis.f \
  PartDrv.f90 \
  ONXParameters.f90 \
  ONXInit.f90 \
  ONXGet.f90 \
  ONXFillOut.f90 \
  Stats.f90 \
  ONXRng.f90 \
  ONX2F77.f \
  ONX2DataType.f90 \
  ONXCtrSclg.f90 \
  ONX2List.f90 \
  ONX2ComputK.f90 \
  ONX2.f90

GONX_SOURCES = \
  QSDis.f \
  PartDrv.f90 \
  ONXParameters.f90 \
  ONXGet.f90 \
  ONX2F77.f \
  ONX2DataType.f90 \
  ONXCtrSclg.f90 \
  GONX2ComputDK.f90 \
  GONX2.f90 \
  ONX2List.f90

ONX_DEPENDENCIES  = \
  $(MondoMods_DEPENDENCY_FILES) \
	../TwoE/GRD/libGRD.a \
	../TwoE/INT/libINT.a \
	../TwoE/VRR/libVRR.a \
	../TwoE/HRR/libHRR.a \
	../TwoE/MSC/libTwoEMisc.a \
	../TwoE/MSC/shellpairstruct.mod \
	../TwoE/MSC/vscratchb.mod \
	../TwoE/MSC/getatompairmod.mod

GONX_DEPENDENCIES = \
  $(MondoMods_DEPENDENCY_FILES) \
	../TwoE/GRD/libGRD.a \
	../TwoE/INT/libINT.a \
	../TwoE/VRR/libVRR.a \
	../TwoE/HRR/libHRR.a \
	../TwoE/MSC/libTwoEMisc.a \
	../TwoE/MSC/shellpairstruct.mod \
	../TwoE/MSC/vscratchb.mod \
	../TwoE/MSC/getatompairmod.mod

SOURCES_DEPCHECK = \
  $(ONX_SOURCES) \
  $(GONX_SOURCES)

EXTRA_DIST = \
  ERIMIC.Inc

CPPMISC = -DONX2_PARALLEL -DGTRESH

EXTRA_INCLUDES = -I../MondoMods -I../TwoE -I../TwoE/MSC -I../TwoE/INT -I../TwoE/GRD

LDADD = $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5) -L../MondoMods -lMondoMods \
  -L../TwoE/GRD -lGRD -L../TwoE/INT -lINT -L../TwoE/VRR -lVRR -L../TwoE/HRR -lHRR -L../TwoE/MSC -lTwoEMisc

ONXInit.o : \
  ONXParameters.o

ONXFillOut.o : \
  ONXParameters.o \
	ONXGet.o

ONXRng.o : \
  ONXParameters.o \
  Stats.o

ONXCtrSclg.o : \
  ONXParameters.o

ONX2List.o : \
  ONX2DataType.o \
	ONXParameters.o

ONX2ComputK.o : \
  ONX2DataType.o \
	ONXParameters.o

ONX2.o : \
  ONX2DataType.o \
	ONX2List.o \
	ONX2ComputK.o \
	ONXInit.o \
	ONXCtrSclg.o \
	ONXRng.o \
	ONXFillOut.o \
	ONXGet.o

GONX2ComputDK.o : \
  ONX2DataType.o \
	ONXParameters.o \
	ONX2List.o

GONX2.o : \
  ONX2DataType.o \
	ONXCtrSclg.o \
	ONX2List.o \
	GONX2ComputDK.o \
	ONXGet.o
