include $(top_srcdir)/make.inc

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
	$(top_srcdir)/TwoE/GRD/libGRD.la \
	$(top_srcdir)/TwoE/INT/libINT.la \
	$(top_srcdir)/TwoE/VRR/libVRR.la \
	$(top_srcdir)/TwoE/HRR/libHRR.la \
	$(top_srcdir)/TwoE/MSC/libTwoEMisc.la \
	$(top_srcdir)/TwoE/MSC/shellpairstruct.mod \
	$(top_srcdir)/TwoE/MSC/vscratchb.mod \
	$(top_srcdir)/TwoE/MSC/getatompairmod.mod

GONX_DEPENDENCIES = \
  $(MondoMods_DEPENDENCY_FILES) \
	$(top_srcdir)/TwoE/GRD/libGRD.la \
	$(top_srcdir)/TwoE/INT/libINT.la \
	$(top_srcdir)/TwoE/VRR/libVRR.la \
	$(top_srcdir)/TwoE/HRR/libHRR.la \
	$(top_srcdir)/TwoE/MSC/libTwoEMisc.la \
	$(top_srcdir)/TwoE/MSC/shellpairstruct.mod \
	$(top_srcdir)/TwoE/MSC/vscratchb.mod \
	$(top_srcdir)/TwoE/MSC/getatompairmod.mod

SOURCES_DEPCHECK = \
  $(ONX_SOURCES) \
  $(GONX_SOURCES)

EXTRA_DIST = \
  ERIMIC.Inc

CPPMISC = -DONX2_PARALLEL -DGTRESH

EXTRA_INCLUDES = -I../MondoMods -I../TwoE -I../TwoE/MSC -I../TwoE/INT -I../TwoE/GRD

LDADD = $(LDADD_MONDOMODS) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5) \
  ../TwoE/GRD/libGRD.la ../TwoE/INT/libINT.la ../TwoE/VRR/libVRR.la ../TwoE/HRR/libHRR.la ../TwoE/MSC/libTwoEMisc.la



GONX2.o: GONX2ComputDK.o ONX2DataType.o ONX2List.o ONXCtrSclg.o ONXGet.o \
	PartDrv.o
GONX2ComputDK.o: ONX2DataType.o ONX2List.o ONXParameters.o
ONX2.o: ONX2ComputK.f90-preprocessed.o ONX2DataType.o ONX2List.o ONXCtrSclg.o \
	ONXFillOut.o ONXGet.o ONXInit.o ONXParameters.o ONXRng.o PartDrv.o
ONX2ComputK.o: ONX2DataType.o ONXGet.o ONXParameters.o
ONX2ComputK.f90-memory.o: ONX2DataType.o ONXGet.o ONXParameters.o
ONX2ComputK.f90-preprocessed.o: ONX2DataType.o ONXGet.o ONXParameters.o
ONX2DataType.o:
ONX2List.o: ONX2DataType.o ONXParameters.o
ONXCtrSclg.o: ONXParameters.o
ONXFillOut.o: ONXGet.o ONXParameters.o
ONXGet.o:
ONXInit.o: ONXParameters.o
ONXParameters.o:
ONXRng.o: ONXParameters.o Stats.o
PartDrv.o:
Stats.o:
