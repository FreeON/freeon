include $(top_srcdir)/make.inc

bin_PROGRAMS = ONX GONX

ONX_SOURCES = \
  QSDis.F \
  PartDrv.F90 \
  ONXParameters.F90 \
  ONXInit.F90 \
  ONXGet.F90 \
  ONXFillOut.F90 \
  Stats.F90 \
  ONXRng.F90 \
  ONX2F77.F \
  ONX2DataType.F90 \
  ONXCtrSclg.F90 \
  ONX2List.F90 \
  ONX2ComputK.F90 \
  ONX2.F90

GONX_SOURCES = \
  QSDis.F \
  PartDrv.F90 \
  ONXParameters.F90 \
  ONXGet.F90 \
  ONX2F77.F \
  ONX2DataType.F90 \
  ONXCtrSclg.F90 \
  GONX2ComputDK.F90 \
  GONX2.F90 \
  ONX2List.F90

ONX_DEPENDENCIES  = \
  $(LIBRARY_DEPENDENCY_FILES) \
  $(top_builddir)/TwoE/GRD/libGRD.la \
  $(top_builddir)/TwoE/INT/libINT.la \
  $(top_builddir)/TwoE/VRR/libVRR.la \
  $(top_builddir)/TwoE/HRR/libHRR.la \
  $(top_builddir)/TwoE/MSC/libTwoEMisc.la \
  $(top_builddir)/TwoE/MSC/shellpairstruct.mod \
  $(top_builddir)/TwoE/MSC/vscratchb.mod \
  $(top_builddir)/TwoE/MSC/getatompairmod.mod

GONX_DEPENDENCIES = \
  $(LIBRARY_DEPENDENCY_FILES) \
  $(top_builddir)/TwoE/GRD/libGRD.la \
  $(top_builddir)/TwoE/INT/libINT.la \
  $(top_builddir)/TwoE/VRR/libVRR.la \
  $(top_builddir)/TwoE/HRR/libHRR.la \
  $(top_builddir)/TwoE/MSC/libTwoEMisc.la \
  $(top_builddir)/TwoE/MSC/shellpairstruct.mod \
  $(top_builddir)/TwoE/MSC/vscratchb.mod \
  $(top_builddir)/TwoE/MSC/getatompairmod.mod

SOURCES_DEPCHECK = \
  $(ONX_SOURCES) \
  $(GONX_SOURCES)

EXTRA_DIST = \
  ERIMIC.Inc

CPPMISC = -DONX2_PARALLEL -DGTRESH

EXTRA_INCLUDES = \
  -I$(top_builddir)/TwoE/INT \
  -I$(top_srcdir)/TwoE/INT \
  -I$(top_builddir)/TwoE/GRD \
  -I$(top_srcdir)/TwoE/GRD \
  -I$(top_builddir)/TwoE/MSC

LDADD = \
  $(LDADD_MODULES) \
  $(builddir)/../TwoE/GRD/libGRD.la \
  $(builddir)/../TwoE/INT/libINT.la \
  $(builddir)/../TwoE/VRR/libVRR.la \
  $(builddir)/../TwoE/HRR/libHRR.la \
  $(builddir)/../TwoE/MSC/libTwoEMisc.la

# In case we have ctags installed, make the necessary tags files.
.PHONY : tags
tags :
	ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{f,f90} ../Modules/*.{c,f,f90}

# Dependencies.
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
  ONXParameters.o \
  ONXGet.o

ONX2.o : \
  ONX2DataType.o \
  ONX2List.o \
  ONX2ComputK.o \
  ONXInit.o \
  ONXCtrSclg.o \
  ONXRng.o \
  ONXFillOut.o \
  ONXGet.o \
  $(top_builddir)/MondoConfig.h

GONX2ComputDK.o : \
  ONX2DataType.o \
  ONXParameters.o \
  ONX2List.o

GONX2.o : \
  ONX2DataType.o \
  ONXCtrSclg.o \
  ONX2List.o \
  GONX2ComputDK.o \
  ONXGet.o \
  $(top_builddir)/MondoConfig.h
