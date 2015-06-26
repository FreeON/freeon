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
  $(top_builddir)/TwoE/libfreeonTwoE.la \
  $(top_builddir)/TwoE/MSC/shellpairstruct.mod \
  $(top_builddir)/TwoE/MSC/vscratchb.mod \
  $(top_builddir)/TwoE/MSC/getatompairmod.mod

GONX_DEPENDENCIES = \
  $(LIBRARY_DEPENDENCY_FILES) \
  $(top_builddir)/TwoE/libfreeonTwoE.la \
  $(top_builddir)/TwoE/MSC/shellpairstruct.mod \
  $(top_builddir)/TwoE/MSC/vscratchb.mod \
  $(top_builddir)/TwoE/MSC/getatompairmod.mod

SOURCES_DEPCHECK = \
  $(ONX_SOURCES) \
  $(GONX_SOURCES)

EXTRA_DIST = ERIMIC.Inc

CPPMISC = -DONX2_PARALLEL -DGTRESH

EXTRA_INCLUDES = \
	-I$(top_builddir)/TwoE/GRD \
	-I$(top_builddir)/TwoE/INT \
	-I$(top_builddir)/TwoE/MSC \
	-I$(top_srcdir)/TwoE/GRD \
	-I$(top_srcdir)/TwoE/INT

LDADD = \
  $(LDADD_MODULES) \
  $(LAPACK_LIBS) \
  $(BLAS_LIBS) \
  $(builddir)/../TwoE/libfreeonTwoE.la

TAGS_FILES = \
	$(top_srcdir)/Modules/*.c
	$(top_srcdir)/Modules/*.F
	$(top_srcdir)/Modules/*.F90

AM_CTAGSFLAGS = --Fortran-kinds=+i+L

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
