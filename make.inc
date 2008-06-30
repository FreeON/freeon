# vim: syntax=make:tw=0:noexpandtab
#
# Some general things everyone needs.

# Set some variables.
CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

if INTERNAL_LAPACK
  LDADD_INTERNAL_LAPACK = ../MondoMods/lapack/lapack/libMondoLAPACK.la \
    ../MondoMods/lapack/blas/libMondoBLAS.la \
    ../MondoMods/lapack/install/libMondoLAPACKEXTRA.la
endif

if INTERNAL_HDF5
  LDADD_INTERNAL_HDF5 = ../MondoMods/hdf5/src/libhdf5.la
endif

LDADD_MONDOMODS = ../MondoMods/libMondoMods.la

# On filesystems which are case-insensitive (which really means HFS+ on OS X)
# we can not preprocess into a filename that differs by the case only, i.e. .C
# -> .c.
PREPROCESSED_SUFFIX = -preprocessed

# For cleanup.
CLEANFILES = *.mod *$(PREPROCESSED_SUFFIX).f *$(PREPROCESSED_SUFFIX).f90 *$(PREPROCESSED_SUFFIX).c *~

# Dependency check (only for debugging).
depcheck :
	for i in $(SOURCES_DEPCHECK); do \
	  echo "********************************** cleaning **********************************"; \
	  $(MAKE) clean; \
	  echo "********************************** verifying `echo $$i | sed -e "s/[.][^.]\+$$//"` **********************************"; \
	  $(MAKE) `echo $$i | sed -e "s/[.][^.]\+$$//"`.lo || exit 1; \
	done

CPP_COMMAND = $(CPP) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(FORTRAN_CPPFLAGS) $(CPPFLAGS) $(CPPMISC)
FC_COMMAND  = $(FC) $(AM_FCFLAGS) $(FCFLAGS) $(OPTIMIZATION_FLAGS) $(EXTRA_INCLUDES)
F77_COMMAND = $(F77) $(AM_FFLAGS) $(FFLAGS) $(OPTIMIZATION_FLAGS) $(EXTRA_INCLUDES)
CC_COMMAND  = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CPPMISC) $(AM_CFLAGS) $(CFLAGS) $(OPTIMIZATION_FLAGS) $(EXTRA_INCLUDES)

# Whether we should keep or remove preprocessed files.
REMOVE_PREPROCESSED_FILES = @REMOVE_PREPROCESSED_FILES@

# Suffix rules.
.f90.o :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $<$(PREPROCESSED_SUFFIX).f90
	$(FC_COMMAND) -c -o $@ $<$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $<$(PREPROCESSED_SUFFIX).f90; fi

.f90.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $<$(PREPROCESSED_SUFFIX).f90
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(FC_COMMAND) -c -o $@ $<$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $<$(PREPROCESSED_SUFFIX).f90; fi

.f.o :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $<$(PREPROCESSED_SUFFIX).f
	$(F77_COMMAND) -c -o $@ $<$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $<$(PREPROCESSED_SUFFIX).f; fi

.f.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $<$(PREPROCESSED_SUFFIX).f
	$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77_COMMAND) -c -o $@ $<$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $<$(PREPROCESSED_SUFFIX).f; fi

# For dependencies of -L and -I.
if INTERNAL_LAPACK
  INTERNAL_LAPACK_DEPENDENCY_FILES = \
    $(top_srcdir)/MondoMods/lapack/blas/libMondoBLAS.la \
    $(top_srcdir)/MondoMods/lapack/lapack/libMondoLAPACK.la \
    $(top_srcdir)/MondoMods/lapack/install/libMondoLAPACKEXTRA.la
endif

if INTERNAL_HDF5
  INTERNAL_HDF5_DEPENDENCY_FILES = $(top_srcdir)/MondoMods/hdf5/src/libhdf5.la
endif

MondoMods_DEPENDENCY_FILES = \
  $(top_srcdir)/MondoMods/libMondoMods.la \
  $(top_srcdir)/MondoMods/ainv.mod \
  $(top_srcdir)/MondoMods/annealmap.mod \
  $(top_srcdir)/MondoMods/atompairs.mod \
  $(top_srcdir)/MondoMods/basissetparameters.mod \
  $(top_srcdir)/MondoMods/boundingbox.mod \
  $(top_srcdir)/MondoMods/brabloks.mod \
  $(top_srcdir)/MondoMods/cellsets.mod \
  $(top_srcdir)/MondoMods/cholfactor.mod \
  $(top_srcdir)/MondoMods/clock.mod \
  $(top_srcdir)/MondoMods/common_debug.mod \
  $(top_srcdir)/MondoMods/derivedtypes.mod \
  $(top_srcdir)/MondoMods/erffunk.mod \
  $(top_srcdir)/MondoMods/fastmatrices.mod \
  $(top_srcdir)/MondoMods/cwrappers.mod \
  $(top_srcdir)/MondoMods/functionals.mod \
  $(top_srcdir)/MondoMods/gammaassymp.mod \
  $(top_srcdir)/MondoMods/gammaf0.mod \
  $(top_srcdir)/MondoMods/gammaf1.mod \
  $(top_srcdir)/MondoMods/gammaf10.mod \
  $(top_srcdir)/MondoMods/gammaf11.mod \
  $(top_srcdir)/MondoMods/gammaf12.mod \
  $(top_srcdir)/MondoMods/gammaf13.mod \
  $(top_srcdir)/MondoMods/gammaf2.mod \
  $(top_srcdir)/MondoMods/gammaf3.mod \
  $(top_srcdir)/MondoMods/gammaf4.mod \
  $(top_srcdir)/MondoMods/gammaf5.mod \
  $(top_srcdir)/MondoMods/gammaf6.mod \
  $(top_srcdir)/MondoMods/gammaf7.mod \
  $(top_srcdir)/MondoMods/gammaf8.mod \
  $(top_srcdir)/MondoMods/gammaf9.mod \
  $(top_srcdir)/MondoMods/gammafunctions.mod \
  $(top_srcdir)/MondoMods/globalcharacters.mod \
  $(top_srcdir)/MondoMods/globalobjects.mod \
  $(top_srcdir)/MondoMods/globalscalars.mod \
  $(top_srcdir)/MondoMods/indexing.mod \
  $(top_srcdir)/MondoMods/inout.mod \
  $(top_srcdir)/MondoMods/invexp.mod \
  $(top_srcdir)/MondoMods/linalg.mod \
  $(top_srcdir)/MondoMods/macros.mod \
  $(top_srcdir)/MondoMods/matfunk.mod \
  $(top_srcdir)/MondoMods/mcmurchie.mod \
  $(top_srcdir)/MondoMods/mechanics.mod \
  $(top_srcdir)/MondoMods/memman.mod \
  $(top_srcdir)/MondoMods/mondologger.mod \
  $(top_srcdir)/MondoMods/mondompi.mod \
  $(top_srcdir)/MondoMods/mpiinclude.mod \
  $(top_srcdir)/MondoMods/opt_trav_band.mod \
  $(top_srcdir)/MondoMods/order.mod \
  $(top_srcdir)/MondoMods/parse.mod \
  $(top_srcdir)/MondoMods/parsingconstants.mod \
  $(top_srcdir)/MondoMods/prettyprint.mod \
  $(top_srcdir)/MondoMods/processcontrol.mod \
  $(top_srcdir)/MondoMods/setxyz.mod \
  $(top_srcdir)/MondoMods/slatec.mod \
  $(top_srcdir)/MondoMods/specfun.mod \
  $(top_srcdir)/MondoMods/specfunmesh.mod \
  $(top_srcdir)/MondoMods/thresholding.mod \
  $(top_srcdir)/MondoMods/utilities.mod \
  $(INTERNAL_LAPACK_DEPENDENCY_FILES) \
  $(INTERNAL_HDF5_DEPENDENCY_FILES)
