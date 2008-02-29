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
	  echo cleaning...; \
	  $(MAKE) clean; \
	  echo verifying `echo $$i | sed -e "s/[.][^.]\+$$//"`; \
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
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	$(FC_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90; fi

.f90.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(FC_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90; fi

.f.o :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	$(F77_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f; fi

.f.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	$(LIBTOOL) --tag=F77 $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=compile $(F77_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f; fi

# For dependencies of -L and -I.
if INTERNAL_LAPACK
  INTERNAL_LAPACK_DEPENDENCY_FILES = \
    ../MondoMods/lapack/blas/libMondoBLAS.la \
    ../MondoMods/lapack/lapack/libMondoLAPACK.la \
    ../MondoMods/lapack/install/libMondoLAPACKEXTRA.la
endif

if INTERNAL_HDF5
  INTERNAL_HDF5_DEPENDENCY_FILES = ../MondoMods/hdf5/src/libhdf5.la
endif

MondoMods_DEPENDENCY_FILES = \
  ../MondoMods/libMondoMods.la \
  ../MondoMods/ainv.mod \
  ../MondoMods/annealmap.mod \
  ../MondoMods/atompairs.mod \
  ../MondoMods/basissetparameters.mod \
  ../MondoMods/boundingbox.mod \
  ../MondoMods/brabloks.mod \
  ../MondoMods/cellsets.mod \
  ../MondoMods/cholfactor.mod \
  ../MondoMods/clock.mod \
  ../MondoMods/common_debug.mod \
  ../MondoMods/derivedtypes.mod \
  ../MondoMods/erffunk.mod \
  ../MondoMods/fastmatrices.mod \
  ../MondoMods/cwrappers.mod \
  ../MondoMods/functionals.mod \
  ../MondoMods/gammaassymp.mod \
  ../MondoMods/gammaf0.mod \
  ../MondoMods/gammaf1.mod \
  ../MondoMods/gammaf10.mod \
  ../MondoMods/gammaf11.mod \
  ../MondoMods/gammaf12.mod \
  ../MondoMods/gammaf13.mod \
  ../MondoMods/gammaf2.mod \
  ../MondoMods/gammaf3.mod \
  ../MondoMods/gammaf4.mod \
  ../MondoMods/gammaf5.mod \
  ../MondoMods/gammaf6.mod \
  ../MondoMods/gammaf7.mod \
  ../MondoMods/gammaf8.mod \
  ../MondoMods/gammaf9.mod \
  ../MondoMods/gammafunctions.mod \
  ../MondoMods/globalcharacters.mod \
  ../MondoMods/globalobjects.mod \
  ../MondoMods/globalscalars.mod \
  ../MondoMods/indexing.mod \
  ../MondoMods/inout.mod \
  ../MondoMods/invexp.mod \
  ../MondoMods/linalg.mod \
  ../MondoMods/macros.mod \
  ../MondoMods/matfunk.mod \
  ../MondoMods/mcmurchie.mod \
  ../MondoMods/mechanics.mod \
  ../MondoMods/memman.mod \
  ../MondoMods/mondologger.mod \
  ../MondoMods/mondompi.mod \
  ../MondoMods/mpiinclude.mod \
  ../MondoMods/opt_trav_band.mod \
  ../MondoMods/order.mod \
  ../MondoMods/parse.mod \
  ../MondoMods/parsingconstants.mod \
  ../MondoMods/prettyprint.mod \
  ../MondoMods/processcontrol.mod \
  ../MondoMods/setxyz.mod \
  ../MondoMods/slatec.mod \
  ../MondoMods/specfun.mod \
  ../MondoMods/specfunmesh.mod \
  ../MondoMods/thresholding.mod \
  ../MondoMods/utilities.mod \
  $(INTERNAL_LAPACK_DEPENDENCY_FILES) \
  $(INTERNAL_HDF5_DEPENDENCY_FILES)
