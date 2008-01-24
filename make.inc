# vim: syntax=make:tw=0:noexpandtab
#
# Some general things everyone needs.

# Set some variables.
CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

if INTERNAL_LAPACK
  LDADD_INTERNAL_LAPACK = -L../MondoMods/lapack/lapack -L../MondoMods/lapack/blas
endif

if INTERNAL_HDF5
  LDADD_INTERNAL_HDF5 = -L../MondoMods/hdf5
endif

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
	  $(MAKE) `echo $$i | sed -e "s/[.][^.]\+$$//"`.o || exit 1; \
	done

CPP_COMMAND = $(CPP) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(FORTRAN_CPPFLAGS) $(CPPFLAGS) $(CPPMISC)
FC_COMMAND  = $(FC) $(AM_FCFLAGS) $(FCFLAGS) $(EXTRA_INCLUDES)
F77_COMMAND = $(F77) $(AM_FFLAGS) $(FFLAGS) $(EXTRA_INCLUDES)
CC_COMMAND  = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CPPMISC) $(AM_CFLAGS) $(CFLAGS) $(EXTRA_INCLUDES)

REMOVE_PREPROCESSED_FILES = @REMOVE_PREPROCESSED_FILES@

# Suffix rules.
.f90.o :
	@$(CPP_COMMAND) $< | grep -v -E "^#" | sed 's/&[ \t]\+/&\n/g' > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	$(FC_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f90; \
	fi

.f.o :
	@$(CPP_COMMAND) $< | grep -v -E "^#" > $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	$(F77_COMMAND) -c -o $@ $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(notdir $(basename $<))$(PREPROCESSED_SUFFIX).f; \
	fi

# For dependencies of -L and -I.
if INTERNAL_HDF5
  INTERNAL_HDF5_DEPENDENCY_FILES = ../MondoMods/hdf5/libMondoHDF5.a
endif

if INTERNAL_LAPACK
  INTERNAL_LAPACK_DEPENDENCY_FILES = \
    ../MondoMods/lapack/blas/libMondoBLAS.a \
    ../MondoMods/lapack/lapack/libMondoLAPACK.a
endif

MondoMods_DEPENDENCY_FILES = \
  ../MondoMods/libMondoMods.a \
  ../MondoMods/globalscalars.mod \
  ../MondoMods/globalcharacters.mod \
  ../MondoMods/derivedtypes.mod \
  ../MondoMods/globalobjects.mod \
  ../MondoMods/basissetparameters.mod \
  ../MondoMods/mpiinclude.mod \
  ../MondoMods/processcontrol.mod \
  ../MondoMods/memman.mod \
  ../MondoMods/clock.mod \
  ../MondoMods/parsingconstants.mod \
  ../MondoMods/parse.mod \
  ../MondoMods/mondompi.mod \
  ../MondoMods/indexing.mod \
  ../MondoMods/setxyz.mod \
  ../MondoMods/inout.mod \
  ../MondoMods/prettyprint.mod \
  ../MondoMods/opt_trav_band.mod \
  ../MondoMods/annealmap.mod \
  ../MondoMods/order.mod \
  ../MondoMods/mcmurchie.mod \
  ../MondoMods/gammaassymp.mod \
  ../MondoMods/specfunmesh.mod \
  ../MondoMods/gammaf0.mod \
  ../MondoMods/gammaf1.mod \
  ../MondoMods/gammaf2.mod \
  ../MondoMods/gammaf3.mod \
  ../MondoMods/gammaf4.mod \
  ../MondoMods/gammaf5.mod \
  ../MondoMods/gammaf6.mod \
  ../MondoMods/gammaf7.mod \
  ../MondoMods/gammaf8.mod \
  ../MondoMods/gammaf9.mod \
  ../MondoMods/gammaf10.mod \
  ../MondoMods/gammaf11.mod \
  ../MondoMods/gammaf12.mod \
  ../MondoMods/gammaf13.mod \
  ../MondoMods/invexp.mod \
  ../MondoMods/erffunk.mod \
  ../MondoMods/gammafunctions.mod \
  ../MondoMods/specfun.mod \
  ../MondoMods/thresholding.mod \
  ../MondoMods/common_debug.mod \
  ../MondoMods/linalg.mod \
  ../MondoMods/matfunk.mod \
  ../MondoMods/ainv.mod \
  ../MondoMods/mechanics.mod \
  ../MondoMods/functionals.mod \
  ../MondoMods/cellsets.mod \
  ../MondoMods/boundingbox.mod \
  ../MondoMods/atompairs.mod \
  ../MondoMods/brabloks.mod \
  ../MondoMods/macros.mod \
  ../MondoMods/fastmatrices.mod \
  ../MondoMods/cholfactor.mod \
  ../MondoMods/slatec.mod \
  $(INTERNAL_LAPACK_DEPENDENCY_FILES) \
  $(INTERNAL_HDF5_DEPENDENCY_FILES)
