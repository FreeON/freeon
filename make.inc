# vim: syntax=make:tw=0:noexpandtab
#
# Some general things everyone needs.

# Set some variables.
CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

# Fix for bug in autoconf-2.59
builddir = .

if INTERNAL_LAPACK
LDADD_INTERNAL_LAPACK = \
  $(top_builddir)/Modules/lapack/lapack/libfreeonlapack.la \
  $(top_builddir)/Modules/lapack/blas/libfreeonblas.la \
  $(top_builddir)/Modules/lapack/install/libfreeonlapackextra.la
endif

if INTERNAL_HDF5
LDADD_INTERNAL_HDF5 = $(top_builddir)/Modules/hdf5-1.8.3/src/libhdf5.la
endif

LDADD_MODULES = $(top_builddir)/Modules/libfreeonmodules.la

# On filesystems which are case-insensitive (which really means HFS+ on OS X)
# we can not preprocess into a filename that differs by the case only, i.e. .C
# -> .c.
PREPROCESSED_SUFFIX = -preprocessed

# For cleanup.
CLEANFILES = \
  *.mod \
  *$(PREPROCESSED_SUFFIX).f \
  *$(PREPROCESSED_SUFFIX).f90 \
  *$(PREPROCESSED_SUFFIX).c \
  *~

# Dependency check (only for debugging).
depcheck-bin :
	for i in $(SOURCES_DEPCHECK); do \
	  echo "********************************** cleaning **********************************"; \
	  $(MAKE) clean; \
	  echo "********************************** verifying `echo $$i | sed -e "s/[.][^.]\+$$//"` **********************************"; \
	  $(MAKE) `echo $$i | sed -e "s/[.][^.]\+$$//"`.o || exit 1; \
	done
	@echo "********************************** done checking, rebuilding **********************************"
	$(MAKE) clean
	$(MAKE)
	@echo "********************************** done checking **********************************"

depcheck-lib :
	for i in $(SOURCES_DEPCHECK); do \
	  echo "********************************** cleaning **********************************"; \
	  $(MAKE) clean; \
	  echo "********************************** verifying `echo $$i | sed -e "s/[.][^.]\+$$//"` **********************************"; \
	  $(MAKE) `echo $$i | sed -e "s/[.][^.]\+$$//"`.lo || exit 1; \
	done
	@echo "********************************** done checking, rebuilding **********************************"
	$(MAKE) clean
	$(MAKE)
	@echo "********************************** done checking **********************************"

CPP_COMMAND = $(CPP) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(FORTRAN_CPPFLAGS) $(CPPFLAGS) $(CPPMISC)

CONFIGURE_CFLAGS = @CONFIGURE_CFLAGS@
CONFIGURE_FFLAGS = @CONFIGURE_FFLAGS@
CONFIGURE_FCFLAGS = @CONFIGURE_FCFLAGS@

AM_CFLAGS  = $(CONFIGURE_CFLAGS) $(EXTRA_INCLUDES)
AM_FFLAGS  = $(CONFIGURE_FFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(DEFAULT_INCLUDES) $(EXTRA_INCLUDES)
AM_FCFLAGS = $(CONFIGURE_FCFLAGS) $(FORTRAN_DEFAULT_INCLUDES) $(DEFAULT_INCLUDES) $(EXTRA_INCLUDES)

if USE_COMPILER_OPTIMIZATIONS
AM_CFLAGS  += $(OPTIMIZATION_FLAGS)
AM_FFLAGS  += $(OPTIMIZATION_FLAGS)
AM_FCFLAGS += $(OPTIMIZATION_FLAGS)
endif

# Whether we should keep or remove preprocessed files.
REMOVE_PREPROCESSED_FILES = @REMOVE_PREPROCESSED_FILES@

# Turn off output of libtool.
AM_LIBTOOLFLAGS = --silent

# Suffix rules.

if USE_COMPILER_OPTIMIZATIONS

.f90.o :
	@$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90 && exit 0; \
	  echo "The following preprocessor command failed:" 1>&2; \
	  echo '$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90' 1>&2; \
	  $(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@echo "compiling $<"
	@$(FCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90 > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(FCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90' 1>&2; \
	  $(FCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f90; \
	fi

.f90.lo :
	@$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90 && exit 0; \
	  echo "The following preprocessor command failed:" 1>&2; \
	  echo '$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90' 1>&2; \
	  $(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@echo "compiling $<"
	@$(LTFCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90 > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(LTFCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90' 1>&2; \
	  $(LTFCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f90; \
	fi

.f.o :
	@$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f && exit 0; \
	  echo "The following preprocessor command failed:" 1>&2; \
	  echo '$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f' 1>&2; \
	  $(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f
	@echo "compiling $<"
	@$(F77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(F77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f' 1>&2; \
	  $(F77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f; \
	fi

.f.lo :
	@$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f && exit 0; \
	  echo "The following preprocessor command failed:" 1>&2; \
	  echo '$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f' 1>&2; \
	  $(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f
	@echo "compiling $<"
	@$(LTF77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f > /dev/null 2>&1 && exit 0; \
	  echo "The following compiler command failed:" 1>&2; \
	  echo '$(LTF77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f' 1>&2; \
	  $(LTF77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f; \
	fi

else

.f90.o :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90
	$(FCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f90; \
	fi

.f90.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f90
	$(LTFCCOMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f90
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f90; \
	fi

.f.o :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f
	@$(F77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f; \
	fi

.f.lo :
	$(CPP_COMMAND) $< | grep -v -E "^#" > $(basename $@)$(PREPROCESSED_SUFFIX).f
	$(LTF77COMPILE) -c -o $@ $(basename $@)$(PREPROCESSED_SUFFIX).f
	@if test "${REMOVE_PREPROCESSED_FILES}" = "yes"; then \
	  rm $(basename $@)$(PREPROCESSED_SUFFIX).f; \
	fi

endif

if USE_DYNAMIC_LINKING
LIBRARY_DEPENDENCY_FILES =
else
LIBRARY_DEPENDENCY_FILES = $(LIBRARY_DEPENDENCY_FILES_LIBS)
endif

LIBRARY_DEPENDENCY_FILES_MODS = \
  $(top_builddir)/Modules/ainv.mod \
  $(top_builddir)/Modules/annealmap.mod \
  $(top_builddir)/Modules/atompairs.mod \
  $(top_builddir)/Modules/basissetparameters.mod \
  $(top_builddir)/Modules/boundingbox.mod \
  $(top_builddir)/Modules/brabloks.mod \
  $(top_builddir)/Modules/pbc.mod \
  $(top_builddir)/Modules/cholfactor.mod \
  $(top_builddir)/Modules/clock.mod \
  $(top_builddir)/Modules/common_debug.mod \
  $(top_builddir)/Modules/derivedtypes.mod \
  $(top_builddir)/Modules/erffunk.mod \
  $(top_builddir)/Modules/fastmatrices.mod \
  $(top_builddir)/Modules/functionals.mod \
  $(top_builddir)/Modules/gammaassymp.mod \
  $(top_builddir)/Modules/gammaf0.mod \
  $(top_builddir)/Modules/gammaf1.mod \
  $(top_builddir)/Modules/gammaf10.mod \
  $(top_builddir)/Modules/gammaf11.mod \
  $(top_builddir)/Modules/gammaf12.mod \
  $(top_builddir)/Modules/gammaf13.mod \
  $(top_builddir)/Modules/gammaf2.mod \
  $(top_builddir)/Modules/gammaf3.mod \
  $(top_builddir)/Modules/gammaf4.mod \
  $(top_builddir)/Modules/gammaf5.mod \
  $(top_builddir)/Modules/gammaf6.mod \
  $(top_builddir)/Modules/gammaf7.mod \
  $(top_builddir)/Modules/gammaf8.mod \
  $(top_builddir)/Modules/gammaf9.mod \
  $(top_builddir)/Modules/gammafunctions.mod \
  $(top_builddir)/Modules/globalcharacters.mod \
  $(top_builddir)/Modules/globalobjects.mod \
  $(top_builddir)/Modules/globalscalars.mod \
  $(top_builddir)/Modules/indexing.mod \
  $(top_builddir)/Modules/inout.mod \
  $(top_builddir)/Modules/invexp.mod \
  $(top_builddir)/Modules/linalg.mod \
  $(top_builddir)/Modules/macros.mod \
  $(top_builddir)/Modules/matfunk.mod \
  $(top_builddir)/Modules/mcmurchie.mod \
  $(top_builddir)/Modules/mechanics.mod \
  $(top_builddir)/Modules/memman.mod \
  $(top_builddir)/Modules/mondologger.mod \
  $(top_builddir)/Modules/mondompi.mod \
  $(top_builddir)/Modules/mpiinclude.mod \
  $(top_builddir)/Modules/opt_trav_band.mod \
  $(top_builddir)/Modules/order.mod \
  $(top_builddir)/Modules/parse.mod \
  $(top_builddir)/Modules/parsingconstants.mod \
  $(top_builddir)/Modules/prettyprint.mod \
  $(top_builddir)/Modules/processcontrol.mod \
  $(top_builddir)/Modules/setxyz.mod \
  $(top_builddir)/Modules/slatec.mod \
  $(top_builddir)/Modules/specfun.mod \
  $(top_builddir)/Modules/specfunmesh.mod \
  $(top_builddir)/Modules/thresholding.mod \
  $(top_builddir)/Modules/utilities.mod

LIBRARY_DEPENDENCY_FILES_LIBS = \
  $(LDADD_INTERNAL_LAPACK) \
  $(LDADD_INTERNAL_HDF5) \
  $(LDADD_MODULES)

$(top_builddir)/MondoLocalVersion.h : $(top_builddir)/config_localversion.h
	@echo "checking MondoLocalVersion.h"
	@grep LOCAL_VERSION $(top_builddir)/config_localversion.h > $(top_builddir)/tempconf
	@if test ! -f $(top_builddir)/MondoLocalVersion.h; then \
	  mv -f $(top_builddir)/tempconf $(top_builddir)/MondoLocalVersion.h; \
	  echo "MondoLocalVersion.h did not exist"; \
	elif (diff -u $(top_builddir)/MondoLocalVersion.h $(top_builddir)/tempconf > $(top_builddir)/MondoLocalVersion.diff); then \
	  rm -f $(top_builddir)/tempconf; \
	  echo "MondoLocalVersion.h unchanged"; \
	else \
	  mv -f $(top_builddir)/tempconf $(top_builddir)/MondoLocalVersion.h; \
	  echo "MondoLocalVersion.h updated"; \
	  cat $(top_builddir)/MondoLocalVersion.diff; \
	fi
