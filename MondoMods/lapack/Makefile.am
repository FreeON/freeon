include $(top_srcdir)/make.inc

if INTERNAL_LAPACK_CHECKS
  ADDITIONAL_DIRS = blas_testing lapack_testing
endif

DIST_SUBDIRS = install blas lapack blas_testing lapack_testing
SUBDIRS      = install blas lapack $(ADDITIONAL_DIRS)

EXTRADIST = COPYING
