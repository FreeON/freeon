include $(top_srcdir)/make.inc

noinst_PROGRAMS = sleeper

sleeper_SOURCES = sleeper.f90
sleeper_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

TESTS = sleepertest

EXTRA_DIST = $(TESTS)
