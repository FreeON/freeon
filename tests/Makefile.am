include $(top_srcdir)/make.inc

noinst_PROGRAMS = sleeper lock

sleeper_SOURCES = sleeper.f90
sleeper_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

lock_SOURCES = lock.f90
lock_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

TESTS = sleepertest locktest

EXTRA_DIST = $(TESTS)
