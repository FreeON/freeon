include $(top_srcdir)/make.inc

noinst_PROGRAMS = \
  sleeper \
  lock \
  sendReceive_DBL_RNK2 \
  sendReceive_HGRho

sleeper_SOURCES = sleeper.f90
sleeper_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

lock_SOURCES = lock.f90
lock_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

sendReceive_DBL_RNK2_SOURCES = sendReceive_DBL_RNK2.f90
sendReceive_DBL_RNK2_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

sendReceive_HGRho_SOURCES = sendReceive_HGRho.f90
sendReceive_HGRho_LDADD = $(LDADD_MODULES) $(LDADD_INTERNAL_LAPACK) $(LDADD_INTERNAL_HDF5)

TESTS = sleepertest locktest

EXTRA_DIST = $(TESTS)
