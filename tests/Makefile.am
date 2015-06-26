include $(top_srcdir)/make.inc

noinst_PROGRAMS = \
  sleeper \
  lock \
  sendReceive_DBL_RNK2 \
  sendReceive_HGRho

sleeper_SOURCES = sleeper.F90
sleeper_LDADD = $(LDADD_MODULES) $(LAPACK_LIBS)

lock_SOURCES = lock.F90
lock_LDADD = $(LDADD_MODULES) $(LAPACK_LIBS)

sendReceive_DBL_RNK2_SOURCES = sendReceive_DBL_RNK2.F90
sendReceive_DBL_RNK2_LDADD = $(LDADD_MODULES) $(LAPACK_LIBS)

sendReceive_HGRho_SOURCES = sendReceive_HGRho.F90
sendReceive_HGRho_LDADD = $(LDADD_MODULES) $(LAPACK_LIBS)

if HAVE_MPIRUN
MPI_TESTS = locktest
endif

TESTS = sleepertest $(MPI_TESTS)

EXTRA_DIST = sleepertest locktest
