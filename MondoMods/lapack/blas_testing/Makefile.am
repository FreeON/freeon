include $(top_srcdir)/make.inc

noinst_PROGRAMS = \
  xblat1s xblat1c xblat1d xblat1z \
  xblat2s xblat2c xblat2d xblat2z \
  xblat3s xblat3c xblat3d xblat3z

TESTS = \
  sblat1_test.sh \
  cblat1_test.sh \
  dblat1_test.sh \
  zblat1_test.sh \
  sblat2_test.sh \
  cblat2_test.sh \
  dblat2_test.sh \
  zblat2_test.sh \
  sblat3_test.sh \
  cblat3_test.sh \
  dblat3_test.sh \
  zblat3_test.sh

noinst_SCRIPTS = $(TESTS)

EXTRA_DIST = \
  $(noinst_SCRIPTS) \
  $(TEST_INPUTS)

TEST_INPUTS = \
  cblat2.in \
  cblat3.in \
  dblat2.in \
  dblat3.in \
  sblat2.in \
  sblat3.in \
  zblat2.in \
  zblat3.in

xblat1s_SOURCES = sblat1.f
xblat1c_SOURCES = cblat1.f 
xblat1d_SOURCES = dblat1.f 
xblat1z_SOURCES = zblat1.f

xblat2s_SOURCES = sblat2.f
xblat2c_SOURCES = cblat2.f 
xblat2d_SOURCES = dblat2.f 
xblat2z_SOURCES = zblat2.f

xblat3s_SOURCES = sblat3.f
xblat3c_SOURCES = cblat3.f 
xblat3d_SOURCES = dblat3.f 
xblat3z_SOURCES = zblat3.f

LDADD = $(LDADD_INTERNAL_LAPACK)

CLEANFILES += \
  blat1.out \
  sblat1.out \
  dblat1.out \
  cblat1.out \
  zblat1.out \
  blat2.out \
  sblat2.out \
  dblat2.out \
  cblat2.out \
  zblat2.out \
  blat3.out \
  sblat3.out \
  dblat3.out \
  cblat3.out \
  zblat3.out
