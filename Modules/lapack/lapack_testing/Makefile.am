include $(top_srcdir)/make.inc

SUBDIRS = matgen lin eig

TESTS = \
  single_lin_test.sh \
  double_lin_test.sh \
  complex_lin_test.sh \
  complex16_lin_test.sh \
  single_double_proto_lin_test.sh \
  complex_complex16_lin_test.sh \
  complex16_eig_test.sh \
  complex_eig_test.sh \
  double_eig_test.sh \
  single_eig_test.sh

noinst_SCRIPTS = \
  $(TESTS) \
	test_functions.sh

EXTRA_DIST = \
  $(noinst_SCRIPTS) \
  $(TEST_INPUTS)

TEST_INPUTS = \
  cbb.in \
  cgbak.in \
  cgg.in \
  ctest.in \
  dbb.in \
  dgbak.in \
  dgg.in \
  dstest.in \
  gqr.in \
  nep.in \
  sbb.in \
  sep.in \
  sgd.in \
  ssg.in \
  zbak.in \
  zctest.in \
  zgbak.in \
  zgg.in \
  ztest.in \
  cbak.in \
  cec.in \
  cgbal.in \
  csb.in \
  dbak.in \
  dec.in \
  dgbal.in \
  dsb.in \
  dtest.in \
  gsv.in \
  sbak.in \
  sec.in \
  sgbak.in \
  sgg.in \
  stest.in \
  zbal.in \
  zec.in \
  zgbal.in \
  zsb.in \
  cbal.in \
  ced.in \
  cgd.in \
  csg.in \
  dbal.in \
  ded.in \
  dgd.in \
  dsg.in \
  glm.in \
  lse.in \
  sbal.in \
  sed.in \
  sgbal.in \
  ssb.in \
  svd.in \
  zbb.in \
  zed.in \
  zgd.in \
  zsg.in

LDADD = $(LDADD_INTERNAL_LAPACK)

CLEANFILES += \
  cnep.out \
  ctest.out \
  dsep.out \
  dtest.out \
  sed.out \
  sgg.out \
  ssep.out \
  stest.out \
  znep.out \
  ztest.out \
  csep.out \
  dnep.out \
  dstest.out \
  sec.out \
  sgd.out \
  snep.out \
  ssvd.out \
  zctest.out \
  zsep.out
