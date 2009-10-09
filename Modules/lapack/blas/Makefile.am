include $(top_srcdir)/make.inc

if USE_DYNAMIC_LINKING
lib_LTLIBRARIES = libfreeonblas.la
else
noinst_LTLIBRARIES = libfreeonblas.la
endif

libfreeonblas_la_SOURCES = \
  $(c_source_files) \
  $(d_source_files) \
  $(s_source_files) \
  $(z_source_files) \
  icamax.f \
  idamax.f \
  isamax.f \
  izamax.f \
  xerbla.f

# This file also exists in install.
#lsame.f

c_source_files = \
  caxpy.f  \
  ccopy.f  \
  cdotc.f  \
  cdotu.f  \
  cgbmv.f  \
  cgemm.f  \
  cgemv.f  \
  cgerc.f  \
  cgeru.f  \
  chbmv.f  \
  chemm.f  \
  chemv.f  \
  cher2.f  \
  cher2k.f \
  cher.f   \
  cherk.f  \
  chpmv.f  \
  chpr2.f  \
  chpr.f   \
  crotg.f  \
  cscal.f  \
  csrot.f  \
  csscal.f \
  cswap.f  \
  csymm.f  \
  csyr2k.f \
  csyrk.f  \
  ctbmv.f  \
  ctbsv.f  \
  ctpmv.f  \
  ctpsv.f  \
  ctrmm.f  \
  ctrmv.f  \
  ctrsm.f  \
  ctrsv.f

d_source_files = \
  dasum.f  \
  daxpy.f  \
  dcabs1.f \
  dcopy.f  \
  ddot.f   \
  dgbmv.f  \
  dgemm.f  \
  dgemv.f  \
  dger.f   \
  dnrm2.f  \
  drot.f   \
  drotg.f  \
  drotm.f  \
  drotmg.f \
  dsbmv.f  \
  dscal.f  \
  dspmv.f  \
  dspr2.f  \
  dspr.f   \
  dswap.f  \
  dsymm.f  \
  dsymv.f  \
  dsyr2.f  \
  dsyr2k.f \
  dsyr.f   \
  dsyrk.f  \
  dtbmv.f  \
  dtbsv.f  \
  dtpmv.f  \
  dtpsv.f  \
  dtrmm.f  \
  dtrmv.f  \
  dtrsm.f  \
  dtrsv.f  \
  dzasum.f \
  dznrm2.f

s_source_files = \
  sasum.f  \
  saxpy.f  \
  scabs1.f \
  scasum.f \
  scnrm2.f \
  scopy.f  \
  sdot.f   \
  sdsdot.f \
  sgbmv.f  \
  sgemm.f  \
  sgemv.f  \
  sger.f   \
  snrm2.f  \
  srot.f   \
  srotg.f  \
  srotm.f  \
  srotmg.f \
  ssbmv.f  \
  sscal.f  \
  sspmv.f  \
  sspr2.f  \
  sspr.f   \
  sswap.f  \
  ssymm.f  \
  ssymv.f  \
  ssyr2.f  \
  ssyr2k.f \
  ssyr.f   \
  ssyrk.f  \
  stbmv.f  \
  stbsv.f  \
  stpmv.f  \
  stpsv.f  \
  strmm.f  \
  strmv.f  \
  strsm.f  \
  strsv.f

z_source_files = \
  zaxpy.f  \
  zcopy.f  \
  zdotc.f  \
  zdotu.f  \
  zdrot.f  \
  zdscal.f \
  zgbmv.f  \
  zgemm.f  \
  zgemv.f  \
  zgerc.f  \
  zgeru.f  \
  zhbmv.f  \
  zhemm.f  \
  zhemv.f  \
  zher2.f  \
  zher2k.f \
  zher.f   \
  zherk.f  \
  zhpmv.f  \
  zhpr2.f  \
  zhpr.f   \
  zrotg.f  \
  zscal.f  \
  zswap.f  \
  zsymm.f  \
  zsyr2k.f \
  zsyrk.f  \
  ztbmv.f  \
  ztbsv.f  \
  ztpmv.f  \
  ztpsv.f  \
  ztrmm.f  \
  ztrmv.f  \
  ztrsm.f  \
  ztrsv.f

EXTRA_DIST = \
  dsdot.f
