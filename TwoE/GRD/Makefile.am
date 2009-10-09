include $(top_srcdir)/make.inc

noinst_LTLIBRARIES = libGRD.la

libGRD_la_SOURCES = \
  BraHRR1010.f90 \
  BraHRR101.f90 \
  BraHRR103.f90 \
  BraHRR106.f90 \
  BraHRR110.f90 \
  BraHRR11.f90 \
  BraHRR13.f90 \
  BraHRR16.f90 \
  BraHRR310.f90 \
  BraHRR31.f90 \
  BraHRR33.f90 \
  BraHRR36.f90 \
  BraHRR610.f90 \
  BraHRR61.f90 \
  BraHRR63.f90 \
  BraHRR66.f90 \
  dIntB10010101.f90 \
  dIntB10010301.f90 \
  dIntB10010303.f90 \
  dIntB10010601.f90 \
  dIntB10010603.f90 \
  dIntB10010606.f90 \
  dIntB10011001.f90 \
  dIntB10030101.f90 \
  dIntB10030301.f90 \
  dIntB10030303.f90 \
  dIntB10030601.f90 \
  dIntB10030603.f90 \
  dIntB10030606.f90 \
  dIntB10031001.f90 \
  dIntB10031003.f90 \
  dIntB10060101.f90 \
  dIntB10060301.f90 \
  dIntB10060303.f90 \
  dIntB10060601.f90 \
  dIntB10060603.f90 \
  dIntB10060606.f90 \
  dIntB10061001.f90 \
  dIntB10061003.f90 \
  dIntB10061006.f90 \
  dIntB10100101.f90 \
  dIntB10100301.f90 \
  dIntB10100303.f90 \
  dIntB10100601.f90 \
  dIntB10100603.f90 \
  dIntB10100606.f90 \
  dIntB10101001.f90 \
  dIntB10101003.f90 \
  dIntB10101006.f90 \
  dIntB10101010.f90 \
  dIntB1010101.f90 \
  dIntB3010101.f90 \
  dIntB3010301.f90 \
  dIntB3030101.f90 \
  dIntB3030301.f90 \
  dIntB3030303.f90 \
  dIntB6010101.f90 \
  dIntB6010301.f90 \
  dIntB6010303.f90 \
  dIntB6010601.f90 \
  dIntB6030101.f90 \
  dIntB6030301.f90 \
  dIntB6030303.f90 \
  dIntB6030601.f90 \
  dIntB6030603.f90 \
  dIntB6060101.f90 \
  dIntB6060301.f90 \
  dIntB6060303.f90 \
  dIntB6060601.f90 \
  dIntB6060603.f90 \
  dIntB6060606.f90

EXTRA_DIST = \
  dERIInterfaceB.Inc

EXTRA_INCLUDES = -I$(builddir)/../MSC

# Some compilers have trouble optimizing some of the sources here. For now we
# will turn the optimization level down...
if USE_COMPILER_OPTIMIZATIONS
OPTIMIZATION_FLAGS = -O1
endif
