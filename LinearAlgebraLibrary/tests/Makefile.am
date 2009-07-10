CLEANFILES = *~

AM_CFLAGS = -I$(top_srcdir)/src

TESTS = $(check_PROGRAMS)

check_PROGRAMS = \
  matrix_allocation \
  get_element \
  set_element \
  matrix_equals \
  matrix_multiply

LDADD = $(top_builddir)/src/liblal.la
