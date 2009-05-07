CLEANFILES = *~

AM_CFLAGS = -I$(top_srcdir)/src

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  matrix_allocation \
  get_element \
  set_element

get_element_SOURCES = get_element.c
get_element_LDADD = $(top_builddir)/src/liblal.la
#get_element_LDFLAGS = -static

set_element_SOURCES = set_element.c
set_element_LDADD = $(top_builddir)/src/liblal.la

matrix_allocation_SOURCES = matrix_allocation.c
matrix_allocation_LDADD = $(top_builddir)/src/liblal.la
