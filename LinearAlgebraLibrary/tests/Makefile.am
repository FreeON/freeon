CLEANFILES = *~

AM_CFLAGS = -I$(top_srcdir)/src

TESTS = $(noinst_PROGRAMS)

noinst_PROGRAMS = \
  get_element \
  matrix_allocation

get_element_SOURCES = get_element.c
get_element_LDADD = $(top_builddir)/src/libmal.la
#get_element_LDFLAGS = -static

matrix_allocation_SOURCES = matrix_allocation.c
matrix_allocation_LDADD = $(top_builddir)/src/libmal.la
