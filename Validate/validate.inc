# vim: syntax=make:tw=0
#
# Common rule to all makefiles here.

TESTS_ENVIRONMENT = \
  FREEON_HOME=$(top_builddir) \
  FREEON_BASISSETS=$(top_srcdir)/BasisSets \
  FREEON_SCRATCH=$(builddir) \
	PATH=$(top_builddir)/FreeON:$(top_builddir)/SCFeqs:$(top_builddir)/OneE:$(top_builddir)/QCTC:$(top_builddir)/HiCu:$(top_builddir)/ONX:${PATH}

DISTCLEANFILES = *.out
