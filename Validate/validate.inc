# vim: syntax=make:noexpandtab:tw=0
#
# Common rule to all makefiles here.

TESTS_ENVIRONMENT = \
  FREEON_HOME=$(top_srcdir) \
  FREEON_BASISSETS=$(top_srcdir)/BasisSets \
	PATH=$(top_builddir)/FreeON:$(top_builddir)/SCFeqs:$(top_builddir)/OneE:$(top_builddir)/QCTC:$(top_builddir)/HiCu:$(top_builddir)/ONX:${PATH}

DISTCLEANFILES = *.out
