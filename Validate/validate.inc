# vim: syntax=make:noexpandtab:tw=0
#
# Common rule to all makefiles here.

TESTS_ENVIRONMENT = \
  MONDO_HOME=$(top_srcdir) \
	PATH=$(top_builddir)/FreeON:$(top_builddir)/SCFeqs:$(top_builddir)/OneE:$(top_builddir)/QCTC:$(top_builddir)/HiCu:$(top_builddir)/ONX:${PATH}
