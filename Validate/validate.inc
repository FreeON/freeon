# vim: syntax=make:noexpandtab:tw=0
#
# Common rule to all makefiles here.

.PHONY : validate

TESTS_ENVIRONMENT = \
  MONDO_HOME=$(top_srcdir) \
	PATH=$(top_builddir)/FreeON:$(top_builddir)/SCFeqs:$(top_builddir)/OneE:$(top_builddir)/QCTC:$(top_builddir)/HiCu:$(top_builddir)/ONX:${PATH}

VALIDATE_TESTS = $(foreach i, $(VALIDATE_INPUTS), $(basename $i).test)
VALIDATE_REFERENCES = $(foreach i, $(VALIDATE_INPUTS), $(basename $i).reference)

EXTRA_DIST = \
	$(VALIDATE_INPUTS) \
  $(VALIDATE_REFERENCES) \
  $(VALIDATE_TESTS)

validate : $(VALIDATE_TESTS)
	@echo "running test $@"
	@for i in $^; do \
	  $(MAKE) $$i; \
	done

%.test :
	@echo "running suffix rule on $@"
	$(top_srcdir)/Validate/runTests.py $@
