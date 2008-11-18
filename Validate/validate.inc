# vim: syntax=make:noexpandtab:tw=0
#
# Common rule to all makefiles here.

.PHONY : validate
validate : $(VALIDATE_TESTS)
	@echo "running test $@"
	for i in $^; do \
	  $(top_srcdir)/Validate/runTests.py $$i; \
	done
