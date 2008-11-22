# vim: syntax=make:noexpandtab:tw=0
#
# Common rule to all makefiles here.

.PHONY : validate validate-recursive

validate : validate-recursive

validate-recursive :
	@for subdir in $(SUBDIRS); do \
		echo "Making validate in $$subdir"; \
		$(MAKE) -C $$subdir validate; \
	done
