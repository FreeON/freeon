SUBDIRS = master slave_1 slave_2 slave_3

TESTS = test_1.sh

TESTS_ENVIRONMENT = PATH="${builddir}/master:${builddir}/slave_1:${builddir}/slave_2:${builddir}/slave_3:${PATH}"

noinst_SCRIPTS = $(TESTS)

EXTRA_DIST = $(TESTS)

install-exec-hook :
	@echo "running hook"
