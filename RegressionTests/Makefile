# The version.
VERSION = 0.1

# The distfile name.
distname = MondoRegressionTests

# The files to distribute.
DISTFILES = \
  COPYING \
  MD_F2_Symplectic_200_STO-2G-SPLIT.inp \
  MD_F2_Symplectic_200_STO-2G-SPLIT.reference \
  MD_F2_Symplectic_200_STO-2G-SPLIT.test \
  MD_F2_Symplectic_50_STO-2G-SPLIT.inp \
  MD_F2_Symplectic_50_STO-2G-SPLIT.reference \
  MD_F2_Symplectic_50_STO-2G-SPLIT.test \
  MD_F2_TRBO_125_STO-2G-SPLIT.inp \
  MD_F2_TRBO_125_STO-2G-SPLIT.reference \
  MD_F2_TRBO_125_STO-2G-SPLIT.test \
  MD_F2_TRBO_50_STO-2G-SPLIT.inp \
  MD_F2_TRBO_50_STO-2G-SPLIT.reference \
  MD_F2_TRBO_50_STO-2G-SPLIT.test \
  MD_H2_TRBO.inp \
  MD_H2_TRBO.reference \
  MD_H2_TRBO.test \
  README \
  build.test \
  example.reference \
  example.test \
  runTests.py

TESTS = \
  build.test \
  MD_H2_TRBO.test \
  MD_F2_Symplectic_200_STO-2G-SPLIT.test \
  MD_F2_Symplectic_50_STO-2G-SPLIT.test \
  MD_F2_TRBO_125_STO-2G-SPLIT.test \
  MD_F2_TRBO_50_STO-2G-SPLIT.test

# Run the tests.
.PHONY : check FORCE
check : $(TESTS)

# Generic rule for the tests.
%.test : FORCE
	@echo "running test $@"
	./runTests.py $@

# Make dist.
dist : clean
	mkdir $(distname)-$(VERSION)
	cp $(DISTFILES) $(distname)-$(VERSION)
	tar cf - $(distname)-$(VERSION) | bzip2 > $(distname)-$(VERSION).tar.bz2
	rm -rf $(distname)-$(VERSION)

# Cryptographically sign the dist.
dist-sign : dist
	rm -f $(distname)-$(VERSION).tar.bz2.sig
	gpg -b $(distname)-$(VERSION).tar.bz2

# Upload dist file to Savannah.
dist-upload : dist-sign
	rsync -v --progress $(distname)-$(VERSION).tar.bz2 $(distname)-$(VERSION).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/

# Clean up.
clean :
	rm -f *~
	rm -f $(distname)-$(VERSION).tar.bz2
	rm -f $(distname)-$(VERSION).tar.bz2.sig
