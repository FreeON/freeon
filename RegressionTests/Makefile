# The distfile name.
distname = MondoRegressionTests

# The files to distribute.
DISTFILES = \
  COPYING \
  MD_H2_TRBO.inp \
  MD_H2_TRBO.test \
  Makefile \
  README \
  build.test \
  example.reference \
  example.test \
  localversion \
  runTests.py

# Fix localversion.
fix-localversion :
	echo "localversion=`git log --pretty=oneline -1 | awk '{print $$1}'`" > localversion

# Make dist.
dist : fix-localversion
	tar cf - $(DISTFILES) | bzip2 > $(distname).tar.bz2

# Cryptographically sign the dist.
dist-sign : dist
	rm -f $(distname).tar.bz2.sig
	gpg -b $(distname).tar.bz2

# Upload dist file to Savannah.
upload-dist : dist-sign
	rsync -v --progress $(distname).tar.bz2 $(distname).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/
