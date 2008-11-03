# The distfile name.
distname = PostProcessing

# The files to distribute.
DISTFILES = \
  COPYING \
  README \
  diffusionCoefficient.py \
  parseMDO.py

# Make dist.
dist :
	tar cf - $(DISTFILES) | bzip2 > $(distname).tar.bz2

# Cryptographically sign the dist.
dist-sign : dist
	rm -f $(distname).tar.bz2.sig
	gpg -b $(distname).tar.bz2

# Upload dist file to Savannah.
upload-dist : dist-sign
	rsync -v --progress $(distname).tar.bz2 $(distname).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/

.PHONY : clean
clean :
	rm *~
