AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = MondoMods LinearAlgebraLibrary FreeON SCFeqs OneE QCTC HiCu TwoE ONX Validate PostProcessing BasisSets

EXTRA_DIST = localversion FreeON.bibtex

CLEANFILES = *~

DISTCLEANFILES = config.h MondoConfig.h MondoLocalVersion.h

include_HEADERS = config.h MondoConfig.h MondoLocalVersion.h

#CONFIG_STATUS_DEPENDENCIES = $(top_srcdir)/localversion

DISTCHECK_CONFIGURE_FLAGS = --disable-optimizations --enable-internal-lapack-checks

# Phony targets.
.PHONY : validate

# Clean out everything except for things that are excrutiatingly long to
# rebuild.
quickclean :
	$(MAKE) -C FreeON clean
	$(MAKE) -C SCFeqs clean
	$(MAKE) -C OneE clean
	$(MAKE) -C HiCu clean
	$(MAKE) -C ONX clean

# Cryptographically sign the dist.
dist-sign : localversion dist
	rm -f $(distdir).tar.bz2.sig
	gpg -b $(distdir).tar.bz2

# Upload dist file to Savannah.
dist-upload : dist-sign
	if test -e freeon-latest.tar.bz2; then \
	  echo "freeon-latest.tar.bz2 already exists"; \
	else \
	  echo "sending dist to savannah"; \
	  ln -s $(distdir).tar.bz2 freeon-latest.tar.bz2; \
	  rsync -av --progress $(distdir).tar.bz2 freeon-latest.tar.bz2 $(distdir).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/; \
	  rm -f freeon-latest.tar.bz2; \
	fi

# Download all files in Downloads area.
download :
	if test ! -d Downloads; then \
	  mkdir Downloads; \
	fi
	rsync -ave ssh --progress --partial --delete nicolasbock@dl.sv.nongnu.org:/releases/freeon/ Downloads/

# Upload all files in Downloads area.
upload :
	if test ! -d Downloads; then \
	  echo "there is no Downloads folder"; \
	  exit 1; \
	fi
	rsync -ave ssh --progress --partial --delete Downloads/ nicolasbock@dl.sv.nongnu.org:/releases/freeon/

# Validate.
validate :
	$(MAKE) -C Validate validate

TESTS = /bin/true

# dependency checks.
depcheck :
	$(MAKE) -C MondoMods depcheck-lib
	$(MAKE) -C FreeON  depcheck-bin
	$(MAKE) -C SCFeqs  depcheck-bin
	$(MAKE) -C OneE  depcheck-bin
	$(MAKE) -C QCTC  depcheck-bin
	$(MAKE) -C HiCu  depcheck-bin
	$(MAKE) -C ONX  depcheck-bin
