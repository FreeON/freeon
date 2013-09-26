ACLOCAL_AMFLAGS = -I m4

CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = \
	Modules \
	tests \
	LinearAlgebraLibrary \
	FreeON \
	SCFeqs \
	OneE \
	QCTC \
	HiCu \
	TwoE \
	ONX \
	Validate \
	PostProcessing \
	BasisSets

EXTRA_DIST = localversion FreeON.bibtex

CLEANFILES = *~

DISTCLEANFILES = MondoConfig.h MondoLocalVersion.h

include_HEADERS = MondoLocalVersion.h

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_srcdir)/localversion

docs :
	doxygen
	$(MAKE) -C latex

#DISTCHECK_CONFIGURE_FLAGS = --disable-optimizations --enable-internal-lapack-checks

# Phony targets.
.PHONY : validate tags

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
	$(MAKE) -C Modules depcheck-lib
	$(MAKE) -C FreeON  depcheck-bin
	$(MAKE) -C SCFeqs  depcheck-bin
	$(MAKE) -C OneE  depcheck-bin
	$(MAKE) -C QCTC  depcheck-bin
	$(MAKE) -C HiCu  depcheck-bin
	$(MAKE) -C ONX  depcheck-bin

$(top_builddir)/MondoLocalVersion.h : $(top_builddir)/config_localversion.h
	@echo "checking MondoLocalVersion.h"
	@grep LOCAL_VERSION $(top_builddir)/config_localversion.h > $(top_builddir)/tempconf
	@if test ! -f $(top_builddir)/MondoLocalVersion.h; then \
	  mv -f $(top_builddir)/tempconf $(top_builddir)/MondoLocalVersion.h; \
	  echo "MondoLocalVersion.h did not exist"; \
	elif (diff -u $(top_builddir)/MondoLocalVersion.h $(top_builddir)/tempconf > $(top_builddir)/MondoLocalVersion.diff); then \
	  rm -f $(top_builddir)/tempconf; \
	  echo "MondoLocalVersion.h unchanged"; \
	else \
	  mv -f $(top_builddir)/tempconf $(top_builddir)/MondoLocalVersion.h; \
	  echo "MondoLocalVersion.h updated"; \
	  cat $(top_builddir)/MondoLocalVersion.diff; \
	fi

# In case we have ctags installed, make the necessary tags files.
tags :
	$(MAKE) -C Modules tags
	$(MAKE) -C FreeON tags
	$(MAKE) -C SCFeqs tags
	$(MAKE) -C OneE tags
	$(MAKE) -C QCTC tags
