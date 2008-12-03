AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = MondoMods FreeON SCFeqs OneE QCTC HiCu TwoE ONX Validate PostProcessing BasisSets

EXTRA_DIST = localversion

CLEANFILES = *~

DISTCLEANFILES = MondoConfig.h MondoLocalVersion.h

CONFIG_STATUS_DEPENDENCIES = $(top_srcdir)/localversion

DISTCHECK_CONFIGURE_FLAGS = --disable-optimizations --enable-internal-lapack-checks

# Phony targets.
.PHONY : always-fix validate

# Fix localversion.
$(top_srcdir)/localversion : always-fix
	@echo "checking localversion"; \
	  branch=`git branch 2>&1 | egrep '^\*' | awk '{print $$2}'`; \
	  localversion=`git log --pretty=oneline -1 2>&1 | awk '{print $$1}'`; \
	  if test -n "$${branch}" -a -n "$${localversion}"; then \
	    echo "branch=$${branch}" > localversion.temp; \
	    echo "localversion=$${localversion}" >> localversion.temp; \
	    if test -f ${top_srcdir}/localversion; then \
	      if test -z "`diff ${top_srcdir}/localversion localversion.temp`"; then \
	        echo "localversion did not change"; \
	        rm localversion.temp; \
	      else \
	        echo "localversion updated"; \
	        diff -u ${top_srcdir}/localversion localversion.temp; \
	        mv -f localversion.temp ${top_srcdir}/localversion; \
	      fi; \
	    else \
	      echo "localversion did not exist"; \
	      mv -f localversion.temp ${top_srcdir}/localversion; \
	    fi; \
	      rm -f localversion.temp; \
	    echo "localversion at $${branch}:$${localversion}"; \
	  fi

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

# Validate.
validate :
	$(MAKE) -C Validate validate

TESTS = /bin/true
