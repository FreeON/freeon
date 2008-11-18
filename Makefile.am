AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = MondoMods FreeON SCFeqs OneE QCTC HiCu TwoE ONX Validate PostProcessing

datadir = ${prefix}/BasisSets

data_DATA = $(BasisSets)

EXTRA_DIST = \
  $(data_DATA) \
  localversion

BasisSets = \
  BasisSets/3-21G-SPLIT.bas \
  BasisSets/3-21G.bas \
  BasisSets/3-21G.bas-SPLIT \
  BasisSets/3-21GSP-SPLIT.bas \
  BasisSets/3-21GSP.bas \
  BasisSets/3-21Gs-SPLIT.bas \
  BasisSets/3-21Gs.bas \
  BasisSets/4-22GSP-SPLIT.bas \
  BasisSets/4-22GSP.bas \
  BasisSets/4-31Gs-SPLIT.bas \
  BasisSets/4-31Gs.bas \
  BasisSets/6-311G-SPLIT.bas \
  BasisSets/6-311G.bas \
  BasisSets/6-311G3df,3pd-SPLIT.bas \
  BasisSets/6-311G3df,3pd.bas \
  BasisSets/6-311Gss-SPLIT.bas \
  BasisSets/6-311Gss.bas \
  BasisSets/6-311ppG3df,3pd-SPLIT.bas \
  BasisSets/6-311ppG3df,3pd.bas \
  BasisSets/6-31G-SPLIT.bas \
  BasisSets/6-31G.bas \
  BasisSets/6-31Gs-SPLIT.bas \
  BasisSets/6-31Gs.bas \
  BasisSets/6-31Gss-SPLIT.bas \
  BasisSets/6-31Gss.bas \
  BasisSets/6-31Gss.bas-SPLIT \
  BasisSets/6-31ppGss-SPLIT.bas \
  BasisSets/6-31ppGss.bas \
  BasisSets/Ahlrichs_TZV-SPLIT.bas \
  BasisSets/Ahlrichs_TZV.bas \
  BasisSets/Ahlrichs_VDZ-SPLIT.bas \
  BasisSets/Ahlrichs_VDZ.bas \
  BasisSets/Ahlrichs_VTZ-SPLIT.bas \
  BasisSets/Ahlrichs_VTZ.bas \
  BasisSets/Ahlrichs_pVDZ-SPLIT.bas \
  BasisSets/Ahlrichs_pVDZ.bas \
  BasisSets/CRENBL_ECP-SPLIT.bas \
  BasisSets/CRENBL_ECP.bas \
  BasisSets/CRENBS_ECP-SPLIT.bas \
  BasisSets/CRENBS_ECP.bas \
  BasisSets/C_C98-SPLIT.bas \
  BasisSets/C_C98.bas \
  BasisSets/Crystal98-SPLIT.bas \
  BasisSets/Crystal98.bas \
  BasisSets/DunningHay_SV-SPLIT.bas \
  BasisSets/DunningHay_SV.bas \
  BasisSets/DunningHay_SVP-SPLIT.bas \
  BasisSets/DunningHay_SVP.bas \
  BasisSets/Dunning_DZ-SPLIT.bas \
  BasisSets/Dunning_DZ.bas \
  BasisSets/Dunning_DZP-SPLIT.bas \
  BasisSets/Dunning_DZP.bas \
  BasisSets/Dunning_TZ-SPLIT.bas \
  BasisSets/Dunning_TZ.bas \
  BasisSets/EvenTmp_PBC-SPLIT.bas \
  BasisSets/EvenTmp_PBC.bas \
  BasisSets/Ghost-SPLIT.bas \
  BasisSets/Ghost.bas \
  BasisSets/LANL2_DZ_ECP-SPLIT.bas \
  BasisSets/LANL2_DZ_ECP.bas \
  BasisSets/MIDI-SPLIT.bas \
  BasisSets/MIDI.bas \
  BasisSets/MINI-SPLIT.bas \
  BasisSets/MINI.bas \
  BasisSets/MgO_C98-SPLIT.bas \
  BasisSets/MgO_C98.bas \
  BasisSets/README \
  BasisSets/SBKJC_VDZ_ECP-SPLIT.bas \
  BasisSets/SBKJC_VDZ_ECP.bas \
  BasisSets/STO-2G-SPLIT.bas \
  BasisSets/STO-2G.bas \
  BasisSets/STO-2G.bas-SPLIT \
  BasisSets/STO-3G-SPLIT.bas \
  BasisSets/STO-3G.bas \
  BasisSets/STO-3G.bas-SPLIT \
  BasisSets/STO-6G-SPLIT.bas \
  BasisSets/STO-6G.bas \
  BasisSets/SrTiOAllE-SPLIT.bas \
  BasisSets/SrTiOAllE.bas \
  BasisSets/Stuttgart_RLC_ECP-SPLIT.bas \
  BasisSets/Stuttgart_RLC_ECP.bas \
  BasisSets/Stuttgart_RSC97_ECP-SPLIT.bas \
  BasisSets/Stuttgart_RSC97_ECP.bas \
  BasisSets/TiO2_C98-SPLIT.bas \
  BasisSets/TiO2_C98.bas \
  BasisSets/User1.bas \
  BasisSets/User2.bas \
  BasisSets/User4.bas \
  BasisSets/User5.bas \
  BasisSets/WTBS-SPLIT.bas \
  BasisSets/WTBS.bas \
  BasisSets/splitbas

CLEANFILES = *~

DISTCLEANFILES = MondoConfig.h MondoLocalVersion.h

CONFIG_STATUS_DEPENDENCIES = $(top_srcdir)/localversion

DISTCHECK_CONFIGURE_FLAGS = --disable-optimizations

# Fix localversion.
.PHONY : always-fix
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
	  rsync -lv --progress $(distdir).tar.bz2 freeon-latest.tar.bz2 $(distdir).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/; \
	  rm -f freeon-latest.tar.bz2; \
	fi
