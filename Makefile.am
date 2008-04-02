AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CPP = @CPP@
AM_CPPFLAGS = @AM_CPPFLAGS@

SUBDIRS = MondoMods MondoSCF SCFeqs OneE QCTC HiCu TwoE ONX

datadir = ${prefix}/BasisSets

data_DATA = $(BasisSets)

include_HEADERS = MondoConfig.h config.h

EXTRA_DIST = \
  $(data_DATA) \
  $(examples) \
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

examples = \
  examples/AcO-.inp \
  examples/AcO-.ref \
  examples/Adenine-NumFreq.frequencies \
  examples/Adenine-NumFreq.out \
  examples/Adenine-Optimize.out \
  examples/Adenine-Optimize.xyz \
  examples/Adenine.inp \
  examples/Al.inp \
  examples/Al.ref \
  examples/BetaCarotene.inp \
  examples/BetaCarotene.ref \
  examples/Biotin.inp \
  examples/Biotin.ref \
  examples/BoronNitride.inp \
  examples/BoronNitride.ref \
  examples/CH2x-6.inp \
  examples/CH2x-6.out \
  examples/CH2x-6.ref \
  examples/CH2x-6.xyz \
  examples/Caffiene.inp \
  examples/Caffiene.ref \
  examples/Chlorophyll.inp \
  examples/Chlorophyll.ref \
  examples/Cytosine-NumFreq.frequencies \
  examples/Cytosine-NumFreq.out \
  examples/Cytosine-Optimize.out \
  examples/Cytosine-Optimize.xyz \
  examples/Cytosine.inp \
  examples/Diamond.inp \
  examples/Diamond.ref \
  examples/Graphite.inp \
  examples/Graphite.ref \
  examples/Guanine-NumFreq.frequencies \
  examples/Guanine-NumFreq.out \
  examples/Guanine-Optimize.out \
  examples/Guanine-Optimize.xyz \
  examples/Guanine.inp \
  examples/H2O.inp \
  examples/H2O.ref \
  examples/Histidine.inp \
  examples/Histidine.ref \
  examples/HybridMC_H2O.inp \
  examples/IceIh.inp \
  examples/IceIh.ref \
  examples/MDO2xyz.py \
  examples/MD_C2F4_Symplectic.inp \
  examples/MD_C2F4_TRBO.inp \
  examples/MD_F2_Symplectic.inp \
  examples/MD_F2_TRBO.inp \
  examples/MD_H2.inp \
  examples/MD_H2_TRBO.inp \
  examples/MD_H2_TRBO_2.inp \
  examples/MD_H2_TRBO_3.inp \
  examples/MD_Water10_T050_1SCF.inp \
  examples/Malonaldahyde.inp \
  examples/NH3.inp \
  examples/NH3.ref \
  examples/NO+.inp \
  examples/NO+.ref \
  examples/O2.inp \
  examples/O2.ref \
  examples/PIce.inp \
  examples/PIce.ref \
  examples/Peroxide.inp \
  examples/Peroxide.ref \
  examples/Quartz.inp \
  examples/Quartz.ref \
  examples/Steam.inp \
  examples/Steam.ref \
  examples/TEMPO.inp \
  examples/TEMPO.ref \
  examples/Vancomycin.inp \
  examples/Vancomycin.ref \
  examples/WaterHexamer.inp \
  examples/WaterHexamer.ref

CLEANFILES = *~

CONFIG_STATUS_DEPENDENCIES = localversion

# Fix localversion.
.PHONY : always-fix
localversion : always-fix
	@echo "checking localversion"; \
	  branch=`git branch 2>&1 | egrep '^\*' | awk '{print $$2}'`; \
	  localversion=`git log --pretty=oneline -1 2>&1 | awk '{print $$1}'`; \
	  if test -n "$${branch}" -a -n "$${localversion}"; then \
	    echo "branch=$${branch}" > localversion.temp; \
	    echo "localversion=$${localversion}" >> localversion.temp; \
	    if test -f localversion -a -z "`diff localversion localversion.temp`"; then \
	      echo "localversion did not change"; \
	      rm localversion.temp; \
	    else \
	      echo "localversion updated"; \
	      mv -f localversion.temp localversion; \
	      touch MondoSCF/PrintParsed.f90; \
	    fi; \
	    echo "localversion at $${branch}:$${localversion}"; \
	  fi

# Cryptographically sign the dist.
dist-sign : localversion dist
	rm -f $(distdir).tar.bz2.sig
	gpg -b $(distdir).tar.bz2

# Upload dist file to Savannah.
dist-upload : dist-sign
	rsync -v --progress $(distdir).tar.bz2 $(distdir).tar.bz2.sig nicolasbock@dl.sv.nongnu.org:/releases/freeon/
