include $(top_srcdir)/make.inc

lib_LIBRARIES = libMondoTMG.a

SCATGEN = slatm1.f slaran.f slarnd.f

SMATGEN = slatms.f slatme.f slatmr.f \
   slagge.f slagsy.f slakf2.f slarge.f slaror.f slarot.f slatm2.f \
   slatm3.f slatm5.f slatm6.f

CMATGEN = clatms.f clatme.f clatmr.f \
   clagge.f claghe.f clagsy.f clakf2.f clarge.f claror.f clarot.f \
   clatm1.f clarnd.f clatm2.f clatm3.f clatm5.f clatm6.f

DZATGEN = dlatm1.f dlaran.f dlarnd.f

DMATGEN = dlatms.f dlatme.f dlatmr.f \
   dlagge.f dlagsy.f dlakf2.f dlarge.f dlaror.f dlarot.f dlatm2.f \
   dlatm3.f dlatm5.f dlatm6.f

ZMATGEN = zlatms.f zlatme.f zlatmr.f \
   zlagge.f zlaghe.f zlagsy.f zlakf2.f zlarge.f zlaror.f zlarot.f \
   zlatm1.f zlarnd.f zlatm2.f zlatm3.f zlatm5.f zlatm6.f

libMondoTMG_a_SOURCES = $(SMATGEN) $(CMATGEN) $(SCATGEN) $(DMATGEN)	$(ZMATGEN) $(DZATGEN)

slaran.o: slaran.f
	@echo "compiling $< with no-optimizations"
	@$(F77) @NOOPT@ -c -o $@ $< && exit 0; \
	  echo "The following comand failed:" 1>&2; \
		echo "$(F77) @NOOPT@ -c -o $@ $<" 1>&2; \
		$(F77) @NOOPT@ -c -o $@ $< > /dev/null 2>&1

dlaran.o: dlaran.f
	@echo "compiling $< with no-optimizations"
	@$(F77) @NOOPT@ -c -o $@ $< && exit 0; \
	  echo "The following comand failed:" 1>&2; \
		echo "$(F77) @NOOPT@ -c -o $@ $<" 1>&2; \
		$(F77) @NOOPT@ -c -o $@ $< > /dev/null 2>&1
