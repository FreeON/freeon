all: normal

normal:
	latex DSEE.tex
	latex DSEE.tex
	xdvi -s 2 DSEE

new: 	clean links figs bib ps

bib:	
	latex DSEE.tex
	bibtex DSEE
	latex DSEE.tex
	bibtex DSEE
	latex DSEE.tex	
	latex DSEE.tex	

links:
	ln -s ../../Bib/mondo.bib ./mondo.bib 
	ln -s ../../RevTex/revtex4.cls ./revtex4.cls
	ln -s ../../RevTex/apsrev.bst  ./apsrev.bst
	ln -s ../../RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s ../../RevTex/revsymb.sty ./revsymb.sty
	ln -s ../../RevTex/aps.rtx ./aps.rtx
	ln -s ../../RevTex/10pt.rtx ./10pt.rtx
	ln -s ../../RevTex/11pt.rtx ./11pt.rtx
	ln -s ../../RevTex/12pt.rtx ./12pt.rtx

figs:
	make -C XE_FIGS
	make -C WATER_FIGS
	make -C FULL_FIGS

ps:
	dvips -o DSEE.ps DSEE

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
