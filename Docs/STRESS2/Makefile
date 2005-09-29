all:	latex #dvi

new: 	clean links dvi

links:
	ln -s $(MONDO_HOME)/Docs/Bib/mondo_new.bib ./mondo_new.bib 
#	ln -s $(REVTEX4_HOME)/revtex4.cls ./revtex4.cls
#	ln -s $(REVTEX4_HOME)/apsrev.bst  ./apsrev.bst
#	ln -s $(REVTEX4_HOME)/apsrmp.bst  ./apsrmp.bst
#	ln -s $(REVTEX4_HOME)/revsymb.sty ./revsymb.sty
#	ln -s $(REVTEX4_HOME)/aps.rtx ./aps.rtx
#	ln -s $(REVTEX4_HOME)/10pt.rtx ./10pt.rtx
#	ln -s $(REVTEX4_HOME)/11pt.rtx ./11pt.rtx
#	ln -s $(REVTEX4_HOME)/12pt.rtx ./12pt.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/revtex4.cls ./revtex4.cls
	ln -s $(MONDO_HOME)/Docs/RevTex/apsrev.bst  ./apsrev.bst
	ln -s $(MONDO_HOME)/Docs/RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s $(MONDO_HOME)/Docs/RevTex/revsymb.sty ./revsymb.sty
	ln -s $(MONDO_HOME)/Docs/RevTex/aps.rtx ./aps.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/10pt.rtx ./10pt.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/11pt.rtx ./11pt.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/12pt.rtx ./12pt.rtx

ps:	latex
#	dvips -t a4 -o Stress2.ps Stress2
	dvips -t letter -o Stress2.ps Stress2

pdf:	ps
	ps2pdf Stress2.ps Stress2.pdf


dvi:	latex
	xdvi -s 3 Stress2 &

latex:
	latex Stress2
	bibtex Stress2
	latex Stress2
	bibtex Stress2
	latex Stress2

spell:
	aspell --mode=tex -c Stress2.tex

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
