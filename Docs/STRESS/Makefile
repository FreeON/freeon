all:	latex #dvi

new: 	clean links dvi

links:
	ln -s $(MONDO_HOME)/Docs/Bib/mondo_new.bib ./mondo_new.bib 
	ln -s $(REVTEX4_HOME)/revtex4.cls ./revtex4.cls
	ln -s $(REVTEX4_HOME)/apsrev.bst  ./apsrev.bst
	ln -s $(REVTEX4_HOME)/apsrmp.bst  ./apsrmp.bst
	ln -s $(REVTEX4_HOME)/revsymb.sty ./revsymb.sty
	ln -s $(REVTEX4_HOME)/aps.rtx ./aps.rtx
	ln -s $(REVTEX4_HOME)/10pt.rtx ./10pt.rtx
	ln -s $(REVTEX4_HOME)/11pt.rtx ./11pt.rtx
	ln -s $(REVTEX4_HOME)/12pt.rtx ./12pt.rtx
#	ln -s $(MONDO_HOME)/Docs/RevTex/revtex4.cls ./revtex4.cls
#	ln -s $(MONDO_HOME)/Docs/RevTex/apsrev.bst  ./apsrev.bst
#	ln -s $(MONDO_HOME)/Docs/RevTex/apsrmp.bst  ./apsrmp.bst
#	ln -s $(MONDO_HOME)/Docs/RevTex/revsymb.sty ./revsymb.sty
#	ln -s $(MONDO_HOME)/Docs/RevTex/aps.rtx ./aps.rtx
#	ln -s $(MONDO_HOME)/Docs/RevTex/10pt.rtx ./10pt.rtx
#	ln -s $(MONDO_HOME)/Docs/RevTex/11pt.rtx ./11pt.rtx
#	ln -s $(MONDO_HOME)/Docs/RevTex/12pt.rtx ./12pt.rtx

ps:	latex
	dvips -t a4 -o Stress.ps Stress

pdf:	ps
	ps2pdf Stress.ps Stress.pdf


dvi:	latex
	xdvi -s 3 Stress &

latex:
	latex Stress
	bibtex Stress
	latex Stress
	bibtex Stress
	latex Stress

spell:
	aspell --mode=tex -c Stress.tex

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
