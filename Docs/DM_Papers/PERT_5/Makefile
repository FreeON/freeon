all:	latex #dvi

new: 	clean links dvi

links:
	ln -s $(MONDO_HOME)/Docs/Bib/mondo_new.bib ./mondo_new.bib 
	ln -s $(MONDO_HOME)/Docs/RevTex/revtex4.cls ./revtex4.cls
	ln -s $(MONDO_HOME)/Docs/RevTex/apsrev.bst  ./apsrev.bst
	ln -s $(MONDO_HOME)/Docs/RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s $(MONDO_HOME)/Docs/RevTex/revsymb.sty ./revsymb.sty
	ln -s $(MONDO_HOME)/Docs/RevTex/aps.rtx ./aps.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/10pt.rtx ./10pt.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/11pt.rtx ./11pt.rtx
	ln -s $(MONDO_HOME)/Docs/RevTex/12pt.rtx ./12pt.rtx

ps:	latex
	dvips  -o DMX.ps DMX

pdf:	pdflatex DMX


dvi:	latex
	xdvi -s 2 DMX &

latex:
	latex DMX
	bibtex DMX
	latex DMX
	bibtex DMX
	latex DMX

spell:
	aspell --mode=tex -c DMX.tex

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
