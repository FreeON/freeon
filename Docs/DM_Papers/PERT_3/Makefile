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
	dvips  -o jcp3.ps Response3

pdf:	ps
	ps2pdf jcp3.ps jcp3.pdf


dvi:	latex
	xdvi  -s 2 Response3 &

latex:
	latex Response3
	bibtex Response3
	latex Response3
	bibtex Response3
	latex Response3

spell:
	aspell --mode=tex -c Response3.tex

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
