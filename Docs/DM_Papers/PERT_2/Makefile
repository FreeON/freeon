all:	latex dvi

new: 	clean links dvi

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

ps:	latex
	dvips -o prl2.ps perturbedprojector

dvi:	latex
	xdvi -s 2 perturbedprojector

latex:
	latex perturbedprojector
	bibtex perturbedprojector
	latex perturbedprojector
	bibtex perturbedprojector
	latex perturbedprojector

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
