all:	bib dvi

latex:
	latex MD_p1.tex
	latex MD_p1.tex

new: 	clean links bib dvi

bib:	
	latex MD_p1.tex
	bibtex MD_p1
	latex MD_p1.tex
	bibtex MD_p1
	latex MD_p1.tex	
#	latex MD_p1.tex	

links:
	ln -s ../Bib/mondo_new.bib ./mondo_new.bib 
	ln -s ../RevTex/revtex4.cls ./revtex4.cls
	ln -s ../RevTex/apsrev.bst  ./apsrev.bst
	ln -s ../RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s ../RevTex/revsymb.sty ./revsymb.sty
	ln -s ../RevTex/aps.rtx ./aps.rtx
	ln -s ../RevTex/10pt.rtx ./10pt.rtx
	ln -s ../RevTex/11pt.rtx ./11pt.rtx
	ln -s ../RevTex/12pt.rtx ./12pt.rtx


ps:
	dvips -o MD_p1.ps MD_p1

dvi:
	xdvi -s 2 MD_p1

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
