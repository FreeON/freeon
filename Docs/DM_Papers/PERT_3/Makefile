all:	latex #dvi

new: 	clean links dvi

links:
	ln -s ../../Bib/mondo_new.bib ./mondo_new.bib 
	ln -s ../../RevTex/revtex4.cls ./revtex4.cls
	ln -s ../../RevTex/apsrev.bst  ./apsrev.bst
	ln -s ../../RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s ../../RevTex/revsymb.sty ./revsymb.sty
	ln -s ../../RevTex/aps.rtx ./aps.rtx
	ln -s ../../RevTex/10pt.rtx ./10pt.rtx
	ln -s ../../RevTex/11pt.rtx ./11pt.rtx
	ln -s ../../RevTex/12pt.rtx ./12pt.rtx
	ln -s ./Gamma_h2o3D_6-31G_6-31Gss_G_T.ps ./figure1.eps
	ln -s ./Gamma_QCTC_T.ps ./figure2.eps
	ln -s Gamma_TC2R_T.ps ./figure3.eps
	ln -s Gamma_ONX_T.ps ./figure4.eps
	ln -s Mix_h2o3D_6-31G_T.ps ./figure5.eps
	ln -s DMix_40pts.ps ./figure6.eps

ps:	latex
	dvips  -o Response3.ps  Response3

pdf:	
	chmod u+x ./pdfconverter
	./pdfconverter
	pdflatex Response3 

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
	rm -rf figure1.eps figure2.eps figure3.eps figure4.eps figure5.eps figure6.eps 
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo.bib mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
