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
	ln -s ./Gamma_h2o3D_6-31G_6-31Gss_G_T.ps ./figure1.ps
	ln -s ./Gamma_QCTC_T.ps ./figure2.ps
	ln -s Gamma_TC2R_T.ps ./figure3.ps
	ln -s Gamma_ONX_T.ps ./figure4.ps
	ln -s Mix_h2o3D_6-31G_T.ps ./figure5.ps
	ln -s DMix_40pts.ps ./figure6.ps

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
