new: 	clean links dvi

latex:
	latex PBCQUICCA.tex
	bibtex PBCQUICCA
	latex PBCQUICCA.tex
	bibtex PBCQUICCA 
	latex PBCQUICCA.tex

links:
	cp mondo_new.bib-save mondo_new.bib
#	ln -s ../../Bib/mondo_new.bib ./mondo_new.bib 
	ln -s ../../RevTex/revtex4.cls ./revtex4.cls
	ln -s ../../RevTex/apsrev.bst  ./apsrev.bst
	ln -s ../../RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s ../../RevTex/revsymb.sty ./revsymb.sty
	ln -s ../../RevTex/aps.rtx ./aps.rtx
	ln -s ../../RevTex/10pt.rtx ./10pt.rtx
	ln -s ../../RevTex/11pt.rtx ./11pt.rtx
	ln -s ../../RevTex/12pt.rtx ./12pt.rtx

ps:	latex
	dvips -t letter -o PBCQUICCA.ps PBCQUICCA

pdf:	clean links
	chmod u+x ./pdfconverter
	./pdfconverter
	pdflatex PBCQUICCA

dvi:	latex
	xdvi -s 2 PBCQUICCA

clean:
	rm -rf *.log *.bak *.bib *.blg *.dvi *~ figure*.pdf
	rm -rf revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
