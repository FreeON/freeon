all:	latex dvi

latex:
	latex PBCQUICCA.tex
	bibtex PBCQUICCA
	latex PBCQUICCA.tex
	bibtex PBCQUICCA 
	latex PBCQUICCA.tex

new: 	clean links dvi

links:
	ln -s /n/srv3/nemeth/MONDO_TEST/Docs/RevTex/revtex4.cls ./revtex4.cls
	ln -s ../../RevTex/apsrev.bst  ./apsrev.bst
	ln -s ../../RevTex/apsrmp.bst  ./apsrmp.bst
	ln -s ../../RevTex/revsymb.sty ./revsymb.sty
	ln -s ../../RevTex/aps.rtx ./aps.rtx
	ln -s ../../RevTex/10pt.rtx ./10pt.rtx
	ln -s ../../RevTex/11pt.rtx ./11pt.rtx
	ln -s ../../RevTex/12pt.rtx ./12pt.rtx

ps:	latex
	dvips -t letter -o prl.ps PBCQUICCA

dvi:	latex
	xdvi -s 2 PBCQUICCA

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
