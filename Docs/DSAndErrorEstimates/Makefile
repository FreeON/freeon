all: normal

normal:
	latex ExactMultipoleError.tex
	latex ExactMultipoleError.tex
	xdvi -s 2 ExactMultipoleError

new: 	clean links figs bib ps

bib:	
	latex ExactMultipoleError.tex
	bibtex ExactMultipoleError

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
	dvips -o ExactMultipoleError.ps ExactMultipoleError

clean:
	rm -rf *.log *.aux *.bbl *.blg *.dvi *~
	rm -rf mondo_new.bib revtex4.cls revsymb.sty aps.rtx 10pt.rtx 11pt.rtx 12pt.rtx apsrev.bst apsrmp.bst
