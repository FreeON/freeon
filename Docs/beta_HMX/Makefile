all: ps

ps: 
	- /bin/rm bhmx.bbl
	latex bhmx.tex
	bibtex bhmx
	latex bhmx.tex
	latex bhmx.tex	
	dvips -t letter -o bhmx.ps bhmx

pdf:
	- /bin/rm bhmx.bbl
	latex bhmx.tex
	bibtex bhmx
	latex bhmx.tex	
	latex bhmx.tex	
	dvips -t letter -P pdf -o bhmx.ps bhmx
	ps2pdf bhmx.ps
	

	
clean:
	- rm *.log *.aux *.dvi *.bbl *.blg
