all: ps

ps: 
	- /bin/rm paraqctc.bbl
	latex paraqctc.tex
	bibtex paraqctc
	latex paraqctc.tex
	latex paraqctc.tex	
	dvips -t letter -o paraqctc.ps paraqctc

pdf:
	- /bin/rm paraqctc.bbl
	latex paraqctc.tex
	bibtex paraqctc
	latex paraqctc.tex	
	latex paraqctc.tex	
	dvips -t letter -P pdf -o paraqctc.ps paraqctc
	ps2pdf paraqctc.ps
	

	
clean:
	- rm *.log *.aux *.dvi
