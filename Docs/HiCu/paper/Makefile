all: normal

normal:
	latex hicu.tex

new: 
	- /bin/rm hicu.bbl
	latex hicu.tex
	bibtex hicu
	latex hicu.tex
	latex hicu.tex

ps:
	dvips -o hicu.ps hicu

clean:
	- rm *.log *.ps *.aux *.dvi
