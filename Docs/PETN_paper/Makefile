all: basic

basic:
	- /bin/rm hicu.bbl
	latex petn.tex
	bibtex petn
	latex petn.tex
	latex petn.tex

ps: basic
	dvips -t letter -o petn.ps petn

clean:
	- /bin/rm *.blg *.aux *.bbl *.dvi *.log
