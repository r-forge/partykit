letter: 
	pdflatex cover_letter.tex
	make clean

paper: 
	pdflatex draft_jrssc.tex
	bibtex draft_jrssc.aux
	pdflatex draft_jrssc.tex
	pdflatex draft_jrssc.tex
	make clean

rebuttal: 
	pdflatex rebuttal.tex
	#bibtex rebuttal
	#pdflatex rebuttal.tex
	pdflatex rebuttal.tex
	make clean


clean:
	rm -f *.out *.aux *.blg *.bbl *.log *.tpt 
