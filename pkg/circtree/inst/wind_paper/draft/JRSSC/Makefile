letter: 
	pdflatex cover_letter.tex
	make clean

paper: 
	pdflatex draft_jrssc.tex
	bibtex draft_jrssc.aux
	pdflatex draft_jrssc.tex
	pdflatex draft_jrssc.tex
	make clean

clean:
	rm *.out *.aux *.blg *.bbl *.log *.tpt 
