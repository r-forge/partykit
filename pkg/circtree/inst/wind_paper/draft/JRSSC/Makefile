letter: 
	pdflatex cover_letter.tex
	make clean

rebuttal: 
	pdflatex rebuttal.tex
	#bibtex rebuttal
	#pdflatex rebuttal.tex
	pdflatex rebuttal.tex
	make clean

paper: 
	pdflatex draft_jrssc.tex
	bibtex draft_jrssc.aux
	pdflatex draft_jrssc.tex
	pdflatex draft_jrssc.tex
	make clean

diff:
	latexdiff -t CTRADITIONAL --exclude-textcmd="section" --exclude-textcmd="subsection" --exclude-textcmd="subsubsection" draft_jrssc_submitted_2020_01_09.tex draft_jrssc.tex > diff.tex
	-pdflatex diff.tex
	-bibtex diff.aux
	-pdflatex diff.tex
	-pdflatex diff.tex
	-rm diff.out diff.aux diff.blg diff.bbl diff.log diff.tex
	-make clean
	-mv diff.pdf diff_draft_jrssc_to_submitted_2020_01_09.pdf

clean:
	rm -f *.out *.aux *.blg *.bbl *.log *.tpt 
