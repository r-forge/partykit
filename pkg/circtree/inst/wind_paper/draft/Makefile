


paper: 
	R CMD Sweave --pdf draft.Rnw
	make clean

clean:
	rm *.out *.aux *.blg *.bbl *.log *.tpt draft.tex

diff:
	R CMD Sweave --pdf draft.Rnw
	latexdiff -t CTRADITIONAL --exclude-textcmd="section" --exclude-textcmd="subsection" --exclude-textcmd="subsubsection" draft_submitted_2020_01_09.tex draft.tex  > diff.tex
	-pdflatex diff.tex
	-bibtex diff.aux
	-pdflatex diff.tex
	-pdflatex diff.tex
	-rm diff.out diff.aux diff.blg diff.bbl diff.log diff.tex
	-make clean
	-mv diff.pdf diff_paper_submitted_2020_01_09.pdf


diff_step1:
	R CMD Sweave --pdf draft.Rnw
	latexdiff -t CTRADITIONAL --exclude-textcmd="section" --exclude-textcmd="subsection" --exclude-textcmd="subsubsection" draft_submitted_2020_01_09.tex draft.tex  > diff.tex


diff_step2:
	-pdflatex diff.tex
	-bibtex diff.aux
	-pdflatex diff.tex
	-pdflatex diff.tex
	-rm diff.out diff.aux diff.blg diff.bbl diff.log diff.tex
	-make clean
	-mv diff.pdf diff_paper_submitted_2020_01_09.pdf
