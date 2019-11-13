


paper: 
	R CMD Sweave --pdf draft.Rnw
	make clean

clean:
	rm *.out *.aux *.blg *.bbl *.log *.tpt draft.tex

diff:
	R CMD Sweave --pdf draft.Rnw
	latexdiff -t CTRADITIONAL --exclude-textcmd="section" --exclude-textcmd="subsection" --exclude-textcmd="subsubsection" draft_2019_11_12.tex draft.tex  > diff.tex
	-pdflatex diff.tex
	-bibtex diff.aux
	-pdflatex diff.tex
	-pdflatex diff.tex
	-rm diff.out diff.aux diff.blg diff.bbl diff.log diff.tex draft.tex
	-mv diff.pdf diff_tease_paper_version_2019_11_12.pdf


