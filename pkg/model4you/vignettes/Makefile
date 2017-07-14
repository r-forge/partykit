model4you.pdf: model4you.tex 
	pdflatex model4you.tex
	bibtex model4you
	texi2pdf --clean model4you.tex

model4you.tex: model4you.Rnw
	Rscript -e "knitr::knit('model4you.Rnw')"
#	R CMD Sweave manuscript_model4you.Rnw

all: model4you.pdf
