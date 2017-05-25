man_model4you.pdf: man_model4you.tex
	texi2pdf --clean man_model4you.tex

man_model4you.tex: man_model4you.Rnw
	Rscript -e "knitr::knit('man_model4you.Rnw')"
#	R CMD Sweave manuscript_model4you.Rnw

all: man_model4you.pdf
