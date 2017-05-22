manuscript_model4you.pdf: manuscript_model4you.tex
	texi2pdf --clean manuscript_model4you.tex

manuscript_model4you.tex: manuscript_model4you.Rnw
	Rscript -e "knitr::knit('manuscript_model4you.Rnw')"
#	R CMD Sweave manuscript_model4you.Rnw

all: manuscript_model4you.pdf
