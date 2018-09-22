# ---------------------------------------------------------------------------
# List of stations to be processed

STATIONS=Axams Lech Zuers SeeimPaznaun Jungholz LadisNeuegg Oetz OchsengartenObergut Ginzling Rotholz Walchsee Koessen Innervillgraten MatreiinOsttirol StJohannimWalde  

# ---------------------------------------------------------------------------

.PHONY: all clean distclean

all: $(STATIONS:=.pdf) clean

# ---------------------------------------------------------------------------
# Generic rules for creating files

%.Rnw: Rain_stationwise_results.Rnw
	echo "x <- readLines('Rain_stationwise_results.Rnw'); x <- gsub('_STATION_', tools::file_path_sans_ext('$@'), x, fixed = TRUE); writeLines(x, '$@')" | R --no-save --no-restore --no-site-file --no-init-file

%.tex: %.Rnw
	echo "Sweave('$<')" | R --no-save --no-restore --no-site-file --no-init-file

%.pdf: %.tex
	texi2pdf $<

# ---------------------------------------------------------------------------
# Clear working space

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb
	rm -f *~
	rm -f *-*.pdf
	rm -f *-*.png
	rm -f Rplots.ps Rplots.pdf
	rm -f $(STATIONS:=.Rnw)
	rm -f $(STATIONS:=.tex)

distclean: clean
	rm -f *.pdf

# ---------------------------------------------------------------------------
