# ---------------------------------------------------------------------------
# List of Files to be processed

RNW_SOURCES=slides_tree slides_forest

# ---------------------------------------------------------------------------

.PHONY: clean distclean

all: slides

slides: $(RNW_SOURCES:=.pdf) clean

scripts: $(RNW_SOURCES:=.R)

handouts: $(RNW_SOURCES:=_1x1.pdf) clean

# ---------------------------------------------------------------------------
# Generic rules for creating files

%.tex: %.Rnw
	echo "Sweave('$<')" | R --no-save --no-restore --no-site-file

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

%.R: %.Rnw
	echo "Stangle('$<')" | R --no-save --no-restore --no-site-file --no-init-file

%_1x1.tex: %.tex
	echo "x <- readLines('$<'); x[1] <- gsub(',t,', ',t,handout,', x[1], fixed = TRUE); writeLines(x, '$@')" | R --no-save --no-restore --no-site-file --no-init-file

%_1x1.pdf: %_1x1.tex
	pdflatex $<
	pdflatex $<

%_2x2.pdf: %_1x1.pdf
	pdfjam --nup 2x2 --landscape --scale 0.95 --outfile $@ $<

# ---------------------------------------------------------------------------
# Clear working space

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb
	rm -f *~
	rm -f Rplots.ps
	rm -f Rplots.pdf
	rm -f *.eps
	rm -f $(RNW_SOURCES:=.tex)
	rm -f $(RNW_SOURCES:=-*.pdf)
	rm -f $(RNW_SOURCES:=-*.png)

distclean: clean
	rm -f $(RNW_SOURCES:=.pdf)
	rm -f $(RNW_SOURCES:=.R)
	rm -f $(RNW_SOURCES:=_1x1.pdf)
	rm -f $(RNW_SOURCES:=_2x2.pdf)
	rm -f *.rda
	rm -f *.rds

# ---------------------------------------------------------------------------
