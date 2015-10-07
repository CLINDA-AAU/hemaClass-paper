##################################################################
# Makefile for LaTeX
##################################################################
# Use:
# make
# make all
# make clean
# options for ps2pdf: http://pages.cs.wisc.edu/~ghost/doc/AFPL/6.50/Ps2pdf.htm

OTHER = *~ *.aux *.dvi *.toc *.bbl *.blg *.out *.thm *.ps *.idx *.ilg *.ind \
 *.tdo *.synctex.gz

paper: pdf clean

all: results pdf clean

results: figures/figure2.pdf figures/figure3.pdf tables/table1.tex
	Rscript --vanilla hemaclass.R

pdf:
	pdflatex --synctex=1 hemaclass.tex
	bibtex hemaclass
	pdflatex --synctex=1 hemaclass.tex
	pdflatex --synctex=1 hemaclass.tex

clean:
	rm -f $(OTHER) $(PS)

