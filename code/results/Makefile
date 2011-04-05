TEX_PLOTS = $(wildcard *.tex)
PDF_PLOTS = ${TEX_PLOTS:.tex=.pdf}

plot: prep $(PDF_PLOTS) post

.SUFFIXES:

.SUFFIXES: .tex .pdf

.tex.pdf:
	pdflatex $<

prep:
	mv *.eps eps
	perl -pi -e 's/results\//eps\//' *.tex

post:
	rm eps/*-converted-to.pdf
	rm *.log *.aux
	mv *.tex tex
	mv *.pdf pdf

clean:
	rm -fr *~ eps/*.eps pdf/*.pdf tex/*.tex
