TEX_PLOTS = $(wildcard *.tex)
PDF_PLOTS = ${TEX_PLOTS:.tex=.pdf}

plot: prep $(PDF_PLOTS) post

.SUFFIXES:

.SUFFIXES: .tex .pdf

.tex.pdf:
	pdflatex $<

prep:
	mkdir -p results
	mv *.eps results

post:
	rm *.log *.aux
	mv results/*.eps eps
	rm -fr results
	mv *.pdf pdf
	mv *.tex tex

clean:
	rm -fr *~ eps/*.eps pdf/*.pdf tex/*.tex
