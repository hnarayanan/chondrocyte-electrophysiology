# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

TEX_PLOTS = $(wildcard *.tex)
PDF_PLOTS = ${TEX_PLOTS:.tex=.pdf}

plots: prep $(PDF_PLOTS) post

.SUFFIXES:

.SUFFIXES: .tex .pdf

.tex.pdf:
	pdflatex -interaction=batchmode $<

prep:
	perl -pi -e 's/..\/results\/epslatex\///g' *.tex

post:
	mv *.pdf ../pdf

clean:
	rm -f *~
	rm -f *.log *.aux

pristine:
	rm -f *.eps
	rm -f *.tex