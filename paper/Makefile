# Modified from the simple build system I used for managing my
# doctoral dissertation.
#
# Copyright (C) 2007--2011 Harish Narayanan

BASENAME = chondrocyte-model

SOURCES = Makefile ${BASENAME}.tex \
	  ${BASENAME}.bib ${BASENAME}.bst \
	  auxiliary/*.tex sections/*.tex \
          ../images/pdf/*.pdf ../results/pdf/*.pdf

all: pristine ${BASENAME}.pdf
	make clean

clean:
	@echo "Cleaning cruft:"
	rm -f *~
	rm -f ${BASENAME}.aux
	rm -f ${BASENAME}.bbl
	rm -f ${BASENAME}.blg
	rm -f ${BASENAME}.brf
	rm -f ${BASENAME}.log
	rm -f ${BASENAME}.out
	rm -f ${BASENAME}.spl
	(cd auxiliary; make clean)
	(cd sections; make clean)

pristine: clean
	@echo "Removing older output files:"
	rm -f *.pdf

${BASENAME}.pdf: ${SOURCES}
	@echo "Creating the PDF file:"
	pdflatex -interaction=batchmode ${BASENAME} && \
	bibtex -terse ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME}
	egrep --color=auto -i '(warn|error|full)' ${BASENAME}.log