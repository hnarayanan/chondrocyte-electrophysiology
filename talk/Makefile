# Modified from the simple build system I used for managing my
# doctoral dissertation.
#
# Copyright (C) 2007--2011  Harish Narayanan

BASENAME = chondrocyte-model-201104

SOURCES = Makefile ${BASENAME}.tex \
	  ../images/pdf/*.pdf ../results/pdf/*.pdf

all: pristine ${BASENAME}.pdf
	make clean

clean:
	@echo "Cleaning cruft:"
	rm -f *~
	rm -f ${BASENAME}.aux
	rm -f ${BASENAME}.log
	rm -f ${BASENAME}.nav
	rm -f ${BASENAME}.out
	rm -f ${BASENAME}.snm
	rm -f ${BASENAME}.toc
	rm -f ${BASENAME}.vrb

pristine: clean
	@echo "Removing older output files:"
	rm -f ${BASENAME}.pdf

${BASENAME}.pdf: ${SOURCES}
	@echo "Creating the PDF file:"
	pdflatex ${BASENAME}
	egrep --color=auto -i '(warn|error|full)' ${BASENAME}.log