# Modified from the simple build system I used for managing my
# doctoral dissertation.
# (c) Harish Narayanan, 2007 -- 2011

# Fundamental variables:

# A versioning scheme to keep track of progress
VERSION = 0.0.1
BASENAME = chondrocyte-model

# Change these to your heart's content.
LATEX  = pdflatex
BIBTEX = bibtex

SOURCES = Makefile ${BASENAME}.tex ${BASENAME}.bib ${BASENAME}.bst \
	  auxiliary/*.tex sections/*.tex \
          ../images/pdf/*.pdf ../results/pdf/*.pdf

# Subtle changes to the command-line flags below can have significant
# impact on the quality of generated documents. Comment-out the ps2pdf
# flags below if on Mac OS X.
LFLAGS  = 					# pdflatex flags
BFLAGS  = -terse				# bibtex flags

# Human-readable targets:

# I find it easiest to use the following names while running make. For
# e.g., running 'make ${BASENAME}' (or just 'make') cleans the source
# distribution of cruft and generates a pristine PDF.

all:	  ${BASENAME}

polish:   ${BASENAME}.pdf

${BASENAME}:   pristine polish
	  make clean

archive:  pristine
	  hg archive -ttbz2 ../older/${BASENAME}-${VERSION}.tar.bz2

# Actual drudgery:

# Shoo, here be dragons.

clean:
	@echo "Cleaning cruft:"
	rm -f *~ *.aux *.bbl *.blg *.brf *.dvi *.loa *.lof *.log *.lot *.out \
	*.ps *.toc *.tmp *.spl
	(cd auxiliary; make clean)
	(cd sections; make clean)
	(cd ../images; make clean)
	(cd ../results; make clean)

pristine: clean
	@echo "Removing older output files:"
	rm -f *.pdf

${BASENAME}.pdf: ${SOURCES}
	@echo "Creating the PDF file:"
	${LATEX}   ${LFLAGS} ${BASENAME} && ${BIBTEX}  ${BFLAGS} ${BASENAME} && \
	${LATEX}   ${LFLAGS} ${BASENAME} && ${LATEX}   ${LFLAGS} ${BASENAME}
	egrep --color=auto -i '(warn|error|full)' ${BASENAME}.log