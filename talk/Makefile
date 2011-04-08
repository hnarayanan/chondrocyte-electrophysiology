file = chondrocyte-201104

all: ${file}.pdf clean

clean:
	rm -f *~
	rm -f ${file}.aux
	rm -f ${file}.log
	rm -f ${file}.nav
	rm -f ${file}.out
	rm -f ${file}.snm
	rm -f ${file}.toc
	rm -f ${file}.vrb

pristine:
	rm -f ${file}.pdf

${file}.pdf: ${file}.tex pristine
	pdflatex ${file}