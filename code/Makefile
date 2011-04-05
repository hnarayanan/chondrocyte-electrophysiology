all:	plots

prep:
	@echo "Preparing calculation"
	mkdir -p results
	mkdir -p results/tex
	mkdir -p results/eps
	mkdir -p results/pdf

calc: prep *.m
	@echo "Running calculation"
	octave driver.m

plots: calc
	@echo "Creating pretty plotss"
	(cd results; make plots)

clean:
	@echo "Cleaning cruft"
	rm -f *~
	(cd results; make clean)