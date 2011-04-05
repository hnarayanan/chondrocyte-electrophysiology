all:	plot

prep:
	@echo "Preparing calculation"
	mkdir -p results
	mkdir -p results/tex
	mkdir -p results/eps
	mkdir -p results/pdf

calc: prep *.m
	@echo "Running calculation"
	octave driver.m

plot: calc
	@echo "Creating pretty plots"
	(cd results; make plot)

clean:
	@echo "Cleaning cruft"
	rm -f *~
	(cd results; make clean)