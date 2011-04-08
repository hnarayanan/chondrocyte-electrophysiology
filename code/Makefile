# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

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