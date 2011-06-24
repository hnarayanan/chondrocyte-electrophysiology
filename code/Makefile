# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

SOURCES = *.m ../data/reference_values/*.data

all: plots
	make clean

calc: ${SOURCES}
	@echo "Running calculation:"
	octave --quiet driver.m

plots: calc
	@echo "Creating pretty plots:"
	(cd ../results; make plots)

clean:
	@echo "Cleaning cruft:"
	rm -f *~
	(cd ../results; make clean)

pristine: clean
	@echo "Removing older output files:"
	(cd ../results; make pristine)

