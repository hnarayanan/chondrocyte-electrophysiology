# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

all: plots paper talk
	make clean

plots:
	@echo "Running calculation and generating result plots"
	(cd code; make)

paper: plots
	@echo "Creating the paper"
	(cd paper; make)

talk: plots
	@echo "Creating the talk"
	(cd talk; make)

clean:
	@echo "Cleaning cruft:"
	rm -f *~
	(cd code; make clean)
	(cd data; make clean)
	(cd images; make clean)
	(cd paper; make clean)
	(cd results; make clean)
	(cd talk; make clean)


pristine:
	@echo "Removing older output files:"
	(cd code; make clean)
	(cd data; make clean)
	(cd images; make clean)
	(cd paper; make pristine)
	(cd results; make pristine)
	(cd talk; make pristine)

