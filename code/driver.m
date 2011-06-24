# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Initialise the environment
initialise;

# Load model parameters
parameters;

# Load utility functions
utility_functions;

# Define different currents
background_currents;
pumps_and_exchangers;
potassium_currents;
other_currents;

# Define the ODE model for the electrophysiology of the chondrocyte
chondrocyte_model;

# Estimate parameters and solve the ODE system
solve_model;

# Extract and postprocess solutions
postprocess_solutions;

# Plot the computed solutions
plot_solutions;