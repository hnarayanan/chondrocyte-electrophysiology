# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the overall model
global model;
x0 = [V_0, K_i_0];
t = linspace(0, t_final, t_final/dt);
theta0 = [0.1];
model.odefcn = @ode_rhs;
model.tplot = t';
model.param = theta0;
model.ic = x0';

# Load experimental measurements from files
table = load('../data/reference_values/generated_I_K_b_small.data');
measure.time = table(:, 1);
measure.data = table(:, 2);
measure.statefcn = @measurefcn;
measure.dstatedx = @measurederiv;

# Define the search space for the parameters
objective.estflag = [1];
objective.paric   = theta0;
objective.parlb   = [0];
objective.parub   = [2];

# Estimate the parameters and corresponding solutions
estimates = parest(model, measure, objective);