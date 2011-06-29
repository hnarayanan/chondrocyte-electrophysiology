# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define initial conditions
x0 = [V_0, Na_i_0, K_i_0, Ca_i_0, H_i_0, a_ur_0, i_ur_0];

# Define the initial values of parameters to be estimated
theta0 = [g_K_b_bar, P_K];

# Define the problem time range
t = linspace(0, t_final, t_final/dt);

# Define the overall model
global model;
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
objective.estflag = [1, 2];
objective.paric   = theta0;
objective.parlb   = [0, 0];
objective.parub   = [1, 1.e-5];

# Estimate the parameters and corresponding solutions
estimates = parest(model, measure, objective);