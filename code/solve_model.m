# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# The solution fields and parameters have the following units:
#
# Voltage: mV
# Current: pa
# Time: s
# Concentration: mM/l
# Conductance: pS
# Capacitance: pF

1;

x0 = [V_0, K_i_0];
t = linspace(0, t_final, t_final/dt);
theta0 = [0.1];

# Define the model
global model;
model.odefcn = @ode_rhs;
model.tplot = t';
model.param = theta0;
model.ic = x0';

# Load measurements from a file
measure.states = [1, 2];
table = load ('../data/reference_values/generated_small.data');
measure.time = table(:, 1);
measure.data = table(:, 1 + measure.states);

# Define the search space for the parameters
objective.estflag = [1];
objective.paric   = theta0;
objective.parlb   = [0];
objective.parub   = [2];

# Estimate the parameters
estimates = parest(model, measure, objective);
disp('Estimated Parameters and Bounding Box')
[estimates.parest estimates.bbox]

# Plot the model fit to the noisy measurements
figure(1);
plot(model.tplot, estimates.x, measure.time, measure.data, 'o');
print -depslatexstandalone "../results/epslatex/potassium_quantities.tex"

table = [model.tplot, estimates.x'];
save ABC.dat table