# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define initial conditions
x0 = [V_0, Na_i_0, K_i_0, Ca_i_0, H_i_0, a_ur_0, i_ur_0];

# Define the initial values of parameters to be estimated
global theta0 = [g_K_b_bar, P_K, Gmax];

# Define the problem time range
t = linspace(0, t_final, t_final/dt);

# Solve the problem
global enable_parest;
if (enable_parest == false) # Without parameter estimation
  x = lsode("ode_rhs", x0, t);
else
  # Define the overall model
  global model;
  model.odefcn = @ode_rhs_parametrized;
  model.tplot = t';
  model.param = theta0;
  model.ic = x0';

  # Load experimental measurements from files
  # table = load('../data/reference_values/generated_I_K_b_small.data');
  table = load('../data/reference_values/Total_Current_Density_vs_Time.data');
  measure.time = table(:, 1);
  measure.data = table(:, 2)*C_m;
  measure.statefcn = @measurefcn;
  measure.dstatedx = @measurederiv;

  # Define the search space for the parameters
  objective.estflag = [1, 2, 3];
  objective.paric   = theta0(objective.estflag);
  objective.parlb   = [0.1, 2.5e-6, 1](objective.estflag);
  objective.parub   = [0.3, 3.5e-6, 3](objective.estflag);

  # Estimate the parameters and corresponding solutions
  estimates = parest(model, measure, objective);

  # Extract solution and parameters
  x = estimates.x';

  # Extract parameters
  disp('Estimated parameters and bounding box')
  [estimates.parest estimates.bbox]

  g_K_b_bar = estimates.parest(1);
  P_K = estimates.parest(2);
  Gmax = estimates.parest(3);
endif