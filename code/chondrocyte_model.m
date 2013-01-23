# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
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

# Define the ODE model for the electrophysiology of the chondrocyte
function xdot = ode_rhs(x, t)

  global theta0;
  xdot = ode_rhs_parametrized(x, t, theta0);

endfunction

function xdot = ode_rhs_parametrized(x, t, theta)

  # Initialize and populate vector of unknowns
  global apply_Vm;
  if (apply_Vm == true)
    V = appliedVoltage(t);
  else
    V = x(1);
  endif
  Na_i = x(2);
  K_i  = x(3);
  Ca_i = x(4);
  H_i  = x(5);
  Cl_i = x(6);
  a_ur = x(7);
  i_ur = x(8);
  global vol_i_0;
  vol_i = vol_i_0; #x(9); #FIXME: Fixing the volume while debugging
  cal   = x(10);

  # Define external concentrations
  K_o = appliedPotassiumConcentration(t);

  # Extract parameters
  g_K_b_bar = theta(1);
  P_K = theta(2);
  Gmax = theta(3);

  # Calculate background currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar);
  I_Cl_b = backgroundChloride(V, Cl_i);

  # Calculate pump and exchanger currents
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i, K_o);
  I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i);
  I_NaH = sodiumHydrogenExchanger(Na_i, H_i);
  I_Ca_ATP = calciumPump(Ca_i);

  # Calculate other potassium currents
  I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur);
  I_K_2pore = twoPorePotassium(V, K_i, K_o, P_K);
  I_K_Ca_act = calciumActivatedPotassium(V, K_i, K_o, Ca_i, Gmax);
  I_K_ATP = potassiumPump(V, K_i, K_o);

  # Calculate other currents
  I_ASIC = voltageActivatedHydrogen();
  I_TRP1 = stretchActivatedTrip(V);
  I_TRP2 = osteoArthriticTrip();
  I_stim = externalStimulation(t);

  # Total ionic contribution
  I_i = I_Na_b + I_K_b + I_Cl_b \
      + I_NaK + I_NaCa + I_Ca_ATP \
      + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP \
      + I_ASIC + I_TRP1 + I_TRP2;

  # Determine incremental changes in evolving quantities
  global F;
  global C_m;

  # Evolve calmodulin
  cal_dot = 200000*Ca_i*(1.0 - cal) - 476.0*cal;

  # Evolve the concentrations
  Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH)/(vol_i*F);
  K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP)/(vol_i*F);
  Ca_i_dot =   (I_NaCa  - I_Ca_ATP)/(vol_i*F) - 0.045*cal_dot;
  H_i_dot =  - (I_NaH)/(vol_i*F);
  Cl_i_dot =   (I_Cl_b)/(vol_i*F);

  [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);

  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur;

  global Na_o;
  global Ca_o;
  global H_o;
  global Cl_o;

  global Na_i_0;
  global K_i_0;
  global Ca_i_0;
  global H_i_0;
  global Cl_i_0;

  osm_i_0 = Na_i_0 + K_i_0 + Ca_i_0 + H_i_0 + Cl_i_0;

  osm_i = Na_i + K_i + Ca_i + H_i + Cl_i;
  osm_o = Na_o + K_o + Ca_o + H_o + Cl_o;

  dosm = osm_i_0 - osm_o;

  P_f = 10.2e-4;
  SA = 6.0^(2.0/3.0)*pi^(1.0/3.0)*vol_i^(2.0/3.0);
  V_W = 18.0;
  vol_i_dot = P_f*SA*V_W*(osm_i - osm_o - dosm);

  xdot = zeros(10, 1);

  global apply_Vm;
  if (apply_Vm == true)
    xdot(1) = 0.0;
  else
    xdot(1) = 1/C_m*(-I_i + I_stim);
  endif

  global clamp_conc;
  if (clamp_conc == true)
    xdot(2) = 0.0;
    xdot(3) = 0.0;
    xdot(4) = 0.0;
    xdot(5) = 0.0;
    xdot(6) = 0.0;
  else
    xdot(2) = Na_i_dot;
    xdot(3) = K_i_dot;
    xdot(4) = Ca_i_dot;
    xdot(5) = H_i_dot;
    xdot(6) = Cl_i_dot;
  endif

  xdot(7) = a_ur_dot;
  xdot(8) = i_ur_dot;

  xdot(9) = vol_i_dot;

  xdot(10) = cal_dot;

endfunction
