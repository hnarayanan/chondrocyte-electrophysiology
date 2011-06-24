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

# Define the ODE model for the electrophysiology of the chondrocyte
function xdot = ode_rhs(x, t, theta)

  # Initialize and populate vector of unknowns
  global apply_Vm;
  if (apply_Vm == true)
    V = appliedVoltage(t);
  else
    V = x(1);
  endif
#  Na_i = x(2);
#  K_i  = x(3);
#  Ca_i = x(4);
#  H_i  = x(5);
#  a_ur = x(6);
#  i_ur = x(7);

  dummy = 0.1;
  Na_i = dummy;
  Ca_i = dummy;
  H_i  = dummy;
  a_ur = dummy;
  i_ur = dummy;
  K_i = x(2);

  g_K_b_bar = theta(1);

  # Calculate background currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i, g_K_b_bar);

  # Calculate pump and exchanger currents
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i);
  I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i);
  I_NaH = sodiumHydrogenExchanger(Na_i, H_i);

  # Calculate other potassium currents
  I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, a_ur, i_ur);
  I_K_2pore = twoPorePotassium(V, K_i);
  I_K_Ca_act = calciumActivatedPotassium(V, K_i, Ca_i);
  I_K_ATP = potassiumPump();

  # Calculate other currents
  I_ASIC = voltageActivatedHydrogen();
  I_TRP1 = stretchActivatedTrip();
  I_TRP2 = osteoArthriticTrip();
  I_stim = externalStimulation(t);

  # Total ionic contribution
  I_i = I_Na_b + I_K_b \
      + I_NaK + I_NaCa + I_NaH \
      + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP \
      + I_ASIC + I_TRP1 + I_TRP2;

  # Determine incremental changes in evolving quantities
  global F, global vol_i;
  global C_m;

  # FIXME: Check signs on the following carefully
  Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH)/(vol_i*F);
  K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP)/(vol_i*F);
  Ca_i_dot =   (I_NaCa)/(vol_i*F);
  H_i_dot =  - (I_NaH)/(vol_i*F);

  [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);

  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur;

# #  xdot = zeros(7, 1);

#   global apply_Vm;
#   if (apply_Vm == true)
#     xdot(1) = 0.0;
#   else
#     xdot(1) = 1/C_m*(-I_i + I_stim);
#   endif

#   global clamp_conc;
#   if (clamp_conc == true)
#     xdot(2) = 0.0;
#     xdot(3) = 0.0;
#     xdot(4) = 0.0;
#     xdot(5) = 0.0;
#   else
#     xdot(2) = Na_i_dot;
#     xdot(3) = K_i_dot;
#     xdot(4) = Ca_i_dot;
#     xdot(5) = H_i_dot;
#   endif

#   xdot(6) = a_ur_dot;
#   xdot(7) = i_ur_dot;

  xdot = zeros(2, 1);

  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = K_i_dot;

endfunction