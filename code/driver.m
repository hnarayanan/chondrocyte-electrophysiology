# The solution fields and parameters have the following units:
#
# Voltage: mV
# Current: pa
# Time: s
# Concentration: mM/l
# Conductance: pS
# Capacitance: pF

# Clear memory
clear all;

# Load model parameters
parameters;

# Load utility functions
utility_functions;

# Define different currents
background_currents;
pumps_and_exchangers;
potassium_currents;
other_currents;

# Define the ODE system
function xdot = f(x, t)

  # Initialize and populate vector of unknowns
  V    = x(1);
  Na_i = x(2);
  K_i  = x(3);
  Ca_i = x(4);
  a_ur = x(5);
  i_ur = x(6);

  # Calculate background currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i);

  # Calculate pump and exchanger currents
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i);
  I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i);
  I_NaH = sodiumHydrogenAntiport();

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
  Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_i*F);
  K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP)/(vol_i*F);
  Ca_i_dot =   (I_NaCa)/(vol_i*F);

  [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);

  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur;

  xdot = zeros(6, 1);
  global clamp_Vm;
  if (clamp_Vm == true)
    global VF, global V0, global t_final;
    xdot(1) = (VF - V0)/(t_final - 0.0);
  else
    xdot(1) = 1/C_m*(-I_i + I_stim);
  endif
  xdot(2) = Na_i_dot;
  xdot(3) = K_i_dot;
  xdot(4) = Ca_i_dot;
  xdot(5) = a_ur_dot;
  xdot(6) = i_ur_dot;

endfunction

# Solve the ODE system for all time t
t = linspace(0, t_final, t_final/dt);
x0 = [V0, Na_i_0, K_i_0, Ca_i_0, a_ur_0, i_ur_0];
x = lsode("f", x0, t);

# Extract and postprocess solutions
postprocess_solutions;

# Plot the computed solutions
plot_solutions;
