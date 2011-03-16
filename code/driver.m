# The solution fields and parameters have the following units:
#
# Voltage: mV
# Current: pa
# Time: s
# Concentration: mM/l
# Conductance: nS
# Capacitance: nF (or pF)

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

  # Load useful constants
  global F, global R, global T;
  global C_m, global vol_i;

  # Initialize and populate vector of unknowns
  V    = x(1);
  Na_i = x(2);
  K_i  = x(3);
  Ca_i = x(4);
  a_ur = x(5);
  I_ur = x(6);

  # Calculate background currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i);

  # Calculate pump and exchanger currents
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i);
  I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i);
  I_NaH = sodiumHydrogenAntiport();

  # Calculate other potassium currents
  I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, a_ur, I_ur);
  I_K_2pore = twoPorePotassium(V, K_i);
  I_K_Ca_act = calciumActivatedPotassium(V, Ca_i);
  I_K_ATP = potassiumPump();

  # Calculate other currents
  I_ASIC = voltageActivatedHydrogen();
  I_stim = externalStimulation(t);

  # Trip channel(s)
  g_TRP      = 1.0; #FIXME: Should be a function of some external agent,
		    #e.g. a steroid or stretch
  I_TRP      = 0.0;#g_TRP*(V - V_Na_b);

  # Total ionic contribution
  I_i = I_Na_b + I_K_b + I_NaK + I_NaCa + I_NaH \
      + I_K_ur + I_K_2pore + I_K_ATP \
      + I_ASIC + I_K_Ca_act;# + I_TRP;

  # Changes in concentration (FIXME: Check these carefully)
  tau_Na = 0.01;
#  Na_i_dot =                      - (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_i*F);
#  Na_c_dot = 0.0;(Na_i - Na_o)/tau_Na + (I_Na_b + 3*I_NaK + 3*I_NaCa +  I_NaH)/(vol_c*F);

  Na_i_dot = -(I_Na_b)/(vol_i*F);
  K_i_dot  = -(I_K_b)/(vol_i*F);
  Ca_i_dot = -(0.0)/(vol_i*F);

  [a_ur_inf, I_ur_inf, tau_a_ur, tau_I_ur] = ultraRapidlyRectifyingPotassiumHelper(V);

  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  I_ur_dot = (I_ur_inf - I_ur)/tau_I_ur;

  alpha_n = 0.01*(V + 10)/(exp((V + 10)/10) - 1);
  beta_n = 0.125*exp(V/80);

  xdot = zeros(6, 1);
  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = Na_i_dot;
  xdot(3) = K_i_dot;
  xdot(4) = Ca_i_dot;
  xdot(5) = a_ur_dot;
  xdot(6) = I_ur_dot;

endfunction

# Solve the ODE system for all time t
t = linspace(0, t_final, t_final/dt);
x0 = [V0, Na_i_0, K_i_0, Ca_i_0, a_ur_0, I_ur_0];
x = lsode("f", x0, t);

# Extract and postprocess solutions
postprocess_solutions;

# Plot the computed solutions
plot_solutions;
