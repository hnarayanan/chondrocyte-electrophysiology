# The solution fields and parameters have the following units:
#
# Voltage: mV
# Current: pa
# Time: s
# Concentration: mM/l
# Conductance: nS
# Capacitance: nF

clear all;

function xdot = f(x, t)

  xdot = zeros(4, 1);

  V    = x(1);
  Na_i = x(2);
  Na_c = x(3);
  n     = x(4);

  # Membrane capacitance
  C_m = 15.0;

  # Cellular volume
  vol_i = 0.005884;
  vol_c = 0.000800224;

  # External stimulation
  t_stim = 0.1;
  t_cycle = 0.5;
  I_stim_bar = 0.0;

  # Background sodium
  g_Na_b_bar = 20;
  V_Na_b = +55;
  I_Na_b = g_Na_b_bar*(V - V_Na_b);

  # Background potassium
  g_K_b_bar = 136;
  V_K_b = -80;
  I_K_b = g_K_b_bar*(V - V_K_b);

  # Sodium-potassium pump
  I_NaK_bar = 68.55;
  K_c = 5.560224;
  K_NaK_K = 1.0;
  K_NaK_Na = 11.0;
  I_NaK = I_NaK_bar*(K_c/(K_c + K_NaK_K))*(Na_i^1.5/(Na_i^1.5 + K_NaK_Na^1.5))*(V + 150.0)/(V + 200.0);

  # Sodium-Calcium exchanger
  K_NaCa = 0.0374842;
  Ca_c = 1.815768;
  gamma_Na = 0.45;
  F = 96487.0;
  R = 8314.0;
  T = 306.15;
  Ca_i = 6.5e-5;
  d_NaCa = 0.0003;
  I_NaCa = K_NaCa*(Na_i^3*Ca_c*exp(gamma_Na*V*F/(R*T)) - Na_c^3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
           / (1.0 + d_NaCa*(Na_c^3*Ca_i + Na_i^3*Ca_c));

  # Sodium-hydrogen exchanger
  I_NaH      = 0.0;

  # Ultra-rapidly rectifying potassium
  I_K_ur     = 0.0;

  # Two-pore potassium channel
  I_K_2pore  = 0.0;

  # Calcium-activated potassium channel
  I_Ca_act_K = 0.0;

  # Trip channel(s)
  I_TRP      = 0.0;

  # Total ionic contribution
  I_i = I_Na_b + I_K_b + I_NaK + I_NaCa + I_NaH \
      + I_K_ur + I_K_2pore + I_Ca_act_K + I_TRP;

  # External stimulation
  I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);

  # Changes in concentration
  tau_Na = 0.01;
  Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_i*F);
  Na_c_dot = (Na_i - Na_c)/tau_Na + (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_c*F);

  alpha_n = 0.01*(V + 10)/(exp((V + 10)/10) - 1);
  beta_n = 0.125*exp(V/80);

  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = Na_i_dot;
  xdot(3) = Na_c_dot;
  xdot(4) = alpha_n*(1 - n) - beta_n*n;

endfunction

% Time stepping information
t_final = 4.0
dt = 0.005

% Initial conditions
V0 = -60
Na_i_0 = 8.516766
Na_c_0 = 130.022096
n0 =  1

t = linspace(0, t_final, t_final/dt);
x0 = [V0, Na_i_0, Na_c_0, n0];
x = lsode("f", x0, t);

clf
figure(1)
plot(t, x(:, 1), "linewidth", 4)
xlabel("Time (s)")
ylabel("Membrane voltage (mV)")
print -depsc2 "membrane_voltage.eps"
print -depslatex "membrane_voltage.tex"