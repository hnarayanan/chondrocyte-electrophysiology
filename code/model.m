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

  xdot = zeros(5, 1);

  V    = x(1);
  Na_i = x(2);
  Na_c = x(3);
  a_ur = x(4);
  I_ur = x(5);

  # Membrane capacitance
  C_m = 15.0;

  # Cellular volume
  vol_i = 0.005884;
  vol_c = 0.000800224;

  # External stimulation
  t_stim = 0.01;
  t_cycle = 1.0;
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
  g_K_ur     = 2.25;
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 6.0)/8.6));
  I_ur_inf = 1.0/(1.0 + exp((V + 7.5)/10.0));
  tau_a_ur   = 0.009/(1.0 + exp((V + 5.0)/12.0)) + 0.0005;
  tau_I_ur   = 0.59/(1.0 + exp((V + 60.0)/10.0)) + 3.05;
  V_K        = -83.042637;
  I_K_ur     = g_K_ur*a_ur*I_ur*(V - V_K);

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
  Na_i_dot =                      - (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_i*F);
  Na_c_dot = (Na_i - Na_c)/tau_Na + (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_c*F);
  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  I_ur_dot = (I_ur_inf - I_ur)/tau_I_ur;

  alpha_n = 0.01*(V + 10)/(exp((V + 10)/10) - 1);
  beta_n = 0.125*exp(V/80);

  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = Na_i_dot;
  xdot(3) = Na_c_dot;
  xdot(4) = a_ur_dot;
  xdot(5) = I_ur_dot;


endfunction

% Time stepping information
t_final = 100.0
dt = 0.005

% Initial conditions
V0 = -62.3
Na_i_0 = 0.516766
Na_c_0 = 130.022096
a_ur_0 = 0.000367
I_ur_0 = 0.967290

t = linspace(0, t_final, t_final/dt);
x0 = [V0, Na_i_0, Na_c_0, a_ur_0, I_ur_0];
x = lsode("f", x0, t);

clf
figure(1)
subplot(3, 2, 1), plot(t, x(:, 1)), legend('Membrane Voltage (mV)')
subplot(3, 2, 2), plot(t, x(:, 2)), legend('[Na^{+}]_{i} (mM/l)')
subplot(3, 2, 3), plot(t, x(:, 3)), legend('[Na^{+}]_{c} (mM/l)')
subplot(3, 2, 4), plot(t, x(:, 4)), legend('a_{{ur}_{0}}')
subplot(3, 2, 5), plot(t, x(:, 5)), legend('I_{{ur}_{0}} (pA)')
xlabel('Time (s)')
hold on

print -depsc2 "output.eps"
#print -depslatex "membrane_voltage.tex"
