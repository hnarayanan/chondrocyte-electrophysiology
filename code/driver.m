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
# potassium_currents;
# etc.

function xdot = f(x, t)

  # Load useful constants
  global F, global R, global T;
  global C_m, global vol_i;
  global t_cycle, global t_stim, global I_stim_bar;

  # Initialize and populate vector of unknowns
  xdot = zeros(5, 1);
  V    = x(1);
  Na_i = x(2);
  K_i  = x(3);
  a_ur = x(4);
  I_ur = x(5);

  # Calculate different currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i);
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i);

  # Sodium-Calcium exchanger
  K_NaCa = 0.0374842;
  Ca_c = 1.815768;
  gamma_Na = 0.45;
  Ca_i = 6.5e-5;
  d_NaCa = 0.0003;
  I_NaCa = 0.0;
#  I_NaCa = K_NaCa*(Na_i^3*Ca_c*exp(gamma_Na*V*F/(R*T)) - Na_c^3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
#           / (1.0 + d_NaCa*(Na_c^3*Ca_i + Na_i^3*Ca_c));

  # Sodium-hydrogen exchanger
  I_NaH      = 0.0;

  # Ultra-rapidly rectifying potassium
  g_K_ur     = 2.25;
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 6.0)/8.6));
  I_ur_inf   = 1.0/(1.0 + exp((V + 7.5)/10.0));
  tau_a_ur   = 0.009/(1.0 + exp((V + 5.0)/12.0)) + 0.0005;
  tau_I_ur   = 0.59/(1.0 + exp((V + 60.0)/10.0)) + 3.05;
  V_K        = -83.042637;
  I_K_ur     = g_K_ur*a_ur*I_ur*(V - V_K);

  # Two-pore potassium channel
  g_2pore    = 1.0; #FIXME: This should be a sane value
  I_K_2pore  = g_2pore*(V - V_K);

  # Calcium-activated (maxi) potassium channel
  # Constants taken straight from H-A's Igor code

  Zj=0.58; Vhj=150; ZL=0.3; L0=1e-6; KDc=11e-6; C=8; D=25; E=2.4; Gmax=1;
  W = [Zj, Vhj, ZL, L0, KDc, C, D, E, Gmax];
  Ca = 11.e-6;
  T = 20;
  kTe = 23.54*((T + 273)/273);
  Lv = W(4)*exp((V*W(3))/kTe);
  Jv = exp(((V - W(2))*W(1))/kTe);
  K = Ca/W(5);
  C = W(6);
  D = W(7);
  E = W(8);
  G_max=W(9);
  N_channel = 1.0;
  V_Ca_act_K = 10.0;

  P0=(Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));

  I_Ca_act_K = N_channel*P0*G_max*(V - V_Ca_act_K);

  # Trip channel(s)
  g_TRP      = 1.0; #FIXME: Should be a function of some external agent,
		    #e.g. a steroid or stretch
  I_TRP      = 0.0;#g_TRP*(V - V_Na_b);

  # Total ionic contribution
  I_i = I_Na_b + I_K_b;# + I_NaK + I_NaCa + I_NaH \
#      + I_K_ur + I_K_2pore + I_Ca_act_K + I_TRP;

  # External stimulation
  I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);

  # Changes in concentration
  tau_Na = 0.01;
#  Na_i_dot =                      - (I_Na_b + 3*I_NaK + 3*I_NaCa + I_NaH)/(vol_i*F);
#  Na_c_dot = 0.0;(Na_i - Na_o)/tau_Na + (I_Na_b + 3*I_NaK + 3*I_NaCa +  I_NaH)/(vol_c*F);
  Na_i_dot = -(I_Na_b)/(vol_i*F);
  K_i_dot  = -(I_K_b)/(vol_i*F);
  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  I_ur_dot = (I_ur_inf - I_ur)/tau_I_ur;

  alpha_n = 0.01*(V + 10)/(exp((V + 10)/10) - 1);
  beta_n = 0.125*exp(V/80);

  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = Na_i_dot;
  xdot(3) = K_i_dot;
  xdot(4) = a_ur_dot;
  xdot(5) = I_ur_dot;

endfunction

% Time stepping information
t_final = 10.0
dt = 0.05

% Initial conditions
V0 = -62.3;
#Na_i_0 = 0.516766;
#Na_c_0 = 130.022096;
a_ur_0 = 0.000367;
I_ur_0 = 0.967290;

t = linspace(0, t_final, t_final/dt);
len_t = size(t, 2);
x0 = [V0, Na_i_0, K_i_0, a_ur_0, I_ur_0];
x = lsode("f", x0, t);

# Extract solution components
V    = x(:, 1);
Na_i = x(:, 2);
K_i  = x(:, 3);

# Compute currents at all times
I_Na_b = zeros(len_t, 1);
for ii = [1:len_t]
  I_Na_b(ii) = backgroundSodium(V(ii), Na_i(ii));
  I_K_b(ii)  = backgroundPotassium(V(ii), K_i(ii));
  I_NaK(ii) = sodiumPotassiumPump(V(ii), Na_i(ii), K_i(ii));
endfor

# Clear the screen and plot solutions
clf;
hold on;

# Just the membrane voltage
figure(1);
plot(t, V); xlabel('t (s)'); ylabel('V_{m} (mV)');
print -depsc2 "output/voltage.eps"

# The different current components
figure(2);
subplot(2, 2, 1), plot(t, I_Na_b), legend('I_{Na_{b}} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_K_b), legend('I_{K_{b}} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/background_currents.eps"

# The different current components
figure(3);
subplot(2, 2, 1), plot(t, I_NaK), legend('I_{NaK} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/pumps_and_exchangers.eps"

# The different concentrations
figure(4);
subplot(2, 2, 1), plot(t, Na_i), legend('[Na^{+}]_{i} (mM/l)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, K_i),  legend('[K^{+}]_{i} (mM/l)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/concentrations.eps"

#print -depslatex "membrane_voltage.tex"