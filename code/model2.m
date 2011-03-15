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

# Universal constants
global R = 8.314472; # Universal gas constant (J K^-1 mol^-1)
global T = 310.15;   # Normal body temperature (K)
global F = 96485.34; # Faraday's constant (C mol^-1)

# Cell parameters
global C_m = 15.0;       # Membrane capacitance
global vol_i = 0.005884; # Internal volume

# Constants related to sodium
global Na_i_0 = 0.516766;    # Initial internal sodium concentration (mM/l)
global Na_o   = 130.022096;  # Clamped external sodium concentration (mM/l)
global z_Na = 1;             # Charge on the sodium ion
global g_Na_b_bar = 20;      # Background sodium leakage conductance (nS)

# Constants related to potassium
global K_i_0 = 129.485991;  # Initial internal potassium concentration (mM/l)
global K_o   = 5.560224;    # Clamped external potassium concentration (mM/l)
global z_K = 1;             # Charge on the potassium ion
global g_K_b_bar = 20;      # Background potassium leakage conductance (nS)

# Potential of an ion X across the membrane (mV).
function E_X = nernstPotential(z, X_i, X_o)
  global R, global T, global F;
  E_X = 1000.0*(R*T)/(z*F)*log(X_o/X_i);
endfunction

# Background sodium current
function I_Na_b = backgroundSodium(V, Na_i)
  global z_Na, global g_Na_b_bar, global Na_o;
  E_Na = nernstPotential(z_Na, Na_i, Na_o)
  I_Na_b = g_Na_b_bar*(V - E_Na);
endfunction

# Background potassium current
function I_K_b = backgroundPotassium(V, K_i)
  global z_K, global g_K_b_bar, global K_o;
  E_K = nernstPotential(z_K, K_i, K_o)
  I_K_b = g_K_b_bar*(V - E_K);
endfunction


function xdot = f(x, t)

  global F, global R, global T;
  global C_m, global vol_i;

  xdot = zeros(5, 1);

  V    = x(1);
  Na_i = x(2);
  K_i  = x(3);
  a_ur = x(4);
  I_ur = x(5);

  # External stimulation
  t_stim = 1.0;
  t_cycle = 1.0;
  I_stim_bar = 0.0;

  # Background sodium
  I_Na_b = backgroundSodium(V, Na_i);

  # Background potassium
  I_K_b = backgroundPotassium(V, K_i);



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
x0 = [V0, Na_i_0, K_i_0, a_ur_0, I_ur_0];
x = lsode("f", x0, t);

clf
figure(1)
plot(t, x(:, 1))
#subplot(3, 2, 1), plot(t, x(:, 1)), legend('Membrane Voltage (mV)')
#subplot(3, 2, 2), plot(t, x(:, 2)), legend('[Na^{+}]_{i} (mM/l)')
#subplot(3, 2, 3), plot(t, x(:, 3)), legend('[Na^{+}]_{c} (mM/l)')
#subplot(3, 2, 4), plot(t, x(:, 4)), legend('a_{{ur}_{0}}')
#subplot(3, 2, 5), plot(t, x(:, 5)), legend('I_{{ur}_{0}} (pA)')
xlabel('Time (s)')
hold on

print -depsc2 "output.eps"
#print -depslatex "membrane_voltage.tex"
