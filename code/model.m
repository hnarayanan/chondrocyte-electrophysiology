# The solution fields and parameters have the following units:
#
# Voltage: xV
# Current: xa
# Time: xs
# Concentration: xM
# Conductance: xS
# Capacitance: xF

function xdot = f (x, t)

  xdot = zeros (4, 1);

  V = x(1);
  m = x(2);
  h = x(3);
  n = x(4);

  # Membrane capacitance
  C_m = 15.0;

  # External stimulation
  t_stim = 0.1;
  t_cycle = 0.5;
  I_stim_bar = 100;

  # Background sodium
  g_Na_b_bar = 120;
  g_Na_b = g_Na_b_bar;
  V_Na_b = -115;
  I_Na_b = g_Na_b*(V - V_Na_b);

  # Background potassium
  g_K_b_bar = 36;
  g_K_b = g_K_b_bar;
  V_K_b = 12;
  I_K_b = g_K_b*(V - V_K_b);

  # Sodium-potassium pump
  I_NaK      = 0.0;

  # Sodium-Calcium exchanger
  I_NaCa     = 0.0;

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

  I_i = I_Na_b + I_K_b + I_NaK + I_NaCa + I_NaH + I_K_ur + I_K_2pore + I_Ca_act_K + I_TRP;
  I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);

  alpha_m = 0.1*(V + 25)/(exp((V + 25)/10) - 1);
  beta_m = 4*exp(V/18);

  alpha_h = 0.07*exp(V/20);
  beta_h = 1/(exp((V + 30)/10) + 1);

  alpha_n = 0.01*(V + 10)/(exp((V + 10)/10) - 1);
  beta_n = 0.125*exp(V/80);

  xdot(1) = 1/C_m*(-I_i + I_stim);
  xdot(2) = alpha_m*(1 - m) - beta_m*m;
  xdot(3) = alpha_h*(1 - h) - beta_h*h;
  xdot(4) = alpha_n*(1 - n) - beta_n*n;

endfunction

% Time stepping information
t_final = 4
dt = 0.01

% Initial conditions
V0 = -86
m0 =  1
h0 =  0
n0 =  1

t = linspace(0, t_final, t_final/dt);
x0 = [V0, m0, h0, n0];
x = lsode("f", x0, t);

clf
figure(1)
plot(t, x(:,1))
%subplot(2,1,1), plot(t, x(:,1)), xlabel('time'), ylabel('Voltage')
%hold on
%subplot(2,1,2), plot(t, x(:,[2:4])), legend('m','h','n')