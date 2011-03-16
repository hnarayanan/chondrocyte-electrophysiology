1;

# Ultra-rapidly rectifying potassium
# FIXME: Insert Maleckar 2009 citation here
function [a_ur_inf, I_ur_inf, tau_a_ur, tau_I_ur] = ultraRapidlyRectifyingPotassiumHelper(V)
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 6.0)/8.6));
  I_ur_inf   = 1.0/(1.0 + exp((V + 7.5)/10.0));
  tau_a_ur   = 0.009/(1.0 + exp((V + 5.0)/12.0)) + 0.0005;
  tau_I_ur   = 0.59/(1.0 + exp((V + 60.0)/10.0)) + 3.05;
endfunction

function I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, a_ur, I_ur)
  global z_K, global g_K_ur, global K_o;
  [a_ur_inf, I_ur_inf, tau_a_ur, tau_I_ur] = ultraRapidlyRectifyingPotassiumHelper(V);
  E_K        = nernstPotential(z_K, K_i, K_o);
  I_K_ur     = g_K_ur*a_ur*I_ur*(V - E_K);
endfunction

# Two-pore potassium current
function I_K_2pore = twoPorePotassium(V, K_i)
  global z_K, global g_K_2pore, global K_o;
  E_K = nernstPotential(z_K, K_i, K_o);
  I_K_2pore = g_K_2pore*(V - E_K);
endfunction

# Calcium-activated potassium current
# FIXME: Clean up the following
# FIXME: Check the following carefully
function I_K_Ca_act = calciumActivatedPotassium(V, Ca_i)
  global T;
  global Zj, global Vhj, global ZL, global L0, global KDc;
  global C, global D, global E;
  global Gmax, global N_channel, global V_K_Ca_act;
  W = [Zj, Vhj, ZL, L0, KDc, C, D, E, Gmax];
  kTe = 23.54*(T/273);
  Lv = W(4)*exp((V*W(3))/kTe);
  Jv = exp(((V - W(2))*W(1))/kTe);
  K = Ca_i/W(5);
  C = W(6);
  D = W(7);
  E = W(8);
  G_max = W(9);
  P0=(Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
  I_K_Ca_act = N_channel*P0*G_max*(V - V_K_Ca_act);
endfunction


# FIXME: Implement the ATP-powered potassium pump
function I_K_ATP = potassiumPump()
  I_K_ATP = 0.0;
endfunction