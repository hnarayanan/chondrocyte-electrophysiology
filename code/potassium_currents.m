1;

# Ultra-rapidly rectifying potassium
# FIXME: Insert Maleckar 2009 citation here
function [a_ur_inf, I_ur_inf, tau_a_ur, tau_I_ur] = ultraRapidlyRectifyingPotassiumHelper(V)
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 6.0)/22.6));
  I_ur_inf   = 1.0/(1.0 + exp((V + 7.5)/25.0));
  tau_a_ur   = 0.009/(1.0 + exp((V + 5.0)/12.0)) + 0.0005;
  tau_I_ur   = 0.59/(1.0 + exp((V + 60.0)/10.0)) + 3.05;
endfunction

function I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, a_ur, I_ur)
  global enable_I_K_ur;
  if (enable_I_K_ur == true)
    global z_K, global g_K_ur, global K_o;
    [a_ur_inf, I_ur_inf, tau_a_ur, tau_I_ur] = ultraRapidlyRectifyingPotassiumHelper(V);
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur     = g_K_ur*a_ur*I_ur*(V - E_K);
  else
    I_K_ur = 0.0;
  endif
endfunction

# Two-pore potassium current
function I_K_2pore = twoPorePotassium(V, K_i)
  global enable_I_K_2pore;
  if (enable_I_K_2pore == true)
    global F, global R, global T;
    global z_K, global P_K, global K_o;
    I_K_2pore = P_K*z_K^2*V*F^2/(R*T)*(K_i - K_o*exp(-z_K*V*F/(R*T)))/(1 - exp(-z_K*V*F/(R*T)));
  else
    I_K_2pore = 0.0;
  endif
endfunction

# Calcium-activated potassium current
# FIXME: Clean up the following
# FIXME: Check the following carefully
function I_K_Ca_act = calciumActivatedPotassium(V, K_i, Ca_i)
  global enable_I_K_Ca_act;
  if (enable_I_K_Ca_act == true)
    global T;
    global Zj, global Vhj, global ZL, global L0, global KDc;
    global C, global D, global E;
    global Gmax, global N_channel, global V_K_Ca_act;
    global z_K, global g_K_b_bar, global K_o;
    W = [Zj, Vhj, ZL, L0, KDc, C, D, E, Gmax];
    kTe = 23.54*(T/273);
    Lv = W(4)*exp((V*W(3))/kTe);
    Jv = exp(((V - W(2))*W(1))/kTe);
    K = Ca_i/W(5);
    C = W(6);
    D = W(7);
    E = W(8);
    G_max = W(9);
    P0 = (Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_Ca_act = N_channel*P0*G_max*(V - E_K);
  else
    I_K_Ca_act = 0.0;
  endif
endfunction


# FIXME: Implement the ATP-powered potassium pump
function I_K_ATP = potassiumPump()
  global enable_I_K_ATP;
  if (enable_I_K_ATP == true)
    I_K_ATP = 0.0;
  else
    I_K_ATP = 0.0;
  endif
endfunction