1;

# Ultra-rapidly rectifying potassium
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

# FIXME: Implement the ATP-powered potassium pump
function I_K_ATP = potassiumPump()
  I_K_ATP = 0.0;
endfunction