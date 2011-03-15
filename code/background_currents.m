1;

# Background sodium current
function I_Na_b = backgroundSodium(V, Na_i)
  global z_Na, global g_Na_b_bar, global Na_o;
  E_Na = nernstPotential(z_Na, Na_i, Na_o);
  I_Na_b = g_Na_b_bar*(V - E_Na);
endfunction

# Background potassium current
function I_K_b = backgroundPotassium(V, K_i)
  global z_K, global g_K_b_bar, global K_o;
  E_K = nernstPotential(z_K, K_i, K_o);
  I_K_b = g_K_b_bar*(V - E_K);
endfunction