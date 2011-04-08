# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Background sodium current
function I_Na_b = backgroundSodium(V, Na_i)
  global enable_I_Na_b;
  if (enable_I_Na_b == true)
    global z_Na, global g_Na_b_bar, global Na_o;
    E_Na = nernstPotential(z_Na, Na_i, Na_o);
    I_Na_b = g_Na_b_bar*(V - E_Na);
  else
    I_Na_b = 0.0;
  endif
endfunction

# Background potassium current
function I_K_b = backgroundPotassium(V, K_i)
  global enable_I_K_b;
  if (enable_I_K_b == true)
    global z_K, global g_K_b_bar, global K_o;
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_b = g_K_b_bar*(V - E_K);
  else
    I_K_b = 0.0;
  endif
endfunction