# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Background sodium current from "Ionic channels of excitable
# membranes," B. Hille.
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

# Background potassium current from "Ionic channels of excitable
# membranes," B. Hille.
function I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar)
  global enable_I_K_b;
  if (enable_I_K_b == true)
    global z_K;
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_b = g_K_b_bar*(V - E_K);
  else
    I_K_b = 0.0;
  endif
endfunction

# Background chloride current from "Ionic channels of excitable
# membranes," B. Hille.
function I_Cl_b = backgroundChloride(V, Cl_i)
  global enable_I_Cl_b;
  if (enable_I_Cl_b == true)
    global z_Cl, global g_Cl_b_bar, global Cl_o;
    E_Cl = nernstPotential(z_Cl, Cl_o, Cl_i);
    I_Cl_b = g_Cl_b_bar*(V - E_Cl);
  else
    I_Cl_b = 0.0;
  endif
endfunction
