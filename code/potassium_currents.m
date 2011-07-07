# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Ultra-rapidly rectifying potassium channel from "Action potential rate
# dependence in the human atrial myocyte," M. M. Maleckar, J. L.
# Greenstein, W. R. Giles and N. A. Trayanova. Am. J. Physiol. Heart.
# Circ. Physiol. 2009; 297; 1398-1410 (Appendix, pp. 1408)

function [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V)
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 30.0)/10.0));
  i_ur_inf   = 1.0/(1.0 + exp((V + 7.5)/10.0));
  tau_a_ur   = 0.005/(1.0 + exp((V + 5.0)/12.0));
  tau_i_ur   = 0.59/(1.0 + exp((V + 10.0)/24.0)) + 0.01;
endfunction

function I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, a_ur, i_ur)
  global enable_I_K_ur;
  if (enable_I_K_ur == true)
    global z_K, global g_K_ur, global K_o;
    [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur     = g_K_ur*a_ur*i_ur*(V - E_K);
  else
    I_K_ur = 0.0;
  endif
endfunction

# Two-pore potassium current
function I_K_2pore = twoPorePotassium(V, K_i, P_K)
  global enable_I_K_2pore;
  if (enable_I_K_2pore == true)
    global F, global R, global T;
    global z_K, global I_K_2pore_0, global K_o;
    I_K_2pore = P_K*z_K^2*V*F^2/(R*T)*(K_i - K_o*exp(-z_K*V*F/(R*T)))/(1 - exp(-z_K*V*F/(R*T))) + I_K_2pore_0;
  else
    I_K_2pore = 0.0;
  endif
endfunction

# Calcium-activated potassium current
# FIXME: Clean up the following
# FIXME: Check the following carefully
function I_K_Ca_act = calciumActivatedPotassium(V, K_i, Ca_i, Gmax)
  global enable_I_K_Ca_act;
  if (enable_I_K_Ca_act == true)
    global T;
    global Zj, global Vhj, global ZL, global L0, global KDc;
    global C, global D, global E;
    global N_channel, global z_K, global K_o;
    kTe = 23.54*(T/273);
    Lv = L0*exp((V*ZL)/kTe);
    Jv = exp(((V - Vhj)*Zj)/kTe);
    K = Ca_i/KDc;
    P0 = (Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_Ca_act = N_channel*P0*Gmax*(V - E_K);
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