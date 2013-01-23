# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Ultra-rapidly rectifying potassium channel from "Action potential rate
# dependence in the human atrial myocyte," M. M. Maleckar, J. L.
# Greenstein, W. R. Giles and N. A. Trayanova. Am. J. Physiol. Heart.
# Circ. Physiol. 2009; 297; 1398-1410 (Appendix, pp. 1408)

function [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V)
  a_ur_inf   = 1.0/(1.0 + exp(-(V + 26.7)/4.1));
  i_ur_inf   = 1.0/(1.0 + exp((V - 30.0)/10.0));
  tau_a_ur   = 0.005/(1.0 + exp((V + 5.0)/12.0));
  tau_i_ur   = 0.59/(1.0 + exp((V + 10.0)/24.0)) + 0.01;
endfunction

function I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur)
  global enable_I_K_ur;
  if (enable_I_K_ur == true)
    global z_K g_K_ur;
    [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur     = g_K_ur*a_ur*(V - E_K);
  else
    I_K_ur = 0.0;
  endif
endfunction

# From Bob Clark et al., J. Physiol. 2011, Figure 4
function I_K_ur_ref = ultrarapidlyRectifyingPotassium_ref(V, K_i, K_o)
    G_K = 28.9;  # pS/pF
    V_h = -26.7; # mV
    S_h = 4.1;   # mV
    global C_m z_K;
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur_ref = G_K*(V - E_K)/(1 + exp(-(V - V_h)/S_h)) * C_m/1000.0;
endfunction

# Two-pore potassium current
function I_K_2pore = twoPorePotassium(V, K_i, K_o, P_K)
  global enable_I_K_2pore;
  if (enable_I_K_2pore == true)
    global F R T z_K I_K_2pore_0;
    I_K_2pore = P_K*z_K^2*V*F^2/(R*T)*(K_i - K_o*exp(-z_K*V*F/(R*T)))/(1 - exp(-z_K*V*F/(R*T))) + I_K_2pore_0;
  else
    I_K_2pore = 0.0;
  endif
endfunction

# Calcium-activated potassium current
# FIXME: Clean up the following
# FIXME: Check the following carefully
function I_K_Ca_act = calciumActivatedPotassium(V, K_i, K_o, Ca_i, Gmax)
  global enable_I_K_Ca_act;
  if (enable_I_K_Ca_act == true)
    global T Zj Vhj ZL L0 KDc C D E N_channel z_K E_K_Ca_act;
    kTe = 23.54*(T/273);
    Lv = L0*exp((V*ZL)/kTe);
    Jv = exp(((V - Vhj)*Zj)/kTe);
    K = Ca_i/KDc;
    P0 = (Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
    # E_K = nernstPotential(z_K, K_i, K_o)
    I_K_Ca_act_temp = N_channel*P0*Gmax*(V - E_K_Ca_act);
    if(I_K_Ca_act_temp < 0.0)
      I_K_Ca_act = 0.0;
    else
      I_K_Ca_act = I_K_Ca_act_temp;
    endif
  else
    I_K_Ca_act = 0.0;
  endif
endfunction

# ATP-powered potassium pump from "Modeling Cardiac Action Potential
# Shortening Driven by Oxidative Stress-Induced Mitochondrial
# Oscillations in Guinea Pig Cardiomyocytes," L. Zhou, S. Cortassa, A-C.
# Wei, M. A. Aon, R. L. Winslow, and B. O'Rourke. Biophys. J. 2009; 97;
# 1843-1852 (p. 1845)
# FIXME: Clean up the following
# FIXME: Check the following carefully

function I_K_ATP = potassiumPump(V, K_i, K_o)
  global enable_I_K_ATP;
  if (enable_I_K_ATP == true)
    sigma   = 0.6;
    g_0     = 30.95/400; # FIXME: Somewhat arbitrary. Scaled this down to match Zhou/Ferrero.
    p_0     = 0.91;
    H_K_ATP = -0.001;
    K_m_ATP = 0.56;
    surf    = 1;

    global V_0;
    ADP_i = 10;
    ATP_i = V - V_0 + ADP_i; # FIXME: Completely arbitrary

    H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i);
    K_m = 35.8 + 17.9*ADP_i^K_m_ATP;
    f_ATP = 1.0/(1.0 + (ATP_i/K_m)^H);

    global z_K;
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_ATP = sigma*g_0*p_0*f_ATP*(V - E_K);
  else
    I_K_ATP = 0.0;
  endif
endfunction
