# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Sodium-potassium pump from "Mathematical Model of an Adult Human
# Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
# Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
# Giles. Circ. Res. 1998; 82; 63-81 (Table 12, pp. 77)

function I_NaK = sodiumPotassiumPump(V, Na_i, K_i, K_o)
  global enable_I_NaK;
  if (enable_I_NaK == true)
    global I_NaK_bar K_NaK_K K_NaK_Na;
    I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) \
        *(Na_i^1.5/(Na_i^1.5 + K_NaK_Na^1.5)) \
        *(V + 150.0)/(V + 200.0);
  else
    I_NaK = 0.0;
  endif
endfunction

# Sodium-calcium exchanger from "Mathematical Model of an Adult Human
# Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
# Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
# Giles. Circ. Res. 1998; 82; 63-81 (Table 13, pp. 77)

function I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i)
  global enable_I_NaCa;
  if (enable_I_NaCa == true)
    global F, global R, global T;
    global Na_o, global Ca_o;
    global K_NaCa, global gamma_Na, global d_NaCa;

    I_NaCa = K_NaCa*(  Na_i^3*Ca_o*exp(gamma_Na*V*F/(R*T)) \
                     - Na_o^3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
             /(1.0 + d_NaCa*(Na_o^3*Ca_i + Na_i^3*Ca_o));
  else
    I_NaCa = 0.0;
  endif
endfunction

# Sodium-hydrogen exchanger from "A Model of Na+/H+ Exchanger and Its
# Central Role in Regulation of pH and Na+ in Cardiac Myocytes," Chae
# Young Cha, Chiaki Oka, Yung E. Earm, Shigeo Wakabayashi, and Akinori
# Noma. Biophysical Journal 2009; 97; 2674-2683 (pp. 2675)

function I_NaH = sodiumHydrogenExchanger(Na_i, H_i)
  global enable_I_NaH;
  if (enable_I_NaH == true)
    global n_H, global K_H_i_mod;
    global k1_p, global k1_m, global k2_p, global k2_m;
    global Na_o, global H_o, global N_NaH_channel;
    global K_Na_o, global K_H_o, global K_Na_i, global K_H_i;

    I_NaH_mod  = 1/(1 + (K_H_i_mod^n_H/H_i^n_H));
    t1 = k1_p*Na_o/K_Na_o / (1 + Na_o/K_Na_o + H_o/K_H_o);
    t2 = k2_p*H_i/K_H_i   / (1 + Na_i/K_Na_i + H_i/K_H_i);
    t3 = k1_m*Na_i/K_Na_i / (1 + Na_i/K_Na_i + H_i/K_H_i);
    t4 = k2_m*H_o/K_H_o   / (1 + Na_o/K_Na_o + H_o/K_H_o);
    I_NaH_exch = (t1*t2 - t3*t4) / (t1 + t2 + t3 + t4);
    I_NaH = N_NaH_channel*I_NaH_mod*I_NaH_exch;
  else
    I_NaH = 0.0;
  endif
endfunction
