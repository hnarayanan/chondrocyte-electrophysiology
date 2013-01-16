# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Extract solution components
len_t = size(t, 2);

global apply_Vm;
if (apply_Vm == true)
  V = zeros(len_t, 1);
  for ii = [1:len_t]
    V(ii) = appliedVoltage(t(ii));
  endfor
else
  V = x(:, 1);
endif
Na_i = x(:, 2);
K_i  = x(:, 3);
Ca_i = x(:, 4);
H_i  = x(:, 5);
Cl_i = x(:, 6);
a_ur = x(:, 7);
I_ur = x(:, 8);
vol_i = x(:, 9);
cal  = x(:, 10);

format long;
[V(end) Na_i(end) K_i(end) Ca_i(end) H_i(end) Cl_i(end) a_ur(end) I_ur(end) cal(end)]'

# Compute currents at all times
K_o        = zeros(len_t, 1);
I_Na_b     = zeros(len_t, 1);
I_K_b      = zeros(len_t, 1);
I_Cl_b     = zeros(len_t, 1);
I_NaK      = zeros(len_t, 1);
I_NaCa     = zeros(len_t, 1);
I_NaH      = zeros(len_t, 1);
I_Ca_ATP   = zeros(len_t, 1);
I_K_ur     = zeros(len_t, 1);
I_K_ur_ref = zeros(len_t, 1);
I_K_2pore  = zeros(len_t, 1);
I_K_Ca_act = zeros(len_t, 1);
I_K_ATP    = zeros(len_t, 1);
I_ASIC     = zeros(len_t, 1);
I_TRP1     = zeros(len_t, 1);
I_TRP2     = zeros(len_t, 1);
I_stim     = zeros(len_t, 1);

for ii = [1:len_t]
  K_o(ii)        = appliedPotassiumConcentration(t(ii));
  I_Na_b(ii)     = backgroundSodium(V(ii), Na_i(ii));
  I_K_b(ii)      = backgroundPotassium(V(ii), K_i(ii), K_o(ii), g_K_b_bar);
  I_Cl_b(ii)     = backgroundChloride(V(ii), Cl_i(ii));
  I_NaK(ii)      = sodiumPotassiumPump(V(ii), Na_i(ii), K_i(ii), K_o(ii));
  I_NaCa(ii)     = sodiumCalciumExchanger(V(ii), Na_i(ii), Ca_i(ii));
  I_NaH(ii)      = sodiumHydrogenExchanger(Na_i(ii), H_i(ii));
  I_Ca_ATP(ii)   = calciumPump(Ca_i(ii));
  I_K_ur(ii)     = ultrarapidlyRectifyingPotassium(V(ii), K_i(ii), K_o(ii), a_ur(ii), I_ur(ii));
  I_K_ur_ref(ii) = ultrarapidlyRectifyingPotassium_ref(V(ii), K_i(ii), K_o(ii));
  I_K_2pore(ii)  = twoPorePotassium(V(ii), K_i(ii), K_o(ii), P_K);
  I_K_Ca_act(ii) = calciumActivatedPotassium(V(ii), K_i(ii), K_o(ii), Ca_i(ii), Gmax);
  I_K_ATP(ii)    = potassiumPump(V(ii), K_i(ii), K_o(ii));
  I_ASIC(ii)     = voltageActivatedHydrogen();
  I_TRP1(ii)     = stretchActivatedTrip(V(ii));
  I_TRP2(ii)     = osteoArthriticTrip();
  I_stim(ii)     = externalStimulation(t(ii));
endfor

# Total ionic currents
I_i = I_Na_b + I_K_b + I_Cl_b \
      + I_NaK + I_NaCa + I_Ca_ATP \
      + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP \
      + I_ASIC + I_TRP1 + I_TRP2;
