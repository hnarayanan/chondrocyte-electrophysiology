# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Extract solution components
len_t = size(t, 2);
x = ones(len_t, 7);
x(:, 1) = estimates.x(1, :);
x(:, 3) = estimates.x(2, :);

# Extract parameters
disp('Estimated parameters and bounding box')
[estimates.parest estimates.bbox]
g_K_b_bar = estimates.parest(1);

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
a_ur = x(:, 6);
I_ur = x(:, 7);

# Compute currents at all times
I_Na_b     = zeros(len_t, 1);
I_K_b      = zeros(len_t, 1);
I_NaK      = zeros(len_t, 1);
I_NaCa     = zeros(len_t, 1);
I_NaH      = zeros(len_t, 1);
I_K_ur     = zeros(len_t, 1);
I_K_2pore  = zeros(len_t, 1);
I_K_Ca_act = zeros(len_t, 1);
I_K_ATP    = zeros(len_t, 1);
I_ASIC     = zeros(len_t, 1);
I_TRP1     = zeros(len_t, 1);
I_TRP2     = zeros(len_t, 1);
I_stim     = zeros(len_t, 1);

for ii = [1:len_t]
  I_Na_b(ii)     = backgroundSodium(V(ii), Na_i(ii));
  I_K_b(ii)      = backgroundPotassium(V(ii), K_i(ii), g_K_b_bar);
  I_NaK(ii)      = sodiumPotassiumPump(V(ii), Na_i(ii), K_i(ii));
  I_NaCa(ii)     = sodiumCalciumExchanger(V(ii), Na_i(ii), Ca_i(ii));
  I_NaH(ii)      = sodiumHydrogenExchanger(Na_i(ii), H_i(ii));
  I_K_ur(ii)     = ultrarapidlyRectifyingPotassium(V(ii), K_i(ii), a_ur(ii), I_ur(ii));
  I_K_2pore(ii)  = twoPorePotassium(V(ii), K_i(ii));
  I_K_Ca_act(ii) = calciumActivatedPotassium(V(ii), K_i(ii), Ca_i(ii));
  I_K_ATP(ii)    = potassiumPump();
  I_ASIC(ii)     = voltageActivatedHydrogen();
  I_TRP1(ii)     = stretchActivatedTrip();
  I_TRP2(ii)     = osteoArthriticTrip();
  I_stim(ii)     = externalStimulation(t(ii));
endfor

# Total ionic currents
I_i =  I_Na_b + I_K_b \
      + I_NaK + I_NaCa + I_NaH \
      + I_K_ur + I_K_2pore + I_K_Ca_act + I_K_ATP \
      + I_ASIC + I_TRP1 + I_TRP2;