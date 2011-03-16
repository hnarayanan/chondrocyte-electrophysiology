1;

# Extract solution components
V    = x(:, 1);
Na_i = x(:, 2);
K_i  = x(:, 3);
Ca_i = x(:, 4);

# Compute currents at all times
I_Na_b = zeros(len_t, 1);
I_K_b  = zeros(len_t, 1);
I_NaK  = zeros(len_t, 1);
I_NaCa = zeros(len_t, 1);
I_NaH  = zeros(len_t, 1);
I_ASIC = zeros(len_t, 1);

for ii = [1:len_t]
  I_Na_b(ii) = backgroundSodium(V(ii), Na_i(ii));
  I_K_b(ii)  = backgroundPotassium(V(ii), K_i(ii));
  I_NaK(ii)  = sodiumPotassiumPump(V(ii), Na_i(ii), K_i(ii));
  I_NaCa(ii) = sodiumCalciumExchanger(V(ii), Na_i(ii), Ca_i(ii));
  I_NaH(ii)  = sodiumHydrogenAntiport();
  I_ASIC(ii) = voltageActivatedHydrogen();
endfor