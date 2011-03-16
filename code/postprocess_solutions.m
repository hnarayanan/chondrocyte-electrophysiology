1;

# Extract solution components
V    = x(:, 1);
Na_i = x(:, 2);
K_i  = x(:, 3);
Ca_i = x(:, 4);
a_ur = x(:, 5);
I_ur = x(:, 6);


# Compute currents at all times
I_Na_b    = zeros(len_t, 1);
I_K_b     = zeros(len_t, 1);
I_NaK     = zeros(len_t, 1);
I_NaCa    = zeros(len_t, 1);
I_NaH     = zeros(len_t, 1);
I_K_ur    = zeros(len_t, 1);
I_K_2pore = zeros(len_t, 1);
I_K_ATP   = zeros(len_t, 1);
I_ASIC    = zeros(len_t, 1);
I_stim    = zeros(len_t, 1);

for ii = [1:len_t]
  I_Na_b(ii)     = backgroundSodium(V(ii), Na_i(ii));
  I_K_b(ii)      = backgroundPotassium(V(ii), K_i(ii));
  I_NaK(ii)      = sodiumPotassiumPump(V(ii), Na_i(ii), K_i(ii));
  I_NaCa(ii)     = sodiumCalciumExchanger(V(ii), Na_i(ii), Ca_i(ii));
  I_NaH(ii)      = sodiumHydrogenAntiport();
  I_K_ATP(ii)    = potassiumPump();
  I_K_ur(ii)     = ultrarapidlyRectifyingPotassium(V(ii), K_i(ii), a_ur(ii), I_ur(ii));
  I_K_2pore(ii)  = twoPorePotassium(V(ii), K_i(ii));
  I_ASIC(ii)     = voltageActivatedHydrogen();
  I_stim(ii)     = externalStimulation(t(ii));
endfor