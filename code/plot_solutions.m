1;

# Clear the screen and plot solutions
clf;
hold on;

# The membrane voltage and total currents
figure(1);
subplot(2, 2, 1), plot(t, V), xlabel('$t (s)$'), legend('$V_{m} (mV)$');
subplot(2, 2, 2), plot(t, I_i), xlabel('$t (s)$'), legend('$I_i (nA)$');
subplot(2, 2, 3), plot(V, I_i), xlabel('$V_{m} (mV)$'), legend('$I_i (nA)$');
print -depslatexstandalone "output/membrane_behaviour.tex"

# The different concentrations
figure(2);
subplot(2, 2, 1), plot(t, Na_i), xlabel('$t (s)$'), legend('$[Na^{+}]_{i} (mM/l)$');
subplot(2, 2, 2), plot(t, K_i), xlabel('$t (s)$'), legend('$[K^{+}]_{i} (mM/l)$');
subplot(2, 2, 3), plot(t, Ca_i), xlabel('$t (s)$'), legend('$[Ca^{2+}]_{i} (mM/l)$');
print -depslatexstandalone "output/concentrations.tex"

# The different background currents
figure(3);
subplot(2, 2, 1), plot(t, I_Na_b), xlabel('$t (s)$'), legend('$I_{Na_{b}} (nA)$');
subplot(2, 2, 2), plot(t, I_K_b), xlabel('$t (s)$'), legend('$I_{K_{b}} (nA)$');
print -depslatexstandalone "output/background_currents.tex"

# The different pump and exchanger currents
figure(4);
subplot(2, 2, 1), plot(t, I_NaK), xlabel('$t (s)$'), legend('$I_{NaK} (nA)$');
subplot(2, 2, 2), plot(t, I_NaCa), xlabel('$t (s)$'), legend('$I_{NaCa} (nA)$');
subplot(2, 2, 3), plot(t, I_NaH), xlabel('$t (s)$'), legend('$I_{NaH} (nA)$');
print -depslatexstandalone "output/pumps_and_exchangers.tex"

# The other potassium currents
figure(5);
subplot(2, 2, 1), plot(t, I_K_ur), xlabel('$t (s)$'), legend('$I_{K_{ur}} (nA)$');
subplot(2, 2, 2), plot(t, I_K_2pore), xlabel('$t (s)$'), legend('$I_{K_{2pore}} (nA)$');
subplot(2, 2, 3), plot(t, I_K_Ca_act), xlabel('$t (s)$'), legend('$I_{K_{Ca-act}} (nA)$');
subplot(2, 2, 4), plot(t, I_K_ATP), xlabel('$t (s)$'), legend('$I_{K_{ATP}} (nA)$');
print -depslatexstandalone "output/potassium_currents.tex"

# The other currents
figure(6);
subplot(2, 2, 1), plot(t, I_stim), xlabel('$t (s)$'), legend('$I_{stim} (nA)$');
subplot(2, 2, 2), plot(t, I_ASIC), xlabel('$t (s)$'), legend('$I_{ASIC} (nA)$');
subplot(2, 2, 3), plot(t, I_TRP1), xlabel('$t (s)$'), legend('$I_{TRP1} (nA)$');
subplot(2, 2, 4), plot(t, I_TRP2), xlabel('$t (s)$'), legend('$I_{TRP2} (nA)$');
print -depslatexstandalone "output/other_currents.tex"