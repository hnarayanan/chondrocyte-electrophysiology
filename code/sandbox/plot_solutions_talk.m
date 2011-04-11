1;

# Clear the screen and plot solutions
clf;
line_width = 3;

# Load reference solution
V_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 1);
I_i_by_C_m_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 2);

# The membrane voltage and total currents
figure(1, 'visible', 'off');
subplot(2, 2, 1), plot(t, V, 'linewidth', line_width), xlabel('$t (s)$'), legend('$V_{m} (mV)$');
subplot(2, 2, 2), plot(t, K_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[K^{+}]_{i} (mM/l)$');
subplot(2, 2, 3), plot(t, I_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_i (pA)$');
subplot(2, 2, 4), plot(V, I_i/C_m, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_i/C_{m} (pA/pF)$');
hold on;
plot(V_ref, I_i_by_C_m_ref, '1', 'linewidth', line_width);
print -depslatexstandalone "../results/epslatex/membrane_behaviour.tex"

# The different concentrations
figure(2, 'visible', 'off');
subplot(2, 2, 1), plot(t, Na_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[Na^{+}]_{i} (mM/l)$');
subplot(2, 2, 2), plot(t, K_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[K^{+}]_{i} (mM/l)$');
subplot(2, 2, 3), plot(t, Ca_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[Ca^{2+}]_{i} (mM/l)$');
print -depslatexstandalone "../results/epslatex/concentrations.tex"

# The different background currents
figure(3, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_Na_b, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{Na_{b}} (pA)$');
subplot(2, 2, 2), plot(t, I_K_b, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{K_{b}} (pA)$');
print -depslatexstandalone "../results/epslatex/background_currents.tex"

# The different pump and exchanger currents
figure(4, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_NaK, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{NaK} (pA)$');
subplot(2, 2, 2), plot(t, I_NaCa, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{NaCa} (pA)$');
subplot(2, 2, 3), plot(t, I_NaH, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{NaH} (pA)$');
print -depslatexstandalone "../results/epslatex/pumps_and_exchangers.tex"

# The other potassium currents
figure(5, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_K_b, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{K_{b}} (pA)$');
subplot(2, 2, 2), plot(t, I_K_Ca_act, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{K_{Ca-act}} (pA)$');
subplot(2, 2, 3), plot(t, I_K_2pore, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{K_{2pore}} (pA)$');
subplot(2, 2, 4), plot(t, I_K_ur, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{K_{ur}} (pA)$');
print -depslatexstandalone "../results/epslatex/potassium_currents.tex"

# The other currents
figure(6, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_stim, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{stim} (pA)$');
subplot(2, 2, 2), plot(t, I_ASIC, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{ASIC} (pA)$');
subplot(2, 2, 3), plot(t, I_TRP1, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{TRP1} (pA)$');
subplot(2, 2, 4), plot(t, I_TRP2, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{TRP2} (pA)$');

print -depslatexstandalone "../results/epslatex/other_currents.tex"