# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Clear the screen and plot solutions
clf;
line_width = 3;

# Load reference solution
V_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 1);
I_i_by_C_m_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 2);

# The membrane voltage and total currents
figure(1, 'visible', 'off');
subplot(2, 2, 1), plot(t, V, 'linewidth', line_width), xlabel('$t (s)$'), legend('$V_{\mathrm{m}} (mV)$');
subplot(2, 2, 2), plot(t, I_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{i}} (pA)$');
subplot(2, 2, 3), plot(V, I_i/C_m, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{i}}/C_{\mathrm{m}} (pA/pF)$');
hold on;
plot(V_ref, I_i_by_C_m_ref, '1', 'linewidth', line_width);
print -depslatexstandalone "../results/epslatex/membrane_behaviour.tex"

# The different concentrations
figure(2, 'visible', 'off');
subplot(2, 2, 1), plot(t, Na_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[Na^{+}]_{\mathrm{i}} (mM/l)$');
subplot(2, 2, 2), plot(t, K_i,  'linewidth', line_width), xlabel('$t (s)$'), legend('$[K^{+}]_{\mathrm{i}} (mM/l)$');
subplot(2, 2, 3), plot(t, Ca_i, 'linewidth', line_width), xlabel('$t (s)$'), legend('$[Ca^{2+}]_{\mathrm{i}} (mM/l)$');
subplot(2, 2, 4), plot(t, H_i,  'linewidth', line_width), xlabel('$t (s)$'), legend('$[H^{+}]_{\mathrm{i}} (mM/l)$');
print -depslatexstandalone "../results/epslatex/concentrations.tex"

# The different background currents (t-I)
figure(3, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_Na_b, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{Na_{b}}} (pA)$');
subplot(2, 2, 2), plot(t, I_K_b,  'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{K_{b}}} (pA)$');
hold on;
plot(measure.time, measure.data, 'o');
print -depslatexstandalone "../results/epslatex/background_currents-ti.tex"

# The different background currents (V-I)
figure(4, 'visible', 'off');
subplot(2, 2, 1), plot(V, I_Na_b, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{Na_{b}}} (pA)$');
subplot(2, 2, 2), plot(V, I_K_b,  'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{K_{b}}} (pA)$');
print -depslatexstandalone "../results/epslatex/background_currents-vi.tex"

# The different pump and exchanger currents (t-I)
figure(5, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_NaK,  'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{NaK}} (pA)$');
subplot(2, 2, 2), plot(t, I_NaCa, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{NaCa}} (pA)$');
subplot(2, 2, 3), plot(t, I_NaH,  'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{NaH}} (pA)$');
print -depslatexstandalone "../results/epslatex/pumps_and_exchangers-ti.tex"

# The different pump and exchanger currents (V-I)
figure(6, 'visible', 'off');
subplot(2, 2, 1), plot(V, I_NaK,  'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{NaK}} (pA)$');
subplot(2, 2, 2), plot(V, I_NaCa, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{NaCa}} (pA)$');
subplot(2, 2, 3), plot(V, I_NaH,  'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{NaH}} (pA)$');
print -depslatexstandalone "../results/epslatex/pumps_and_exchangers-vi.tex"

# The other potassium currents (t-I)
figure(7, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_K_ur,     'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{K_{ur}}} (pA)$');
subplot(2, 2, 2), plot(t, I_K_2pore,  'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{K_{2pore}}} (pA)$');
subplot(2, 2, 3), plot(t, I_K_Ca_act, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{K_{Ca-act}}} (pA)$');
subplot(2, 2, 4), plot(t, I_K_ATP,    'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm{K_{ATP}}} (pA)$');
print -depslatexstandalone "../results/epslatex/potassium_currents-ti.tex"

# The other potassium currents (V-I)
figure(8, 'visible', 'off');
subplot(2, 2, 1), plot(V, I_K_ur,     'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{K_{ur}}} (pA)$');
subplot(2, 2, 2), plot(V, I_K_2pore,  'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{K_{2pore}}} (pA)$');
subplot(2, 2, 3), plot(V, I_K_Ca_act, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{K_{Ca-act}}} (pA)$');
subplot(2, 2, 4), plot(V, I_K_ATP,    'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm{K_{ATP}}} (pA)$');
print -depslatexstandalone "../results/epslatex/potassium_currents-vi.tex"

# The other currents (t-I)
figure(9, 'visible', 'off');
subplot(2, 2, 1), plot(t, I_stim, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm stim} (pA)$');
subplot(2, 2, 2), plot(t, I_ASIC, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm ASIC} (pA)$');
subplot(2, 2, 3), plot(t, I_TRP1, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm TRP1} (pA)$');
subplot(2, 2, 4), plot(t, I_TRP2, 'linewidth', line_width), xlabel('$t (s)$'), legend('$I_{\mathrm TRP2} (pA)$');
print -depslatexstandalone "../results/epslatex/other_currents-ti.tex"

# The other currents (V-I)
figure(10, 'visible', 'off');
subplot(2, 2, 1), plot(V, I_stim, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm stim} (pA)$');
subplot(2, 2, 2), plot(V, I_ASIC, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm ASIC} (pA)$');
subplot(2, 2, 3), plot(V, I_TRP1, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm TRP1} (pA)$');
subplot(2, 2, 4), plot(V, I_TRP2, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), legend('$I_{\mathrm TRP2} (pA)$');
print -depslatexstandalone "../results/epslatex/other_currents-vi.tex"

# Some special plots of interest
figure(11, 'visible', 'off');
subplot(2, 2, 1), plot(K_i, V, 'linewidth', line_width), xlabel('$[K^{+}]_{\mathrm{i}} (mM/l)$'), legend('$V_{\mathrm{m}} (mV)$');
print -depslatexstandalone "../results/epslatex/special_plots.tex"
