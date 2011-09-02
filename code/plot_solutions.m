# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Clear the screen and define plotting parameters
clf;
line_width = 5;
blue = [0.00, 0.61, 0.99];
red =  [1.00, 0.17, 0.00];

# Plot directly to file
figure(1, 'visible', 'off');

# Load reference solutions
V_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 1);
I_i_by_C_m_ref = csvread('../data/reference_values/Total_Current_Density_vs_Membrane_Voltage.data')(:, 2);

V_ref_without_TEA = csvread('../data/reference_values/Total_Current_vs_Voltage_without_TEA.data')(:, 1);
I_ref_without_TEA = csvread('../data/reference_values/Total_Current_vs_Voltage_without_TEA.data')(:, 2);
V_ref_with_TEA = csvread('../data/reference_values/Total_Current_vs_Voltage_with_TEA.data')(:, 1);
I_ref_with_TEA = csvread('../data/reference_values/Total_Current_vs_Voltage_with_TEA.data')(:, 2);
I_ref_without_TEA_int = interp1(V_ref_without_TEA, I_ref_without_TEA, V);
I_ref_with_TEA_int = interp1(V_ref_with_TEA, I_ref_with_TEA, V);

V_ref_without_BUP = csvread('../data/reference_values/Total_Current_vs_Voltage_without_BUP.data')(:, 1);
I_ref_without_BUP = csvread('../data/reference_values/Total_Current_vs_Voltage_without_BUP.data')(:, 2);
V_ref_with_BUP = csvread('../data/reference_values/Total_Current_vs_Voltage_with_BUP.data')(:, 1);
I_ref_with_BUP = csvread('../data/reference_values/Total_Current_vs_Voltage_with_BUP.data')(:, 2);
I_ref_without_BUP_int = interp1(V_ref_without_BUP, I_ref_without_BUP, V);
I_ref_with_BUP_int = interp1(V_ref_with_BUP, I_ref_with_BUP, V);

V_ref_without_DTX = csvread('../data/reference_values/Total_Current_Density_vs_Voltage_without_DTX.data')(:, 1);
I_ref_density_without_DTX = csvread('../data/reference_values/Total_Current_Density_vs_Voltage_without_DTX.data')(:, 2);
V_ref_with_DTX = csvread('../data/reference_values/Total_Current_Density_vs_Voltage_with_DTX.data')(:, 1);
I_ref_density_with_DTX = csvread('../data/reference_values/Total_Current_Density_vs_Voltage_with_DTX.data')(:, 2);

I_ref_without_DTX_int = interp1(V_ref_without_DTX, I_ref_density_without_DTX, V)*C_m;
I_ref_with_DTX_int = interp1(V_ref_with_DTX, I_ref_density_with_DTX, V)*C_m;

# Plot the membrane voltage and total currents
plot(t, V, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $V_{\mathrm{m}}\,(mV)$');
print -depslatexstandalone "../results/epslatex/t-V.tex"
plot(t, I_i, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{i}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_i.tex"
plot(V, I_i/C_m, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{i}}/C_{\mathrm{m}}\,(pA/pF)$');
hold on;
plot(V_ref, I_i_by_C_m_ref, '1', 'linewidth', line_width, 'color', red);
hold off;
print -depslatexstandalone "../results/epslatex/t-I_i_by_Cm.tex"

# Plot the different concentrations
plot(t, Na_i, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[Na^{+}]_{\mathrm{i}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-Na_i.tex"
plot(t, K_i,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[K^{+}]_{\mathrm{i}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-K_i.tex"
plot(t, Ca_i, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[Ca^{2+}]_{\mathrm{i}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-Ca_i.tex"
plot(t, H_i,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[H^{+}]_{\mathrm{i}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-H_i.tex"
plot(t, Cl_i,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[Cl^{-}]_{\mathrm{i}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-Cl_i.tex"

# Plot the different background currents (t-I)
plot(t, I_Na_b, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{Na_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_Na_b.tex"
plot(t, I_K_b,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_K_b.tex"
plot(t, I_Cl_b, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{Cl_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_Cl_b.tex"

# Plot the different background currents (V-I)
plot(V, I_Na_b, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{Na_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_Na_b.tex"
plot(V, I_K_b,  'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_K_b.tex"
plot(V, I_Cl_b, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{Cl_{b}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_Cl_b.tex"

# Plot the different pump and exchanger currents (t-I)
plot(t, I_NaK,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{NaK}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_NaK.tex"
plot(t, I_NaCa, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{NaCa}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_NaCa.tex"
plot(t, I_NaH,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{NaH}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_NaH.tex"

# Plot the different pump and exchanger currents (V-I)
plot(V, I_NaK,  'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{NaK}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_NaK.tex"
plot(V, I_NaCa, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{NaCa}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_NaCa.tex"
plot(V, I_NaH,  'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{NaH}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_NaH.tex"

# Plot the other potassium currents (t-I)
plot(t, I_K_ur,     'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{ur}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_K_ur.tex"
plot(t, I_K_2pore,  'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{2pore}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_K_2pore.tex"
plot(t, I_K_Ca_act, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{Ca-act}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_K_Ca_act.tex"
plot(t, I_K_ATP,    'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{ATP}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_K_ATP.tex"

# Plot the other potassium currents (V-I)
plot(V, I_K_ur,     'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{ur}}}\,(pA)$');
hold on;
plot(V, I_ref_without_DTX_int - I_ref_with_DTX_int, '1', 'linewidth', line_width, 'color', red);
hold off;
print -depslatexstandalone "../results/epslatex/V-I_K_ur.tex"
plot(V, I_K_2pore,  'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{2pore}}}\,(pA)$');
hold on;
plot(V, I_ref_without_BUP_int - I_ref_with_BUP_int, '1', 'linewidth', line_width, 'color', red);
hold off;
print -depslatexstandalone "../results/epslatex/V-I_K_2pore.tex"
plot(V, I_K_Ca_act, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{Ca-act}}}\,(pA)$');
hold on;
plot(V, I_ref_without_TEA_int - I_ref_with_TEA_int, '1', 'linewidth', line_width, 'color', red);
hold off;
print -depslatexstandalone "../results/epslatex/V-I_K_Ca_act.tex"
plot(V, I_K_ATP,    'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm{K_{ATP}}}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_K_ATP.tex"

# Plot the other currents (t-I)
plot(t, I_stim, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm stim}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_stim.tex"
plot(t, I_ASIC, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm ASIC}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_ASIC.tex"
plot(t, I_TRP1, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm TRP1}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_TRP1.tex"
plot(t, I_TRP2, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $I_{\mathrm TRP2}\,(pA)$');
print -depslatexstandalone "../results/epslatex/t-I_TRP2.tex"

# Plot the other currents (V-I)
plot(V, I_stim, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm stim}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_stim.tex"
plot(V, I_ASIC, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm ASIC}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_ASIC.tex"
plot(V, I_TRP1, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm TRP1}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_TRP1.tex"
plot(V, I_TRP2, 'linewidth', line_width, 'color', blue), xlabel('\Huge $V_{m}\,(mV)$'), grid(), ylabel('\Huge $I_{\mathrm TRP2}\,(pA)$');
print -depslatexstandalone "../results/epslatex/V-I_TRP2.tex"

# Some special plots of interest
plot(K_o, V, 'linewidth', line_width, 'color', blue), xlabel('\Huge $[K^{+}]_{\mathrm{o}}\,(mM/l)$'), grid(), ylabel('\Huge $V_{\mathrm{m}}\,(mV)$');
print -depslatexstandalone "../results/epslatex/K_o-V.tex"
plot(t, K_o, 'linewidth', line_width, 'color', blue), xlabel('\Huge $t\,(s)$'), grid(), ylabel('\Huge $[K^{+}]_{\mathrm{o}}\,(mM/l)$');
print -depslatexstandalone "../results/epslatex/t-K_o.tex"