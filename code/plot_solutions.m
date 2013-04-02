# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Clear the screen and define plotting parameters
clf;
line_width = 5;
blue = [0.00, 0.61, 0.99];
red =  [1.00, 0.17, 0.00];

# Plot directly to file
h = figure(1, 'visible', 'off');
set (h,'papertype', '<custom>');
set (h,'paperunits','inches');
set (h, "paperorientation", "landscape");
papersize = [5.905, 4.134];
set (h, "papersize", papersize);
set (h, "paperposition", [0.25 0.25, papersize - 0.5]);
set (0,'defaultaxesfontsize', 14);

I_i_by_Cm_ref = csvread('../data/reference_values/I_i.data');
I_K_2pore_ref = csvread('../data/reference_values/I_K_2pore.data');
I_K_Ca_act_ref = csvread('../data/reference_values/I_K_Ca_act.data');
I_TRPv4_ref = csvread('../data/reference_values/I_TRPv4.data');

# Plot the membrane voltage and total currents
plot(t, V, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-V.tex"
plot(t, I_i, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_i.tex"
h1 = errorbar(I_i_by_Cm_ref(:, 1), I_i_by_Cm_ref(:, 2), I_i_by_Cm_ref(:, 3));
set(h1(1), "color", red)
set(h1(1), "linewidth", line_width)
hold on;
plot(V(5:end), I_i(5:end)/C_m, 'linewidth', line_width, 'color', blue);
hold off;
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_i_by_Cm.tex"
V_comp = V(5:end);
I_i_by_C_m = I_i(5:end)/C_m;
if (enable_I_TRP1 == false)
  save computation_no_TRP V_comp I_i_by_C_m;
else
  save computation_TRP V_comp I_i_by_C_m;
endif
# Plot the different concentrations
plot(t, Na_i, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-Na_i.tex"
plot(t, K_i,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-K_i.tex"
plot(t, Ca_i*1.e6, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-Ca_i.tex"
plot(t, H_i*1.e10,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-H_i.tex"
plot(t, Cl_i,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-Cl_i.tex"
plot(t, cal,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-cal.tex"

# Plot the different background currents (t-I)
plot(t, I_Na_b, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_Na_b.tex"
plot(t, I_K_b,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_K_b.tex"
plot(t, I_Cl_b, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_Cl_b.tex"

# Plot the different background currents (V-I)
plot(V, I_Na_b, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_Na_b.tex"
plot(V, I_K_b,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_K_b.tex"
plot(V, I_Cl_b, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_Cl_b.tex"

# Plot the different pump and exchanger currents (t-I)
plot(t, I_NaK,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "bottom");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_NaK.tex"
plot(t, I_NaCa, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_NaCa.tex"
plot(t, I_NaH,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_NaH.tex"
plot(t, I_Ca_ATP,    'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_Ca_ATP.tex"

# Plot the different pump and exchanger currents (V-I)
plot(V, I_NaK,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "bottom");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_NaK.tex"
plot(V(5:end), I_NaCa(5:end), 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_NaCa.tex"
plot(V, I_NaH,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_NaH.tex"
plot(V, I_Ca_ATP,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_Ca_ATP.tex"

# Plot the other potassium currents (t-I)
plot(t, I_K_ur,     'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_K_ur.tex"
plot(t, I_K_2pore,  'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_K_2pore.tex"
plot(t, I_K_Ca_act, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_K_Ca_act.tex"
plot(t, I_K_ATP,    'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_K_ATP.tex"

# Plot the other potassium currents (V-I)
plot(V(10:end), I_K_ur(10:end),     'linewidth', line_width, 'color', blue);
hold on;
plot(V, I_K_ur_ref, '1', 'linewidth', line_width, 'color', red);
hold off;
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_K_ur.tex"
plot(I_K_2pore_ref(:, 1), I_K_2pore_ref(:, 2), '.', 'linewidth', line_width, 'color', red);
hold on;
plot(V, I_K_2pore,  'linewidth', line_width, 'color', blue);
hold off;
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_K_2pore.tex"
plot(I_K_Ca_act_ref(:, 1), I_K_Ca_act_ref(:, 2), '.', 'linewidth', line_width, 'color', red);
hold on;
plot(V, I_K_Ca_act, 'linewidth', line_width, 'color', blue);
hold off;
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_K_Ca_act.tex"
plot(V, I_K_ATP,    'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_K_ATP.tex"

# Plot the other currents (t-I)
plot(t, I_stim, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_stim.tex"
plot(t, I_ASIC, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_ASIC.tex"
plot(t, I_TRP1, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_TRP1.tex"
plot(t, I_TRP2, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-I_TRP2.tex"

# Plot the other currents (V-I)
plot(V, I_stim, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_stim.tex"
plot(V, I_ASIC, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_ASIC.tex"
plot(I_TRPv4_ref(:, 1), I_TRPv4_ref(:, 2), '.', 'linewidth', line_width, 'color', red);
hold on;
plot(V, I_TRP1, 'linewidth', line_width, 'color', blue);
hold off;
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_TRP1.tex"
plot(V, I_TRP2, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-I_TRP2.tex"

# Some special plots of interest
plot(K_o, V, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/K_o-V.tex"
plot(t, K_o, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-K_o.tex"
plot(t, vol_i, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/t-vol_i.tex"
plot(V, a_ur, 'linewidth', line_width, 'color', blue);
# axis ([-150 100 0 1])
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-a_dr.tex"
plot(V, tau_a_ur, 'linewidth', line_width, 'color', blue);
set (gca, "xaxislocation", "zero");
set (gca, "yaxislocation", "zero");
box off;
print -depslatexstandalone "../results/epslatex/V-tau_a_dr.tex"
