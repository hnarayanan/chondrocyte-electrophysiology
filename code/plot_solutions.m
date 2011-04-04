1;

# Clear the screen and plot solutions
clf;
line_width = 3;

# Load reference solution
V_ref = csvread('references/Total_Current_Density_vs_Membrane_Voltage.data')(:, 1);
I_i_by_C_m_ref = csvread('references/Total_Current_Density_vs_Membrane_Voltage.data')(:, 2);

# Hide plots
figure(1, 'visible', 'off')

# Plot membrane voltage vs time
plot(t, V, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$V_{m} (mV)$');
print -depslatexstandalone "output/t_V.tex"

# Plot total ionic current vs time
plot(t, I_i, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_i (pA)$');
print -depslatexstandalone "output/t_I_i.tex"

# Plot total ionic current density vs time (including reference values)
plot(V, I_i/C_m, 'linewidth', line_width), xlabel('$V_{m} (mV)$'), ylabel('$I_i/C_{m} (pA/pF)$');
hold on;
plot(V_ref, I_i_by_C_m_ref, '1', 'linewidth', line_width);
hold off;
print -depslatexstandalone "output/V_I_i.tex"

# Plot interior sodium concentration vs time
plot(t, Na_i, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$[Na^{+}]_{i} (mM/l)$');
print -depslatexstandalone "output/t_Na_i.tex"

# Plot interior potassium concentration vs time
plot(t, K_i, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$[K^{+}]_{i} (mM/l)$');
print -depslatexstandalone "output/t_K_i.tex"

# Plot interior calcium concentration vs time
plot(t, Ca_i, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$[Ca^{2+}]_{i} (mM/l)$');
print -depslatexstandalone "output/t_Ca_i.tex"

# Plot background sodium current vs time
plot(t, I_Na_b, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{Na_{b}} (pA)$');
print -depslatexstandalone "output/t_I_Na_b.tex"

# Plot background potassium current vs time
plot(t, I_K_b, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{K_{b}} (pA)$');
print -depslatexstandalone "output/t_I_K_b.tex"

# Plot sodium-potassium pump current vs time
plot(t, I_NaK, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{NaK} (pA)$');
print -depslatexstandalone "output/t_I_NaK.tex"

# Plot sodium-calcium exchanger current vs time
plot(t, I_NaCa, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{NaCa} (pA)$');
print -depslatexstandalone "output/t_I_NaCa.tex"

# Plot sodium-hydrogen antiport current vs time
plot(t, I_NaH, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{NaH} (pA)$');
print -depslatexstandalone "output/t_I_NaH.tex"

# Plot ultra-rapidly rectifying potassium current vs time
plot(t, I_K_ur, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{K_{ur}} (pA)$');
print -depslatexstandalone "output/t_I_K_ur.tex"

# Plot two-pore potassium current vs time
plot(t, I_K_2pore, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{K_{2pore}} (pA)$');
print -depslatexstandalone "output/t_I_K_2pore.tex"

# Plot calcium-activated potassium current vs time
plot(t, I_K_Ca_act, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{K_{Ca-act}} (pA)$');
print -depslatexstandalone "output/t_I_K_Ca_act.tex"

# Plot ATP-powered potassium pump current vs time
plot(t, I_K_ATP, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{K_{ATP}} (pA)$');
print -depslatexstandalone "output/t_I_K_ATP.tex"

# Plot external stimulation current vs time
plot(t, I_stim, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{stim} (pA)$');
print -depslatexstandalone "output/t_I_stim.tex"

# Plot voltage-activated hydrogen current vs time
plot(t, I_ASIC, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{ASIC} (pA)$');
print -depslatexstandalone "output/t_I_ASIC.tex"

# Plot stretch-activated trip current vs time
plot(t, I_TRP1, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{TRP1} (pA)$');
print -depslatexstandalone "output/t_I_TRP1.tex"

# Plot osteo-arthritic trip current vs time
plot(t, I_TRP2, 'linewidth', line_width), xlabel('$t (s)$'), ylabel('$I_{TRP2} (pA)$');
print -depslatexstandalone "output/t_I_TRP2.tex"