1;

# Clear the screen and plot solutions
clf;
hold on;

# The membrane voltage and total currents
figure(1);
subplot(2, 2, 1), plot(t, V), legend('V_{m} (mV)'); xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_i), legend('I_i (nA)'); xlabel('t (s)');
subplot(2, 2, 3), plot(V, I_i), legend('I_i (nA)'); xlabel('V_{m} (mV)');
print -depsc2 "output/membrane_behaviour.eps"

# The different concentrations
figure(2);
subplot(2, 2, 1), plot(t, Na_i), legend('[Na^{+}]_{i} (mM/l)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, K_i),  legend('[K^{+}]_{i} (mM/l)'), xlabel('t (s)');
subplot(2, 2, 3), plot(t, Ca_i),  legend('[Ca^{2+}]_{i} (mM/l)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/concentrations.eps"

# The different background currents
figure(3);
subplot(2, 2, 1), plot(t, I_Na_b), legend('I_{Na_{b}} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_K_b), legend('I_{K_{b}} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/background_currents.eps"

# The different pump and exchanger currents
figure(4);
subplot(2, 2, 1), plot(t, I_NaK), legend('I_{NaK} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_NaCa), legend('I_{NaCa} (nA)'), xlabel('t (s)');
subplot(2, 2, 3), plot(t, I_NaH), legend('I_{NaH} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/pumps_and_exchangers.eps"

# The other potassium currents
figure(5);
subplot(2, 2, 1), plot(t, I_K_ur), legend('I_{K_{ur}} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_K_2pore), legend('I_{K_{2pore}} (nA)'), xlabel('t (s)');
subplot(2, 2, 3), plot(t, I_K_Ca_act), legend('I_{K_{Ca-act}} (nA)'), xlabel('t (s)');
subplot(2, 2, 4), plot(t, I_K_ATP), legend('I_{K_{ATP}} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/potassium_currents.eps"

# The other currents
figure(6);
subplot(2, 2, 1), plot(t, I_stim), legend('I_{stim} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_ASIC), legend('I_{ASIC} (nA)'), xlabel('t (s)');
subplot(2, 2, 3), plot(t, I_TRP1), legend('I_{TRP1} (nA)'), xlabel('t (s)');
subplot(2, 2, 4), plot(t, I_TRP2), legend('I_{TRP2} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/other_currents.eps"

#print -depslatex "membrane_voltage.tex"