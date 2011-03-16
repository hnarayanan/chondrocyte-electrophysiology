1;

# Clear the screen and plot solutions
clf;
hold on;

# Just the membrane voltage
figure(1);
plot(t, V); xlabel('t (s)'); ylabel('V_{m} (mV)');
print -depsc2 "output/voltage.eps"

# The different current components
figure(2);
subplot(2, 2, 1), plot(t, I_Na_b), legend('I_{Na_{b}} (nA)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, I_K_b), legend('I_{K_{b}} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/background_currents.eps"

# The different current components
figure(3);
subplot(2, 2, 1), plot(t, I_NaK), legend('I_{NaK} (nA)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/pumps_and_exchangers.eps"

# The different concentrations
figure(4);
subplot(2, 2, 1), plot(t, Na_i), legend('[Na^{+}]_{i} (mM/l)'), xlabel('t (s)');
subplot(2, 2, 2), plot(t, K_i),  legend('[K^{+}]_{i} (mM/l)'), xlabel('t (s)');
xlabel('t (s)');
print -depsc2 "output/concentrations.eps"

#print -depslatex "membrane_voltage.tex"