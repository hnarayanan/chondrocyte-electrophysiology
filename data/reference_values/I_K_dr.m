G_K = 28.9;  # pS/pF
E_K = -83;   # mV
V_h = -26.7; # mV
S_h = 4.1;   # mV
C_m = 8.5;

V_m = [-100:1.0:75];
I = zeros(size(V_m));

for i = 1:size(V_m, 2)
    I(i) = G_K*(V_m(i) - E_K) / (1 + exp(-(V_m(i) - V_h)/S_h)) * C_m/1000.0;
endfor

plot(V_m, I)
