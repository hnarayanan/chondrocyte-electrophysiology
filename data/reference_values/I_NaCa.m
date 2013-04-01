F =  96485.34;
R = 8314.472;
T = 310.15;

Na_i  =    1.22582880260390e+00;
Na_o = 140;
Ca_i  = 5.e-8;#   1.20992489429946e-07;
Ca_o = 2;
K_NaCa = 0.02*8;
gamma_Na = 0.45;
d_NaCa = 0.0003;
V_m = [-150:1.0:90];
I = zeros(size(V_m));

for i = 1:size(V_m, 2)
    I(i) = K_NaCa*(  Na_i^3*Ca_o*exp(gamma_Na*V_m(i)*F/(R*T)) \
                     - Na_o^3*Ca_i*exp((gamma_Na - 1.0)*V_m(i)*F/(R*T))) \
             /(1.0 + d_NaCa*(Na_o^3*Ca_i + Na_i^3*Ca_o));
endfor

set (0,'defaultaxesfontsize', 12)
plot(V_m, I)
