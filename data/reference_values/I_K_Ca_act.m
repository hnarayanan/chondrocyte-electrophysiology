Zj = 0.70;
Vhj = 250;
ZL = 0.1;
L0 = 12e-6;
KDc = 3e-6;
C = 8;
D = 25;
E = 2.4;
Gmax = 3.8*1.6;
N_channel = 1.0;
E_K_Ca_act = 42;
T = 310.15;

Ca_i_range = [5.e-8, 1.e-7, 5e-7, 1e-6, 5e-6]

V = [-150:1.0:90];

axis([0 100 0 40]);
hold

for j = 1:5
    Ca_i = Ca_i_range(j)
  I = zeros(size(V));

  for i = 1:size(V, 2)
    kTe = 23.54*(T/273);
    Lv = L0*exp((V(i)*ZL)/kTe);
    Jv = exp(((V(i) - Vhj)*Zj)/kTe);
    K = Ca_i/KDc;
    P0 = (Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
    I(i) = N_channel*P0*Gmax*(V(i) - E_K_Ca_act);
  endfor

  set (0,'defaultaxesfontsize', 12)
  plot(V, I)
endfor


print -depslatexstandalone "temporary.tex"
