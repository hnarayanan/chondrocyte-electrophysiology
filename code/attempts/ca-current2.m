function I_Ca_act_K = current(V)

# Straight from the Igor code
Zj=0.58; Vhj=150; ZL=0.3; L0=1e-6; KDc=11e-6; C=8; D=25; E=2.4; Gmax=1;
W = [Zj, Vhj, ZL, L0, KDc, C, D, E, Gmax];
Ca = 11.e-6;
T = 20;
kTe = 23.54*((T + 273)/273);
Lv = W(4)*exp((V*W(3))/kTe);
Jv = exp(((V - W(2))*W(1))/kTe);
K = Ca/W(5);
C = W(6);
D = W(7);
E = W(8);
Gmax=W(9);

P0=(Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
I_Ca_act_K = P0;

endfunction

currents = [];
V_range = [-100:0.1:300];
for V = V_range;
  currents = [currents current(V)];
endfor

plot(V_range, currents);
xlabel("Voltage (mV)");
ylabel("Channel open probability");

print -depsc2 "output.eps"