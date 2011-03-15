function I_Ca_act_K = current(V)

# Constants that correspond to the paper and make sense
e_ = 1.60217646e-19; # C
T = 306.15;          # K
k = 1.3806503e-23;   # m^2 kg s^-2 K^-1
L0 = 1.e-6;          # Slightly different from the paper
Z_L = 0.3*e_;        # C
Vh_J = 150;          # mV
Z_J = 0.58*e_;       # C
K_D = 1.1e-5;        # M
C = 8;
D = 25;
E = 2.4;

# Constants and relationships that need to be checked
Ca = 1.e-6;                 # Calcium concentration should vary
J0 = exp(Vh_J*Z_J/(k*T));   # Derived from the paper, check sign

# Relationships from the paper
J_  = J0*exp(-Z_J*V/(k*T));
K   = Ca/K_D;
L   = L0*exp(Z_L*V/(k*T))

P0_nr = L*(1 + K*C + J_*D + J_*K*C*D*E)^4;
P0  = P0_nr/(P0_nr + (1 + J_ + K + J_*K*E)^4); # Open probability

G_max = 1; # Maximum single channel conductance
# N_channel = 1.e4;

#V_Ca_act_K = 10;
#I_Ca_act_K = N_channel*P0*G_max*(V - V_Ca_act_K);

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