# Universal constants
global R = 8314.472; # Universal gas constant (mJ K^-1 mol^-1)
global T = 310.15;   # Normal body temperature (K)
global F = 96485.34; # Faraday's constant (C mol^-1)

# Cell parameters
global C_m = 15.0;       # Membrane capacitance
global vol_i = 0.005884; # Internal volume

# Time-stepping information
t_final = 10.0           # Final time(s)
dt = 0.05                # Time increment(s)

# Initial conditions
global V0 = -62.3;          # Initial membrane potential (mV)
global Na_i_0 = 0.516766;   # Initial internal sodium concentration (mM/l)
global K_i_0 = 129.485991;  # Initial internal potassium concentration (mM/l)
global Ca_i_0 = 6.5e-5;     # Initial internal calcium concentration (mM/l)
global a_ur_0 = 0.000367;
global I_ur_0 = 0.967290;

# Constants related to external stimulation
global t_cycle = 1.0;    # Total cycle time (s)
global t_stim = 0.1;     # Stimulation time/cycle (s)
global I_stim_bar = 0.0; # Stimulation current magnitude ()

# Constants related to sodium
global Na_o = 130.022096;   # Clamped external sodium concentration (mM/l)
global z_Na = 1;            # Charge on the sodium ion

# Constants related to potassium
global K_o = 5.560224;      # Clamped external potassium concentration (mM/l)
global z_K = 1;             # Charge on the potassium ion

# Constants related to calcium
global Ca_o = 1.815768;     # Clamped external calcium concentration (mM/l)
global z_Ca = 2;            # Charge on the calcium ion

# Background conductances
global g_Na_b_bar = 20;    # Background sodium leakage conductance (nS)
global g_K_b_bar = 20;     # Background potassium leakage conductance (nS)

# Constants related to the sodium-potassium pump
global I_NaK_bar = 68.55;
global K_NaK_K = 1.0;
global K_NaK_Na = 11.0;

# Constants related to the sodium-calcium exchanger
global K_NaCa = 0.0374842;
global gamma_Na = 0.45;
global d_NaCa = 0.0003;

# Constants related to the ultra-rapidly rectifying potassium channel
global g_K_ur = 2.25;

# Constants related to the two-pore potassium channel
# FIXME: Find out the value of the conductance
global g_K_2pore = 0.0;

# Constants related to the calcium-activated potassium channel
global Zj = 0.58;
global Vhj = 150;
global ZL = 0.3;
global L0 = 1e-6;
global KDc = 11e-6;
global C = 8;
global D = 25;
global E = 2.4;
global Gmax = 1;
global N_channel = 1.0;
global V_K_Ca_act = 10.0;
