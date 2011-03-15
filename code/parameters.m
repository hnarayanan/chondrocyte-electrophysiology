# Universal constants
global R = 8.314472; # Universal gas constant (J K^-1 mol^-1)
global T = 310.15;   # Normal body temperature (K)
global F = 96485.34; # Faraday's constant (C mol^-1)

# Cell parameters
global C_m = 15.0;       # Membrane capacitance
global vol_i = 0.005884; # Internal volume

# Constants related to external stimulation
global t_cycle = 1.0;    # Total cycle time (s)
global t_stim = 0.1;     # Stimulation time/cycle (s)
global I_stim_bar = 0.0; # Stimulation current magnitude ()

# Constants related to sodium
global Na_i_0 = 0.516766;   # Initial internal sodium concentration (mM/l)
global Na_o   = 130.022096; # Clamped external sodium concentration (mM/l)
global z_Na   = 1;          # Charge on the sodium ion
global g_Na_b_bar = 20;     # Background sodium leakage conductance (nS)

# Constants related to potassium
global K_i_0 = 129.485991; # Initial internal potassium concentration (mM/l)
global K_o   = 5.560224;   # Clamped external potassium concentration (mM/l)
global z_K   = 1;          # Charge on the potassium ion
global g_K_b_bar = 20;     # Background potassium leakage conductance (nS)