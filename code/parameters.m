# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# The solution fields and parameters have the following units:
#
# Voltage: mV
# Current: pa
# Time: s
# Concentration: mM/l
# Conductance: pS
# Capacitance: pF

# Toggle estimation of parameters
global enable_parest = false;

# Toggle clamping of the internal concentrations
global clamp_conc = false;

# Toggle setting of the membrane voltage
global apply_Vm = true;

# If apply_Vm is true, then one of clamp_Vm, ramp_Vm and step_Vm must be
# true to define what voltage is to be applied
global clamp_Vm = false;
global step_Vm = false;
global ramp_Vm = true;
global V_final = 90.0;  # Final value of membrane voltage when ramped (mV)

# Toggle individual currents
global enable_I_Na_b = false;
global enable_I_K_b = true;
global enable_I_NaK = false;
global enable_I_NaCa = false;
global enable_I_NaH = false;
global enable_I_K_ur = true;
global enable_I_K_2pore = true;
global enable_I_K_Ca_act = true;
global enable_I_K_ATP = false;
global enable_I_ASIC = false;
global enable_I_TRP1 = false;
global enable_I_TRP2 = false;
global enable_I_stim = false;

# Time-stepping information
global t_final = 0.4;       # Final time (s)
global dt = t_final/1000;   # Time increment (s)

# External concentrations
global Na_o = 320;       # Clamped external sodium concentration (mM/l)
global K_o  = 10;        # Clamped external potassium concentration (mM/l)
global Ca_o = 16;        # Clamped external calcium concentration (mM/l)
global H_o  = 0.125893;  # Clamped external hydrogen concentration (mM/l)

# Initial conditions
global V_0    = -130.0;    # Initial membrane potential (mV)
global Na_i_0 = 0.516766;  # Initial internal sodium concentration (mM/l)
global K_i_0  = 110.00;    # Initial internal potassium concentration (mM/l)
global Ca_i_0 = 0.001;     # Initial internal calcium concentration (mM/l)
global H_i_0  = 0.5;       # Initial internal hydrogen concentration (mM/l)
global a_ur_0 = 0.000367;
global i_ur_0 = 0.967290;

# Universal constants
global R = 8314.472; # Universal gas constant (mJ K^-1 mol^-1)
global T = 310.15;   # Normal body temperature (K)
global F = 96485.34; # Faraday's constant (C mol^-1)

# Charges on each of the ions
global z_Na = 1;         # Charge on the sodium ion
global z_K  = 1;         # Charge on the potassium ion
global z_Ca = 2;         # Charge on the calcium ion
global z_H  = 1;         # Charge on the calcium ion

# Cell parameters
global C_m = 15.0;       # Membrane capacitance
global vol_i = 0.005884; # Internal volume

# Constants related to external stimulation
global t_cycle = 1.0;    # Total cycle time (s)
global t_stim = 0.1;     # Stimulation time/cycle (s)
global I_stim_bar = 0.0; # Stimulation current magnitude ()

# Background conductances
global g_Na_b_bar = 0.2;   # Background sodium leakage conductance (pS)
global g_K_b_bar = 0.2;    # Background potassium leakage conductance (pS)

# Constants related to the sodium-potassium pump
global I_NaK_bar = 68.55;
global K_NaK_K = 1.0;
global K_NaK_Na = 11.0;

# Constants related to the sodium-calcium exchanger
global K_NaCa = 0.001;
global gamma_Na = 0.45;
global d_NaCa = 0.0003;

# Constants related to the sodium-hydrogen exchanger
global n_H = 1;
global m_H = 3;
global K_H_i_mod = 3.07e-5;
global K_H_o_mod = 4.8e-7;
global k1_p = 10.5;
global k1_m = 0.201;
global k2_p = 15.8;
global k2_m = 183;
global K_Na_i = 16.2
global K_Na_o = 195
global K_H_i = 6.05e-4;
global K_H_o = 1.62e-3;
global N_NaH_channel = 4899;

# Constants related to the ultra-rapidly rectifying potassium channel
global g_K_ur = 1.0;

# Constants related to the two-pore potassium channel
global P_K = 3.2e-6;

# Constants related to the calcium-activated potassium channel
global Zj = 1.10;
global Vhj = 180;
global ZL = 0.3;
global L0 = 9e-6;
global KDc = 3e-6;
global C = 8;
global D = 25;
global E = 2.4;
global Gmax = 2.1;
global N_channel = 1.0;
