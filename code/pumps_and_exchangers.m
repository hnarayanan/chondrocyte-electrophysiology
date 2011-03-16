1;

# Sodium-potassium pump from "Mathematical Model of an Adult Human
# Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
# Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
# Giles. Circ. Res. 1998; 82; 63-81 (Table 12, pp. 77)

function I_NaK = sodiumPotassiumPump(V, Na_i, K_i)
  global I_NaK_bar, global K_o;
  global K_NaK_K, global K_NaK_Na;
  I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) \
                   *(Na_i^1.5/(Na_i^1.5 + K_NaK_Na^1.5)) \
                   *(V + 150.0)/(V + 200.0);
endfunction

# Sodium-calcium exchanger from "Mathematical Model of an Adult Human
# Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
# Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
# Giles. Circ. Res. 1998; 82; 63-81 (Table 13, pp. 77)

function I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i)
  global F, global R, global T;
  global Na_o, global Ca_o;
  global K_NaCa, global gamma_Na, global d_NaCa;

  I_NaCa = K_NaCa*(  Na_i^3*Ca_o*exp(gamma_Na*V*F/(R*T)) \
                   - Na_o^3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
      / (1.0 + d_NaCa*(Na_o^3*Ca_i + Na_i^3*Ca_o));
endfunction

# FIXME: Implement the sodium-hydrogen antiport

function I_NaH = sodiumHydrogenAntiport()
  I_NaH = 0.0;
endfunction

# FIXME: Implement the voltage-activated hydrogen channel

function I_ASIC = voltageActivatedHydrogen()
  I_ASIC = 0.0;
endfunction