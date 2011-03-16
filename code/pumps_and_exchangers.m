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