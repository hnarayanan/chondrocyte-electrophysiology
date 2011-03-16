1;

# Potential of an ion X across the membrane (mV).
function E_X = nernstPotential(z, X_i, X_o)
  global R, global T, global F;
  E_X = (R*T)/(z*F)*log(X_o/X_i);
endfunction