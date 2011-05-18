# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

1;

# Potential of an ion X across the membrane (mV).
function E_X = nernstPotential(z, X_i, X_o)
  global R, global T, global F;
  E_X = (R*T)/(z*F)*log(X_o/X_i);
endfunction

# Applied voltage (mV)
function V = appliedVoltage(t)
  global clamp_Vm, global ramp_Vm, global step_Vm;
  if (clamp_Vm == true)
    global V_0;
    V = V_0;
  elseif (ramp_Vm == true)
    global V_0, global V_final, global t_final;
    V = V_0 + (V_final - V_0)*t/t_final;
  elseif (step_Vm == true)
    if (t <= 10)
      V = -60;
    elseif (t > 10 & t <= 20)
      V = -40;
    elseif (t > 20 & t <= 30)
      V = -20;
    elseif (t > 30 & t <= 40)
      V = -1;
    else
      V = -60;
    endif
  endif