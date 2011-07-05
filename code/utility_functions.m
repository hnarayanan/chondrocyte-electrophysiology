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
    global t_cycle t_stim;
    V = (ceil((t - 30)/t_cycle).*square((t - 30)*2*pi/t_cycle, t_stim/t_cycle) + ceil((t - 30)/t_cycle))/2*10 - 90;
    if(V == 0) V = 0.01; endif
  endif