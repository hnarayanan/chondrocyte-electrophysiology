1;

# External stimulation
function I_stim = externalStimulation(t)
  global enable_I_stim;
  if (enable_I_stim == true)
    global t_cycle, global t_stim, global I_stim_bar;
    I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);
  else
    I_stim = 0.0;
  endif
endfunction

# FIXME: Implement the voltage-activated hydrogen channel
function I_ASIC = voltageActivatedHydrogen()
  global enable_I_ASIC;
  if (enable_I_ASIC == true)
    I_ASIC = 0.0;
  else
    I_ASIC = 0.0;
  endif
endfunction

# FIXME: Implement the stretch-activated trip channel
function I_TRP1 = stretchActivatedTrip()
  global enable_I_TRP1;
  if (enable_I_TRP1 == true)
    I_TRP1 = 0.0;
  else
    I_TRP1 = 0.0;
  endif
endfunction

# FIXME: Implement the osteo-arthritic trip channel
function I_TRP2 = osteoArthriticTrip()
  global enable_I_TRP2;
  if (enable_I_TRP2 == true)
    I_TRP2 = 0.0;
  else
    I_TRP2 = 0.0;
  endif
endfunction