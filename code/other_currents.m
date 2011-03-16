1;

# External stimulation
function I_stim = externalStimulation(t)
  global t_cycle, global t_stim, global I_stim_bar;
  I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);
endfunction

# FIXME: Implement the voltage-activated hydrogen channel
function I_ASIC = voltageActivatedHydrogen()
  I_ASIC = 0.0;
endfunction