# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the measured function(s)
function retval = measurefcn(x)
    global g_K_b_bar, global P_K, global Gmax;
    V = x(1);
    Na_i = x(2);
    K_i  = x(3);
    Ca_i = x(4);
    H_i  = x(5);
    a_ur = x(6);
    I_ur = x(7);
    retval = backgroundPotassium(V, K_i, g_K_b_bar) \
	   + twoPorePotassium(V, K_i, P_K) \
	   + calciumActivatedPotassium(V, K_i, Ca_i, Gmax);
endfunction