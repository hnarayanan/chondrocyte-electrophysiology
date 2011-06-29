# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the measured function(s)
function retval = measurefcn(x)
    global g_K_b_bar;
    V = x(1);
    K_i = x(3);
    retval = backgroundPotassium(V, K_i, g_K_b_bar);
endfunction