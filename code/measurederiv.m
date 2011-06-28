# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the derivative of the measured function(s) with respect to the
# different parameters
function retval = measurederiv(x)
    global z_K, global K_o, global g_K_b_bar;
    global R, global T, global F, global z_K;
    K_i = x(2);
    dstatedx = zeros(size(x'));
    dstatedx(1) = g_K_b_bar;
    dstatedx(2) = g_K_b_bar*R*T/(z_K*F)/K_i;
    retval = dstatedx;
endfunction