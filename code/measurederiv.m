# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2011  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the derivative of the measured function(s) with respect to the
# different parameters
function retval = measurederiv(x)
    global z_K, global K_o;
    global R, global T, global F, global z_K;
    global g_K_b_bar, global P_K;

    V = x(1);
    Na_i = x(2);
    K_i  = x(3);
    Ca_i = x(4);
    H_i  = x(5);
    a_ur = x(6);
    I_ur = x(7);

    dstatedx = zeros(size(x'));

    global enable_I_K_b;
    if (enable_I_K_b == true)
      dstatedx(1) += g_K_b_bar;
      dstatedx(3) += g_K_b_bar*R*T/(z_K*F)/K_i;
    endif

    global enable_I_K_2pore;
    if (enable_I_K_2pore == true)
      dstatedx(1) += -exp(-F*V*T^(-1)*z_K*R^(-1))*P_K*F^3*V*T^(-2)*z_K^3*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-2)*K_o+exp(-F*V*T^(-1)*z_K*R^(-1))*P_K*F^3*V*T^(-2)*(exp(-F*V*T^(-1)*z_K*R^(-1))*K_o-K_i)*z_K^3*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-2)*R^(-2)+P_K*F^2*T^(-1)*(exp(-F*V*T^(-1)*z_K*R^(-1))*K_o-K_i)*z_K^2*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-1);
      dstatedx(3) += -P_K*F^2*V*T^(-1)*z_K^2*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-1);
    endif

    retval = dstatedx;
endfunction