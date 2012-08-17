# This file is part of the chondrocyte modelling project at Simula
# Research Laboratory, Norway. Refer to the files README and COPYING for
# more information about the project as well as terms of distribution.
#
# Copyright (C) 2010--2012  Harish Narayanan
# Licensed under the GNU GPL Version 3

# Define the derivative of the measured function(s) with respect to the
# different parameters
function retval = measurederiv(x)
    global z_K, global K_o;
    global R, global T, global F, global z_K;

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
      global g_K_b_bar;
      dstatedx(1) += g_K_b_bar;
      dstatedx(3) += g_K_b_bar*R*T/(z_K*F)/K_i;
    endif

    global enable_I_K_2pore;
    if (enable_I_K_2pore == true)
      global P_K;
      dstatedx(1) += -exp(-F*V*T^(-1)*z_K*R^(-1))*P_K*F^3*V*T^(-2)*z_K^3*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-2)*K_o+exp(-F*V*T^(-1)*z_K*R^(-1))*P_K*F^3*V*T^(-2)*(exp(-F*V*T^(-1)*z_K*R^(-1))*K_o-K_i)*z_K^3*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-2)*R^(-2)+P_K*F^2*T^(-1)*(exp(-F*V*T^(-1)*z_K*R^(-1))*K_o-K_i)*z_K^2*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-1);
      dstatedx(3) += -P_K*F^2*V*T^(-1)*z_K^2*(-1+exp(-F*V*T^(-1)*z_K*R^(-1)))^(-1)*R^(-1);
    endif

  global enable_I_K_Ca_act;
  if (enable_I_K_Ca_act == true)
    global N_channel, global Gmax;
    global Zj, global Vhj, global ZL, global L0, global KDc;
    global C, global D, global E;
    dstatedx(1) += ((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-1)*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel+(11.597281223449447748)*((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-1)*T^(-1)*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*(V-R*T*F^(-1)*log(K_i^(-1)*K_o)*z_K^(-1))*ZL*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel-((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-2)*(4*((11.597281223449447748)*Ca_i*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*Zj*KDc^(-1)+(11.597281223449447748)*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*Zj)*(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^3+(11.597281223449447748)*T^(-1)*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*ZL*exp((11.597281223449447748)*V*T^(-1)*ZL)+4*((11.597281223449447748)*Ca_i*C*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*Zj*KDc^(-1)*D+(11.597281223449447748)*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*Zj*D)*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^3*L0*exp((11.597281223449447748)*V*T^(-1)*ZL))*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*(V-R*T*F^(-1)*log(K_i^(-1)*K_o)*z_K^(-1))*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel+4*((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-1)*((11.597281223449447748)*Ca_i*C*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*Zj*KDc^(-1)*D+(11.597281223449447748)*T^(-1)*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*Zj*D)*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^3*L0*(V-R*T*F^(-1)*log(K_i^(-1)*K_o)*z_K^(-1))*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel;
    dstatedx(3) += R*((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-1)*T*F^(-1)*K_i^(-1)*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*z_K^(-1)*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel;
    dstatedx(4) += -4*((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-2)*((exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)+KDc^(-1))*(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^3+(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^3*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)*(C*KDc^(-1)+C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D))*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*(V-R*T*F^(-1)*log(K_i^(-1)*K_o)*z_K^(-1))*exp((11.597281223449447748)*V*T^(-1)*ZL)*N_channel+4*((1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^4*L0*exp((11.597281223449447748)*V*T^(-1)*ZL)+(1+Ca_i*KDc^(-1)+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)+Ca_i*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1))^4)^(-1)*Gmax*(1+exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*D+Ca_i*C*KDc^(-1)+Ca_i*C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)^3*L0*(V-R*T*F^(-1)*log(K_i^(-1)*K_o)*z_K^(-1))*exp((11.597281223449447748)*V*T^(-1)*ZL)*(C*KDc^(-1)+C*exp((11.597281223449447748)*T^(-1)*(V-Vhj)*Zj)*E*KDc^(-1)*D)*N_channel;
  endif

  retval = dstatedx;
endfunction
