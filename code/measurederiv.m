function retval = measurederiv(x)
    global z_K, global K_o, global g_K_b_bar;
    global R, global T, global F, global z_K;
    K_i = x(2);
    dstatedx = zeros(size(x'));
    dstatedx(1) = g_K_b_bar;
    dstatedx(2) = g_K_b_bar*R*T/(z_K*F)/K_i;
    retval = dstatedx;
endfunction