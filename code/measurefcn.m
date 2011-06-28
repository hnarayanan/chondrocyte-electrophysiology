function retval = measurefcn(x)
    global g_K_b_bar;
    V = x(1);
    K_i = x(2);
    retval = backgroundPotassium(V, K_i, g_K_b_bar);
endfunction