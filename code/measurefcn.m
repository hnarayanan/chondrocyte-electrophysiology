function retval = measurefcn(x)
    global z_K, global K_o, global g_K_b_bar;
    V = x(1);
    K_i = x(2);
    E_K = nernstPotential(z_K, K_i, K_o);
    retval = g_K_b_bar*(V - E_K);
endfunction