V = [-100:1.0:60];
I = zeros(size(V));

for i = 1:size(V, 2)
  sigma   = 0.6;
  g_0     = 4;
  p_0     = 0.91;
  E_K_ATP = -94.02;
  H_K_ATP = -0.001;
  K_m_ATP = 0.56;
  surf    = 1;

  ATP_i = V(i) + 101;
  ADP_i = 10;

  H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i);
  K_m = 35.8 + 17.9*ADP_i^K_m_ATP;
  f_ATP = 1.0/(1.0 + (ATP_i/K_m)^H);

  I(i) = sigma*surf*g_0*p_0*f_ATP*(V(i) - E_K_ATP);
endfor

plot(V, I)
