a = 100;
b = -1000;

V_m = [-100:1.0:100];
I = zeros(size(V_m));

for i = 1:size(V_m, 2)
    if(V_m(i) < 0)
      I(i) = 1e-6*(
		  b*V_m(i) + (1-b)*a*(1-(1-(V_m(i)/a))*(1-(V_m(i)/a))*(1-(V_m(i)/a))));
    else
      I(i) = 2.e-6*V_m(i)**3;
    endif
endfor

plot(V_m, I)
