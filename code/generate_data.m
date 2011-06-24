parameters;
utility_functions;
background_currents;

states = [1, 2];
table = load('../data/reference_values/generated_small.data');
t = table(:, 1); len_t = size(t, 1);
data = table(:, 1 + states);

V = data(:, 1);
K_i = data(:, 2);
I_K_b = zeros(len_t, 1);

for ii = [1:len_t]
  I_K_b(ii) = backgroundPotassium(V(ii), K_i(ii), g_K_b_bar);
endfor

table = [t, I_K_b]
save('../data/reference_values/generated_I_K_b_small.data', 'table');