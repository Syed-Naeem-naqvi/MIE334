function [g_new] = naqvi_1006980073_g_values(td, n)

% For i = 1
g_new = zeros(1,n-1)';
gi = td(1,2);
fi = td(1,1);
gi_prime = gi/fi;
g_new(1) = gi_prime;
gi_prev = gi_prime;


for i = 2:n-1
    g_new(i) = td(i,i+1)/(td(i,i)-td(i,i-1)*gi_prev);
    gi_prev = g_new(i);
end

