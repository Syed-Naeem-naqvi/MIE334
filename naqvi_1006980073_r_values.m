function [r_new] = naqvi_1006980073_r_values(td, n, r, g_new)
r_new = zeros(n,1);
r_new(1) = r(1)/td(1,1);
r_prev = r_new(1);
for i = 2:n
    r_new(i) = (r(i) - td(i,i-1)*r_prev)/(td(i,i) - td(i,i-1)*g_new(i-1));
    r_prev = r_new(i);
end
