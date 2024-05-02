function[T] = naqvi_1006980073_TDS_solver(td, r)
s = size(td);
n = s(1);

T = zeros(n,1);
T(n) = r(n);

for i = n-1:-1:1
    T(i) = r(i) - td(i,i+1)*T(i+1);
end
