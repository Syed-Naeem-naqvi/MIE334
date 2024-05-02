function [f, x_, r2] = naqvi_1006980073_fit_poly_curve(x, r, m, ni)

% Takes in the x and r values, the order and the number of points to
% interpolate and returns the corresponding polynomial, its x values (x_) and its coeffeficint
% of determination.

% Define an empty mxm matrix 

A = zeros(m+1,m+1);    % Define Empty mxm matrix
b = zeros(1,m+1)';     % Define Empty 1xm vector

s = size(x);
n = s(2);

% Fill out the matrix

for row = 1:m+1
    for col = 1:m+1
        A(row,col) = sum(x.^((row-1)+(col-1)));
    b(row) = sum(x.^(row-1).*r);
    end
end

A(1,1) = n;

% Solve by finding the inverse
I = inv(A);
coeffs = I*b;

% Evaluate the polynomial for the resulting coefficents at the linspaced
% points

xmax = max(x);
x_ = linspace(0,xmax,ni);
length_coeffs = length(coeffs);
f = [];
for i = 1:ni
    f_total = 0;
    for j = 1:length_coeffs
        f_total = f_total + coeffs(j)*x_(i).^(j-1);
    end
    f(i) = f_total;
end


% Error calculation
% Evaluate the function at the original points

f_at_x_pts = [];
for i = 1:n
    f_total = 0;
    for j = 1:length(coeffs)
        f_total = f_total + coeffs(j)*x(i).^(j-1);
    end
    f_at_x_pts(i) = f_total;
end

S_r = sum((r-f_at_x_pts).^2);
r_bar = sum(r)/n;
S_t = sum((r-r_bar).^2);

r2 = (S_t-S_r)/S_t;