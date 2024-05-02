function [x, T] = naqvi_1006980073_project1(properties, dimensions, n)
% Transpose inputs
properties = properties';
dimensions = dimensions';

% Define variables
qb = properties(1); % W/m^2
Tinf = properties(2); % K
hb = properties(3); % W/m^2K
ht = properties(4); % W/m^2K
k = properties(5); % W/mK
L = dimensions(1); % m
b = dimensions(2); % m
% W = 1m

% Define helper stuff

A = zeros(n,n+1);           % Empty matrix to store equations
x = linspace(0,L,n);        % x locations
dx = L/(n-1);               % segment of fin length, L
h_x = (ht - hb).*x/L + hb;  % Heat transfer Coef. as a func. of x location

lambda = sqrt(L^2+(b/2)^2);  % Hypotenuse of the fin
dlambda = lambda/(n-1);      % Segment of the Hypotenuse of the fin
Ac_half = dlambda/2;         % Area at the 1st and nth (last) convection site
Ac_full = dlambda;           % Area of convection sites 2, 3, .. n-1
A_cond = b;                  % Areas at each of the conduction sites
dy = sqrt(dlambda^2 - dx^2); % Segment of fin height, b
b_middle = b-dy:-2*dy:0;     % heights of conduction areas
A_cond = [A_cond,b_middle];  % Conduction areas
dy = sqrt(dlambda^2 - dx^2); % Segment of fin height, b


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the First equation

f1 = (-k*A_cond(2)/dx - 2*h_x(1)*Ac_half);
g1 = (k*A_cond(2)/dx);
r1 =  -qb*A_cond(1)-2*h_x(1)*Ac_half*Tinf;
A(1,1) = f1;
A(1,2) = g1;
A(1, n+1) =  -qb*A_cond(1)-2*h_x(1)*Ac_half*Tinf;

% Add the Last equation

en = (k*A_cond(end)/dx);
fn = (-k*A_cond(end)/dx - 2*h_x(end)*Ac_half);
rn = -2*h_x(end)*Ac_half*Tinf;
A(n, end-2) = en;
A(n, end-1) = fn;
A(n, end) = rn;

% Create the system of equations

position_counter = 0;
for i = 2:n-1
    ei = (k*A_cond(i))/dx;
    fi = (-k*A_cond(i)/dx - k*A_cond(i+1)/dx - 2*h_x(i)*Ac_full);
    gi = (k*A_cond(i+1))/dx;
    ri = -2*h_x(i)*Ac_full*Tinf;
    A(i, 1+position_counter) = ei;
    A(i, 2+position_counter) = fi;
    A(i, 3+position_counter) = gi;
    A(i, end) = ri;
    position_counter = position_counter + 1;
end

% extract td and r
td = A(:,1:end-1);
r = A(:,end);

% Calculate g_new, r_new
g_new = naqvi_1006980073_g_values(td, n);
r_new = naqvi_1006980073_r_values(td, n, r, g_new);

% Update r and g
r = r_new;
for i = 1:n-1
    td(i,i+1) = g_new(i);
end

% Solve for temperature values
T = naqvi_1006980073_TDS_solver(td, r);
x = x';

% Graph the results
figure;
plot(x*100, T)
xlabel('x (cm)')
ylabel('Temperature (K)')
title('Naqvi')

