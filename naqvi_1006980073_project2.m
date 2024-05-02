%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIE334 Project 2, Naeem Naqvi, 1006980073, April
function [CA, V, r2] = naqvi_1006980073_project2(x, r, m, ni)
% Takes in the experimental x-values and r-values, the order of the
% polynomial m and the number of interpolated polynomial points ni and
% returns the contact angle CA, volume V and and the coefficinet of
% determination r2
%%% NOTE: The code associated with the ORIGINAL DATA POINTS has been commented out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINDING THE POLYNOMIAL AND GRAPHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[func, func_x_vals, coeff_of_det] = naqvi_1006980073_fit_poly_curve(x, r, m, ni);

f = func;
x_ = func_x_vals;

% Graph the results: Plot both the original data points and the fitted
% polynomial function

figure;
grid on;
hold on;
scatter(r,x, 'filled');
plot(f,x_)
hold off;
title('x versus r = f(x)');
ylabel('X [mm]');
xlabel('r [mm]');
pbaspect([1 1 1])
txt = ['R^2 = ', num2str(coeff_of_det)];
text(0.5,0.5,txt)
legend('Data Points', 'Fitted Function')

r2 = coeff_of_det; % Return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTACT ANGLE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A. Using fitted polynmial

x1 = x_(1);
x2 = x_(2);
r1 = f(1);
r2_ = f(2);
del_x = x2 - x1;
del_r = r2_ - r1;
del_h = sqrt(del_x.^2 + del_r.^2);
theta = acos(del_x/del_h);
theta_deg = rad2deg(theta);
theta_total = theta_deg + 90;

CA = theta_total; % Return

% B. Using the Original Points of data
% 
% x1p = x(1);
% x2p = x(2);
% r1p = r(1);
% r2p = r(2);
% del_xp = x2p - x1p;
% del_rp = r2p - r1p;
% del_hp = sqrt(del_xp.^2 + del_rp.^2);
% thetap = acos(del_xp/del_hp);
% thetap_deg = rad2deg(thetap);
% thetap_total = thetap_deg + 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING VOLUME USING TRAPEZOIDAL RULE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. Using Polynomial

h = (x_(end)-x_(1))/ni;
sum_vol = f(1).^2;

for i = 1:ni-1
    sum_vol = sum_vol + 2*(f(i)).^2;
end
sum_vol = sum_vol + f(end).^2;
I_vol_using_poly = pi*h/2*sum_vol


% B. Using original data points
% 
% I = [];
% for i = 1:n-1
%     h = x(i+1)-x(i);
%     term = r(i+1).^2 + r(i).^2;
%     I(i) = h*term/2;
% end
% I_vol_using_pts = pi*sum(I)

V = I_vol_using_poly; % Return