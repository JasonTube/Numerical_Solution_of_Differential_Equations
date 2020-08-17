% Lax-Wendroff scheme for the one-way wave equation:
% u_t + a*u_x = 0, u(x,0) = ufunc(x,0)
clear;clc;close all
uexact = @(x,t)(x-3*t).^2; % exact solution for the Cauchy problem 

a = 3;
x_final = 1; t_final = 1; % solve the PDE at x = 1, t = 1

h = 0.05;
tau = 0.01;

k = t_final/tau; % how many iterations are needed.

% dependent interval
x = x_final(1)-k*h : h : x_final(end)+k*h;

miu = a*tau/h;

% initial data of dependent interval
U = u0(x);
n = length(U);

for t = tau:tau:t_final
    % Lax-Wendroff Update
    U = U(2:n-1) + (miu^2 * (U(1:n-2)-2*U(2:n-1)+U(3:n))...
        - miu * (U(3:n)-U(1:n-2)))/2;
    n = length(U);
end

disp(['Computed Solution at x = 1, t = 1: ',num2str(U)])
disp(['Exact Solution at x = 1, t = 1   : ',num2str(uexact(1,1))])
disp(['h = ',num2str(h),' tau = ',num2str(tau),' Error: ',num2str(U-uexact(1,1))])