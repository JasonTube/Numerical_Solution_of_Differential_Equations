% Main function for Solving non-linear ODE
% theta'' + K*sin(theta) = 0
% theta(0) = theta1, theta(2*pi) = theta2
clear;clc;
%% input
K = 20;
theta = [0,-0.05];
%% solve
n = 40;% grid
h = 2*pi/n;
t = h:h:(2*pi-h);t=t';
x_initial = -0.1*cos(sqrt(K).*t+pi/2);% initial guess
% call the Newton solver
x_nonlin_numerical=...
    Newton(x_initial,@f_nonlin,@df_nonlin,K,theta,n,1.e-6);
x_lin_numerical=...
    Newton(x_initial,@f_lin,@df_lin,K,theta,n,1.e-6);
%% the exact solution of linear problem: theta''+K*theta=0
lin_exact = dsolve('D2u+K*u=0','u(0)=a','u(2*pi)=b');
lin_exact = matlabFunction(lin_exact);
x_lin_exact = lin_exact(K,theta(1),theta(2),t);
%% plot
figure();
plot([0;t;2*pi],[theta(1);x_nonlin_numerical;theta(2)],'-*','lineWidth',1.2);hold on
plot([0;t;2*pi],[theta(1);x_lin_numerical;theta(2)],'-+','lineWidth',1.2);hold on
plot([0;t;2*pi],[theta(1);x_lin_exact;theta(2)],'-','lineWidth',1.2);
legend('Numerical for $\frac{{{d^2}\theta }}{{d{t^2}}} + K\sin \theta  = 0$',...
    'Numerical for $\frac{{{d^2}\theta }}{{d{t^2}}} + K\theta  = 0$',...
    'Exact for $\frac{{{d^2}\theta }}{{d{t^2}}} + K\theta  = 0$',...
    'Interpreter','latex','FontSize',12)
grid minor;xlabel('t');ylabel('\theta')
title(['Solution plots with n = ',num2str(n)],'Interpreter','latex')
%% functions for Discretizations and their Jacobian Matrix
function f = f_nonlin(x,K,theta,n)
f = zeros(n-1,1);
h = 2*pi/n; h2 = h^2;
f(2:n-2) = (x(1:n-3)-2*x(2:n-2)+x(3:n-1))/h2 + K*sin(x(2:n-2));
f(1) = (theta(1)-2*x(1)+x(2))/h2 + K*sin(x(1));
f(n-1) = (x(n-2)-2*x(n-1)+theta(2))/h2 + K*sin(x(n-1));
end

function df = df_nonlin(x,K,n)
df = zeros(n-1,n-1);
h = 2*pi/n; h2 = h^2;
df = df + diag(-2/h2 + K*cos(x));
df(1:n-2,2:n-1) = df(1:n-2,2:n-1) + eye(n-2)/h2;
df(2:n-1,1:n-2) = df(2:n-1,1:n-2) + eye(n-2)/h2;
end

function f = f_lin(x,K,theta,n)
f = zeros(n-1,1);
h = 2*pi/n; h2 = h^2;
f(2:n-2) = (x(1:n-3)-2*x(2:n-2)+x(3:n-1))/h2 + K*x(2:n-2);
f(1) = (theta(1)-2*x(1)+x(2))/h2 + K*x(1);
f(n-1) = (x(n-2)-2*x(n-1)+theta(2))/h2 + K*x(n-1);
end

function df = df_lin(x,K,n)
df = zeros(n-1,n-1);
h = 2*pi/n; h2 = h^2;
df = df + (-2/h2 + K)*eye(n-1);
df(1:n-2,2:n-1) = df(1:n-2,2:n-1) + eye(n-2)/h2;
df(2:n-1,1:n-2) = df(2:n-1,1:n-2) + eye(n-2)/h2;
end
