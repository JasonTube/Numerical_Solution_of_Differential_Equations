% Main function for Solving the elliptic problem
% u_xx + p*u_yy + r*u = f a<x<b c<y<d 
% u(a,y) = 0, u(x,c) = 0, u(x,d) = 0, u_x(b,y) = -pi*sin(pi*y)
clear;clc;
%% Input
a = 0; b = 1; c = 0; d = 1;
n = 16;% grid
omega = 1.85;% the best omega found by testing
tol = 1e-5;
maxIter = 1e4;
%%
hx = (b-a)/n; hx2=hx*hx;
hy = (d-c)/n; hy2=hy*hy;
x = linspace(a,b,n+1);
y = linspace(c,d,n+1);

P = zeros(n+1,n+1);% generate P_ij
R = zeros(n+1,n+1);% generate R_ij
F = zeros(n+1,n+1);% generate F_ij
B = zeros(1,n+1);
for j = 1:n+1
    for i = 1:n+1
        P(i,j) = p_func(x(i),y(j));
        R(i,j) = r_func(x(i),y(j));
        F(i,j) = f_func(x(i),y(j));
    end
    B(j) = b_func(y(j));
end
u_old = ones(n+1,n+1);
u_new = zeros(n+1,n+1);
for i=1:n+1
    u_new(1,i) = u_func(a,y(i));
    u_new(i,1) = u_func(x(i),c);
    u_new(i,n+1) = u_func(x(i),d);
end
%% SOR with the best omega
k = 0;
tic;
while max(max(abs(u_new-u_old))) > tol
    k = k + 1;
    u_old = u_new;
    for i = 2:n+1
        for j = 2:n
            if i==n+1
                u_new(i,j) = (1-omega)*u_old(i,j) + omega*(...
                    ( F(i,j) - 2*B(j)/hx - 2*u_new(i-1,j)/hx2 ...
                    - P(i,j).*(u_new(i,j+1)+u_new(i,j-1))/hy2 )...
                    ./(-2/hx2 + -2*P(i,j)/hy2 + R(i,j)) );
            else
                u_new(i,j) = (1-omega)*u_old(i,j) + omega*(...
                    ( F(i,j) - (u_new(i+1,j)+u_new(i-1,j))/hx2 ...
                    - P(i,j).*(u_new(i,j+1)+u_new(i,j-1))/hy2 ) ...
                    ./(-2/hx2 + -2*P(i,j)/hy2 + R(i,j)));
            end
        end
    end
end
toc;
display(['Iterations:' num2str(k)]);
%% exact solution
u_exact = zeros(n+1,n+1);
for i = 1:n+1
    for j = 1:n+1
        u_exact(i,j) = u_func(x(i),y(j));
    end
end
%%
e = max(max(abs(u_new-u_exact)));
display(['Error:' num2str(e)]);
figure();
subplot(1,2,1); mesh(x,y,u_new'); title('The computed solution')
xlabel('x');ylabel('y');zlabel('u');grid minor
subplot(1,2,2); mesh(x,y,u_new'-u_exact'); title('The error plot')
xlabel('x');ylabel('y');zlabel('error');grid minor
suptitle('Results of SOR Method with the best \omega = 1.85 and n = 32')
%% test functions
function u = u_func(x,y)
u = sin(pi*x)*sin(pi*y);
end
function p = p_func(x,y)
p = 1 + x^2 + y^2;
end
function r = r_func(x,y)
r = -x*y;
end
function f = f_func(x,y)
p = 1 + x^2 + y^2;
r = -x*y;
u = sin(pi*x)*sin(pi*y);
u_xx = -(pi^2)*sin(pi*x)*sin(pi*y);
u_yy = -(pi^2)*sin(pi*x)*sin(pi*y);
f = u_xx + p*u_yy + r*u;
end
function b = b_func(y)
b = -pi*sin(pi*y);
end
