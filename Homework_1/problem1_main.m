% Main function for Solving the self_adjoint BVP
% (beta(x)u')'- gamma(x)u(x) = f(x) in (a,b) 
% u(lb) = ua, au(ub) + bu'(ub) = c
clear;close all;
%% Input
lb = 0;
ub = 1;
ua = ufunc(lb);
a = 2;
b = -3;
c = a*ufunc(ub)+b*dufunc(ub);
%% solve
n = 80;% grid
[x,U] = problem1_solver(lb,ub,ua,a,b,c,@func,@beta,@gamma,n);
for i = 1:n
u(i) = ufunc(x(i));
end
%% plot
figure()
subplot(1,2,1)
plot(x,u,'o');hold on
plot(x,U,'-');
grid minor;title('Solution plots');xlabel('x');ylabel('u')
legend('exact solution', 'numerical solution')
subplot(1,2,2)
plot(x,U-u');
grid minor;title('Error plot');xlabel('x');ylabel('error')
suptitle('Numerical Results with n = 80')
%% grid redinement analysis
n = 4;
for i = 1:17
H(i) = 1/n;
[x,U] = problem1_solver(lb,ub,ua,a,b,c,@func,@beta,@gamma,n);
u = [];
for j = 1:n
u(j) = ufunc(x(j));
end
error(i) = max(abs(U-u'));
n = 2*n;
end
p = log(error(1:16)./error(2:17))/log(2);
figure();
subplot(1,2,1);loglog(H,error);grid minor;
axis('equal'); axis('square'); xlabel('h');ylabel('error')
gtext('Slope of the method is 2')
subplot(1,2,2);plot(log10(H(2:17)),p);grid minor;
axis('equal'); axis('square'); xlabel('log_{10}h');ylabel('Order p')
suptitle('Grid refinement analysis')
%% functions
function b = beta(x)
b = 1 + x^2;
end

function g = gamma(x)
g = x;
end

function u = ufunc(x)
u = exp(-x)*((x-1)^2);
end

function du = dufunc(x)
du = exp(-x)*(2*x - 2) - exp(-x)*(x - 1)^2;
end

function d2u = d2ufunc(x)
d2u = 2*exp(-x) - 2*exp(-x)*(2*x - 2) + exp(-x)*(x - 1)^2;
end

function f=func(x)
% (beta(x)u')'- gamma(x)u(x) = f(x)
b = 1 + x^2;
db = 2*x;
g = x;
u = exp(-x)*((x-1)^2);
du =  exp(-x)*(2*x - 2) - exp(-x)*(x - 1)^2;
d2u =  2*exp(-x) - 2*exp(-x)*(2*x - 2) + exp(-x)*(x - 1)^2;
f = (db*du)+(b*d2u)-g*u;
end
