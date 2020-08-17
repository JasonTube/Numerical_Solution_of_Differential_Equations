function [x,U] = problem1_solver(lb,ub,ua,a,b,c,f,beta,gamma,n)
% solving the self-adjoint boundary value problem:
% (beta(x)u')'- gamma(x)u(x) = f(x) in (a,b) 
% u(lb) = ua  au(ub) + bu'(ub) = c
% xi=lb+ih
h = (ub-lb)/n;
A = sparse(n,n);
F = zeros(n,1);
h2 = h^2;
for i = 2:n-1
    beta1 = feval(beta,lb+(i-1/2)*h);
    beta2 = feval(beta,lb+(i+1/2)*h);
    A(i,i-1) = beta1/h2;
    A(i,i) = -(beta1+beta2)/h2 - feval(gamma,lb+i*h);
    A(i,i+1) = beta2/h2;
end
beta1 = feval(beta,lb+(1/2)*h);
beta2 = feval(beta,lb+(3/2)*h);
A(1,1) = -(beta1+beta2)/h2 - feval(gamma,lb+h);
A(1,2) = beta2/h2;

beta1 = feval(beta,lb+(n-1/2)*h);
beta2 = feval(beta,lb+(n+1/2)*h);
A(n,n-1) = (beta1+beta2)/h2;
A(n,n) = - beta1/h2 - (1+2*h*a/b)*beta2/h2 - feval(gamma,lb+n*h);

for i=1:n
    x(i)=lb+i*h;
    F(i)=feval(f,x(i));
end

beta1 = feval(beta,lb+(1/2)*h);
F(1) = F(1) - beta1*ua/h2;

beta2 = feval(beta,lb+(n+1/2)*h);
F(n) = F(n) - beta2*2*h*c/(h2*b);

if b==0 % deal with b = 0
    A(n,n)=1;
    F(n)=c;
end

U = A\F;
return
