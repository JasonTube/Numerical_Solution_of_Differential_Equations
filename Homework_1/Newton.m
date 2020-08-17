function x = Newton(x,func,dfunc,K,theta,n,tol)
% Newton Iteration Method
k = 0; b = ones(n-1,1);
while norm(b) > tol
    b = feval(func,x,K,theta,n);
    A = feval(dfunc,x,K,n);
    deta = linsolve(A,-b);
    x = x+deta;
    k = k + 1;
    if norm(deta)<= tol^2 % MinStep
        disp('Newton: Change in X too small.');return;
    elseif(k > 1000) % MaxIter
        disp('Newton£ºToo many function iterations.');return;
    end
end
disp('Newton converged to a root!');
end
