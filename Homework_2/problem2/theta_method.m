function [U,x]=theta_method(a,b,h,t,tau,u0func,uafunc,ubfunc,theta)
% implicit solver for the initial-boundary value problem:
% u_t = u_xx, u(x,0)=u0, u(a,t)=ua, u(b,t)=ub
% using the Theta-Method

% theta = 1: Forward Euler Method
% theta = 0: Backward Euler Method
% theta =0.5: Crank-Nicolson Method

n = (b-a)/h;
x = a+h:h:b-h; x = x';% generate a grid
U = feval(u0func,x);

try
    A = (2*(1-theta)/(h^2)+1/tau) * eye(n-1);
    A(2:n-1,1:n-2) = A(2:n-1,1:n-2)+((theta-1)/(h^2)) * eye(n-2);
    A(1:n-2,2:n-1) = A(1:n-2,2:n-1)+((theta-1)/(h^2)) * eye(n-2);
catch Exception
    % if memory is not enough
    disp(Exception)
    disp('Switch to sparse form')
    A = sparse(n-1,n-1);
    tmp1 = 2*(1-theta)/(h^2)+1/tau;
    tmp2 = (theta-1)/(h^2);
    for i = 1:n-2
        A(i,i) = tmp1; A(i,i+1) = tmp2; A(i+1,i) = tmp2;
    end
    A(n-1,n-1) = tmp1;
end

bar = waitbar(0,'Solving the problem','name',['Compute for h = ',num2str(h),' & tau = ',num2str(tau)]);
for timer = tau:tau:t
    tmp = [feval(uafunc,timer-tau);U(1:n-2)] - 2*U + [U(2:n-1);feval(uafunc,timer-tau)];
    % tmp = U_k_(i-1) - 2*U_k_i + U_k_(i+1)
    F = U/tau + theta * tmp / (h^2);
    F(1) = F(1) + ((1-theta)/(h^2)) * feval(uafunc,timer);
    F(n-1) = F(n-1) + ((1-theta)/(h^2)) * feval(ubfunc,timer);
    
    U = A\F;
    waitbar(timer/t,bar,['Solving the problem: ' num2str(roundn((timer/t*100),-2)) '%']);
end
close(bar)

% add boundary values
U = [feval(uafunc,t);U;feval(uafunc,t)];
x = [a;x;b];
end
