% Lax-Wendroff scheme for the one-way wave equation:
% u_t + a*u_x = 0, u(x,0) = u0(x)
clear; clc;close all

a = 3;
t_final = 1;

N = [1,2,4,8,16];% prepare for grid refinement analysis

for r = 1:length(N)
    
    h = 0.05/N(r);
    tau = 0.01/N(r);
    
    x = 0:h:2;% using the data on [0, 2]
    
    n = length(x)-1;
    miu = a*tau/h;
    
    U = u0(x); % initial data
    
    figure(1);plot(x,U);hold on
    xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
    title('Dynamic presentation (plot every 10 steps)','Interpreter','latex');grid minor;
    
    bar = waitbar(0,'Solving the problem','name',['Compute for h = ',num2str(h),' & tau = ',num2str(tau)]);
    for t = tau:tau:t_final
        % Lax-Wendroff Update
        U(2:n) = (miu^2+miu)*U(1:n-1)/2 + (-miu^2+1)*U(2:n) + ...
            (miu^2-miu)*U(3:n+1)/2;
        
        U(1) = bc(t); % update left bound
        
        % update right bound
        % use previous time levels for the interpolation
        U(n+1) = U(n)*(x(n+1)-x(n-1))/(x(n)-x(n-1)) + ...
            U(n-1)*(x(n+1)-x(n))/(x(n-1)-x(n));
        
        waitbar(t/t_final,bar,['Solving the problem: ' num2str(roundn((t/t_final*100),-2)) '%']);
        
        % plot
        if mod(t,10*tau)==0 % plot every 10 steps
            plot(x,U);hold on
        end
    end
    close(bar)
    pause(1)
    close(figure(1))
    
    uexact = u(x,t_final);
    figure(2)
    subplot(1,2,1)
    plot(x,uexact);hold on 
    plot(x,U,'o');
    xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
    title('Compare the final solution','Interpreter','latex');grid minor;
    L = legend('Exact solution','Computed solution'); 
    set(L,'Interpreter','latex','Location','Southeast')
    subplot(1,2,2)
    plot(x,U - uexact);
    xlabel('x','Interpreter','latex'); ylabel('error','Interpreter','latex')
    title('Error plot','Interpreter','latex');grid minor;
    pause(2);close(figure(2))
    
    % error
    E(r) = max(abs(U - uexact));
    display(['h = ' num2str(h) ' tau = ' num2str(tau) ' Error:' num2str(E(r))]);
end

% Compute the order
p = log(E(1:length(N)-1)./E(2:length(N)))/log(2);
display('======== Grid refinement analysis ========')
for r = 1:length(N)-1
    h = 0.05/N(r);
    tau = 0.01/N(r);
    display(['h1 = ' num2str(h) ' tau1 = ' num2str(tau)...
        ' h2 = ' num2str(h/2) ' tau2 = ' num2str(tau/2)...
        ' p = ' num2str(p(r))])
end
