%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Box scheme for one-way wave equation:                        %
%                 u_t + a(x,t)u_x = f(x,t)                       %
%                                                                %
%    Test problem:                                               %
%      u(x,t) = sin(x + t)                                       %
%      a(x,t) = 1 + xt                                           %
%      f(x,t) = (2+ xt)cos(x+t)                                  %
%    Boundary condition:                                         %
%      u(x,0) = sin(x)   IC                                      %
%      u(0,t) = sin(t)   BC                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

c = 0; d = 1;
t_final = 1;

N = [8,16,32,64,128];

for r = 1:length(N)
    n_space = N(r); n_time = N(r);
    
    h = (d-c)/n_space;
    tau = t_final/n_time;
    
    x = c:h:d; % generate grid
    
    u1 = u0(x);         % save u_k
    u2 = zeros(1,length(x));  % save u_{k+1}
    
    figure(1); % dynamic presentation
    plot(x,u1);
    grid minor;title('Dynamic presentation','Interpreter','latex')
    xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
    
    bar = waitbar(0,'Solving the problem','name',['Compute for h = ',num2str(h),' & tau = ',num2str(tau)]);
    for t = tau:tau:t_final
        % update left boundary
        u2(1) = bc(t);
        % update rest nodes
        for i = 2:length(x)
            u2(i) = ( f((x(i-1)+x(i))/2,t-tau/2) ...
                + (-u2(i-1)+u1(i-1)+u1(i))/(2*tau) ...
                + a((x(i-1)+x(i))/2,t-tau/2)*(u2(i-1)-u1(i)+u1(i-1)/(2*h)) ) ...
                / (1/(2*tau)+a((x(i-1)+x(i))/2,t-tau/2)/(2*h));
        end
        u1 = u2;
        plot(x,u1);hold on
        waitbar(t/t_final,bar,['Solving the problem: ' num2str(roundn((t/t_final*100),-2)) '%']);
    end
    close(bar)
    pause(0.5)
    close(figure(1))
    
    % exact solution
    uexact = u(x,t_final);
    
    % error
    error = abs(u1 - uexact);
    E(r) = max(abs(u1 - uexact));
    display(['h = ' num2str(h) ' tau = ' num2str(tau)...
        ' Error:' num2str(E(r))]);
end

% Compute the order
p = log(E(1:length(N)-1)./E(2:length(N)))/log(2);
display('======== Grid refinement analysis ========')
for r = 1:length(N)-1
    n_space = N(r); n_time = N(r);
    display(['h1 = ' num2str((d-c)/n_space) ' tau1 = ' num2str(t_final/n_time)...
        ' h2 = ' num2str((d-c)/(2*n_space)) ' tau2 = ' num2str(t_final/(2*n_time))...
        ' p = ' num2str(p(r))])
end
