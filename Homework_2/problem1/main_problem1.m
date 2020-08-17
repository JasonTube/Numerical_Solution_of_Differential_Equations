%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Crank-Nicolson Method for 2D heat equation:                       %
%          u_t = beta(x,y,t) * (u_{xx} + u_{yy}) + f(x,t)              %
%                                                                      %
%    Test problem:                                                     %
%      u(t,x,y) = exp(-t)*sin(pi*x)*sin(pi*y)                          %
%      beta(x,y,t) = 1                                                 %
%      f(t,x,y) = exp(-t)*sin(pi*x)*sin(pi*y)*(2pi^2-1)                %
%      Dirichlet boundary condition                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all

% Input
N = [8 16 32 64];

a = 0; b = 1; c = 0; d = 1;
t_final = 1;
tol = 1e-8;
maxIter = 1e4;

for r = 1:length(N) % Loop for grid refinement
    n_space = N(r);n_time = N(r);

    % Generate matrix
    hx = (b-a)/n_space; hx2=hx*hx;
    hy = (d-c)/n_space; hy2=hy*hy;
    x = linspace(a,b,n_space+1);
    y = linspace(c,d,n_space+1);
    
    tau = t_final/n_time;
    t = linspace(0,t_final,n_time+1);
    
    [X,Y] = meshgrid(x,y);
    
    u_k = u(X,Y,0);            % u_k and initialize
    
    f_k0 = f(X,Y,0);           % f_k
    f_k1 = f(X,Y,tau);         % f_{k+1}
    
    beta_k0 = beta(X,Y,0);     % beta_k
    beta_k1 = beta(X,Y,tau);   % beta_{k+1}
    
    % Solve
    bar = waitbar(0,'Solving the problem','name',['Compute for n = ',num2str(N(r))]);
    for t = tau:tau:t_final % time loop
        
        % Iterative Method Start
        u_old = ones(n_space+1,n_space+1);
        u_new = u_k;
        q = 0;
        
        while max(max(abs(u_new-u_old))) > tol && q < maxIter
            q = q + 1;
            u_old = u_new;
            
            for i = 2:n_space
                for j = 2:n_space
                    u_new(i,j) =  ( ...
                        beta_k0(i,j)*( (u_k(i-1,j)-2*u_k(i,j)+u_k(i+1,j))/(2*hx2) ...
                        + (u_k(i,j-1)-2*u_k(i,j)+u_k(i,j+1))/(2*hy2) ) ...
                        + u_k(i,j)/tau + (f_k0(i,j)+f_k1(i,j))/2 ...
                        + beta_k1(i,j)*( (u_new(i-1,j)+u_new(i+1,j))/(2*hx2) ...
                        + (u_new(i,j-1)+u_new(i,j+1))/(2*hy2) ) ...
                        ) / (1/tau+1/hx2+1/hy2);
                end
            end
        end
        % Iterative Method end
        
        u_k = u_new; % save u_{k+1}
        
        % update boundary
        u_k(1,:) = u(a,y,t);
        u_k(n_space+1,:) = u(b,y,t);
        u_k(:,1) = u(x',c,t);
        u_k(:,n_space+1) = u(x',d,t);
        
        % update beta and f
        f_k0 = f_k1;
        f_k1 = f(X,Y,t+tau);
        
        beta_k0 = beta_k1;
        beta_k1 = beta(X,Y,t+tau);
        
        waitbar(t/t_final,bar,['Solving the problem: ' num2str(roundn((t/t_final*100),-2)) '%']);
    end
    close(bar)
    % Exact solution at t_final
    u_exact = u(X,Y,t_final);
    
    % Error
    E(r) = max(max(abs(u_k-u_exact)));
    display(['h = ' num2str(hx) ' tau = ' num2str(t_final/n_time)...
        ' Error:' num2str(E(r))]);
    
    % Plots
    figure(r);
    subplot(1,2,1); mesh(x,y,u_k); title('The computed solution');
    xlabel('x');ylabel('y');zlabel('u');grid minor
    subplot(1,2,2); mesh(x,y,u_k-u_exact); title('The error plot');
    xlabel('x');ylabel('y');zlabel('error');grid minor
    suptitle(['Results of Crank-Nicolson with n = ' num2str(n_space)])
    pause(1.5)
    close(figure(r))
    
end

% Compute the order
p = log(E(1:length(N)-1)./E(2:length(N)))/log(2);
display('======== Grid refinement analysis ========')
for r = 1:length(N)-1
    n_space = N(r);n_time = N(r);

    display(['h1 = ' num2str((b-a)/n_space) ' tau1 = ' num2str(t_final/n_time)...
        ' h2 = ' num2str((b-a)/(2*n_space)) ' tau2 = ' num2str(t_final/(2*n_time))...
        ' p = ' num2str(p(r))])
end







