%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Linear finite element method for Sturm-Liouville Problem:   %
%           - (pu')' + qu = f, 0 < x < 2pi                     %
%  With Dirichlet boundary condition u(0) = u(2pi) = 0.        %
%  Test functions:                                             %
%         p(x) = 1 + x, q(x) = x,                              %
%         f(x) = (2x+1)sinx - cosx,                            %
%  Exact solution: u(x) = sinx.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all

a = 0; b = 2*pi; % solution interval (a,b)

N = [32,64,128,256,512,1024]; % prepare for grid refinement analysis

for r = 1:length(N) % r + 1 -> grid refinement once
    n = N(r); % number of subintervals (x_{i-1},x_i)
    h = (b-a)/n;
    
    X = linspace(a,b,n+1); X = X'; % generate grid
    
    % generate matrix
    A = zeros(n-1,n-1);   % stiffness matrix
    F = zeros(n-1,1);     % load vector
    
    % integral default setting: absTol = 1e-10, RelTol = 1e-6
    % The integral function wouldn't influence our order.
    A(1,1) = integral(@(x)( p(x).*dhat(x,X(1:3)).^2 ...
        + q(x).*hat(x,X(1:3)).^2 ),X(1),X(2));
    A(n-1,n-1) = integral(@(x)( p(x).*dhat(x,X(n-1:n+1)).^2 ...
        + q(x).*hat(x,X(n-1:n+1)).^2 ),X(n),X(n+1));
    
    F(1) = integral(@(x)( f(x).*hat(x,X(1:3)) ),X(1),X(2));
    F(n-1) = integral(@(x)( f(x).*hat(x,X(n-1:n+1)) ),X(n),X(n+1));
    
    % assemble stiffness matrix and load vector element by element
    bar = waitbar(0,'Assembling Matrices','name',['Computing for n = ',num2str(N(r))]);
    for i = 1:n-2
        A(i,i) = A(i,i) + integral(@(x)( p(x).*dhat(x,X(i:i+2)).^2 ...
            + q(x).*hat(x,X(i:i+2)).^2 ),X(i+1),X(i+2));
        A(i,i+1) = A(i,i+1) + integral(@(x)( p(x).*dhat(x,X(i:i+2)).*dhat(x,X(i+1:i+3)) ...
            + q(x).*hat(x,X(i:i+2)).*hat(x,X(i+1:i+3)) ),X(i+1),X(i+2));
        A(i+1,i) = A(i+1,i) + integral(@(x)( p(x).*dhat(x,X(i+1:i+3)).*dhat(x,X(i:i+2)) ...
            + q(x).*hat(x,X(i+1:i+3)).*hat(x,X(i:i+2)) ),X(i+1),X(i+2));
        A(i+1,i+1) = A(i+1,i+1) + integral(@(x)( p(x).*dhat(x,X(i+1:i+3)).^2 ...
            + q(x).*hat(x,X(i+1:i+3)).^2 ),X(i+1),X(i+2));
        
        F(i) = F(i) + integral(@(x)( f(x).*hat(x,X(i:i+2)) ),X(i+1),X(i+2));
        F(i+1) = F(i+1) + integral(@(x)( f(x).*hat(x,X(i+1:i+3)) ),X(i+1),X(i+2));
        waitbar(i/(n-2),bar,['Assembling Matrices: ' num2str(roundn((i/(n-2)*100),-2)) '%']);
    end
    close(bar);
    U = A\F; U = [0; U; 0]; % FEM solution
    
    % exact solution and error
    uexact = u(X);
    error = U - uexact;
    E(r) = max(abs(error));
    display( sprintf('%-12s %-18s', ['h = 2pi/' num2str(N(r))], ...
        [' Error = ' num2str(E(r))]) );
    
    % plot
    figure(1)
    subplot(1,2,1) % Comparation Plot
    plot(X,U,'o'); hold on;
    plot(X,uexact,'lineWidth',1); grid minor
    xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
    L = legend('FEM solution','Exact solution');
    set(L,'Interpreter','latex')
    title('Comparation Plot','Interpreter','latex')
    subplot(1,2,2) % Error Plot
    plot(X,error,'lineWidth',1); grid minor
    xlabel('x','Interpreter','latex'); ylabel('error','Interpreter','latex')
    title('Error Plot','Interpreter','latex')
    T = suptitle(['Solution Plots $h = 2\pi/' num2str(N(r)) '$']);
    set(T,'Interpreter','latex')
    pause(2)
    close(figure(1))
    
end

% Results of grid refinement analysis
p = log(E(1:length(N)-1)./E(2:length(N)))/log(2);
display( sprintf('\n========== Grid refinement analysis ==========') )
display( sprintf('%12s %12s %12s','h1', 'h2', 'Order') )
for r = 1:length(N)-1
    display( sprintf('%12s %12s %12.4f', ['2pi/' num2str(N(r))], ...
        ['2pi/' num2str(N(r+1))], p(r)) )
end