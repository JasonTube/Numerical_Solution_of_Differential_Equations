clear; clc; close all
% Theta-Method for the initial-boundary value problem:
% u_t = u_xx, u(x,0)=u0, u(a,t)=ua, u(b,t)=ub
%% Input
a = 0; b = 1;                   % left and right boundaries
t = 0.1;                        % final time
h = [0.1 0.05 0.025 0.0125];    % space step
tau = h/2;                      % time step
%% Solve and Plot
figure();

subplot(1,2,1);% solution plot
subplot(1,2,2);% error plot

for i=1:length(h)
    % call the Theta-Method solver: 
    % theta = 0.5: Crank-Nicolson Scheme, theta = 0: Backward-Euler Method
    [U,x]=theta_method(a,b,h(i),t,tau(i),@u0func,@uafunc,@ubfunc,0.5);
    
    subplot(1,2,1);plot(x,U,'-+','lineWidth',1);hold on
    
    % Error
    %E(i) = max(max(abs(U-ufunc(x,t,5))));% maximum error
    E(i) = sum(abs(U-ufunc(x,t,5)))/length(U);% average error
    display(['h = ' num2str(h(i)) ' tau = ' num2str(tau(i))...
        ' Error:' num2str(E(i))]);
    
    subplot(1,2,2);plot(x,U-ufunc(x,t,5),'lineWidth',1);hold on
end

subplot(1,2,1);
plot(x,ufunc(x,t,5),'lineWidth',2);% plot of exact solution

grid minor;title('Solution plots','Interpreter','latex')
xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
L = legend('Computed solution h = 0.1','Computed solution h = 0.05',...
    'Computed solution h = 0.025','Computed solution h = 0.0125',...
    'Exact solution'); set(L,'Interpreter','latex','Location','South')

subplot(1,2,2);
grid minor;title('Error plots','Interpreter','latex')
xlabel('x','Interpreter','latex'); ylabel('u','Interpreter','latex')
L = legend('Computed solution h = 0.1','Computed solution h = 0.05',...
    'Computed solution h = 0.025','Computed solution h = 0.0125'); 
set(L,'Interpreter','latex')

%% Compute the order
p = log(E(1:length(h)-1)./E(2:length(h)))/log(2);
display('======== Grid refinement analysis ========')
for i = 1:length(h)-1
    display(['h1 = ' num2str(h(i)) ' tau1 = ' num2str(tau(i))...
        ' h2 = ' num2str(h(i+1)) ' tau2 = ' num2str(tau(i+1))...
        ' p = ' num2str(p(i))])
end