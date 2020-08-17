% solve the steady state temperature distribution
clear;clc;
%% Input
a = 0; b = 1; c = 0; d = 1;
n = 36;% grid
tol = 1e-5;
maxIter = 1e4;
f= @(x,y,n) n^2*(x==0.5)*(y==0.5);% deta function
%%
hx = (b-a)/n; hx2=hx*hx;
hy = (d-c)/n; hy2=hy*hy;
x = linspace(a,b,n+1);
y = linspace(c,d,n+1);

P = ones(n+1,n+1);% generate P_ij
R = zeros(n+1,n+1);% generate R_ij
F = zeros(n+1,n+1);% generate F_ij
B = -ones(1,n+1);% generate Neumann boundary
for j = 1:n+1
    for i = 1:n+1
        F(i,j) = f(x(i),y(j),n);
    end
end
u_old = ones(n+1,n+1);
u_new = zeros(n+1,n+1);
for i=1:n+1
    u_new(1,i) = 0;
    u_new(i,1) = 0;
    u_new(i,n+1) = 0;
end
%% Gauss_Seidel
k = 0;
tic;
while max(max(abs(u_new-u_old))) > tol
    k = k + 1;
    u_old = u_new;
    for i = 2:n+1
        for j = 2:n
            if i==n+1
                u_new(i,j) = ( F(i,j) - 2*B(j)/hx - 2*u_new(i-1,j)/hx2 ...
                    - P(i,j).*(u_new(i,j+1)+u_new(i,j-1))/hy2 )...
                    ./(-2/hx2 + -2*P(i,j)/hy2 + R(i,j));
            else
                u_new(i,j) = ( F(i,j) - (u_new(i+1,j)+u_new(i-1,j))/hx2 ...
                    - P(i,j).*(u_new(i,j+1)+u_new(i,j-1))/hy2 ) ...
                    ./(-2/hx2 + -2*P(i,j)/hy2 + R(i,j));
            end
        end
    end
end
toc;
display(['Iterations:' num2str(k)]);
%% plot
figure();
subplot(1,2,1)
mesh(x,y,u_new'); grid minor; title('3D Mesh Plot')
xlabel('x');ylabel('y');zlabel('u')
subplot(1,2,2)
contour(x,y,u_new',30); grid minor; title('Contour Plot')
xlabel('x');ylabel('y');
suptitle('The computed solution with n = 36')
