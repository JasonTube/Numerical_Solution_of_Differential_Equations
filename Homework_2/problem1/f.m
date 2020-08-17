function f = f(x,y,t) 
f = exp(-t).*sin(pi*x).*sin(pi*y).*(2*pi^2 .* beta(x,y,t) -1);
end