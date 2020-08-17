function u=ufunc(x,t,k)
u = zeros(length(x),1);
for i = 1:k
    u = u + (1/i^2)*exp(-i^2*pi^2*t)*sin(i*pi/2)*sin(i*pi*x);
end
u = u * (8/pi^2);
end