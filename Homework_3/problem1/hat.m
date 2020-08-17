function y = hat(x,X)
y = (x-X(1))/(X(2)-X(1)).*(X(1)<=x).*(x<X(2)) ...
    + (X(3)-x)/(X(3)-X(2)).*(X(2)<=x).*(x<X(3));
end
