function y = dhat(x,X)
y = 1/(X(2)-X(1)).*(X(1)<=x).*(x<X(2)) ...
    - 1/(X(3)-X(2)).*(X(2)<=x).*(x<X(3));
end