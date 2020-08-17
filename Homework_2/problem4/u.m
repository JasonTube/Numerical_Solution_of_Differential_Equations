function y = u(x,t)
% exact solution
y = (x-3*t).^2 .* (x >= 3*t) + cos(pi*(t-x/3)-pi/2) .* (x < 3*t);

%if x >= 3*t
%    y = (x-3*t).^2;
%else
%    y = cos(pi*(t-x/3)-pi/2);
%end
end
