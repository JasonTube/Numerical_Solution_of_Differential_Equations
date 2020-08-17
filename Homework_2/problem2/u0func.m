function u = u0func(x)
u = 2.*x.*(x<=0.5)+2.*(1-x).*(x>0.5);
end