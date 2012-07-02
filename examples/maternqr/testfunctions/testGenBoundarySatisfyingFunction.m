f = exp(x);
beta = 5;
lbound = 0;
ubound = 1;
lbval = 0;
ubval = 1;
g = genBoundarySatisfyingFunction( f, beta, lbound, ubound, lbval, ubval )

range = [lbound,ubound];
clf
hold on
ezplot(g,range);
for i = 2.*(1:beta)
    ezplot(diff(g,i),range);
end
hold off