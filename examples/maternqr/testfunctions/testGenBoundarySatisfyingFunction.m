syms x
f = sym(franke(x,0.5));
lbsat = 0;
rbsat = 2;
lbound = 0;
rbound = 1;
lbval = 0;
rbval = 0;
g = genBoundarySatisfyingFunction( f, lbsat, rbsat, lbound, rbound, lbval, rbval )

range = [lbound,rbound];
clf

i = 0;
d = symfun(diff(g,2*i),x);
ezplot(d,range);
eval(d(lbound))
eval(d(rbound))

hold off