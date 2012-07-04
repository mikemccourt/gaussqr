syms x
f = sym(franke(sym('x'),0.5));

i = 5; % 2*i is the order of the derivative to examine
lbsat = 5;
rbsat = 5;
lbound = 0;
rbound = 1;
lbval = 1e10;
rbval = 1e10;
range = [lbound,rbound];
test = 0; % 0 = add, 1 = mult

if test == 0 % test additive version
    g = forceBCsatADD( f, lbsat, rbsat, lbound, rbound, lbval, rbval )
    g = g / 10^(2*i) % normalization factor
    d = symfun(diff(g,2*i),x);
    figure;
    ezplot(d,range);
    disp(['L boundary val: ',num2str(eval(d(lbound)))]);
    disp(['R boundary val: ',num2str(eval(d(rbound)))]);
    
elseif test == 1 % test multiplicative version
    g = forceBCsatMULT( f, lbsat, rbsat, lbound, rbound )
    d = symfun(diff(g,2*i),x);
    figure;
    ezplot(d,range);
    disp(['L boundary val: ',num2str(eval(d(lbound)))]);
    disp(['R boundary val: ',num2str(eval(d(rbound)))]);
    
end