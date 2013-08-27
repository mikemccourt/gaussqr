%ex19a_gqr.m
%This example will help us continue our comparison that we saw in
%ex19_gqr.m, except that we will assign a value for epsilon and explore the
%plots for different values of N
global GAUSSQR_PARAMETERS

Nvec = logspace(10, 1000, 31);

epsilon = 10^(-2);
NN = 200;
% x = pickpoints(-1,1,
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
yf = @(x) x.^3-3*x.^2+2*x+1;
fstring = 'y(x) = x^3-3x^2+2x+1';
yf = @(x) 4*tan(2*x+6);
fstring  = 'y(x) = 4tan(2x+6)';
fstring = sprintf('%s, \epsilon = %d',fstring,epsilon);

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;