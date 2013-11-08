% ex12.m
% This example shows the effect of ill-conditioning on derivatives
% This considers a few different methods of approximation
%
% Specifically, the functions we are considering here are
%   u(x) = x^2+x
%   u(x) = x*sin(x)
%   u(x) = 1/(1+x^2)
%   u(x) = sinh(x)/(1+cosh(x))
% And we consider different methods of approximating the first two
% derivatives on the domain [-1,1]
%
% The methods we consider are:
%   Trefethen: Cheb
%   Collocation
%   GaussQR
%   GaussQRr

rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4; % Use absolute error
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

N = 25;
fopt = 3;

switch fopt
    case 1
        uf = @(x) x.^2+x;
        ufd = @(x) 2*x+1;
        ufd2 = @(x) 2*ones(size(x));
    case 2
        uf = @(x) x.*sin(x);
        ufd = @(x) sin(x) + x.*cos(x);
        ufd2 = @(x) 2*cos(x)-x.*sin(x);
    case 3
        uf = @(x) 1./(1+x.^2);
        ufd = @(x) -2*x./(1+x.^2).^2;
        ufd2 = @(x) (6*x.^2-2)./(1+x.^2).^3;
    case 4
        uf = @(x) sinh(x)./(1+cosh(x));
        ufd = @(x) 1./(1+cosh(x));
        ufd2 = @(x) -sinh(x)./(1+cosh(x)).^2;
end

rbf = @(e,r) exp(-(e*r).^2);
drbf = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

[D,x] = cheb(N);
u = uf(x);
ud = ufd(x);
ud2 = ufd2(x);
D2 = D^2;
err_Trefethen_1d = errcompute(D*u,ud);
err_Trefethen_2d = errcompute(D2*u,ud2);

epvec = logspace(-2,1,60);
r = DistanceMatrix(x,x);
dx = DifferenceMatrix(x,x);
errvec1d = [];
errvec2d = [];
k = 1;
for ep=epvec
  Amat = rbf(ep,r);
  Dmat = drbf(ep,r,dx);
  D2mat = d2rbf(ep,r);
  b = Amat\u;
  errvec1d(k) = errcompute(Dmat*b,ud);
  errvec2d(k) = errcompute(D2mat*b,ud2);
  k = k + 1;
end

alpha = 3;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .25;

errvecQ1d = [];
errvecQ2d = [];
errvecR2d = [];
errvecR1d = [];
k = 1;
for ep=epvec
    GQR = gqr_solve(x,u,ep,alpha);
    errvecQ1d(k) = errcompute(gqr_eval(GQR,x,1),ud);
    errvecQ2d(k) = errcompute(gqr_eval(GQR,x,2),ud2);
    
    GQR = gqr_rsolve(x,u,ep,alpha);
   errvecR1d(k) = errcompute(gqr_eval(GQR,x,1),ud);
    errvecR2d(k) = errcompute(gqr_eval(GQR,x,2),ud2);
  
    k = k + 1;
end

loglog(epvec,errvec1d,'--b','linewidth',2),hold on
loglog(epvec,errvec2d,'--r','linewidth',2)
loglog(epvec,errvecQ1d,'b','linewidth',3)
loglog(epvec,errvecQ2d,'r','linewidth',3)
loglog(epvec,errvecR1d,'-^b','linewidth',3)
loglog(epvec,errvecR2d,'-^r','linewidth',3)
loglog(epvec,err_Trefethen_1d*ones(size(epvec)),':b')
loglog(epvec,err_Trefethen_2d*ones(size(epvec)),':r')
hold off
