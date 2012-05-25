function y = gqr_eval(GQR,x,deriv)
% function y = gqr_eval(GQR,x,deriv)
% Here, the GQR object is created by gqr_solve or gqr_rsolve
% which is storage for the various items needed to do GaussQR.
%
% The x should be the locations you want the interpolant computed at.
%  optional: deriv is if you want to evaluate the derivative
%            it is passed the same as in gqr_phi

ep    = GQR.ep;
alpha = GQR.alpha;
coef  = GQR.coef;
Marr  = GQR.Marr;
N     = GQR.N;
reg   = GQR.reg;

if nargin==2
    deriv = zeros(1,size(Marr,1));
end

if reg
    phiEval = gqr_phi(Marr,x,ep,alpha,deriv);
    y = phiEval*coef;
else
    Rbar = GQR.Rbar;
    phiEval1 = gqr_phi(Marr(:,1:N),x,ep,alpha,deriv);
    phiEval2 = gqr_phi(Marr(:,N+1:end),x,ep,alpha,deriv);
    y = phiEval1*coef+phiEval2*Rbar*coef;
end
