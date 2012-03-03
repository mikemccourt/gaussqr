function y = rbfqr_eval(rbfqrOBJ,x,deriv)
% function y = rbfqr_eval(rbfqrOBJ,x,deriv)
% Here, the rbfqrOBJ object is created by rbfqr_solve or rbfqrr_solve
% which is storage for the various items needed to do RBR-QR.
%
% The x should be the locations you want the interpolant computed at.
%  optional: deriv is if you want to evaluate the derivative
%            it is passed the same as in rbfphi

ep    = rbfqrOBJ.ep;
alpha = rbfqrOBJ.alpha;
coef  = rbfqrOBJ.coef;
Marr  = rbfqrOBJ.Marr;
N     = rbfqrOBJ.N;
reg   = rbfqrOBJ.reg;

if nargin==2
    deriv = [0,0];
end

if reg
    phiEval = rbfphi(Marr,x,ep,alpha,deriv);
    y = phiEval*coef;
else
    Rbar = rbfqrOBJ.Rbar;
    phiEval1 = rbfphi(Marr(:,1:N),x,ep,alpha,deriv);
    phiEval2 = rbfphi(Marr(:,N+1:end),x,ep,alpha,deriv);
    y = phiEval1*coef+phiEval2*Rbar*coef;
end
