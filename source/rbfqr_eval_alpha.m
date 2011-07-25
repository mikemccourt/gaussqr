function y = rbfqr_eval_alpha(rbfqrOBJ,x)
% function y = rbfqr_eval_alpha(rbfqrOBJ,x)
% Here, the rbfqrOBJ object is created by rbfqr_solve_alpha which is
% storage for the various items needed to do RBR-QR.
%
% The x should be the locations you want the interpolant computed at.

ep    = rbfqrOBJ.ep;
alpha = rbfqrOBJ.alpha;
coef  = rbfqrOBJ.coef;
Marr  = rbfqrOBJ.Marr;
N     = rbfqrOBJ.N;
reg   = rbfqrOBJ.reg;

if reg
    phiEval = rbfphialpha(Marr,x,ep,alpha);
    y = phiEval*coef;
else
    Rbar = rbfqrOBJ.Rbar;
    phiEval1 = rbfphialpha(Marr(:,1:N),x,ep,alpha);
    phiEval2 = rbfphialpha(Marr(:,N+1:end),x,ep,alpha);
    y = phiEval1*coef+phiEval2*Rbar*coef;
end