function y = rbfqr_eval(rbfqrOBJ,x)
% function y = rbfqr_eval(rbfqrOBJ,x)
% Here, the rbfqrOBJ object is created by rbfqr_solve which is storage
% for the various items needed to do RBR-QR.
%
% The x should be the locations you want the interpolant computed at.

ep   = rbfqrOBJ.ep;
a    = rbfqrOBJ.a;
beta = rbfqrOBJ.beta;
Marr = rbfqrOBJ.Marr;
N    = rbfqrOBJ.N;
reg  = rbfqrOBJ.reg;

if reg
    phiEval = rbfphi(Marr,x,ep,a);
    y = phiEval*beta;
else
    Rbar = rbfqrOBJ.Rbar;
    phiEval1 = rbfphi(Marr(:,1:N),x,ep,a);
    phiEval2 = rbfphi(Marr(:,N+1:end),x,ep,a);
    y = phiEval1*beta+phiEval2*Rbar*beta;
end