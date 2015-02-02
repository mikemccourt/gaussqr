function gw = CGW_eval(x,y,rbf,eparr,xeval)
% function gw = CGW_eval(x,y,rbf,eparr,xeval)
% Inputs: x     - data sites
%         y     - data values
%         rbf   - radial kernel in anisotropic form
%         eparr - anisotropic shape parameter
%         xeval - data sites at which to test the power function
% Outputs: gw   - the C_GW criterion (Kriging variance)
K = rbf(DistanceMatrix(x,x,eparr));
Keval = rbf(DistanceMatrix(xeval,x,eparr));
L = chol(K,'lower');
Linvy = L\y;
hsnorm = norm(Linvy);
pow = norm(sqrt(1 - sum((Keval/L').^2,2)),'inf');
gw = hsnorm*pow;