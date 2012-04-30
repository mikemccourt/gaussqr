function y = mqr_eval(MQR,x,deriv)
% function y = mqr_eval(MQR,x,deriv)
% Here, the MQR object is created by mqr_solve
% which is storage for the various items needed to do MaternQR.
%
% Inputs: MQR - the MaternQR object created by mqr_solve
%         x - column vector of values at which to evaluate interpolant
%         deriv - <default=0> interpolant derivative requested
% Outputs: y - interpolant evaluated at x

ep   = MQR.ep;
L    = MQR.L;
beta = MQR.beta;
N    = MQR.N;
M    = MQR.Mmax;
Rbar = MQR.Rbar;
coef = MQR.coef;

if nargin==2
    deriv = 0;
end

phiEval1 = mqr_phi(1:N,x,L,deriv);
phiEval2 = mqr_phi(N+1:M,x,L,deriv);
y = phiEval1*coef+phiEval2*Rbar*coef;
