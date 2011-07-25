function X = ranksolve(U,VT,B)
% This solves a system (eye(n)+U*VT)*X=B
% where [n m]=size(U) and [m n]=size(VT)
% using the identity inv(eye(n)+U*VT)=eye(n)-U*inv(eye(m)+VT*U)*VT
% 
% Depending on the size of U and VT it may be better to directly compute
% the product U*VT and solve I+U*VT directly.  The value
% GAUSSQR_PARAMETERS.RANKSOLVE_PROPORTION set in rbfsetup can choose the
% transition point.

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
r = GAUSSQR_PARAMETERS.RANKSOLVE_PROPORTION;

[n m] = size(U);
[mV nV] = size(VT);
[nB mB] = size(B);

if n~=nV || m~=mV || n~=nB
    error('Unacceptable sizes in ranksolve, remember transpose')
end

if m>=r*n
    X = (eye(n)+U*VT)\B;
else
    X = B-U*((eye(m)+VT*U)\(VT*B));
end