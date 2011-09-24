function [alpha,k] = rbfalphasearch(ep,a,b)
% function alpha = rbfalphasearch(ep,a,b)
% ep - shape parameter for the Gaussians
% a - minimum bound(s) of domain
% b - maximum bound(s) of domain
%   In multiple dimensions these should be a=[low_1 ... low_d]
%                                          b=[high_1 ... high_d]
%
% This function tries to return the smallest alpha for which orthonormality
% is maintained throughout the domain until a set eigenfunction.
%
% You can set the top eigenfunction for which orthonornmality should be
% maintained in rbfsetup with the value
%      GAUSSQR_PARAMETERS.DEFAULT_ORTH_REQUESTED
%
% If this function cannot find an alpha up to what was requested, it will
% try to get as high as possible.  The index up to which orthogonality was
% maintained is returned as k.

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
orthmax = GAUSSQR_PARAMETERS.DEFAULT_ORTH_REQUESTED;

k = orthmax;
alpha = fminbnd(@(alpha)-rbforthintegral(k,ep,alpha,a,b),1e-3,1e3);

end

function integral = rbforthintegral(k,ep,alpha,a,b)
% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
tol = GAUSSQR_PARAMETERS.DEFAULT_ORTH_TOLERANCE;

weight = @(a,x) a/sqrt(pi)*exp(-(a*x).^2);
if exist('quadgk')
%    quadgkEXISTS = true;
    intappx = quadgk(@(x)rbfphialpha(k,x',ep,alpha)'.^2.*weight(alpha,x),a,b);
else
    intappx = quadl(@(x)rbfphialpha(k,x',ep,alpha)'.^2.*weight(alpha,x),a,b);
end
integral = (abs(1-intappx)<tol)*(1/alpha);

end