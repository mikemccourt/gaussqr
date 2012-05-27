function err = errcompute(x,y,n)
% function err = errcompute(x,y,n)
% This function takes two vectors and computes their difference
% It may use multiple different ways to compute difference
%
% Input arguments
%   x - computed vector
%   y - true vector
% Optional arguments
%   n - length of vector x and y
%       to be used if some values of x and y are identical by design

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
normtype = GAUSSQR_PARAMETERS.NORM_TYPE;
errstyle = GAUSSQR_PARAMETERS.ERROR_STYLE;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

if any(size(x)~=size(y))
    error('Input vectors have different sizes')
elseif size(x,2)~=1 || size(y,2)~=1
    if alertuser
        warning('GAUSSQR:nonVectorError','Non-vectors passed, computing 2 norm of (:) form')
    end
    x = x(:);
    y = y(:);
end

if not(exist('n'))
    n = size(x,1);
end

if ~(normtype==1 || normtype==2 || normtype==inf)
    error('GP.NORM_TYPE = %g is unacceptable',normtype)
end

switch errstyle
    case 1
        err = (norm((x-y)./(abs(y)+eps),normtype)+eps)/sqrt(n);
    case 2
        err = norm(x-y,normtype);
    case 3
        err = norm(x-y,normtype)/norm(y,normtype);
    case 4
        err = norm(x-y,normtype)/norm(y,normtype)/sqrt(n);
    otherwise
        error('GP.ERROR_STYLE = %d, unacceptable',errstyle)
end
