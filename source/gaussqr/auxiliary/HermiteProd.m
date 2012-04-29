function [Hx,negInd] = HermiteProd(m,x,logoption)
% This function takes a vector m and a matrix x
% [mr mc] = size(m)
% [xr xc] = size(x)
% We need xc==mr since that is the number of dimensions
% We need mc==1 since we can only evaluate one product at a time
% xr is the number of points we need to evaluate at
% If logoption is activated, the value returned will be log of Hx
%   Note that with logoption on you will have complex Hx
%   You should see real(exp(HermiteProd(m,x,1)))=HermiteProd(m,x) to near machine precision 
%
% At some point it will probably be cost effective to exploit the
% recurrence relation when the log option isn't needed
% To do so I'll have to change the structure from
%   Compute each column individually
% to
%   Compute 3 or more columns at a time, and use H_n+1 = 2xH_n - 2nH_n-1
% Doing so wouldn't be appropriate for the log option, I don't think

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GPHAMI = GAUSSQR_PARAMETERS.HERMITE_ASYMPTOTIC_MIN_INDEX;
GPHC = GAUSSQR_PARAMETERS.HERMITE_COEFFICIENTS;

[mr mc] = size(m);
[xr xc] = size(x);

if xc~=mr
    error(sprintf('Dimension mismatch: xc=%d, mr=%d',xc,mr))
end
if mc~=1
    error(sprintf('m must be a column vector: mr=%d, mc=%d',mr,mc))
end

if nargin==3 && (logoption~=0 || logoption~=false)
    logoption=1;
    Hx = zeros(xr,1);
else
    logoption=0;
    Hx = ones(xr,1);
end
negInd = zeros(xr,1);

for k=1:length(m)
    if m(k)<GPHAMI
        if m(k)<=length(GPHC)-1 % To allow for changes on the fly
            H = polyval(GPHC{m(k)+1},x(:,k));
        else
            H = polyval(HermitePoly(m(k)),x(:,k));
        end
        if logoption
            nI = H<0;
            H = log(abs(H)+eps);
        end
    else
        [H,nI] = HermiteAppx(m(k),x(:,k),logoption);
    end
    
    if logoption
        Hx = Hx + H;
        negInd = xor(nI,negInd);
    else
        Hx = Hx.*H;
    end
end