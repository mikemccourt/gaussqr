function s = ppsplinekernel(x,z,L,beta,deriv)
% function s = ppsplinekernel(x,z,L,beta,deriv)
% This function evaluates the kernel associated with the ep->0 limit of the
% Compact Matern kernel.  It has the value
%    (2L)^(2beta-1)/(2beta!) (B_2beta((x-z)/(2L)) - B_2beta((x+z)/(2L)))
%          for x>z, and
%    (2L)^(2beta-1)/(2beta!) (B_2beta((z-x)/(2L)) - B_2beta((x+z)/(2L)))
%          for x<z
%
% Inputs : x - column vector of evaluation points
%          z - column vector of kernel centers
%            * IF z is a single value, this will use that as all centers
%          L - size of domain of BVP
%          beta - <optional, default = 1> order of the differential operator
%          deriv - <optional, default = 0> order of derivative
% Outputs : s - kernel function evaluation

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end

% This allows for passing just one center if it is the same for all the
% functions that we want to evaluate
if length(z)==1
    z = z*ones(size(x));
end

% This replaces defaults if none are passed
if nargin<5
    deriv = 0;
    if nargin<4
        beta = 1;
    end
end

[rx,cx] = size(x);
[rz,cz] = size(z);

if rx~=rz
    error('dimension mismatch; size(x,1)=%d, size(z,1)=%d',rx,rz)
elseif cx~=1 | cz~=1
    error('only column vectors for x and z allowed; size(x,2)=%d, size(z,2)=%d',cx,cz)
elseif L<=0 | length(L)~=1 | imag(L)~=0
    error('unacceptable length value L=%g',L)
end

if min(x)<0 | max(x)>L
    error('This function can only be evaluated on [0,%g]; min(x)=%g, max(x)=%g',L,min(x),max(x))
elseif min(z)<0 | max(z)>L
    error('This function must have centers in [0,%g]; min(z)=%g, max(z)=%g',L,min(z),max(z))
end

if beta<1 | length(beta)~=1 | imag(beta)~=0 | floor(beta)~=beta
    warning('unacceptable order beta=%g, reset to 1',beta)
end
if deriv<0 | deriv>4 | length(deriv)~=1 | imag(deriv)~=0 | floor(deriv)~=deriv
    warning('unacceptable derivative deriv=%g, reset to 0',deriv)
    deriv = 0;
end

if deriv>0 & beta==1
    error('First order kernel is not differentiable; deriv=%d, beta=%d',deriv,beta)
elseif deriv>2 & beta==2
    error('Second order kernel can have at most 2 derivatives; deriv=%d, beta=%d',deriv,beta)
end

bp = BernoulliPoly(2*beta-deriv,(x+z)/(2*L));
b1 = BernoulliPoly(2*beta-deriv,(x-z)/(2*L));
b2 = BernoulliPoly(2*beta-deriv,(z-x)/(2*L));
s1 = (b1 - bp).*(x>=z);
s2 = ((-1)^deriv*b2 - bp).*(x<z);
s = (-1)^(beta+1)*exp((2*beta-1-deriv)*log(2*L)-gammaln(2*beta+1-deriv))*(s1+s2);