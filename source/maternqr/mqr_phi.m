function p = mqr_phi(Marr,x,L,deriv)
% function p = mqr_phi(Marr,x,L,deriv)
% This function computes the eigenfunctions for the Compact Matern kernel
%
% Inputs: Marr - row vector of function indices (usually just 1:M)
%         x - column vector of values to evaluate the function at
%         L - length of the eigenfunction domain
%         deriv - <default=0> derivatives of the eigenfunctions
%                 this must be 0<=deriv<=4
%
% Note: if you pass x values beyond [0,L], this will error out
%
% Note: minimum M value is 1
%
% We may want to reconsider the L value - since we primarily work with L=1
% and nothing in the book uses L!=1

[Mr Mc] = size(Marr);
[n xc] = size(x);

if Mr~=1
    error('array of function indices Marr must be row vector, size(Marr)=[%d,%d]',Mr,Mc)
elseif xc~=1
    error('data locations x must be column vector, size(x)=[%d,%d]',n,xc)
end

if nargin<3
    error('Insufficient parameters passed')
elseif nargin==3
    deriv = 0; % Default to no derivatives being used
else
    [dr dc] = size(deriv);
    if dc>1 | dr>1
        error('deriv should be a scalar, size(deriv)=[%d,%d]',dr,dc)
    elseif deriv<0
        error('negative derivative requested, deriv=%g',deriv)
    elseif deriv>4 % Right now you only have 4 derivatives
        error('only up to 4 derivatives allowed, deriv=%g',deriv)
    elseif deriv~=floor(real(deriv))
        error('deriv must be an integer, deriv=%g',deriv)
    end
end

if any(Marr<0)
    error('Unacceptable eigenfunction index, min(Marr)=%g',min(Marr))
elseif any(Marr~=floor(real(Marr)))
    error('Marr must be integers')
end

if L~=real(abs(L)) | L==0
    error('L=%g unacceptable; must be real and positive',L)
elseif min(x)<0
    error('Values for MaternQR must be in [0,L], min(x)=%g',min(x))
elseif max(x)>L
    error('Values for MaternQR must be in [0,L], max(x)=%g',max(x))
end

switch deriv
    case 0
        p = sqrt(2/L)*sin(pi*x*Marr/L);
    case 1
        p = (sqrt(2/L)*ones(size(x))*(pi*Marr/L)).*cos(pi*x*Marr/L);
    case 2
        p = (-sqrt(2/L)*ones(size(x))*(pi*Marr/L).^2).*sin(pi*x*Marr/L);
    case 3
        p = (-sqrt(2/L)*ones(size(x))*(pi*Marr/L).^3).*cos(pi*x*Marr/L);
    case 4
        p = (sqrt(2/L)*ones(size(x))*(pi*Marr/L).^4).*sin(pi*x*Marr/L);
    otherwise
        error('Unacceptable derivative called ... somehow ..., deriv=%d',deriv)
end
