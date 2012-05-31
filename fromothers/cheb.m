% function [D,t] = cheb(N,a,b)
% This creates a 1D differentiation matrix, and returns also
% the Chebyshev nodes associate with it.
% Input : N - Number of points requested
%         a,b - <default=[-1,1]> Bounds of [a,b] domain
% Output : D - The differentiation matrix
%          t - The chebyshev nodes for the differentiation
%
% Borrowed from Trefethen's book: Spectral Methods in Matlab
% Modified to allow general domains, not just [-1,1]
function [D,t] = cheb(N,a,b)

if N~=abs(real(N)) | N==0
    error('N must be a positive integer, N=%g',N)
end

switch nargin
    case 1
        a = -1;
        b = 1;
    case 2
        error('Must pass both upper and lower bounds')
    case 3
        if a~=real(a) | b~=real(b)
            error('Cannot pass complex domains')
        elseif a>=b
            error('Must have a<b, a=%g, b=%g',a,b)
        end
    otherwise
        error('Unacceptable parameters have been passed')
end

t = pickpoints(a,b,N,'cheb');

c = [2; ones(N-2,1); 2].*(-1).^(0:N-1)';
T = repmat(t,1,N);
dT = T-T';
D  = (c*(1./c)')./(dT+(eye(N)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
