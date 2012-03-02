% CHEB  compute D = differentiation matrix, t = Chebyshev grid
% Borrowed from Trefethen's book: Spectral Methods in Matlab

% This creates a 1D differentiation matrix, and returns also
% the Chebyshev nodes associate with it.

  function [D,t] = cheb(N)
  if N==0, D=0; t=1; return, end
  t = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  T = repmat(t,1,N+1);
  dT = T-T';                  
  D  = (c*(1./c)')./(dT+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
