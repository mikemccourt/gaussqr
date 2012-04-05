clear all;
clc

format short e

N = 40;
M = 160;
sigma = 1;
L = 1;

spaceopt = 'cheb';
fopt = 'sinh';

[yf,fstr] = pickfunc(fopt,1);

aa = 0; bb = L;

x = pickpoints(aa,bb,N,spaceopt);

n =1:M;
S = sqrt(2/L)*sin(pi*x*n/L); % size(S)
lambda = 1./((pi*n/L).^2+sigma^2);
D = diag(lambda); % size(D)
% S*D*S'

% KERNEL:
MINVAL = min(repmat(x,1,N),repmat(x',N,1));
MAXVAL = max(repmat(x,1,N),repmat(x',N,1));
K = sinh(sigma*MINVAL).*sinh(sigma*(L-MAXVAL))/(sigma*sinh(L*sigma));

error = norm(K-S*D*S')