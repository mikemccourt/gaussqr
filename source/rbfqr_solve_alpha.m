function rbfqrOBJ = rbfqr_solve(x,y,ep,alpha,M)
% function rbfqrOBJ = rbfqr_solve(x,y,ep,alpha,M)
% This function accepts required inputs x, y and ep
% x is a Nxd vector of input data points of the form
%     x = [x1';x2';...;xN'], x1 is a column vector
% y is a Nx1 vector of data values at the x locations
% ep is the traditional RBF shape parameter
%
% There is also a recommended input alpha which is the
% global scale parameter which determines the
% orthogonality enjoyed by the eigenfunctions
% The default is set in rbfsetup, which can be used by passing []
% For the time being alpha is only a single value
% In the future we may allow for different alpha in each dimension
%
% Note that only positive values of ep and a are used
%
% Optional input is M>N which is the length of the
% expansion.  If you don't pass M, this will choose an
% M based on the eigenvalues
%
% What is returned from this function is rbfqrOBJ which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
alphaDefault = GAUSSQR_PARAMETERS.ALPHA_DEFAULT;

if nargin<3
    error('Insufficient inputs')
end

if sum(size(x)~=size(y))
    error('Different sized x (input) and y (output) vectors')
end

if not(exist('alpha'))
    alpha = alphaDefault;
elseif length(alpha)==0
    alpha = alphaDefault;
end

a = abs(a);
ep = abs(ep);

N = size(y,1);
nu = (2*ep/alpha)^2;
lam = nu/(2+nu+2*sqrt(1+nu));

if not(exist('M'))
    M = ceil(min(N+8/max(.1,-log10(ep)),N+30));
else
    if length(M)==0
        error('Empty M passed')
    elseif M<N
        warning(sprintf('Input length: %d<%d unacceptable',M,N))
        M = ceil(min(N+8/max(.1,-log10(ep)),N+20));
    end
end

b = ep^2;
c = sqrt(a^2+2*a*b);

Marr = rbfformMarr(M+1);
phiMat = rbfphi(Marr,x,ep,a);

[Q,R] = qr(phiMat);
R1 = R(:,1:N);
R2 = R(:,N+1:end);
iRdiag = diag(1./diag(R1));
R1s = iRdiag*R1;
opts.UT = true;
Rhat = linsolve(R1s,iRdiag*R2,opts);

lam = b/(a+b+c);
D = lam.^(toeplitz(sum(Marr(N+1:end),1),sum(Marr(N+1:-1:2),1)));
Rbar = D.*Rhat';
beta = ranksolve(Rhat,Rbar,linsolve(R1s,iRdiag*(Q'*y),opts));

rbfqrOBJ.reg  = false;
rbfqrOBJ.ep   = ep;
rbfqrOBJ.a    = a;
rbfqrOBJ.b    = b;
rbfqrOBJ.c    = c;
rbfqrOBJ.N    = N;
rbfqrOBJ.beta = beta;
rbfqrOBJ.Rbar = Rbar;
rbfqrOBJ.Marr = Marr;