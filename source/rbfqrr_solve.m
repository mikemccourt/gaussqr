function rbfqrOBJ = rbfqrr_solve(x,y,ep,a,M)
% function rbfqrOBJ = rbfqrr_solve(x,y,ep,a,M)
% This function accepts required inputs x, y and ep
% x is a Nxd vector of input data points of the form
%     x = [x1';x2';...;xN'], x1 is a column vector
% y is a Nx1 vector of data values at the x locations
% ep is the traditional RBF shape parameter
%
% There is also a recommended input a which is the
% global scale parameter which determines the
% orthogonality enjoyed by the eigenfunctions
% NOTE: if a==0, this program picks a for you.
%
% Optional input is M which is the max length of the
% regression.  If you don't pass M, this will choose an
% M which is kinda random
%
% Note that only positive values of ep and a are used
%
% What is returned from this function is rbfqrOBJ which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.
if nargin<3
    error('Insufficient inputs')
end
if not(exist('a'))
    a = 1;
end
N = size(y,1);
if not(exist('M'))
    M = .4*N;
end
if M>N
    warning(sprintf('Input length: %d>%d unacceptable',M,N))
    M = .4*N;
end
if size(x,1)~=size(y,1)
    error('Different numbers of inputs and outputs')
end

a = abs(a);
ep = abs(ep);
b = ep^2;

Marr = rbfformMarr(M);
if a==0
    [a,kappa] = fminbnd(@(a) rbfsvdF(Marr,x,ep,a,b,sqrt(a^2+2*a*b)),.01,10);
end
c = sqrt(a^2+2*a*b);
phiMat = rbfphi(Marr,x,ep,a);
beta = phiMat\y;

rbfqrOBJ.reg   = true;
rbfqrOBJ.ep    = ep;
rbfqrOBJ.a     = a;
rbfqrOBJ.b     = b;
rbfqrOBJ.c     = c;
rbfqrOBJ.N     = N;
rbfqrOBJ.beta  = beta;
rbfqrOBJ.kappa = kappa;
rbfqrOBJ.Marr  = Marr;