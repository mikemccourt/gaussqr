% ex5
% This is an example showing how to use rbfqr for BVP
% This is a 1D example:
%   u''(x) = cosh(x)
%   u(0) = 1
%   u(1) = cosh(1)
% The answer is u(x)=cosh(x) so we can test the method
%
% Eventually I'll try to encapsulate some of this into functions, but at
% the moment I don't have a smart way to work with general linear
% operators.  I'll think on it.
rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end

% This might get carried over from other examples
clear opts

epvec = logspace(-1,1,30);
N = 25;
NN = 200;

spaceopt = 'cheb';
[x,spacestr] = pickpoints(0,1,N,spaceopt);

truesol = @(x)cosh(x);
f = @(x)cosh(x);
y = truesol(x);
xx = pickpoints(0,1,NN);
yy = truesol(xx);

errvec = zeros(size(epvec));
errvecREG = zeros(size(epvec));
errvecINT = zeros(size(epvec));
alpha = 2;

ie = 1;
for ep=epvec
    % Compute the 2-pt BVP solution
    GQR = gqr_solveprep(0,x,ep,alpha);
    phiMat = gqr_phi(GQR,x([1,end]));
    phiMatD2 = gqr_phi(GQR,x([2:end-1]),2);
    A = [phiMat;phiMatD2]*[eye(N);GQR.Rbar];
    rhs = [truesol(x([1,end]));f(x([2:end-1]))];
    
    warning off MATLAB:nearlySingularMatrix
    coef = A\rhs;
    warning on MATLAB:nearlySingularMatrix
    GQR.coef  = coef;
    
    yp = gqr_eval(GQR,xx);
    errvec(ie) = errcompute(yp,yy);
    
    % Consider the regression solution for comparison
    GQRreg = gqr_solveprep(1,x,ep,alpha);
    phiMat = gqr_phi(GQRreg,x([1,end]));
    phiMatD2 = gqr_phi(GQRreg,x(2:end-1),2);
    A = [phiMat;phiMatD2];
    rhs = [truesol(x([1,end]));f(x([2:end-1]))];
    
    warning off MATLAB:nearlySingularMatrix
    coef = A\rhs;
    warning on MATLAB:nearlySingularMatrix
    GQRreg.coef  = coef;
    
    yp = gqr_eval(GQRreg,xx);
    errvecREG(ie) = errcompute(yp,yy);
    
    ie = ie + 1;
end

clf reset
loglog(epvec,[errvec;errvecREG],'linewidth',3)
title(sprintf('Collocation for u_{xx}=cosh(x), Dirichlet BC, N=%d',N))
ylabel('Relative Error')
xlabel('\epsilon')
legend('RBF-QR','RBF-QRr','Location','NorthWest')
