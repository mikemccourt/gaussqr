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
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mfactor = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;

epvec = logspace(-1,1,50);
N = 25;
NN = 200;

spaceopt = 'cheb';
[x,spacestr] = pickpoints(0,1,N,spaceopt);

truesol = @(x)cosh(x);
y = truesol(x);
xx = pickpoints(0,1,NN);
yy = truesol(xx);

errvec = zeros(size(epvec));
errvecREG = zeros(size(epvec));
errvecINT = zeros(size(epvec));
alpha = 2;

ie = 1;
for ep=epvec
    nu = (2*ep/alpha)^2;
    lam = nu/(2+nu+2*sqrt(1+nu));
    M = ceil(N+log(eps)/log(lam));
    if Mextramax~=0
        M = min(M,abs(Mextramax));
    end
    Marr = rbfformMarr(M)+1;
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMatD2 = rbfphi(Marr,x(2:end-1),ep,alpha,2); % Only interior derivatives needed

    [Q,R] = qr(phiMat);
    R1 = R(:,1:N);
    R2 = R(:,N+1:end);
    iRdiag = diag(1./diag(R1));
    R1s = iRdiag*R1;
    opts.UT = true;
    Rhat = linsolve(R1s,iRdiag*R2,opts);
    D = lam.^(toeplitz(sum(Marr(N+1:end),1),sum(Marr(N+1:-1:2),1)));
    Rbar = D.*Rhat';
    
    A = [phiMat(1,:);phiMatD2;phiMat(end,:)]*[eye(N);Rbar];
    rhs = [1;cosh(x(2:end-1));cosh(1)];
    warning off MATLAB:nearlySingularMatrix
    coef = A\rhs;
    warning on MATLAB:nearlySingularMatrix
    
    rbfqrOBJ.reg   = false;
    rbfqrOBJ.ep    = ep;
    rbfqrOBJ.alpha = alpha;
    rbfqrOBJ.N     = N;
    rbfqrOBJ.coef  = coef;
    rbfqrOBJ.Rbar  = Rbar;
    rbfqrOBJ.Marr  = Marr;
    
    yp = rbfqr_eval(rbfqrOBJ,xx);
    errvec(ie) = norm((yy-yp)./(abs(yy)+eps));
    
    M = Mfactor*N;
    Marr = rbfformMarr(M)+1;
    phiMat = rbfphi(Marr,x,ep,alpha);
    phiMatD2 = rbfphi(Marr,x(2:end-1),ep,alpha,2); % Only interior derivatives needed
    
    A = [phiMat(1,:);phiMatD2;phiMat(end,:)];
    rhs = [1;cosh(x(2:end-1));cosh(1)];
    warning off MATLAB:nearlySingularMatrix
    coef = A\rhs;
    warning on MATLAB:nearlySingularMatrix
    
    rbfqrOBJ.reg   = true;
    rbfqrOBJ.ep    = ep;
    rbfqrOBJ.alpha = alpha;
    rbfqrOBJ.N     = N;
    rbfqrOBJ.coef  = coef;
    rbfqrOBJ.Marr  = Marr;
    
    yp = rbfqr_eval(rbfqrOBJ,xx);
    errvecREG(ie) = norm((yy-yp)./(abs(yy)+eps));
    
    ie = ie + 1;
end

clf reset
loglog(epvec,[errvec;errvecREG])
title(sprintf('Collocation for u_{xx}=cosh(x), Dirichlet BC, N=%d',N))
ylabel('Relative Error')
xlabel('\epsilon')
legend('RBF-QR','RBF-QRr','Location','NorthWest')