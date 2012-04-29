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
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mfactor = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;

epvec = logspace(-1,1,30);
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
    Marr = gqr_formMarr(M);
    phiMat = gqr_phi(Marr,x,ep,alpha);
    phiMatD2 = gqr_phi(Marr,x(2:end-1),ep,alpha,2); % Only interior derivatives needed

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
    
    GQR.reg   = false;
    GQR.ep    = ep;
    GQR.alpha = alpha;
    GQR.N     = N;
    GQR.coef  = coef;
    GQR.Rbar  = Rbar;
    GQR.Marr  = Marr;
    
    yp = gqr_eval(GQR,xx);
    errvec(ie) = errcompute(yp,yy);
    
    M = Mfactor*N;
    Marr = gqr_formMarr(M);
    phiMat = gqr_phi(Marr,x,ep,alpha);
    phiMatD2 = gqr_phi(Marr,x(2:end-1),ep,alpha,2); % Only interior derivatives needed
    
    A = [phiMat(1,:);phiMatD2;phiMat(end,:)];
    rhs = [1;cosh(x(2:end-1));cosh(1)];
    warning off MATLAB:nearlySingularMatrix
    coef = A\rhs;
    warning on MATLAB:nearlySingularMatrix
    
    GQR.reg   = true;
    GQR.ep    = ep;
    GQR.alpha = alpha;
    GQR.N     = N;
    GQR.coef  = coef;
    GQR.Marr  = Marr;
    
    yp = gqr_eval(GQR,xx);
    errvecREG(ie) = errcompute(yp,yy);
    
    ie = ie + 1;
end

clf reset
loglog(epvec,[errvec;errvecREG])
title(sprintf('Collocation for u_{xx}=cosh(x), Dirichlet BC, N=%d',N))
ylabel('Relative Error')
xlabel('\epsilon')
legend('RBF-QR','RBF-QRr','Location','NorthWest')
