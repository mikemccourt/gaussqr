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
errvecINT = zeros(size(epvec));
alpha = 2;

ie = 1;
for ep=epvec
    nu = (2*ep/alpha)^2;
    lam = nu/(2+nu+2*sqrt(1+nu));
    M = ceil(N+log(eps)/log(lam));
    if Mextramax~=0
        M = min(M,Mextramax);
    end
    Marr = rbfformMarr(M)+1;
    phiMat = rbfphialpha(Marr,x,ep,alpha);
    phiMatD2 = rbfphialpha(Marr,x(2:end-1),ep,alpha,2); % Only interior derivatives needed

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
    
    yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
    errvec(ie) = norm((yy-yp)./(abs(yy)+eps));
    
    rbfqrOBJ = rbfqr_solve_alpha(x,y,ep,alpha);
    yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
    errvecINT(ie) = norm((yy-yp)./(abs(yy)+eps));
    
    ie = ie + 1;
end

clf reset
loglog(epvec,[errvec;errvecINT])
title(sprintf('u_{xx}=cosh(x), Dirichlet BC, N=%d, alpha=%g',N,alpha))
ylabel('Error: ||(y-yans)/yans||')
xlabel('\epsilon')
legend('Collocation','Interpolation','Location','NorthWest')