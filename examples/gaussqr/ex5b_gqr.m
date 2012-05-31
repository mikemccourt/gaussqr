% ex5b.m
% Comparing various methods of 2-pt BVP solvers
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  GaussQR interpolation
%  GaussQR regression

% The problem we are interested in solving is
%    u_xx = exp(4x)  u(-1)=u(1)=0

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error
Mextramax = GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC;
Mfactor = GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC;

N = 16;
NN = 200;
xx = pickpoints(-1,1,NN);
%f = @(x) exp(4*x);
%exact = @(x) 1/16*(exp(4*x)-sinh(4)*x-cosh(4));
%f = @(x) cosh(x).*x+2*sinh(x);
%exact = @(x) cosh(x).*x-x-1;
f = @(x) -sinh(x)./(1+cosh(x)).^2;
exact = @(x) sinh(x)./(1+cosh(x));
rbf = @(e,r) exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

warning off MATLAB:nearlySingularMatrix

% Trefethen's method
[D,t] = cheb(N);
D2 = D^2;
D2([1,end],:) = [1,zeros(1,N-1);zeros(1,N-1),1];
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
u = D2\rhs;
[pcoef,S,mu] = polyfit(t,u,N-1);
err_Trefethen = errcompute(polyval(pcoef,xx,S,mu),exact(xx))

% Kansa's method
epvec = logspace(-1,1,41);
k = 1;
err_Kansa = [];
r = DistanceMatrix(t,t);
reval = DistanceMatrix(xx,t);
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
for ep=epvec
  A = rbf(ep,r);
  D2 = d2rbf(ep,r);
  Aeval = rbf(ep,reval);
  D2([1,end],:) = A([1,end],:);
  coef = D2\rhs;
  err_Kansa(k) = errcompute(Aeval*coef,exact(xx));
  k = k+1;
end

% GaussQR interpolation
err_GQR = [];
alpha = 1;
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
k = 1;
for ep = epvec
  nu = (2*ep/alpha)^2;
  lam = nu/(2+nu+2*sqrt(1+nu));
  M = ceil(N+log(eps)/log(lam));
  if Mextramax~=0
      M = min(M,abs(Mextramax));
  end
  Marr = gqr_formMarr(M);
  phiMat = gqr_phi(Marr,t,ep,alpha);
  phiMatD2 = gqr_phi(Marr,t(2:end-1),ep,alpha,2);

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
  coef = A\rhs;

  GQR.reg   = false;
  GQR.ep    = ep;
  GQR.alpha = alpha;
  GQR.N     = N;
  GQR.coef  = coef;
  GQR.Rbar  = Rbar;
  GQR.Marr  = Marr;

  err_GQR(k) = errcompute(gqr_eval(GQR,xx),exact(xx));
  k = k+1;
end

% GaussQR regression
err_GQRr = [];
Mfactor = .5;
M = Mfactor*N;
Marr = gqr_formMarr(M);
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
k = 1;
for ep=epvec
%  alpha = gqr_alphasearch(ep,-1,1); % bounds of the problem: [-1,1]
  alpha = 1;
  phiMat = gqr_phi(Marr,t,ep,alpha);
  phiMatD2 = gqr_phi(Marr,t(2:end-1),ep,alpha,2);
  A = [phiMat(1,:);phiMatD2;phiMat(end,:)];
  coef = A\rhs;

  GQR.reg   = true;
  GQR.ep    = ep;
  GQR.alpha = alpha;
  GQR.N     = N;
  GQR.coef  = coef;
  GQR.Marr  = Marr;

  err_GQRr(k) = errcompute(gqr_eval(GQR,xx),exact(xx));
  k = k+1;
end

warning on MATLAB:nearlySingularMatrix

% Plot the results
loglog(epvec,err_Kansa,'b','LineWidth',3),hold on
loglog(epvec,err_GQR,'r','LineWidth',3)
loglog(epvec,err_GQRr,'g','LineWidth',3)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('absolute error')
title(sprintf('Collocation for a 2-pt BVP, N=%d',N))
legend('Kansa','GaussQR','GaussQRr (M=.5N)','Direct','Location','East')
