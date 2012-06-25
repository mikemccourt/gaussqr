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
rhs = [exact(t([1,end]));f(t(2:end-1))];
k = 1;
for ep = epvec
    GQR = gqr_solveprep(0,t,ep,alpha);
    phi = gqr_phi(GQR,t([1,end]));
    phi2 = gqr_phi(GQR,t(2:end-1),2);
    A = [phi;phi2]*[eye(N);GQR.Rbar];

    coef = A\rhs;
    GQR.coef  = coef;

    err_GQR(k) = errcompute(gqr_eval(GQR,xx),exact(xx));
    k = k+1;
end

% GaussQR regression
err_GQRr = [];
alpha = 1;
rhs = [exact(t([1,end]));f(t(2:end-1))];
k = 1;
for ep=epvec
    GQR = gqr_solveprep(1,t,ep,alpha);
    phi = gqr_phi(GQR,t([1,end]));
    phi2 = gqr_phi(GQR,t(2:end-1),2);
    A = [phi;phi2];

    coef = A\rhs;
    GQR.coef  = coef;

    err_GQRr(k) = errcompute(gqr_eval(GQR,xx),exact(xx));
    k = k+1;
end

warning on MATLAB:nearlySingularMatrix

% Plot the results
loglog(epvec,err_GQR,'r','LineWidth',3),hold on
loglog(epvec,err_Kansa,'b','LineWidth',2)
%loglog(epvec,err_GQRr,'g','LineWidth',3)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('absolute 2-norm error')
title(sprintf('Collocation for a 2-pt BVP, N=%d',N))
% legend('Kansa','GaussQR','GaussQRr (M=.5N)','Direct','Location','East')
legend('True Gauss Solution','Direct Gauss Collocation','Polynomial Collocation','Location','SouthEast')
