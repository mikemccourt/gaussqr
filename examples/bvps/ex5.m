% ex5b
% Comparing various methods of 2-pt BVP solvers
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  HS-SVD basis
% The problem we are interested in solving is
%    u_xx = -sinh(x)/(1+cosh(x))^2,  u(-1) = u(1) = 0

f = @(x) -sinh(x)./(1+cosh(x)).^2;
exact = @(x) sinh(x)./(1+cosh(x));

N = 16;
NN = 200;
xx = pickpoints(-1,1,NN);
uu = exact(xx);
epvec = logspace(-1,1,41);

% Define the Gaussian kernel and its second derivative
rbf = @(e,r) exp(-(e*r).^2);
d2rbf = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

% Trefethen's method for polynomial collocation
[D,t] = cheb(N);
D2 = D^2;
D2([1,end],:) = [1,zeros(1,N-1);zeros(1,N-1),1];
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
u = D2\rhs;
[pcoef,S,mu] = polyfit(t,u,N-1);
err_Trefethen = errcompute(polyval(pcoef,xx,S,mu),uu);

% Kansa's method
err_Kansa = zeros(size(epvec));
r = DistanceMatrix(t,t);
reval = DistanceMatrix(xx,t);
rhs = [exact(t(1));f(t(2:end-1));exact(t(end))];
k = 1;
for ep=epvec
  A = rbf(ep,r);
  D2 = d2rbf(ep,r);
  Aeval = rbf(ep,reval);
  D2([1,end],:) = A([1,end],:);
  warning('off','MATLAB:nearlySingularMatrix')
  coef = D2\rhs;
  warning('on','MATLAB:nearlySingularMatrix')
  err_Kansa(k) = errcompute(Aeval*coef,uu);
  k = k+1;
end

% GaussQR interpolation
err_GQR = zeros(size(epvec));
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

    err_GQR(k) = errcompute(gqr_eval(GQR,xx),uu);
    k = k+1;
end

% Plot the results
loglog(epvec,err_GQR,'r','LineWidth',3),hold on
loglog(epvec,err_Kansa,'b','LineWidth',2)
loglog(epvec,err_Trefethen*ones(size(epvec)),'--k','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('absolute 2-norm error')
title(sprintf('Collocation for a 2-pt BVP, N=%d',N))
legend('True Gauss Solution','Direct Gauss Collocation','Polynomial Collocation','Location','North')
