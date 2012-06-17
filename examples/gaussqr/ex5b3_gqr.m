% ex5b3.m
% Comparing various methods of 2-pt BVP solvers
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  GaussQR regression (M<N)

% The problem we are interested in solving is
%    Lap(u) = -36*J0(6r),  u=(-1,1)^2
%         u = J0(6r),      boundary
% where r = sqrt(x^2_y^2)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Use absolute error
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;

f = @(x) -100*besselj(0,10*sqrt(x(:,1).^2+x(:,2).^2));
exact = @(x) besselj(0,10*sqrt(x(:,1).^2+x(:,2).^2));
f = @(x) -36*besselj(0,6*sqrt(x(:,1).^2+x(:,2).^2));
exact = @(x) besselj(0,6*sqrt(x(:,1).^2+x(:,2).^2));
% r24 = @(x) sqrt(x(:,1).^2+x(:,2).^4+eps);
% f = @(x) -36*(x(:,1).^2+4*x(:,2).^6)./r24(x).^2.*besselj(0,6*r24(x)) - 6*(6*x(:,1).^2.*x(:,2).^2-x(:,1).^2-2*x(:,2).^6+x(:,2).^4)./r24(x).^3.*besselj(1,6*r24(x));
% exact = @(x) besselj(0,6*r24(x));

rbf = @(e,r) exp(-(e*r).^2);
Lrbf = @(e,r) 4*e^2*((e*r).^2-1).*exp(-(e*r).^2);

N = 17;
x = pick2Dpoints([-1,-1],[1,1],N,'cheb');
u = exact(x);
NN = 30;
xx = pick2Dpoints([-1,-1],[1,1],NN);
uu = exact(xx);

% Identify the boundary points and interior points
ib = find((abs(x(:,1))==1) + (abs(x(:,2))==1));
xb = x(ib,:);
ii = setdiff(1:size(x,1),ib);
xi = x(ii,:);

alpha = 1;
epvec = logspace(-2,1,25);

warning off MATLAB:nearlySingularMatrix

errvecR2D = [];
errvecI2D = [];
m = 1;
for ep=epvec
    GQR = gqr_rsolve(x,u,ep,alpha);
    [uI,GQR] = gqr_eval(GQR,xx); % Store phi_eval to save flops later
    errvecI2D(m) = errcompute(uI,uu);
    
    [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
    phiMat = gqr_phi(Marr,xb,ep,alpha);
    phiMat2d = gqr_phi(Marr,xi,ep,alpha,[2,0])+gqr_phi(Marr,xi,ep,alpha,[0,2]);
    A = [phiMat;phiMat2d];
    rhs = [exact(xb);f(xi)];
    GQR.coef = A\rhs;
    errvecR2D(m) = errcompute(gqr_eval(GQR,xx),uu);

    fprintf('%d\n',m)
    m = m+1;
end

warning on MATLAB:nearlySingularMatrix

% Plot the results
loglog(epvec,errvecI2D,'r','LineWidth',3),hold on
loglog(epvec,errvecR2D,'b','LineWidth',3),hold off
xlabel('\epsilon')
ylabel('absolute 2-norm error')
% legend(sprintf('GaussQRr (M=%d)',size(GQR.Marr,2)),'Direct Gauss Collocation','Polynomial Collocation','Location','NorthEast')