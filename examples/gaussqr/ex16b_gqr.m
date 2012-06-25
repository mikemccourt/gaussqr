% ex16b_gqr.m
% Comparing various methods of 2D BVP solvers
% These methods all employ a tensor product structure
% differentiation matrix
%
% The methods we consider are
%  Trefethen's cheb matrix
%  Kansa's nonsymmetric
%  GaussQR collocation
%
% The problem we are interested in solving is the Helmholtz equation:
%    Lap(u)+k^2 u = f
%    (x,y)=[-1,1]^2     Dirichlet BC
%
% Unlike ex16_gqr, this considers a range of ep and alpha values and that
% allows you to see the effect of both on the accuracy

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

alphavec = logspace(-2,1,17);
epvec = logspace(-2,1,16);


errmatR2D = zeros(length(alphavec),length(epvec));
errmatI2D = zeros(length(alphavec),length(epvec));
m = 1;
for ep=epvec
    k = 1;
    for alpha=alphavec
        GQR = gqr_rsolve(x,u,ep,alpha);
        [uI,GQR] = gqr_eval(GQR,xx);
        errmatI2D(k,m) = errcompute(uI,uu);
    
        GQRr = gqr_solveprep(1,x,ep,alpha);
        phiMat = gqr_phi(GQRr,xb);
        phiMat2d = gqr_phi(GQRr,xi,[2,0])+gqr_phi(GQRr,xi,[0,2]);
        A = [phiMat;phiMat2d];
        rhs = [exact(xb);f(xi)];
warning off MATLAB:rankDeficientMatrix
        GQRr.coef = A\rhs;
warning on MATLAB:rankDeficientMatrix
        errmatR2D(k,m) = errcompute(gqr_eval(GQRr,xx),uu);
        
        fprintf('%d ',k)
        k = k+1;
    end
    fprintf('%d\n',m)
    m = m+1;
end


% Plot the results
[A,E] = meshgrid(alphavec,epvec);
% surf(A,E,log10(min(errmatI2D',1e5)));
surf(A,E,log10(min(errmatR2D',1e5)));
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('\alpha')
ylabel('\epsilon')