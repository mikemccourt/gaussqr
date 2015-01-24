% HermiteConvergence
% This example presents a study of global differentiation matrices in a 2D
% domain and their convergence behavior.  We sample data at the 2D
% Chebyshev tensor points and evaluate on a uniform grid.

% Choose to use the 2-norm relative error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Define the some kernels & Laplacians
% Note that we only allow for isotropic kernels here
rbfM4 = @(e,r) (1+(e*r)+(e*r).^2/3).*exp(-e*r);
rbfM4L = @(e,r) 1/3*e^2*exp(-e*r).*((e*r).^2-2*e*r-2);
rbfM6 = @(e,r) (1+(e*r)+2/5*(e*r).^2+1/15*(e*r).^3).*exp(-e*r);
rbfM6L = @(e,r) 1/15*e^2*exp(-e*r).*(-6-6*(e*r)-(e*r).^2+(e*r).^3);

% Choose the kernel and shape parameter
rbf = rbfM4;
rbfL = rbfM4L;
ep = 4;

% Define a function and its Laplacian
f = @(x,y) (x.^2+y.^4).^(7/2) + x.*y;
fL = @(x,y) 7*(x.^2+y.^4).^(3/2).* ...
           (6*x.^2.*y.^2+6*x.^2+26*y.^6+y.^4);
% This is an analytic function if that is preferred
% f = @(x,y) exp(x.*y);
% fL = @(x,y) (x.^2+y.^2).*exp(x.*y);

% Define some points for testing
xeval = pick2Dpoints([-1,-1],[1 1],24);

% Evaluate at the chosen locations
ueval = f(xeval(:,1),xeval(:,2));
uLeval = fL(xeval(:,1),xeval(:,2));

% Define a range of N values to test convergence over
Nvec = [6,10,17,30,50,100];
hvec = zeros(size(Nvec));

% Loop over the required N values and compute error
errvec = zeros(size(Nvec));
errvecL = zeros(size(Nvec));
k = 1;
for N = Nvec
    % Create the data with the right amount of points
    x = pick2Dpoints([-1,-1],[1 1],N,'cheb');
    u = f(x(:,1),x(:,2));
    
    % Interpolate and evaluate the interpolant
    DM = DistanceMatrix(x,x);
    DMeval = DistanceMatrix(xeval,x);
    K = rbf(ep,DM);
    Keval = rbf(ep,DMeval);
    KLeval = rbfL(ep,DMeval);
    intcoef = K\u;
    seval = Keval*intcoef;
    sLeval = KLeval*intcoef;
    
    % Record the errors and the fill distance
    errvec(k) = errcompute(seval,ueval);
    errvecL(k) = errcompute(sLeval,uLeval);
    hvec(k) = max(min(DM+2*eye(N^2)));
    k = k + 1;
end

% Compute the lines of best fit
hvplot = logspace(-2,0,3);
p = polyfit(log10(hvec),log10(errvec),1);
pL = polyfit(log10(hvec),log10(errvecL),1);

% Plot the results for interpolation and derivatives
h = figure;
loglog(hvec,errvec,'-^','linewidth',2)
hold on
loglog(hvec,errvecL,'linewidth',2)
loglog(hvplot,10^p(2)*hvplot.^p(1),'--^k')
loglog(hvplot,10^pL(2)*hvplot.^pL(1),'--k')
hold off
xlabel('number of data locations')
ylabel('relative RMS 2-norm error')
legend('Interpolant','Laplacian',...
       sprintf('O(h^{%2.1f})',p(1)),sprintf('O(h^{%2.1f})',pL(1)),...
       'location','southeast')