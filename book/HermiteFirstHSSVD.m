% HermiteFirstHSSVD
% I don't think this example is going to stick around, but I need to test
% some things

% Choose to use the absolute max-norm error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Choose a simple functions to approximation
uf = @(x) cos(2*pi*x);
uxf = @(x) -2*pi*sin(2*pi*x);

% Choose a small number of points to work on
N = 30;
x = pickpoints(-1,1,N);
xeval = pickpoints(-1,1,500);

% Create data values and test values
u = uf(x);
ueval = uf(xeval);
uxeval = uxf(xeval);

% Define the HS-SVD basis
ep = 1;alpha = 1;
GQR = gqr_solve(x,u,ep,alpha);

% Evaluate the HS-SVD basis solution
uhatHSeval = gqr_eval(GQR,xeval);
uhatHSxeval = gqr_eval(GQR,xeval,1);

% Create the standard Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);
rbfx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);

% Solve for the interpolant in the standard basis
coef = rbf(ep,DistanceMatrix(x,x))\u;
uhateval = rbf(ep,DistanceMatrix(xeval,x))*coef;
uhatxeval = rbfx(ep,DistanceMatrix(xeval,x),DifferenceMatrix(xeval,x))*coef;

% Plot the results of the two solution methods
subplot(2,2,1)
plot(xeval,abs(uhateval-ueval),'linewidth',2)
ylabel('$|\hat{u}-u|$','interpreter','latex')
title('Standard basis, function values')
subplot(2,2,2)
plot(xeval,abs(uhatxeval-uxeval),'linewidth',2)
ylabel('$|\hat{u}-u|$','interpreter','latex')
title('Standard basis, derivative values')
subplot(2,2,3)
plot(xeval,abs(uhatHSeval-ueval),'linewidth',2)
ylabel('$|\hat{u}-u|$','interpreter','latex')
title('HS-SVD basis, function values')
subplot(2,2,4)
plot(xeval,abs(uhatHSxeval-uxeval),'linewidth',2)
ylabel('$|\hat{u}-u|$','interpreter','latex')
title('HS-SVD basis, derivative values')