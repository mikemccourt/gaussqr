% HermiteNonsymmetric
% This example demonstrates nonsymmetric Hermite interpolation in Matlab

global GAUSSQR_PARAMETERS;
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Define data locations and evaluation locations
x = pickpoints(0,1,10,'rand');
xeval = pickpoints(0,1,500);

% Choose a function and evaluate it at those locations
uf = @(x) cos(2*pi*x);
uxf = @(x) -2*pi*sin(2*pi*x);
ux = uxf(x);
ueval = uf(xeval);

% Consider the C2 Matern kernel
rbf = @(e,r) (1+(e*r)).*exp(-(e*r));
rbfx = @(e,r,dx) -e^2*dx.*exp(-(e*r));
ep = 4;

% Form the necessary distance and difference matrices
DM = DistanceMatrix(x,x);
DiffM = DifferenceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);

% Compute the Hermite interpolation matrix and solve
Kx = rbfx(ep,DM,DiffM);
c = Kx\ux;

% Evaluate the interpolant at the desired locations
Keval = rbf(ep,DMeval);
uhateval = Keval*c;

% Create a function to evaluate the interpolation for any scalar epsilon
% and at any locations
uhatepf = @(ep,xeval) rbf(ep,DistanceMatrix(xeval,x))*(rbfx(ep,DM,DiffM)\ux);

% Use this function to find the epsilon so that u(0)=1
fzero(@(ep)uhatepf(ep,0)-1,1)