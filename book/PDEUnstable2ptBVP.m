% PDEUnstable2ptBVP
% This example shows that standard collocation schemes can be effective but
% are susceptible to the same ill-conditioning as interpolation problems
% The problem of interest here is
%             u'' = -sin(x)
%         u(0) = 0, u(pi) = 0
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define the inverse multiquadric and its 2nd derivative
rbf = @(e,r) 1./sqrt(1+(e*r).^2);
rbfxx = @(e,r) e^2*(2*(e*r).^2-1)./(1+(e*r).^2).^(5/2);
ep = .045;

% Define functions for the interior and boundary
fint = @(x) -sin(x);
fbc = @(x) zeros(size(x,1),1);
uf = @(x) sin(x);

% Choose a set of collocation points
% Separate them into a interior and boundary region
N = 25;
xall = pickpoints(0,pi,N);
xint = xall(xall~=0 & xall~=pi);
xbc = xall(xall==0 | xall==pi);
x = [xint;xbc];

% Choose some evaluation points
Neval = 200;
xeval = pickpoints(0,pi,Neval);

% Define the necessary distance matrices
DMint = DistanceMatrix(xint,x);
DMbc = DistanceMatrix(xbc,x);
DMeval = DistanceMatrix(xeval,x);

% Form the right hand side, in the same order as the points
rhsint = fint(xint);
rhsbc = fbc(xbc);
rhs = [rhsint;rhsbc];

% Form the kernel matrices, both collocation and evaluation
Kxxint = rbfxx(ep,DMint);
Kbc = rbf(ep,DMbc);
A = [Kxxint;Kbc];
Keval = rbf(ep,DMeval);

% Find the solution coefficients and evaluate the solution
coef = A\rhs;
ucoll = Keval*coef;
plot(xeval,ucoll,'b','linewidth',2)
xlim([0,pi])
ylim([0,1.2])