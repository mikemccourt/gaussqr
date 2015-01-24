% HermiteEpProfile
% This script studies the effect of epsilon on the quality of approximating
% a derivative.  It studies a Gaussian approximation to the mixed
% derivative of a function and shows that the stable basis can be a useful
% mechanism for approximating derivatives with small epsilon

% Choose to use the absolute max-norm error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 200;

% Choose the Gaussian RBFs for HS-SVD comparison
rbf = @(e,r) exp(-(e*r).^2);
rbfxy = @(e,r,dx,dy) 4*e^4*dx.*dy.*exp(-(e*r).^2);

% Choose a range of shape parameters to consider
epvec = logspace(-1,1,30);

% Choose a function for testing
uf = @(x,y) cos(pi*(x.^2+y));
ufx = @(x,y) -2*pi*x.*sin(pi*(x.^2+y));
ufxy = @(x,y) -2*pi^2*x.*cos(pi*(x.^2+y));

% Create the data and evaluation points
N = 8;Neval = 16;
x = pick2Dpoints([-1 -1],[1 1],N,'halt');
xeval = pick2Dpoints([-1 -1],[1 1],Neval);

% Evaluate the distance and difference matrices
DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);
DiffMxeval = DifferenceMatrix(xeval(:,1),x(:,1));
DiffMyeval = DifferenceMatrix(xeval(:,2),x(:,2));

% Evaluate the functions
u = uf(x(:,1),x(:,2));
ueval = uf(xeval(:,1),xeval(:,2));
uxeval = ufx(xeval(:,1),xeval(:,2));
uxyeval = ufxy(xeval(:,1),xeval(:,2));

% Define a function to compute the standard basis derivatives
uhatxyeval = @(ep) rbfxy(ep,DMeval,DiffMxeval,DiffMyeval)* ...
                   (rbf(ep,DM)\u);
               
% Define a function to compute the HS-SVD derivatives
uhatHSxyeval = @(ep) gqr_eval(gqr_solve(x,u,ep,1),xeval,[1 1]);
uhatHSxeval = @(ep) gqr_eval(gqr_solve(x,u,ep,1),xeval,[1 0]);

% Compute all the errors at once
warning('off','MATLAB:nearlySingularMatrix')
errvec = arrayfun(@(ep)errcompute(uhatxyeval(ep),uxyeval),epvec);
warning('on','MATLAB:nearlySingularMatrix')
% errvecHS = arrayfun(@(ep)errcompute(uhatHSxyeval(ep),uxyeval),epvec);
errvecHS = arrayfun(@(ep)errcompute(uhatHSxeval(ep),uxeval),epvec);
% errvec = arrayfun(@(ep)errcompute(rbf(ep,DMeval)*(rbf(ep,DM)\u),ueval),epvec);
% errvecHS = arrayfun(@(ep)errcompute(gqr_eval(gqr_solve(x,u,ep,1),xeval),ueval),epvec);

% Plot the errors together
h = figure;
loglog(epvec,errvec,'--','linewidth',2)
hold on
loglog(epvec,errvecHS,'linewidth',2)
hold off