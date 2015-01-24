% HermiteEpProfile
% This script studies the effect of epsilon on the quality of approximating
% a derivative.  It studies a Gaussian approximation to the mixed
% derivative of a function and shows that the stable basis can be a useful
% mechanism for approximating derivatives with small epsilon

% Choose to use the absolute max-norm error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 100;

% Choose the Gaussian RBFs for HS-SVD comparison
rbf = @(e,r) exp(-(e*r).^2);
rbfx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
rbfxx = @(e,r) 2*e^2*exp(-(e*r).^2).*(2*(e*r).^2-1);

% Choose a range of shape parameters to consider
epvec = logspace(-1,1,30);

% Choose a function for testing
uf = @(x) 1./(1+x.^2);
ufx = @(x) -2*x./(1+x.^2).^2;
ufxx = @(x) (6*x.^2-2)./(1+x.^2).^3;

% Create the data and evaluation points
N = 30;
Neval = 50;
x = pickpoints(-1,1,N,'halt');
xeval = pickpoints(-1,1,Neval);

% Evaluate the distance and difference matrices
DM = DistanceMatrix(x,x);
DMeval = DistanceMatrix(xeval,x);
DiffMxeval = DifferenceMatrix(xeval,x);

% Evaluate the functions
u = uf(x);
ueval = uf(xeval);
uxeval = ufx(xeval);
uxxeval = ufxx(xeval);

% Initialize the error containers
errvec = zeros(size(epvec));errvecx = errvec;errvecxx = errvec;
errvecHS = errvec;errvecHSx = errvec;errvecHSxx = errvec;

% Loop through the epsilon values
k = 1;
for ep=epvec
    % Solve and evaluate with the standard basis
    c = rbf(ep,DM)\u;
    uhateval = rbf(ep,DMeval)*c;
    uhatxeval = rbfx(ep,DMeval,DiffMxeval)*c;
    uhatxxeval = rbfxx(ep,DMeval)*c;
    
    % Solve and evaluate with the HS-SVD basis
    gqr_alpha = 3;
    GQR = gqr_solve(x,u,ep,gqr_alpha);
    uhatHSeval = gqr_eval(GQR,xeval);
    uhatHSxeval = gqr_eval(GQR,xeval,1);
    uhatHSxxeval = gqr_eval(GQR,xeval,2);
    
    % Compute the errors and store them
    errvec(k) = errcompute(ueval,uhateval);
    errvecx(k) = errcompute(uxeval,uhatxeval);
    errvecxx(k) = errcompute(uxxeval,uhatxxeval);
    errvecHS(k) = errcompute(ueval,uhatHSeval);
    errvecHSx(k) = errcompute(uxeval,uhatHSxeval);
    errvecHSxx(k) = errcompute(uxxeval,uhatHSxxeval);
    k = k + 1;
end

% Plot the errors together
h_plot = figure;
h = loglog(epvec,errvec,'b--','linewidth',1);
hold on
hHS = loglog(epvec,errvecHS,'b','linewidth',2);
hx = loglog(epvec,errvecx,'r--^','linewidth',1);
hHSx = loglog(epvec,errvecHSx,'r-^','linewidth',2);
hxx = loglog(epvec,errvecxx,'k--+','linewidth',1);
hHSxx = loglog(epvec,errvecHSxx,'k-+','linewidth',2);
hold off
xlabel('$\varepsilon$','interpreter','latex')
ylabel('absolute max norm error')
legend([hHS,hHSx,hHSxx],{'values','1 deriv','2 deriv'},...
       'location','southeast')