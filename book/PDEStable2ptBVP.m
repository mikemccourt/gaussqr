% PDEStable2ptBVP
% This example demonstrates how the HS-SVD basis can be used to evaluate
% the solution to PDEs stably
% We use the Gaussian here because we have access to its eigenvalues
% The problem of interest here is
%             u'' = -sin(x)
%         u(0) = 0, u(pi) = 0
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define the Gaussian and its 2nd derivative
rbf = @(e,r) exp(-(e*r).^2);
rbfxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

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

% Evaluate the solution where it is needed
ueval = uf(xeval);

% Define a range of epsilon parameters to study
epvec = logspace(-1,1,20);
errvec = zeros(size(epvec));
errvecH = zeros(size(epvec));

k = 1;
for ep=epvec
    % Solve with the stable basis
    ustandard = rbf(ep,DMeval)*([rbfxx(ep,DMint);rbf(ep,DMbc)]\rhs);
    
    % Use the HS-SVD basis
    GQR = gqr_solveprep(0,x,ep,1);
    IdenCorr = [eye(N);GQR.CbarT];
    Psixxint = gqr_phi(GQR,xint,2)*IdenCorr;
    Psibc = gqr_phi(GQR,xbc)*IdenCorr;
    A = [Psixxint;Psibc];
    GQR.coef = A\rhs;
    uhssvd = gqr_eval(GQR,xeval);
    
    % Evaluate the errors
    errvec(k) = errcompute(ustandard,ueval);
    errvecH(k) = errcompute(uhssvd,ueval);
    k = k + 1;
end

% Plot the results
h = figure;
loglog(epvec,errvec,'--','linewidth',2)
hold on
loglog(epvec,errvecH,'linewidth',2)
hold off
xlabel('$\varepsilon$','interpreter','latex')
ylabel('relative rms error')
legend('standard basis','hssvd basis','location','north')