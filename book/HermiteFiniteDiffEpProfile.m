% HermiteFiniteDiffConv
% This example studies the quality of the finite difference approximation
% with Gaussians as a function of epsilon.  A single stencil is
% considered and the approximation is studied with a range of shape
% parameters for the central point.  Then the approximation is considered
% for one of the points off center

% Choose to use the absolute max-norm error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 200;

% Choose a function to test
uf = @(x) exp(-2*x(:,1).*x(:,2));
ufxy = @(x) (4*x(:,1).*x(:,2)-2).*exp(-2*x(:,1).*x(:,2));

% Define an epsilon range to study
% Top one for squished stencil, bottom for round
epvec = [logspace(-3,-2.2,6), ...
         logspace(-2.18,-1.8,30), ...
         logspace(-1.75,-1.06,6), ...
         logspace(-1,-.6,30), ...
         logspace(-.55,0,30)];
% epvec = [logspace(-3,-2.5,30), ...
%          logspace(-2.48,-.9,7), ...
%          logspace(-.87,-.2,55), ...
%          logspace(-.18,0,3)];

% Create a stencil to be used for differencing
% Change xtest size for different shapes: 10 for squished, 12 for round
xeval = [0 0;.5 0];
Nx = 80;
xtest = pick2Dpoints(-1,1,10,'halt');
x = xtest(knnsearch(xtest,xeval(1,:),'K',Nx),:);
u = uf(x);
uxyeval = ufxy(xeval);

h_waitbar = waitbar(0,'Initializing');
errvec = zeros(2,length(epvec));
k = 1;
alpha = 1;
for ep=epvec
    waitbar((k-1)/length(epvec),h_waitbar,sprintf('Computing \\epsilon=%g',ep));
    % Find the HS-SVD differentiation matrix
    GQR = gqr_solveprep(0,x,ep,alpha);
    CbarT = GQR.Rbar;
    Psi = GQR.stored_phi1 + GQR.stored_phi2*CbarT;
    Phixy = gqr_phi(GQR,xeval,[1 1]);
    Psixy = Phixy*[eye(Nx);CbarT];
    FDmat = Psixy/Psi;
    uhatxyeval = FDmat*u;
    
    % Compute the errors
    errvec(:,k) = abs(uhatxyeval-uxyeval);
    k = k + 1;
end

waitbar(1,h_waitbar,sprintf('Plotting'));
% Plot the stencil
h_stencil = figure;
plot(x(:,1),x(:,2),'.r','linewidth',2,'markersize',12)
hold on
h_points(1) = plot(xeval(1,1),xeval(1,2),'x','linewidth',2,'markersize',12);
h_points(2) = plot(xeval(2,1),xeval(2,2),'^','linewidth',2,'markersize',12);
xlim([-1 1]),ylim([-1 1])
hold off
legend(h_points,'centered','skewed','location','northwest')

% Plot the results
h_errs = figure;
loglog(epvec,errvec(1,:),'linewidth',2);
hold on
loglog(epvec,errvec(2,:),'--','linewidth',2);
hold off
ylim([1e-8 1])
xlabel('$\varepsilon$','interpreter','latex')
ylabel('absolute max norm error')
legend({'Centered','Skewed'},'location','northwest')

close(h_waitbar)