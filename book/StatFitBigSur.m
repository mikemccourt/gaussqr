% StatFitBigSur
% This example loads the Big Sur data set from 
%           Interpolation of track data with radial basis methods
%           R. E. Carlson and T. A. Foley
%           Comput. Math. Appl., 24:27-34, 1992.
% We can access this from the GaussQR data directory
% To demonstrate the use of the HS-SVD, we compute with Gaussians
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 250;

h_waitbar = waitbar(0,'Initializing Big Sur data');

% Load the data in from the repository as needed
% This data comes with sizes Nx and Ny where length(x) = Nx*Ny
gqr_downloaddata('BigSur_data.mat')
load BigSur_data
[x,t] = rescale_data(locations,temperatures);
N = size(x,1);

% Choose to work with the Gaussians in the HS-SVD basis
gqr_alpha = 1;

% Choose some points at which to evaluate the power function
xeval = pick2Dpoints(-1,1,4,'halt');

% Consider a range of epsilon values
epvec = logspace(0,1.1,15);

mplevec = zeros(size(epvec));
gwvec = zeros(size(epvec));
cvvec = zeros(size(epvec));
k = 1;
for ep=epvec
    waitbar(k/length(epvec),h_waitbar,sprintf('\\epsilon=%g',ep));
    GQR = gqr_solveprep(0,x,ep,gqr_alpha);
    
    % Compute the log(det(K))
    Phi1 = GQR.stored_phi1;
    lamvec1 = GQR.eig(GQR.Marr(:,1:N));
    lamvec2 = GQR.eig(GQR.Marr(:,N+1:end));
    CbarT = GQR.CbarT;
    Psi = Phi1 + GQR.stored_phi2*CbarT;
    logdetK = log(abs(det(Phi1))) + sum(log(lamvec1)) + log(abs(det(Psi)));
    
    % Solve the interpolation problem
    b = Psi\t;
    
    % Compute the Hilbert space norm
    hsnorm = sqrt(norm(bsxfun(@ldivide,sqrt(lamvec1)',b)).^2 + ...
                  norm(bsxfun(@ldivide,sqrt(lamvec2)',CbarT*b)).^2);
    
	% Compute the power function at the sample points
    poweval = zeros(size(xeval,1),1);
    j = 1;
    for xp=xeval'
        GQRp = gqr_solveprep(0,[x;xp'],ep,gqr_alpha);
        Phi1 = GQRp.stored_phi1;
        lamvec1 = GQRp.eig(GQRp.Marr(:,1:N));
        lamvec2 = GQRp.eig(GQRp.Marr(:,N+1:end));
        Psi = Phi1 + GQRp.stored_phi2*GQRp.CbarT;
        logdetKp = log(abs(det(Phi1))) + sum(log(lamvec1)) + log(abs(det(Psi)));
        poweval(j) = sqrt(exp(logdetKp - logdetK));
    end
    pow = norm(poweval,'inf');
    
    % Compute the LOOCV residual
    cvvec(k) = crossval('mse',x,t,...
         'Predfun',@(x,y,xeval) gqr_eval(gqr_solve(x,y,ep,gqr_alpha),xeval),...
                       'leaveout',1);
    
    gwvec(k) = hsnorm*pow;
    mplevec(k) = 2*N*log(hsnorm) + logdetK;
    k = k + 1;
end

waitbar(1,h_waitbar,'Plotting')

h_ep = figure;
[AX,h1,h2] = plotyy(epvec,[mplevec;cvvec],epvec,gwvec,'semilogx','loglog');
set(AX(1),'ylim',[0 1000])
set(AX(1),'ytick',[0 500 1000])
set(AX,'xlim',[1 12])
set(get(AX(1),'xlabel'),'string','$\varepsilon$','interpreter','latex')
set(AX,{'ycolor'},{'k';'k'})
set(h1(2),'linestyle','--')
set(h2,'linestyle','-.')
set([h1;h2],'linewidth',2)
legend([h1;h2]',{'C_{MPLE}','C_{CV}','C_{GW}'},'location','northeast')

h_int = figure;
epplot = 8;
Nplot = 35;
xplot = pick2Dpoints(-1,1,Nplot);
splot = gqr_eval(gqr_solve(x,t,epplot,gqr_alpha),xplot);
X = reshape(xplot(:,1),Nplot,Nplot);
Y = reshape(xplot(:,2),Nplot,Nplot);
S = reshape(splot,Nplot,Nplot);
h_surf = surf(X,Y,S);
hold on
plot3(x(:,1),x(:,2),t,'or','linewidth',2)
hold off
xlabel('x')
ylabel('y')
zlabel('temperature')
view([-34.5 42])

close(h_waitbar)