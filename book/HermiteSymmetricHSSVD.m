% HermiteSymmetricHSSVD
% This example demonstrates how to conduct Hermite interpolation in the
% HS-SVD basis to avoid ill-conditioning.
% The function of interest is
%            f(x,y) = tanh(5(x-y));
% Data will be sampled on a uniform grid, and the data will consist of
% function values and gradient values
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 100;

% Create the necessary functions
uf = @(x) tanh(5*(x(:,1)-x(:,2)));
uxf = @(x) 5*sech(5*(x(:,1)-x(:,2))).^2;
uyf = @(x) -5*sech(5*(x(:,1)-x(:,2))).^2;

% Define the isotropic Gaussian kernel, so we can use the HS-SVD
rbf = @(e,r) exp(-(e*r).^2);
rbfx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
rbfy = @(e,r,dy) -2*e^2*dy.*exp(-(e*r).^2);
rbfxx = @(e,r,dx) 2*e^2*(2*e^2*dx.^2-1).*exp(-(e*r).^2);
rbfyy = @(e,r,dy) 2*e^2*(2*e^2*dy.^2-1).*exp(-(e*r).^2);
rbfxy = @(e,r,dx,dy) 4*e^4*dx.*dy.*exp(-(e*r).^2);

% Choose some data locations
x = pick2Dpoints(-1,1,7);
N = size(x,1);
Neval = 30;
xeval = pick2Dpoints(-1,1,Neval);

% Choose a range of epsilon values to study
epvec = [logspace(-1,.5,10),3];
gqr_alpha = 3;

% Start the waitbar
h_waitbar = waitbar(0,'Initializing','visible','on');

% Create the data to be used for interpolation
u = uf(x);
ux = uxf(x);
uy = uyf(x);
udata = [u;ux;uy];
ueval = uf(xeval);

% Form the necessary distance matrices
DM = DistanceMatrix(x,x);
DiffMx = DifferenceMatrix(x(:,1),x(:,1));
DiffMy = DifferenceMatrix(x(:,2),x(:,2));
DMeval = DistanceMatrix(xeval,x);
DiffMxeval = DifferenceMatrix(xeval(:,1),x(:,1));
DiffMyeval = DifferenceMatrix(xeval(:,2),x(:,2));

errvec = zeros(size(epvec));
errvecHS = zeros(size(epvec));
k = 1;
for ep=epvec
    h_waitbar = waitbar((k-1)/length(epvec),h_waitbar,sprintf('Computing ep=%g',ep));
    % Evaluate the appropriate kernel matrices
    K = rbf(ep,DM);
    Kx = rbfx(ep,DM,DiffMx);
    Ky = rbfy(ep,DM,DiffMy);
    Kxx = rbfxx(ep,DM,DiffMx);
    Kyy = rbfyy(ep,DM,DiffMy);
    Kxy = rbfxy(ep,DM,DiffMx,DiffMy);
    
    % Form the Hermite interpolation matrix
    H = [ K  -Kx  -Ky;
         Kx -Kxx -Kxy;
         Ky -Kxy -Kyy];
    
    % Evaluate the solution coefficients
    c = H\udata;
    
    % Evaluate the solution at the locations
    Keval = rbf(ep,DMeval);
    Kxeval = rbfx(ep,DMeval,DiffMxeval);
    Kyeval = rbfy(ep,DMeval,DiffMyeval);
    uhateval = [Keval -Kxeval -Kyeval]*c;
    
    % Create the GQR object
    GQR = gqr_solveprep(0,x,ep,gqr_alpha);
    GQR_x = GQR;
    GQR_y = GQR;
    Marrall = GQR.Marr;
    Phi = gqr_phi(GQR,x);
    Phix = gqr_phi(GQR,x,[1 0]);
    Phiy = gqr_phi(GQR,x,[0 1]);
    
    % Pivot out unnecessary columns
    Min = 1:size(Marrall,2);
    [~,R,e] = qr(Phi(:,Min(1:N)),'vector');
    badind = e(abs(diag(R))<1e-10);
    Min = Min(setdiff(1:length(Min),badind));
    while not(isempty(badind))
        [~,R,e] = qr(Phi(:,Min(1:N)),'vector');
        badind = e(abs(diag(R))<1e-10);
        Min = Min(setdiff(1:length(Min),badind));
    end
    GQR.Marr = Marrall(:,Min);
    Phi = gqr_phi(GQR,x);
    Phi1 = Phi(:,1:N);
    Phi2 = Phi(:,N+1:end);
    
    % Pivot out unnecessary columns
    Min = 1:size(Marrall,2);
    [~,R,e] = qr(Phix(:,Min(1:N)),'vector');
    badind = e(abs(diag(R))<1e-10);
    Min = Min(setdiff(1:length(Min),badind));
    while not(isempty(badind))
        [~,R,e] = qr(Phix(:,Min(1:N)),'vector');
        log10(abs(diag(R)));
        badind = e(abs(diag(R))<1e-10);
        Min = Min(setdiff(1:length(Min),badind));
    end
    GQR_x.Marr = Marrall(:,Min);
    Phix_x = gqr_phi(GQR_x,x,[1 0]);
    Phix_x1 = Phix_x(:,1:N);
    Phix_x2 = Phix_x(:,N+1:end);
    
    % Pivot out unnecessary columns
    Min = 1:size(Marrall,2);
    [~,R,e] = qr(Phiy(:,Min(1:N)),'vector');
    badind = e(abs(diag(R))<1e-10);
    Min = Min(setdiff(1:length(Min),badind));
    while not(isempty(badind))
        [~,R,e] = qr(Phiy(:,Min(1:N)),'vector');
        badind = e(abs(diag(R))<1e-10);
        Min = Min(setdiff(1:length(Min),badind));
    end
    GQR_y.Marr = Marrall(:,Min);
    Phiy_y = gqr_phi(GQR_y,x,[0 1]);
    Phiy_y1 = Phiy_y(:,1:N);
    Phiy_y2 = Phiy_y(:,N+1:end);
    
    % Create the matrix that operates as Lam2*...*inv(Lam1)
    lamvec = GQR.eig(GQR.Marr);
    lamFull = bsxfun(@rdivide,lamvec(N+1:end)',lamvec(1:N));
    lamvec_x = GQR_x.eig(GQR_x.Marr);
    lamFull_x = bsxfun(@rdivide,lamvec_x(N+1:end)',lamvec_x(1:N));
    lamvec_y = GQR_y.eig(GQR_y.Marr);
    lamFull_y = bsxfun(@rdivide,lamvec_y(N+1:end)',lamvec_y(1:N));
    
    % Create the 3 CbarT matrices for the HS-SVD
    CbarT   =   lamFull.*(Phi2'/Phi1');
    CbarT_x = lamFull_x.*(Phix_x2'/Phix_x1');
    CbarT_y = lamFull_y.*(Phiy_y2'/Phiy_y1');
    GQR.CbarT = CbarT;
    GQR_x.CbarT = CbarT_x;
    GQR_y.CbarT = CbarT_y;
    
    % Evaluate the remaining necessary Phi matrices
    Phix = gqr_phi(GQR,x,[1 0]);
    Phiy = gqr_phi(GQR,x,[0 1]);
    Phi_x = gqr_phi(GQR_x,x);
    Phiy_x = gqr_phi(GQR_x,x,[0 1]);
    Phi_y = gqr_phi(GQR_y,x);
    Phix_y = gqr_phi(GQR_y,x,[1 0]);
    
    % Form the 9 Psi matrices
    Psi = Phi1 + Phi2*CbarT;
    Psix = Phix*[eye(N);CbarT];
    Psiy = Phiy*[eye(N);CbarT];
    Psi_x = Phi_x*[eye(N);CbarT_x];
    Psix_x = Phix_x1 + Phix_x2*CbarT_x;
    Psiy_x = Phiy_x*[eye(N);CbarT_x];
    Psi_y = Phi_y*[eye(N);CbarT_y];
    Psix_y = Phix_y*[eye(N);CbarT_y];
    Psiy_y = Phiy_y1 + Phiy_y2*CbarT_y;
    
    % Form the Hermite interpolation matrix
    bigPsi = [Psi   Psi_x  Psi_y;
              Psix Psix_x Psix_y;
              Psiy Psiy_x Psiy_y];
    
    % Form the Psi matrix for evaluating the interpolant
    Phieval = gqr_phi(GQR,xeval);
    Phieval_x = gqr_phi(GQR_x,xeval);
    Phieval_y = gqr_phi(GQR_y,xeval);
    Psieval = [Phieval*[eye(N);CbarT],...
               Phieval_x*[eye(N);CbarT_x],...
               Phieval_y*[eye(N);CbarT_y]];
    
    % Evaluate the HS-SVD basis interpolant
    uhatHSeval = Psieval*(bigPsi\udata);
    
    errvec(k) = errcompute(uhateval,ueval);
    errvecHS(k) = errcompute(uhatHSeval,ueval);
    k = k + 1;
end

h_waitbar = waitbar(1,h_waitbar,'Plotting');

h_errs = figure;
loglog(epvec(1:end-1),[errvec(1:end-1);errvecHS(1:end-1)],'linewidth',2)
xlabel('$\varepsilon$','interpreter','latex')
xlabel('absolute max norm error','interpreter','latex')

% Create a surface plot of the data
h_surf = figure;
X = reshape(xeval(:,1),Neval,Neval);
Y = reshape(xeval(:,2),Neval,Neval);
Uhat = reshape(uhatHSeval,Neval,Neval);
surf(X,Y,Uhat)

close(h_waitbar)