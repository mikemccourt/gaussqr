% Test studying the kriging variance
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = -300;

rbf = @(e,r) exp(-(e*r).^2);
yf = @(x) cos(pi*x);
%yf = @(x) besselj(0,abs(4*x));

spaceopt = 'even';
N = 10;
x = pickpoints(-1,1,N,spaceopt);
y = yf(x);

NN = 100;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

alpha = 1;
I = eye(N);
psi_len = N;

epvec = logspace(-2,1,200);
kvvec = zeros(size(epvec));
servec = zeros(size(epvec));
errvec = zeros(size(epvec));
DM_INT = DistanceMatrix(x,x);
DM_EVAL = DistanceMatrix(xx,x);
DM_KV = DistanceMatrix(x,xx);
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
%     yp = gqr_eval(GQR,xx);
%     yp = rbf(ep,DM_EVAL)*(rbf(ep,DM_INT)\y);
%     errvec(k) = errcompute(yp,yy);
    
    Rbar = GQR.Rbar;
    Psi = gqr_phi(GQR,x)*[I;Rbar];
    Lam = diag(GQR.eig(GQR.Marr(1:N)));
    
    kx = rbf(ep,DM_KV);
    psix = gqr_phi(GQR,xx)*[I;Rbar];
    kvvec(k) = sqrt(errcompute(rbf(ep,0)-sum((psix/Psi)'.*kx)));
    errvec(k) = sqrt(errcompute(rbf(ep,0)-sum((rbf(ep,DM_EVAL)/rbf(ep,DM_INT))'.*kx)));
%     servec(k) = sqrt(errcompute(rbf(ep,0)-sum((psix*Lam).*psix,2)));
    servec(k) = sqrt(errcompute(rbf(ep,0)-sum((psix(:,1:psi_len)*Lam(1:psi_len,1:psi_len)).*psix(:,1:psi_len),2)));
    
    k = k + 1;
end

h = figure;
loglog(epvec,errvec.^2,'b','linewidth',3)
hold on
% loglog(epvec,servec,'k','linewidth',3)
loglog(epvec,kvvec.^2,'--r','linewidth',3)
hold off
xlabel('\epsilon')
% legend('Power function','Interp error','location','northwest')
legend('Standard Basis','HS-SVD','location','northwest')