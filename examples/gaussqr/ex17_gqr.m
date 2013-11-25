% Test studying the kriging variance
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = -300;

rbf = @(e,r) exp(-(e*r).^2);
yf = @(x) cos(pi*x);

N = 10;
x = pickpoints(-1,1,N);
y = yf(x);

NN = 100;
xx = pickpoints(-1,1,NN);
yy = yf(xx);

alpha = 1;
I = eye(N);

epvec = logspace(-2,1,50);
kvvec = [];
errvec = [];
k = 1;
for ep=epvec
    GQR = gqr_solveprep(0,x,ep,alpha);
    Rbar = GQR.Rbar;
    Phi = gqr_phi(GQR,x);
    Phi1 = Phi(:,1:N);
    Psi = Phi*[I;Rbar];
    invPsi = pinv(Psi);
    
    kx = rbf(ep,DistanceMatrix(x,xx));
    psix = gqr_phi(GQR,xx)*[I;Rbar];
    kvvec(k) = sqrt(errcompute(rbf(ep,0)-sum(psix'.*(invPsi*kx))));
    
    GQR.coef = invPsi*y;
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    k = k + 1;
    fprintf('%g\t%g\t%g\n',k,ep,cond(Psi))
end

loglog(epvec,kvvec,'linewidth',3)
hold on
loglog(epvec,errvec,'r','linewidth',3)
hold off
xlabel('\epsilon')
legend('Power function','Interp error','location','northwest')