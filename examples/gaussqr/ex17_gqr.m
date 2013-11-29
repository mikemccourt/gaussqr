% Test studying the kriging variance
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
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

epvec = [logspace(-2,-1,5),logspace(-1,1,500)];
kvvec = [];
errvec = [];
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Rbar = GQR.Rbar;
    Psi = gqr_phi(GQR,x)*[I;Rbar];
    
    kx = rbf(ep,DistanceMatrix(x,xx));
    psix = gqr_phi(GQR,xx)*[I;Rbar];
    kvvec(k) = sqrt(errcompute(rbf(ep,0)-sum((psix/Psi)'.*kx)));
    
    k = k + 1;
end

loglog(epvec,kvvec,'linewidth',3)
hold on
loglog(epvec,errvec,'r','linewidth',3)
hold off
xlabel('\epsilon')
legend('Power function','Interp error','location','northwest')