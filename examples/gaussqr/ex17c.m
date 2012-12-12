% MLE using truncated Hilbert-Schmidt SVD

% Test studying the leave half out cross-validation
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

rbf = @(e,r) exp(-(e*r).^2);
 f = @(x) sin(2*pi*x);
%f = @(x) 1-x.^2;

N = 24;
x = pickpoints(-1,1,N);
y = f(x);

NN = 100;
xx = pickpoints(-1,1,NN);
yy = f(xx);

alpha = 2;

THStol = 1e-12;

DM = DistanceMatrix(x,x);

epvec = [logspace(-2,-1,5),logspace(-1,.3,33),logspace(.3,1,6)];
GQRvec = [];
THSvec = [];
MLEvec = [];
errvec = [];
k = 1;
for ep=epvec
    GQR = gqr_solveprep(0,x,ep,alpha);
    phi = gqr_phi(GQR,x);
    phi1 = phi(:,1:N);
    psi = phi*[eye(N);GQR.Rbar];

    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    laminv = 1./lamvec;
    
    
    GQR_ip = log((phi1\y)'*((psi\y).*laminv));
    GQR_det = log(det(psi)) + log(det(phi1)) + sum(log(lamvec));
    GQRvec(k) = GQR_ip + 1/N*GQR_det;
    
    
    lam_drop = find(lamvec<THStol);
    lamvec(lam_drop) = ones(size(lam_drop));
    laminv(lam_drop) = zeros(size(lam_drop));
    GQR_ip = log((phi1\y)'*((psi\y).*laminv));
    GQR_det = log(det(psi)) + log(det(phi1)) + sum(log(lamvec));
    THSvec(k) = GQR_ip + 1/N*GQR_det;
    
    
    A = rbf(ep,DM);
    [U,S,V] = svd(A);
    MLEvec(k) = log((V\y)'*((1./diag(S)).*(U\y))) + 1/N*sum(log(diag(S)));
%     MLEvec(k) = log(y'*(A\y)) + 1/N*log(det(A));
    
    
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    
    fprintf('k=%d \t ep=%g \n',k,ep)
    k = k + 1;
end

% loglog(epvec,MLEvec,'--k','linewidth',3)
% loglog(epvec,GQRvec,'r','linewidth',3), hold on
% loglog(epvec,THSvec,'g','linewidth',3)
semilogx(epvec,GQRvec,'r','linewidth',3), hold on
semilogx(epvec,THSvec,'g','linewidth',3), hold off
%loglog(epvec,errvec,'linewidth',3), hold off
xlabel('\epsilon')
ylabel('error')
legend('MLE GaussQR','MLE Truncated','True solution')