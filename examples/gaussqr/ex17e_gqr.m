% ex17e_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error

epvec = logspace(-2,1,20);

N = 10;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
%yf = @(x) x.^3-3*x.^2+2*x+1;
%fstring = 'y(x) = x^3-3x^2+2x+1';

y = yf(x);
xx = pickpoints(-1,1,NN,'cheb');
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;

errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];

rbf = @(e,r) exp(-(e*r).^2);
DM = DistanceMatrix(x,x);
EM = DistanceMatrix(xx,x);

k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Phi1 = GQR.stored_phi1;
    [U,S,V] = svd(Phi1);
    yPhi = V*((1./diag(S)).*(U'*y));
    logdetPhi = sum(log(diag(S)));
    
    Psi = Phi1 + GQR.stored_phi2*GQR.Rbar;
    [U,S,V] = svd(Psi);
    yPsi = V*((1./diag(S)).*(U'*y));
    logdetPsi = sum(log(diag(S)));

    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    lamvec = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    
    logdetK = logdetPsi + logdetPhi + sum(log(lamvec));
    
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);
    
    mahaldist = yPhi'*(lamsave.*yPsi);
    
    mvec(k) = log(abs(mahaldist));
    detvec(k) = 1/N*logdetK;
    lvec(k) = log(abs(mahaldist)) + 1/N*logdetK;

    A = rbf(ep,DM);
    kbasis = rbf(ep,EM);
    yp = kbasis*(A\y);
    derrvec(k) = errcompute(yp,yy);
    
%    [U,S,V] = svd(A);
%    MLEvec(k) = real(log((V\y)'*((1./diag(S)).*(U\y))) + 1/N*sum(log(diag(S))));
    dmvec(k) = log(abs(y'*(A\y)));
    ddetvec(k) = 1/N*log(det(A));
    dlvec(k) =  dmvec(k) + ddetvec(k);

    k = k + 1;
end

loglog(epvec,errvec,'color',[0 .5 0],'linewidth',3), hold on
loglog(epvec,exp(lvec)/exp(lvec(end)),'r','linewidth',3)
loglog(epvec,derrvec,'color',[0 .5 0],'linestyle','--','linewidth',3), hold on
loglog(epvec,exp(dlvec)/exp(dlvec(end)),'--r','linewidth',3)
legend('error HS-SVD','HS-SVD MLE','error direct','MLE direct')
xlabel('\epsilon')
ylabel('Error')
title(fstring), hold off
figure
semilogx(epvec,lvec,'r','linewidth',3), hold on
semilogx(epvec,dlvec,'--r','linewidth',3)
semilogx(epvec,mvec,'b','linewidth',3)
semilogx(epvec,dmvec,'--b','linewidth',3)
semilogx(epvec,detvec,'k','linewidth',3)
semilogx(epvec,ddetvec,'--k','linewidth',3)
legend('- log-like HS-SVD','- log-like direct','log(H_K-norm) HS','log(H_K-norm) direct','logdet(K) HS','logdet(K) direct')
xlabel('\epsilon')
ylabel('log-like function ish')
title(fstring), hold off