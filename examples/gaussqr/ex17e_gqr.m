% ex17e_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,1,31);

N = 15;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x+1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
yf = @(x) x.^3-3*x.^2+2*x+1;
fstring = 'y(x) = x^3-3x^2+2x+1';

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

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

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Phi1 = GQR.stored_phi1;
    [U,S,V] = svd(Phi1);
    dS = diag(S);
    yPhi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    logdetPhi = sum(log(dS));
    
    Psi = Phi1 + GQR.stored_phi2*GQR.Rbar;
    [U,S,V] = svd(Psi);
    dS = diag(S);
    yPsi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    logdetPsi = sum(log(dS));

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
    warning off
    yp = kbasis*(A\y);
    derrvec(k) = errcompute(yp,yy);
    
%    [U,S,V] = svd(A);
%    MLEvec(k) = real(log((V\y)'*((1./diag(S)).*(U\y))) + 1/N*sum(log(diag(S))));
    dmvec(k) = log(abs(y'*(A\y)));
    ddetvec(k) = 1/N*log(det(A));
    dlvec(k) =  dmvec(k) + ddetvec(k);
    warning on

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