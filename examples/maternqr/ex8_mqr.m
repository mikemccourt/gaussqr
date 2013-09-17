% ex8_mqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(0,2,31);
%epvec = fliplr(epvec);

beta = 5; L = 1;
N = 20;
NN = 200;
x = pickpoints(0,L,N+2); x = x(2:end-1);
HMN = 2; gamval = .0567; CC = 1/(.5-gamval)^2; 
yf = @(x) CC^(HMN).*(max(0,x-gamval)).^HMN.*(max(0,(1-gamval)-x)).^HMN.*exp(-36*(x-0.4).^2);
fstring = 'Hubbert-Müller';

y = yf(x);
xx = pickpoints(0,1,NN);
yy = yf(xx);
lamratio = 1e-12;
%lamratio = 0;
pinvtol = 1e-12;

A = zeros(N);
kbasis = zeros(NN,N);
errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];
lambdas = zeros(N,length(epvec));
lambdainv = zeros(N,length(epvec));
lambdasave = zeros(N,length(epvec));

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    MQR = mqr_solve(x,y,L,ep,beta);
    yp = mqr_eval(MQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    Phi1 = MQR.stored_phi1;
    %condPhi1(k) = cond(Phi1);
    [U,S,V] = svd(Phi1);
    dS = diag(S);
    yPhi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    %yPhi = yPhi.*(yPhi>1e-4);
    logdetPhi = sum(log(dS));
    
    Psi = Phi1 + MQR.stored_phi2*MQR.Rbar;
    %condPsi(k) = cond(Psi);
    [U,S,V] = svd(Psi);
    dS = diag(S);
    yPsi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    %yPsi = yPsi.*(yPsi>1e-4);
    logdetPsi = sum(log(dS));

    diffy(k) = norm((yPhi-yPsi)./yPhi);
    diffPsi(k) = norm(Phi1-Psi);
    
    lamvec = (pi^2/L^2*(1:N).^2+ep^2).^(-beta)';
    %lambdas(:,k) = lamvec;
 
    logdetK = logdetPsi + logdetPhi + sum(log(lamvec));
    
    laminv = 1./lamvec;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);
    %lambdainv(:,k) = laminv;
    %lambdasave(:,k) = lamsave;
    
    mahaldist = yPhi'*(lamsave.*yPsi);
    mdist1(k) = mahaldist;
    mdist2(k) = yPhi'*(laminv.*yPsi);
    %mahaldist = mdist2(k);
    
    mvec(k) = log(abs(mahaldist));
    detvec(k) = 1/N*logdetK;
    lvec(k) = log(abs(mahaldist)) + 1/N*logdetK;

    for j=1:N
        A(:,j) = cmatern(x,x(j),L,ep,beta);
        kbasis(:,j) = cmatern(xx,x(j),L,ep,beta);
    end
    warning off
    b = linsolve(A,y);
    yp = kbasis*b;
    %yp = kbasis*(A\y);

    derrvec(k) = errcompute(yp,yy);
    
%    [U,S,V] = svd(A);
%    MLEvec(k) = real(log((V\y)'*((1./diag(S)).*(U\y))) + 1/N*sum(log(diag(S))));
    %dmvec(k) = log(abs(y'*(A\y)));
    dmvec(k) = log(abs(y'*b));
    ddetvec(k) = 1/N*log(det(A));
    dlvec(k) =  dmvec(k) + ddetvec(k);
    warning on

    k = k + 1;
end

figure
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
legend('- log-like HS-SVD','- log-like direct','log(H_K-norm) HS','log(H_K-norm) direct','logdet(K)/N HS','logdet(K)/N direct')
xlabel('\epsilon')
ylabel('log-like function')
title(fstring), hold off