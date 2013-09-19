% ex17e_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,1,31);
%epvec = logspace(-10,0,31);
%epvec = fliplr(epvec);

N = 15;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
%yf = @(x) x + 1./(1+x.^2);
%fstring = 'y(x) = x + 1/(1+x^2)';
yf = @(x) x.^3-3*x.^2+2*x+1 + 1e-10*cos(10*x);
fstring = 'y(x) = x^3-3x^2+2x+1 + 10^{-10}cos(10x)';
%yf = @(x) x;% + 0.001*cos(10*x);
%fstring = 'y(x) = x';% + cos(10x)/1000';
%yf = @(x) 0.75*exp(-((9*(x+1)/2-2).^2)/4)+0.75*exp(-((9*(x+1)/2+1).^2/49))+0.5*exp(-((9*(x+1)/2-7).^2)/4)-0.2*exp(-((9*(x+1)/2-4).^2));
%fstring = '"Franke"';
%yf = @(x) tanh(9*(x-1))+1;
%fstring = 'tanh(9(x-1))+1';

% Why does this artificial "noise" make things work?

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
%lamratio = 1e-12;
lamratio = 0;
pinvtol = 1e-12;

errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
bvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];
lambdas = zeros(N,length(epvec));
lambdainv = zeros(N,length(epvec));
lambdasave = zeros(N,length(epvec));

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
    Phi2 = GQR.stored_phi2;
    %condPhi1(k) = cond(Phi1);
    [U,S,V] = svd(Phi1);
    dS = diag(S);
    yPhi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    %yPhi = yPhi.*(yPhi>1e-4);
    logdetPhi = sum(log(dS));
    
    Psi = Phi1 + Phi2*GQR.Rbar;
    %condPsi(k) = cond(Psi);
    [U,S,V] = svd(Psi);
    dS = diag(S);
    yPsi = V*((1./dS.*(dS/max(dS)>pinvtol)).*(U'*y));
    %yPsi = yPsi.*(yPsi>1e-4);
    logdetPsi = sum(log(dS));

    diffy(k) = norm((yPhi-yPsi)./yPhi);
    diffPsi(k) = norm(Phi1-Psi);
    
    beta = (1+(2*ep/alpha)^2)^.25;
    delta2 = alpha^2/2*(beta^2-1);
    ead = ep^2 + alpha^2 + delta2;
    Lambda1 = sqrt(alpha^2/ead)*(ep^2/ead).^(0:N-1)';
    Lambda2 = sqrt(alpha^2/ead)*(ep^2/ead).^(N:size(GQR.Marr,2)-1)';
     
    logdetK = logdetPsi + logdetPhi + sum(log(Lambda1));
    
    laminv = 1./Lambda1;
    lamsave = laminv.*(laminv/laminv(end)>lamratio);
    
    % Mahaldist first version
    mahaldist = yPhi'*(lamsave.*yPsi);
    mdist1(k) = mahaldist;

    % Mahaldist second version
    mdist2(k) = yPhi'*(laminv.*yPsi);
    %mahaldist = mdist2(k);
    
    % Mahaldist third version
    b = Psi\y;
    %bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b))'*((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
    bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
    bvec(k) = bvector'*bvector;
    mahaldist = b'*(lamsave.*b)+ bvec(k);
    mdist3(k) = mahaldist;
    
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
    S = svd(A);
    ddetvec(k) = 1/N*sum(log(S));
    %ddetvec(k) = 1/N*log(det(A));
    dlvec(k) =  dmvec(k) + ddetvec(k);
    warning on

    k = k + 1;
end

figure
loglog(epvec,errvec,'color',[0 .5 0],'linewidth',3), hold on
loglog(epvec,exp(lvec)/exp(lvec(end)),'r','linewidth',3)
loglog(epvec,derrvec,'color',[0 .5 0],'linestyle','--','linewidth',3), hold on
loglog(epvec,exp(dlvec)/exp(dlvec(end)),'--r','linewidth',3)
legend('error HS-SVD','MLE HS-SVD','error direct','MLE direct')
xlabel('\epsilon')
ylabel('Error')
title([fstring,', N = ',num2str(N)]), hold off
figure
loglog(epvec,exp(lvec),'r','linewidth',3), hold on
loglog(epvec,exp(dlvec),'--r','linewidth',3)
loglog(epvec,exp(mvec),'b','linewidth',3)
loglog(epvec,exp(dmvec),'--b','linewidth',3)
loglog(epvec,exp(detvec),'k','linewidth',3)
loglog(epvec,exp(ddetvec),'--k','linewidth',3)
legend('MLE HS-SVD','MLE direct','H_K-norm HS','H_K-norm direct','det(K) HS','det(K) direct')
xlabel('\epsilon')
ylabel('log-like function')
title([fstring,', N = ',num2str(N)]), hold off
figure
loglog(epvec,mdist3,'color',[0 .5 0],'linewidth',3), hold on
loglog(epvec,mdist3-bvec,'--c','linewidth',3)
loglog(epvec,bvec,'b','linewidth',3)
legend('H_K-norm HS','lower bound (L1inv)','lower bound (L2)')
xlabel('\epsilon')
ylabel('log-like function')
title([fstring,', N = ',num2str(N)]), hold off