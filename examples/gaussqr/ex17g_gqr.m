% ex17g_gqr.m
% This should compute a hybrid likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,1,31);
%epvec = fliplr(epvec);

N = 50;
NN = 200;
x = pickpoints(-1,1,N,'cheb');
yf = @(x) x + 1./(1+x.^2);
fstring = 'y(x) = x + 1/(1+x^2)';
%yf = @(x) x.^3-3*x.^2+2*x+1;% + 0.001*cos(10*x);
%yf = @(x) 0*x+1;% + 0.001*cos(10*x);
% Why does this artificial "noise" make things work?
%fstring = 'y(x) = x^3-3x^2+2x+1';

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
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

    A = rbf(ep,DM);
    kbasis = rbf(ep,EM);
    c = A\y;
    warning off
    yp = kbasis*c;
    derrvec(k) = errcompute(yp,yy);
    
    if ep < 1
        % Mahaldist first version
        %mahaldist = yPhi'*(lamsave.*yPsi);
        
        % Mahaldist second version
        %mahaldist = yPhi'*(laminv.*yPsi);
    
        % Mahaldist third version
        b = Psi\y;
        bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(lamsave.*b));
        mahaldist = b'*(lamsave.*b)+ bvector'*bvector;
    else
        %    [U,S,V] = svd(A);
        %    MLEvec(k) = real(log((V\y)'*((1./diag(S)).*(U\y))) + 1/N*sum(log(diag(S))));
        mahaldist = y'*c;
        %logdetK = log(det(A));
        logdetK = sum(log(svd(A)));
    end
    warning on
    
    mvec(k) = log(abs(mahaldist));
    detvec(k) = 1/N*logdetK;
    lvec(k) = log(abs(mahaldist)) + 1/N*logdetK;

    k = k + 1;
end

figure
loglog(epvec,errvec,'color',[0 .5 0],'linewidth',3), hold on
loglog(epvec,exp(lvec)/exp(lvec(end)),'r','linewidth',3)
loglog(epvec,derrvec,'color',[0 .5 0],'linestyle','--','linewidth',3), hold on
legend('error HS-SVD','hybrid MLE','error direct')
xlabel('\epsilon')
ylabel('Error')
title(fstring), hold off
figure
semilogx(epvec,lvec,'r','linewidth',3), hold on
semilogx(epvec,mvec,'b','linewidth',3)
semilogx(epvec,detvec,'k','linewidth',3)
legend('- log-like hybrid','hybrid log(H_K-norm)','logdet(K)/N hybrid')
xlabel('\epsilon')
ylabel('log-like function')
title(fstring), hold off