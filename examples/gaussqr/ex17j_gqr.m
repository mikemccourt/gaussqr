% ex2
% This example compares RBF-Direct to RBF-QRr
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 1;

epvecd = logspace(-2,1,20);
epvecr = logspace(-2,1,20);
Nvec = [400,800,1600,3200];
%epvecd = logspace(-2,1,10);
%epvecr = logspace(-2,1,10);
%Nvec = [400,500,600];
NN = 4000;

spaceopt = 'halton';
rbf = @(ep,x) exp(-(ep*x).^2);

yf = @(x) 1+(x(:,1)+x(:,2)+x(:,3)).^2.*(x(:,4)-x(:,5)).^2.*(x(:,1)+x(:,4));
yf = @(x) sin(.2*sum(x,2))

xx = 2*haltonseq(NN,5)-1;
yy = yf(xx);
errvecr = zeros(size(Nvec,2),length(epvecr));
errvecd = zeros(size(Nvec,2),length(epvecd));
lvecr = zeros(size(Nvec,2),length(epvecr));
lvecd = zeros(size(Nvec,2),length(epvecd));

tic
status = 'Finding alpha values'
k = 1;
for ep=epvecr
    alphavals(k) = gqr_alphasearch(ep,[-1,-1,-1,-1,-1],[1,1,1,1,1]);
    k = k+1;
end
Total_time = toc

tic
status = 'Performing RBF-QRr'
j = 1;
for N=Nvec
    x = 2*haltonseq(N,5)-1;
    y = yf(x);
    k = 1;
    for ep=epvecr
        GQR = gqr_solve(x,y,ep,alphavals(k));
        yp = gqr_eval(GQR,xx);
        errvecr(j,k) = errcompute(yp,yy,NN-N);
        fprintf(' %g ',GQR.alpha)
        % Now, MLE computation
        Phi1 = GQR.stored_phi1;
        Phi2 = GQR.stored_phi2;
        S = svd(Phi1); logdetPhi = sum(log(S));
        Psi = Phi1 + Phi2*GQR.Rbar;
        S = svd(Psi); logdetPsi = sum(log(S));
        beta = (1+(2*ep/GQR.alpha)^2)^.25;
        delta2 = GQR.alpha^2/2*(beta^2-1);
        ead = ep^2 + GQR.alpha^2 + delta2;
        Lambda1 = sqrt(GQR.alpha^2/ead)*(ep^2/ead).^(0:N-1)';
        Lambda2 = sqrt(GQR.alpha^2/ead)*(ep^2/ead).^(N:size(GQR.Marr,2)-1)';
        logdetK = logdetPsi + logdetPhi + sum(log(Lambda1));
        laminv = 1./Lambda1;
        % Mahaldist
        warning off
        b = Psi\y;
        bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(laminv.*b));
        warning on
        bvec = bvector'*bvector;
        mahaldist = b'*(laminv.*b) + bvec;
        % Log-likelihood
        lvecr(j,k) = log(abs(mahaldist)) + 1/N*logdetK;
%         fprintf(' %d ',k)
        k = k+1;
    end
    fprintf(' %d \n',N)
    j = j+1;
end
Total_time = toc

tic
status = 'Performing RBF-Direct'
j = 1;
for N=Nvec
    x = 2*haltonseq(N,5)-1;
    y = yf(x);
    k = 1;
    for ep=epvecd
        DM_DATA = DistanceMatrix(x,x);
        IM = rbf(ep,DM_DATA);
        warning off % I know it's bad
        beta = IM\y;
        S = svd(IM); % for MLE computation in next line
        lvecd(j,k) =  log(abs(y'*beta)) + 1/N*sum(log(S));
        warning on
        DM_EVAL = DistanceMatrix(xx,x);
        EM = rbf(ep,DM_EVAL);
        yp = EM*beta;
        errvecd(j,k) = errcompute(yp,yy,NN-N);
        fprintf(' %d ',k)
        k = k+1;
    end
    fprintf(' %d \n',N)
    j = j+1;
end
Total_time = toc

loglog(epvecd,errvecd,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecd),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('average error')
title('sin(mean(x)), Halton points, Direct')
legend('N=400 (error)','N=900 (error)','N=1600 (error)','N=2500 (error)',...
       'N=400 (MLE)',    'N=900 (MLE)',    'N=1600 (MLE)',    'N=2500 (MLE)','Location','SouthEast')

figure
loglog(epvecr,errvecr,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecr),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('average error')
title('sin(mean(x)), Halton points, GaussQR')
legend('N=400 (error)','N=900 (error)','N=1600 (error)','N=2500 (error)',...
       'N=400 (MLE)',    'N=900 (MLE)',    'N=1600 (MLE)',    'N=2500 (MLE)','Location','SouthEast')
