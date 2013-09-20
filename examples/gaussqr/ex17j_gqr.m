% ex17j_gqr
% This example plots "optimal"-epsilon curves based on errors relative to
% the exact solution and via MLE. We use direct computation and HS-SVD in
% 5D.
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 1;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 200; % Set to, e.g., -50 for 1.5N
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

epvecd = logspace(-2,1,20);
epvecr = logspace(-2,1,20);
Nvec = [400,800,1600,3200];
% epvecd = logspace(-2,1,10);
% epvecr = logspace(-2,1,10);
% Nvec = [400,500,600];
NN = 4000;

% If you want to keep track of how long these computations take
active_timer = 1;

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
logdetK = zeros(size(Nvec,2),length(epvecr));
logdetKd = zeros(size(Nvec,2),length(epvecd));
mahaldist = zeros(size(Nvec,2),length(epvecr));
mahaldistd = zeros(size(Nvec,2),length(epvecd));
bvec = zeros(size(Nvec,2),length(epvecr));

if active_timer tic, end
status = 'Finding alpha values'
k = 1;
for ep=epvecr
    alphavals(k) = gqr_alphasearch(ep,[-1,-1,-1,-1,-1],[1,1,1,1,1]);
    k = k+1;
end
if active_timer Total_time = toc, end

if active_timer tic, end
status = 'Performing RBF-QR'
j = 1;
for N=Nvec
    x = 2*haltonseq(N,5)-1;
    y = yf(x);
    k = 1;
    for ep=epvecr
        alpha = alphavals(k);
        GQR = gqr_solve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvecr(j,k) = errcompute(yp,yy,NN-N);
        fprintf(' %g ',alpha)
        % Now, MLE computation
        Phi1 = GQR.stored_phi1;
        Phi2 = GQR.stored_phi2;
        Marr = GQR.Marr;
        S = svd(Phi1); logdetPhi = sum(log(S));
        Psi = Phi1 + Phi2*GQR.Rbar;
        S = svd(Psi); logdetPsi = sum(log(S));
        beta = (1+(2*ep/alpha)^2)^.25;
        delta2 = alpha^2/2*(beta^2-1);
        ead = ep^2 + alpha^2 + delta2;
        Lambda1 = sqrt(alpha^2/ead)*(ep^2/ead).^sum(Marr(:,1:N),1)';
        Lambda2 = sqrt(alpha^2/ead)*(ep^2/ead).^sum(Marr(:,N+1:end),1)';
        logdetK(j,k) = (logdetPsi + logdetPhi + sum(log(Lambda1)))/N;
        laminv = 1./Lambda1;
        % Mahaldist
        warning off
        b = Psi\y;
        bvector = ((Lambda2.^(.5))'*(Phi2')/(Phi1')*(laminv.*b));
        warning on
        bvec(j,k) = bvector'*bvector;
        mahaldist(j,k) = b'*(laminv.*b) + bvec(j,k);
        % Log-likelihood
        lvecr(j,k) = log(abs(mahaldist(j,k))) + logdetK(j,k);
%         fprintf(' %d ',k)
        k = k+1;
    end
    fprintf(' %d \n',N)
    j = j+1;
end
if active_timer Total_time = toc, end

if active_timer tic, end
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
        logdetKd(j,k) = 1/N*sum(log(S));
        mahaldistd(j,k) = y'*beta;
        lvecd(j,k) = log(abs(mahaldistd(j,k))) + logdetKd(j,k);
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
if active_timer Total_time = toc, end

legvals = {};
%for k=1:size(Nvec,2)
%    legvals{k} = sprintf('N=%d (error)',Nvec(k));
%end
%for k=1:size(Nvec,2)
%    legvals{k} = sprintf('N=%d (MLE)',Nvec(k));
%end
for k=1:size(Nvec,2)
    legvals{k} = sprintf('N=%d',Nvec(k));
end

figure
loglog(epvecd,errvecd,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecd),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('error')
title('MLE vs Error (direct)')
legend(legvals,'Location','SouthEast')


figure
loglog(epvecr,errvecr,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecr),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('error')
title('MLE vs Error (HS-SVD)')
legend(legvals,'Location','SouthEast')

figure
loglog(epvecr,errvecr,'LineWidth',3)
hold on
loglog(epvecd,errvecd,'-.','LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('error')
title('Error (direct vs HS-SVD)')
legend(legvals,'Location','SouthEast')

figure
loglog(epvecd,exp(lvecr),'LineWidth',3)
hold on
loglog(epvecd,exp(lvecd),'-.','LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('error')
title('MLE (direct vs HS-SVD)')
legend(legvals,'Location','SouthEast')

figure
loglog(epvecr,exp(logdetK),'linewidth',3), hold on
loglog(epvecd,exp(logdetKd),'-.','linewidth',3)
legend(legvals,'Location','SouthEast')
xlabel('\epsilon')
ylabel('logdet(K)/N')
title('5D: y(x)=sin(mean(x)), direct vs HS-SVD'), hold off

figure
loglog(epvecr,mahaldist,'linewidth',3), hold on
loglog(epvecd,abs(mahaldistd),'-.','linewidth',3)
legend(legvals,'Location','NorthEast')
xlabel('\epsilon')
ylabel('Native space norm')
title('5D: y(x)=sin(mean(x)), direct vs HS-SVD'), hold off

figure
loglog(epvecr,mahaldist(2,:),'color',[0 .5 0],'linewidth',3), hold on
loglog(epvecr,mahaldist(2,:)-bvec(2,:),'--c','linewidth',3)
loglog(epvecr,bvec(2,:),'b','linewidth',3)
legend('H_K-norm HS','lower bound (L1inv)','lower bound (L2)')
xlabel('\epsilon')
ylabel('log-like function')
title(['5D: y(x)=sin(mean(x)), N = ',num2str(Nvec(2))]), hold off