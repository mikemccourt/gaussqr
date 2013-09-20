% ex17j_gqr
% This example compares RBF-QR to RBF-Direct in 5D
% The goal of this example is to 
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
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
if active_timer Total_time = toc, end

legvals = {};
for k=1:size(Nvec,2)
    legvals{k} = sprintf('N=%d (error)',Nvec(k));
end
for k=1:size(Nvec,2)
    legvals{k} = sprintf('N=%d (MLE)',Nvec(k));
end

loglog(epvecd,errvecd,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecd),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('average error')
title('sin(mean(x)), Halton points, Direct')
legend(legvals,'Location','SouthEast')

figure
loglog(epvecr,errvecr,'LineWidth',3)
hold on
loglog(epvecd,exp(lvecr),'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('average error')
title('sin(mean(x)), Halton points, GaussQR')
legend(legvals,'Location','SouthEast')
