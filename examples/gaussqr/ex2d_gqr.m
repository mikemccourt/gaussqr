% ex2
% This example compares RBF-Direct to RBF-QRr
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 1;

epvecd = logspace(-2,1,20);
epvecr = logspace(-2,1,20);
Nvec = [400,900,1600,2500];
NN = 4000;

spaceopt = 'halton';
rbf = @(ep,x) exp(-(ep*x).^2);

yf = @(x) 1+(x(:,1)+x(:,2)+x(:,3)).^2.*(x(:,4)-x(:,5)).^2.*(x(:,1)+x(:,4));
yf = @(x) sin(.2*sum(x,2))

xx = 2*haltonseq(NN,5)-1;
yy = yf(xx);
errvecr = zeros(size(Nvec,2),length(epvecr));
errvecd = zeros(size(Nvec,2),length(epvecd));

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
        GQR = gqr_rsolve(x,y,ep,alphavals(k));
        yp = gqr_eval(GQR,xx);
        errvecr(j,k) = errcompute(yp,yy,NN-N);
        fprintf(' %g ',GQR.alpha)
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

loglog(epvecd,errvecd,'-.','LineWidth',2)
hold on
loglog(epvecr,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
title('sin(mean(x)), Halton points')
legend('N=400 (Direct)','N=900 (Direct)','N=1600 (Direct)','N=2500 (Direct)',...
       'N=400 (QR)',    'N=900 (QR)',    'N=1600 (QR)',    'N=2500 (QR)','Location','SouthEast')
