% ex2
% This example compares RBF-Direct to RBF-QR and RBF-QRr
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;
GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 4;

epvecd = logspace(-1,1,20);
epvecr = logspace(-1,1,20);
Nvec = [11,15,19;11,15,19];
NN = [50;50];

spaceopt = 'even';
fopt = 'sin';
rbf = @(ep,x) exp(-(ep*x).^2);

[yf,fstr] = pickfunc(fopt,2);

aa = [-1 -1];bb = [1 1];
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);
errvecr = zeros(size(Nvec,2),length(epvecr));
errvecd = zeros(size(Nvec,2),length(epvecd));

status = 'Performing RBF-QRr'
j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    k = 1;
    for ep=epvecr
        rbfqrOBJ = rbfqrr_solve(x,y,ep);
        yp = rbfqr_eval(rbfqrOBJ,xx);
        errvecr(j,k) = errcompute(yp,yy);
        k = k+1;
    end
    j = j+1;
end

status = 'Performing RBF-Direct'
j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    k = 1;
    for ep=epvecd
        [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
        y = yf(x);
        DM_DATA = DistanceMatrix(x,x);
        IM = rbf(ep,DM_DATA);
        warning off % I know it's bad
        beta = IM\y;
        warning on
        DM_EVAL = DistanceMatrix(xx,x);
        EM = rbf(ep,DM_EVAL);
        yp = EM*beta;
        errvecd(j,k) = errcompute(yp,yy);
        k = k+1;
    end
    j = j+1;
end

loglog(epvecd,errvecd,'-.','LineWidth',2)
hold on
loglog(epvecr,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([1e-16,1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=11^2 (Direct)','N=15^2 (Direct)','N=19^2 (Direct)','N=11^2 (QR)','N=15^2 (QR)','N=19^2 (QR)','Location','SouthEast')
