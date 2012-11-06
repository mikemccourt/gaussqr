% ex2
% This example compares RBF-Direct to RBF-QR
rbfsetup
global GAUSSQR_PARAMETERS
%GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 30;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

epvecd = logspace(-2,1,21);
epvecr = logspace(-2,1,61);
Nvec = [3,5,9;3,5,9];
NN = [30;30];

spaceopt = 'halton';
fopt = 'sin';
rbf = @(ep,x) exp(-(ep*x).^2);

[yf,fstr] = pickfunc(fopt,2);
fstr = 'f(x)=sinc((x+1)/2)*sinc((y+1)/2)';
yf = @(x) sinc((x(:,1)+1)/2).*sinc((x(:,2)+1)/2);

aa = [-1 -1];bb = [1 1];
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);
errvecr = zeros(size(Nvec,2),length(epvecr));
errvecd = zeros(size(Nvec,2),length(epvecd));

alpha = 1;
tic
status = 'Finding alpha values'
k = 1;
for ep=epvecr
%    alphavals(k) = gqr_alphasearch(ep,[-1,-1],[1,1]);
    k = k+1;
end
Total_time = toc

status = 'Performing RBF-QR'
j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    k = 1;
    for ep=epvecr
        GQR = gqr_solve(x,y,ep,alpha);
%        GQR = gqr_solve(x,y,ep,alphavals(k));
        yp = gqr_eval(GQR,xx);
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

figure
loglog(epvecr,errvecr,'LineWidth',3)
hold on
loglog(epvecd,errvecd,'-.','LineWidth',2)
hold off
xlabel('\epsilon')
ylabel('Error')
ylim([1e-10,10])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=9','N=25','N=81','Location','SouthEast')
