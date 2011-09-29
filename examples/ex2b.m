% ex2
% This example compares RBF-Direct to RBF-QR and RBF-QRr
rbfsetup

epvecd = logspace(-2,1,20);
epvecr = logspace(-2,1,20);
N = [7;7];
NN = [50;50];

spaceopt = 'cheb';
fopt = 'tanh';
rbf = @(ep,x) exp(-(ep*x).^2);

[yf,fstr] = pickfunc(fopt,2);

aa = [-1 -1];bb = [1 1];
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);
errvecr = zeros(size(epvecr));
errvecd = zeros(size(epvecd));

status = 'Performing RBF-QRr'
k = 1;
for ep=epvecr
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    rbfqrOBJ = rbfqrr_solve_alpha(x,y,ep);
    yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
    errvecr(k) = norm((yy-yp)./(abs(yy)+eps))/prod(NN);
    fprintf(' %d ',k)
    k = k+1;
end

status = '\nPerforming RBF-Direct'
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
    errvecd(k) = norm((yy-yp)./(abs(yy)+eps))/prod(NN);
    fprintf(' %d ',k)
    k = k+1;
end

loglog(epvecd,errvecd,'-.')
hold on
loglog(epvecr,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([10^-16 1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('Direct','QR','Location','SouthWest')