% ex2b
% This example compares low rank eigenfunction approximate interpolants to
% the standard basis for a 2D problem
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .5;

epvec = logspace(-1,1,20);
Nvec = [11,15,19];
NN = 25;

rbf = @(ep,r) exp(-(ep*r).^2);

spaceopt = 'even';
fopt = 'sin';

[yf,fstr] = pickfunc(fopt,2);

aa = -1;bb = 1;
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);
errvecr = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvec));

alpha = 1;

j = 1;
for N=Nvec
    [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt);
    y = yf(x);
    k = 1;
    for ep=epvec
        GQR = gqr_rsolve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvecr(j,k) = errcompute(yp,yy);
        
        K = rbf(ep,DistanceMatrix(x,x));
        warning('off','MATLAB:nearlySingularMatrix')
        beta = K\y;
        warning('on','MATLAB:nearlySingularMatrix')
        Keval = rbf(ep,DistanceMatrix(xx,x));
        yp = Keval*beta;
        errvecd(j,k) = errcompute(yp,yy);
        k = k + 1;
    end
    j = j + 1;
end

loglog(epvec,errvecd,'-.','LineWidth',2)
hold on
loglog(epvec,errvecr,'LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([1e-18,.1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=11^2 (Direct)','N=15^2 (Direct)','N=19^2 (Direct)','N=11^2 (QR)','N=15^2 (QR)','N=19^2 (QR)','Location','SouthEast')
