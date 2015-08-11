% ex2
% This example compares the convergence behavior of the low-rank
% eigenfunction approximate basis for different
% eigenfunction series lengths.
global GAUSSQR_PARAMETERS

Nvec = 2:4:30;
NN = 15;
ep = .1;
alpha = 1;
Mfracvec = .1:.2:.9;

[yf,fstr] = pickfunc('franke_centered');

xx = pick2Dpoints(-1,1,NN);
yy = yf(xx);
errvecr = zeros(length(Mfracvec),size(Nvec,2));

k = 1;
for N=Nvec
    x = pick2Dpoints(-1,1,N,'cheb');
    y = yf(x);
    
    i = 1;
    for Mfrac=Mfracvec
        GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = Mfrac;
        GQR = gqr_rsolve(x,y,ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvecr(i,k) = errcompute(yp,yy);
        i = i + 1;
    end
    
    k = k + 1;
end

loglog(Nvec,errvecr,'linewidth',2);
xlabel('N');
ylabel('absolute error');
ptsstr=strcat(', x\in[-1,1],');
title(strcat(fstr,ptsstr,'cheb'))
legend('M=.1N','M=.3N','M=.5N','M=.7N','M=.9N','Location','SouthWest')
