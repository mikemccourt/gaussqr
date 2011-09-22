% Tests RBF-QR for a 2D example
rbfsetup

epvecd = logspace(-2,1,20);
epvec = logspace(-2,0.3,30);
Nvec = [3,7;3,7];
NN = [30;30];

spaceopt = 'cheb';
fopt = 'sin';
rbf = @(ep,x) exp(-(ep*x).^2);

[yf,fstr] = pickfunc(fopt,2);

aa = [-3 -3];bb = [3 3];
xx = pick2Dpoints(aa,bb,NN);
yy = yf(xx);
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvecd));
alpha = 2;

n = 1;
for N=Nvec
    k = 1;
    for ep=epvec
        [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
        y = yf(x);
        rbfqrOBJ = rbfqr_solve_alpha(x,y,ep,alpha);
        yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
        errvec(n,k) = norm((yy-yp)./(abs(yy)+eps))/prod(NN);
        k = k+1;
        fprintf(' %d ',k-1)
    end

    k = 1;
    for ep=epvecd
        [x,spacestr] = pick2Dpoints(aa,bb,N,spaceopt,ep);
        y = yf(x);
        DM_data = DistanceMatrix(x,x);
        IM = rbf(ep,DM_data);
        warning off
        beta = IM\y;
        warning on
        DM_eval = DistanceMatrix(xx,x);
        EM = rbf(ep,DM_eval);
        yp = EM*beta;
        errvecd(n,k) = norm((yy-yp)./(abs(yy)+eps))/prod(NN);
        k = k+1;
    end
    n = n+1;
    fprintf(' %d \n',n-1)
end

figure
loglog(epvecd,errvecd(1,:)/prod(NN),'-bx'), hold on
loglog(epvecd,errvecd(2,:)/prod(NN),'-g+')
%loglog(epvecd,errvecd(3,:)/prod(NN),'-r^')
loglog(epvec,errvec(1,:)/prod(NN),'b','LineWidth',3)
loglog(epvec,errvec(2,:)/prod(NN),'g','LineWidth',3)
%loglog(epvec,errvec(3,:)/prod(NN),'r','LineWidth',3)
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([10^-17 1])
ptsstr=', x\in[-3,3]^2,';
%title(strcat(fstr,ptsstr,spacestr))
title(sprintf('alpha=%g',alpha))
legend('N=7x7 (Direct)','N=11x11 (Direct)','N=15x15 (Direct)','N=7x7 (QR)','N=11x11 (QR)','N=15x15 (QR)','Location','SouthEast')