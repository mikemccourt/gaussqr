% ex2
% This example compares RBF-Direct to RBF-QR and RBF-QRr
rbfsetup

epvecd = logspace(-2,1,35);
epvecr = logspace(-2,1,35);
epvec = logspace(-2,0,22);
N = 150;
NN = 1000;

spaceopt = 'cheb';
fopt = 'runge';

[yf,fstr] = pickfunc(fopt,1);

aa = -4;bb = 4;
xx = linspace(aa,bb,NN)';
yy = yf(xx);
errvec = zeros(size(epvec));
errvecr = zeros(size(epvec));
errvecd = zeros(size(epvecd));

status = 'Performing RBF-QRr'
k = 1;
for ep=epvecr
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    rbfqrOBJ = rbfqrr_solve_alpha(x,y,ep);
    yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
    errvecr(k) = norm((yy-yp)./(abs(yy)+eps))/NN;
    fprintf(' %d ',k)
    k = k+1;
end

status = '\nPerforming RBF-QR'
k = 1;
for ep=epvec
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    rbfqrOBJ = rbfqr_solve_alpha(x,y,ep);
    yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
    errvec(k) = norm((yy-yp)./(abs(yy)+eps))/NN;
    fprintf(' %d ',k)
    k = k+1;
end

status = '\nPerforming RBF-Direct'
k = 1;
for ep=epvecd
    K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
    warning off % I know it's bad
    beta = K\y;
    warning on
    yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
    errvecd(k) = norm((yy-yp)./(abs(yy)+eps))/NN;
    fprintf(' %d ',k)
    k = k+1;
end

status = '\nPolynomial computation'
[x,spacestr] = pickpoints(aa,bb,N,spaceopt);
y = yf(x);
warning off
[ppoly,spoly,mupoly] = polyfit(x,y,90);
yp = polyval(ppoly,xx,spoly,mupoly);
warning on
errpoly = norm((yy-yp)./(abs(yy)+eps))/NN;

loglog(epvecd,errvecd,'--')
hold on
loglog(epvec,errvec,'-.','LineWidth',2)
loglog(epvecr,errvecr,'LineWidth',3)
loglog(epvecd,errpoly*ones(size(epvecd)))
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([10^-16 1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('Direct','QR','Regression','polyfit(90)','Location','SouthWest')