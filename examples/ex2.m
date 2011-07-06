% ex2
% This example compares RBF-Direct to RBF-QR and RBF-QRr
rbfsetup

epvecd = logspace(-2,1,29);
epvecr = logspace(-2,1,36);
epvec = logspace(-2,0,21);
N = 100;
NN = 1000;

warning off

spaceopt = 'cheb';
fopt = 'runge';

[yf,fstr] = pickfunc(fopt,1);

aa = -3;bb = 3;
xx = linspace(aa,bb,NN)';
yy = yf(xx);
errvec = zeros(size(epvec));
errvecr = zeros(size(epvec));
errvecd = zeros(size(epvecd));

status = 'Performing RBF-QR'
k = 1;
for ep=epvec
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    a = 100;
    rbfqrOBJ = rbfqr_solve(x,y,ep,a);
    yp = rbfqr_eval(rbfqrOBJ,xx);
    errvec(k) = norm((yy-yp)./(abs(yy)+eps));
    k = k+1;
end

status = 'Performing RBF-QRr'
k = 1;
for ep=epvecr
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
    y = yf(x);
    M = 45 + (log10(ep)+2)*7; % Chosen sort of randomly
    rbfqrOBJ = rbfqrr_solve(x,y,ep,0,M);
    yp = rbfqr_eval(rbfqrOBJ,xx);
    errvecr(k) = norm((yy-yp)./(abs(yy)+eps));
    k = k+1;
end

status = 'Performing RBF-Direct'
k = 1;
for ep=epvecd
    K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
    beta = K\y;
    yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
    errvecd(k) = norm((yy-yp)./(abs(yy)+eps));
    k = k+1;
end

[x,spacestr] = pickpoints(aa,bb,N,spaceopt);
y = yf(x);
yp = polyval(polyfit(x,y,44),xx);
errpoly = norm((yy-yp)./(abs(yy)+eps));

a = 1;

loglog(epvecd,errvecd/NN,'--')
hold on
loglog(epvec,errvec/NN,'-.','LineWidth',2)
loglog(epvecr,errvecr/NN,'LineWidth',3)
loglog(epvecd,errpoly*ones(size(epvecd)))
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([10^-13 1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],',' a=',num2str(a),',');
title(strcat(fstr,ptsstr,spacestr))
legend('Direct','QR','Regression','polyfit(44)','Location','SouthWest')