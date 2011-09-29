% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
rbfsetup

epvecd = logspace(-2,1,40);
epvec = logspace(-2,0.3,60);
Nvec = [11,31,51];
NN = 100;

spaceopt = 'cheb';
fopt = 'exp';

[yf,fstr] = pickfunc(fopt,1);

aa = -3;bb = 3;
xx = linspace(aa,bb,NN)';
yy = yf(xx);
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvecd));
alpha = 3;

n = 1;
for N=Nvec
    k = 1;
    for ep=epvec
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
        rbfqrOBJ = rbfqr_solve_alpha(x,yf(x),ep,alpha);
        yp = rbfqr_eval_alpha(rbfqrOBJ,xx);
        errvec(n,k) = norm((yy-yp)./(abs(yy)+eps));
        k = k+1;
    end

    k = 1;
    for ep=epvecd
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
        y = yf(x);
        K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
        beta = K\y;
        yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
        errvecd(n,k) = norm((yy-yp)./(abs(yy)+eps));
        k = k+1;
    end
    n = n+1;
end

[x,spacestr] = pickpoints(aa,bb,Nvec(1),spaceopt);
y = yf(x);
yp = polyval(polyfit(x,y,Nvec(1)-1),xx);
errpoly1 = norm((yy-yp)./(abs(yy)+eps))/NN;

[x,spacestr] = pickpoints(aa,bb,Nvec(2),spaceopt);
y = yf(x);
yp = polyval(polyfit(x,y,Nvec(2)-1),xx);
errpoly2 = norm((yy-yp)./(abs(yy)+eps))/NN;

[x,spacestr] = pickpoints(aa,bb,Nvec(3),spaceopt);
y = yf(x);
yp = polyval(polyfit(x,y,Nvec(3)-1),xx);
errpoly3 = norm((yy-yp)./(abs(yy)+eps))/NN;

loglog(epvecd,errvecd(1,:)/NN,'-bx')
hold on
loglog(epvecd,errvecd(2,:)/NN,'-g+')
loglog(epvecd,errvecd(3,:)/NN,'-r^')
loglog(epvec,errvec(1,:)/NN,'b','LineWidth',3)
loglog(epvec,errvec(2,:)/NN,'g','LineWidth',3)
loglog(epvec,errvec(3,:)/NN,'r','LineWidth',3)
loglog(epvecd,errpoly1*ones(size(epvecd)),'--b')
loglog(epvecd,errpoly2*ones(size(epvecd)),'--g')
loglog(epvecd,errpoly3*ones(size(epvecd)),'--r')
hold off
xlabel('\epsilon')
ylabel('average error')
ylim([10^-17 1])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=11 (Direct)','N=31 (Direct)','N=51 (Direct)','N=11 (QR)','N=31 (QR)','N=51 (QR)','Location','SouthEast')