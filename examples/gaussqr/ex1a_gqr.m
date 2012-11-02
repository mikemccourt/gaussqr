% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

epvecd = logspace(-2,1.5,41);
epvec = logspace(-2,1.5,161); % Change to 81 if taking too long
Nvec = [3,5,9,17];
NN = 100;

spaceopt = 'even';
fopt = 'sin';

[yf,fstr] = pickfunc(fopt,1);

fstr = 'y(x)=sinc(x)';
yf = @(x) sinc((x+1)/2);

aa = -1;bb = 1;
xx = pickpoints(aa,bb,NN);
yy = yf(xx);
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvecd));
alpha = 1; % For larger epsilon, may need to be varied

n = 1;
for N=Nvec
    k = 1;
    for ep=epvec
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
        GQR = gqr_solve(x,yf(x),ep,alpha);
        yp = gqr_eval(GQR,xx);
        errvec(n,k) = errcompute(yp,yy);
        k = k+1;
    end

    k = 1;
    for ep=epvecd
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
        y = yf(x);
        K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
        warning off % We know it's bad
        beta = K\y;
        warning on
        yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
        errvecd(n,k) = errcompute(yp,yy);
        k = k+1;
    end
    n = n+1;
end

[x,spacestr] = pickpoints(aa,bb,Nvec(1),spaceopt);
y = yf(x);
warning off
[p,S,mu] = polyfit(x,y,Nvec(1)-1);
warning on
yp = polyval(p,xx,S,mu);
errpoly1 = errcompute(yp,yy);
errpoly1 = max(abs(yp-yy));

[x,spacestr] = pickpoints(aa,bb,Nvec(2),spaceopt);
y = yf(x);
warning off
[p,S,mu] = polyfit(x,y,Nvec(2)-1);
warning on
yp = polyval(p,xx,S,mu);
errpoly2 = errcompute(yp,yy);

[x,spacestr] = pickpoints(aa,bb,Nvec(3),spaceopt);
y = yf(x);
warning off
[p,S,mu] = polyfit(x,y,Nvec(3)-1);
warning on
yp = polyval(p,xx,S,mu);
errpoly3 = errcompute(yp,yy);

[x,spacestr] = pickpoints(aa,bb,Nvec(4),spaceopt);
y = yf(x);
warning off
[p,S,mu] = polyfit(x,y,Nvec(4)-1);
warning on
yp = polyval(p,xx,S,mu);
errpoly4 = errcompute(yp,yy);

loglog(epvec,errvec(1,:),'r','LineWidth',3)
hold on
loglog(epvec,errvec(2,:),'g','LineWidth',3)
loglog(epvec,errvec(3,:),'b','LineWidth',3)
loglog(epvec,errvec(4,:),'c','LineWidth',3)
loglog(epvecd,errvecd(1,:),'-rx')
loglog(epvecd,errvecd(2,:),'-g+')
loglog(epvecd,errvecd(3,:),'-b^')
loglog(epvecd,errvecd(4,:),'-c^')
loglog(epvecd,errpoly1*ones(size(epvecd)),'--r')
loglog(epvecd,errpoly2*ones(size(epvecd)),'--g')
loglog(epvecd,errpoly3*ones(size(epvecd)),'--b')
loglog(epvecd,errpoly4*ones(size(epvecd)),'--c')
hold off
xlabel('\epsilon')
ylabel('Error')
ylim([10^-15 10])
xlim([min([epvec,epvecd]),max([epvec,epvecd])])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
%legend('N=3 (Direct)','N=5 (Direct)','N=9 (Direct)','N=17 (Direct)','N=3 (QR)','N=5 (QR)','N=9 (QR)','N=17 (QR)','Location','SouthEast')
legend('N=3','N=5','N=9','N=17','Location','SouthEast')
