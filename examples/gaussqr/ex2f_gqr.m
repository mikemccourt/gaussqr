% ex2
% This example compares the convergence behavior of RBF-QR and RBF-QRr
rbfsetup
global GAUSSQR_PARAMETERS;
GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED = 30;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = 90;

Nvec = [8,16,32,64,128,256,512];
NN = 1000;
ep=.1;
alpha=1;
alpha = gqr_alphasearch(ep,-1,1);
Mfrac = .1:.2:.9;
    
%yf = @(x) x + 1./(1+x.^2);
%fstr = 'y(x) = x + 1/(1+x^2)';
%yf = @(x) x.^3-3*x.^2+2*x+1 + 1e-10*cos(10*x);
%fstr = 'y(x) = x^3-3x^2+2x+1 + 10^{-10}cos(10x)';
%yf = @(x) x;% + 0.001*cos(10*x);
%fstr = 'y(x) = x';% + cos(10x)/1000';
yf = @(x) 0.75*exp(-((9*(x+1)/2-2).^2)/4)+0.75*exp(-((9*(x+1)/2+1).^2/49))+0.5*exp(-((9*(x+1)/2-7).^2)/4)-0.2*exp(-((9*(x+1)/2-4).^2));
fstr = '"Franke"';
%yf = @(x) tanh(9*(x-1))+1;
%fstr = 'tanh(9(x-1))+1';

xx = pickpoints(-1,1,NN);
yy = yf(xx);
errvec = zeros(size(Nvec));
errvecr = zeros(length(Mfrac),size(Nvec,2));
errvecd = zeros(size(Nvec));
errpoly = zeros(size(Nvec));

status = 'Performing RBF-QRr'
for i=1:length(Mfrac)
    k = 1;
    for N=Nvec
        x = pickpoints(-1,1,N,'cheb');
        y = yf(x);
        M = round(N*Mfrac(i));
        GQR = gqr_rsolve(x,y,ep,alpha,M);
        yp = gqr_eval(GQR,xx);
        errvecr(i,k) = errcompute(yp,yy);
        k = k+1;
    end
end

status = 'Performing RBF-QR'
k = 1;
for N=Nvec
    x = pickpoints(-1,1,N,'cheb');
    y = yf(x);
    GQR = gqr_solve(x,y,ep);
    GQR.Marr(end)
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    k = k+1;
end

% status = 'Performing RBF-Direct'
% k = 1;
% for N=Nvec
%     x = pickpoints(-1,1,N,'cheb');
%     y = yf(x);
%     K = exp(-ep^2*(repmat(x,1,N)-repmat(x',N,1)).^2);
%     warning off % I know it's bad
%     beta = K\y;
%     warning on
%     yp = exp(-ep^2*(repmat(x',NN,1)-repmat(xx,1,N)).^2)*beta;
%     errvecd(k) = errcompute(yp,yy);
%     k = k+1;
% end

% status = 'Polynomial computation'
% k = 1;
% for N=Nvec
%     x = pickpoints(-1,1,N,'cheb');
%     y = yf(x);
%     warning off
%     [ppoly,spoly,mupoly] = polyfit(x,y,N);
%     yp = polyval(ppoly,xx,spoly,mupoly);
%     warning on
%     errpoly(k) = errcompute(yp,yy);
%     k = k+1;
% end

% loglog(Nvec,errvecd,'r','linewidth',2);
% hold on
% loglog(Nvec,errvec,'g','linewidth',2);
% loglog(Nvec,errvecr,'b','linewidth',2);
% loglog(Nvec,errpoly,'c','linewidth',2);
% xlabel('N');
% ylabel('absolute error');
% hold off
% ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
% title(strcat(fstr,ptsstr,spacestr))
% legend('Direct','QR','Regression','polyfit','Location','SouthWest')

loglog(Nvec,errvec,'g','linewidth',2);
hold on
loglog(Nvec,errvecr,'linewidth',2);
xlabel('N');
ylabel('absolute error');
hold off
ptsstr=strcat(', x\in[-1,1],');
title(strcat(fstr,ptsstr,'cheb'))
legend('QR','M=.1N','M=.3N','M=.5N','M=.7N','M=.9N','Location','SouthWest')
