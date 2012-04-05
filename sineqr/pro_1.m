clear all;
clc

sigmavecd = logspace(-3,1,40);  %epvecd = logspace(-2,1,40);
                                %sigmavec = logspace(-2,0.3,100);  %epvec = logspace(-2,0.3,100);
Nvec = [10,30,90];
NN = 100;

spaceopt = 'cheb';
fopt = 'sinh';

[yf,fstr] = pickfunc(fopt,1);

aa = -3;bb = 3;
xx = linspace(aa,bb,NN)';
yy = yf(xx);
%errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(sigmavecd));
                 %alpha = .65;

 n = 1;
for N=Nvec
%     k = 1;
%     for ep=epvec
%         [x,spacestr] = pickpoints(aa,bb,N,spaceopt,ep);
%         rbfqrOBJ = rbfqr_solve(x,yf(x),ep,alpha);
%         yp = rbfqr_eval(rbfqrOBJ,xx);
%         errvec(n,k) = errcompute(yp,yy);
%         k = k+1;
%     end

    k = 1;
    for sigma=sigmavecd
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,sigma);
        y = yf(x);
        K = exp(-sigma*abs(repmat(x,1,N)-repmat(x',N,1))); % C0
        %K = (1+sigma*abs(repmat(x,1,N)-repmat(x',N,1))).*exp(-sigma*abs(repmat(x,1,N)-repmat(x',N,1))); % C2
        %K = (3+3*sigma*abs(repmat(x,1,N)-repmat(x',N,1))+...
        %sigma^2*abs(repmat(x,1,N)-repmat(x',N,1)).^2).*exp(-sigma*abs(repmat(x,1,N)-repmat(x',N,1))); % C4
        c(n,k)=cond(K); % conditioning
        warning off % We know it's bad
%         beta = pinv(K)*y;
        beta = K\y;
        warning on
        yp = exp(-sigma*abs(repmat(x',NN,1)-repmat(xx,1,N)))*beta; % C0
        %yp =((1+sigma*abs(repmat(x',NN,1)-repmat(xx,1,N))).*exp(-sigma*abs(repmat(x',NN,1)-repmat(xx,1,N))))*beta; % C2
        %yp = ((3+3*sigma*abs(repmat(x',NN,1)-repmat(xx,1,N))+...
        %     sigma^2*abs(repmat(x',NN,1)-repmat(xx,1,N)).^2).*exp(-sigma*abs(repmat(x',NN,1)-repmat(xx,1,N))))*beta; % C4
        errvecd(n,k) = abs(max(yy-yp)); %errvecd(n,k) = errcompute(yp,yy);
        k = k+1;
    end
    n = n+1;
end

% [x,spacestr] = pickpoints(aa,bb,Nvec(1),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(1)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly1 = max(abs(yy-yp))   %errpoly1 = errcompute(yp,yy);
% 
% [x,spacestr] = pickpoints(aa,bb,Nvec(2),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(2)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly2 = max(abs(yy-yp))   %errcompute(yp,yy);
% 
% [x,spacestr] = pickpoints(aa,bb,Nvec(3),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(3)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly3 = max(abs(yy-yp))   %errcompute(yp,yy);

loglog(sigmavecd,errvecd(1,:),'-bx')
hold on
loglog(sigmavecd,errvecd(2,:),'-g+')
loglog(sigmavecd,errvecd(3,:),'-r^')
% loglog(sigmavec,errvec(1,:),'b','LineWidth',3)
% loglog(sigmavec,errvec(2,:),'g','LineWidth',3)
% loglog(sigmavec,errvec(3,:),'r','LineWidth',3)
% loglog(sigmavecd,errpoly1*ones(size(sigmavecd)),'--b')
% loglog(sigmavecd,errpoly2*ones(size(sigmavecd)),'--g')
% loglog(sigmavecd,errpoly3*ones(size(sigmavecd)),'--r')
hold off
xlabel('\sigma')
ylabel('max abs error') % ylabel('average error')
%ylim([10^-15 10])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=10 (Direct)','N=30 (Direct)','N=90 (Direct)','Location','Best') %legend('N=10 (Direct)','N=20 (Direct)','N=30 (Direct)','N=10 (QR)','N=20 (QR)','N=30 (QR)','Location','SouthEast')
