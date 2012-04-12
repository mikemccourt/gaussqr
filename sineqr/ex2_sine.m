% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
rbfsetup

sigmavecd = logspace(-1,2,40);
sigmavec =  logspace(-1,2,100); % [1e-3, 1/3, 1, 3, 100];
Nvec = [10,20,30];
NN = 100;

spaceopt = 'halton';
fopt = 'sin';

L = 1;
% [yf,fstr] = pickfunc(fopt,1);
yf = @(x) x.*(L-x);
%yf = @(x) x.*(L-x).*sqrt(x);

aa = 0; bb = L;
xx = linspace(aa,bb,NN)';
yy = yf(xx);
errvec = zeros(length(Nvec),length(sigmavec));
errvecd = zeros(length(Nvec),length(sigmavecd));
%alpha = .65;

i = 1;
for N=Nvec
       k = 1;
       for sigma=sigmavec
           [x,spacestr] = pickpoints(aa,bb,N,spaceopt,sigma);
           y = yf(x);
           M = ceil(2.1*N);
           n = 1:M;
          S = sqrt(2/L)*sin(pi*x*n/L); % size(S)  K = Se*D*Sn'  psi = Se*(I;D2Sn2'Sn1^-TD1^-1)
          [Q,R] = qr(S);
             R1 = R(:,1:N);
             R2 = R(:,N+1:end);
          lambda = ((pi*n/L).^2+sigma^2).^(-1);
          D = diag(lambda); % size(D)
          D1 = diag(lambda(1:N));
          D2 = diag(lambda(N+1:end));
%           b = (Q*R1)\y;
%           beta = (eye(N)+inv(R1)*R2*D2*R2'*inv(R1')*inv(D1))\b;
          Rbar = D2*(R2'/R1')/D1;
          beta = (S*[eye(N); Rbar])\y;
          SS = sqrt(2/L)*sin(pi*xx*n/L);                        % <----- ??
          yp = (SS*[eye(N); Rbar])*beta;     % <----- ??
%          rbfqrOBJ = rbfqr_solve(x,yf(x),sigma,alpha);
%          yp = rbfqr_eval(rbfqrOBJ,xx);
          errvec(i,k) = errcompute(yp,yy);
         k = k+1;
       end
 
     k = 1;
    for sigma=sigmavecd
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt,sigma);
        y = yf(x);
        MINVAL = min(repmat(x,1,N),repmat(x',N,1));
        MAXVAL = max(repmat(x,1,N),repmat(x',N,1));
        K = sinh(sigma*MINVAL).*sinh(sigma*(L-MAXVAL))/(sigma*sinh(L*sigma));
        warning off % We know it's bad
        beta = K\y;
        warning on
        MINGR = min(repmat(x',NN,1),repmat(xx,1,N));
        MAXGR = max(repmat(x',NN,1),repmat(xx,1,N));
        yp = (sinh(sigma*MINGR).*sinh(sigma*(L-MAXGR))/(sigma*sinh(L*sigma)))*beta;
        errvecd(i,k) = abs(max(yy-yp));%errcompute(yp,yy);
        errvecd(i,k) = errcompute(yp,yy);
        k = k+1;
    end
    i = i+1;
end

% [x,spacestr] = pickpoints(aa,bb,Nvec(1),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(1)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly1 = errcompute(yp,yy);
% 
% [x,spacestr] = pickpoints(aa,bb,Nvec(2),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(2)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly2 = errcompute(yp,yy);
% 
% [x,spacestr] = pickpoints(aa,bb,Nvec(3),spaceopt);
% y = yf(x);
% warning off
% [p,S,mu] = polyfit(x,y,Nvec(3)-1);
% warning on
% yp = polyval(p,xx,S,mu);
% errpoly3 = errcompute(yp,yy);

loglog(sigmavecd,errvecd(1,:),'-bx')
hold on
loglog(sigmavecd,errvecd(2,:),'-g+')
loglog(sigmavecd,errvecd(3,:),'-r^')
loglog(sigmavec,errvec(1,:),'b','LineWidth',3)
loglog(sigmavec,errvec(2,:),'g','LineWidth',3)
loglog(sigmavec,errvec(3,:),'r','LineWidth',3)
% loglog(sigmavecd,errpoly1*ones(size(sigmavecd)),'--b')
% loglog(sigmavecd,errpoly2*ones(size(sigmavecd)),'--g')
% loglog(sigmavecd,errpoly3*ones(size(sigmavecd)),'--r')
hold off
xlabel('\sigma')
ylabel('average error')
%ylim([10^-15 10])
ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(fstr,ptsstr,spacestr))
legend('N=10 (Direct)','N=20 (Direct)','N=30 (Direct)','N=10 (QR)','N=20 (QR)','N=30 (QR)', 'Location', 'Best');
