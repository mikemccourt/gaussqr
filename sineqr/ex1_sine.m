clear all;
%clc

format short e

%N = 40;
N = 100;
%Mvec = [40, 80, 160, 320];
Mvec = [10, 100, 200, 400, 800];
sigmavec = [1e-3 1/3, 1, 3, 100];
L = 1;

spaceopt = 'cheb';
% rbf = @(sigma,x,z) x*ones(size(z))'; % x,z are column vectors

%r = DistanceMatrix(x,y);

aa = 0; bb = L;

[x,spacestr] = pickpoints(aa,bb,N,spaceopt);

m = 1;
for M = Mvec
    k = 1;
    for sigma = sigmavec
        n = 1:M;
        S = sqrt(2/L)*sin(pi*x*n/L); % size(S)
        lambda = ((pi*n/L).^2+sigma^2).^(-1);
        D = diag(lambda); % size(D)
        % S*D*S'

        % KERNEL:
%        DM = DistanceMatrix(x,x);
%        K = rbf(sigma,r);
        MINVAL = min(repmat(x,1,N),repmat(x',N,1));
        MAXVAL = max(repmat(x,1,N),repmat(x',N,1));
        K = sinh(sigma*MINVAL).*sinh(sigma*(L-MAXVAL))/(sigma*sinh(L*sigma));

        errvec(m,k) = norm(K-S*D*S');
        k = k+1;
    end
    m = m+1;
end
ERROR = errvec

loglog(sigmavec,errvec(1,:),'-bx')
hold on
loglog(sigmavec,errvec(2,:),'-g+')
loglog(sigmavec,errvec(3,:),'-r^')
loglog(sigmavec,errvec(4,:),':kd')
loglog(sigmavec,errvec(5,:),'-.mo')
hold off
xlabel('\sigma')
ylabel('norm error') % ylabel('average error')
%ylim([10^-15 10])
ptsstr = strcat('x\in[',num2str(aa),',',num2str(bb),'],');
title(strcat(ptsstr,spacestr))
legend('M=10','M=100','M=200','M=400','M=800','Location','Best')
