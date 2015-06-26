
%test for cheb fun function
clear N L ptspace K_F H_mat x z Int_Kh errc errgk l i p f fc Truth phi phicc phigk;
N = 40;
L = 1;
K = 1;
B = 3;
l = 1;
errcc = [];
errgk = [];
phi = HSeigsolve(N,K,B);
phicc = HSeigsolveqd(N,K,B,3);
phigk = HSeigsolveqd(N,K,B,2);
A = phi.A;
Acc = phicc.K;
Agk = phigk.K;
errcc =abs((A(l,:) - Acc(l,:))./A(l,:));
errgk =abs((A(l,:) - Agk(l,:))./A(l,:));
semilogy((1:N),errcc,'r*');
hold on;
semilogy((1:N),errgk,'bs');
hold on;
legend('errcc','errgk');
title(sprintf('text_chebfun r=%d N = %d',l,N));
xlabel('basis function');
ylabel('bas error');

