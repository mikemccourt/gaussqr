basis = 2;
N = 50;
kernel = 3;
M = 3;
figure()
phigk = HSeigsolveqd(N,kernel,basis,M);
%phicheb = HSeigsolveqd(N,kernel,basis,M);
semilogy((1:N),phigk.eigvals,'bs')

hold on
semilogy((1:N),phicheb.eigvals,'r*')
title('eigenvals')
legend('qdgk','chebfun')
