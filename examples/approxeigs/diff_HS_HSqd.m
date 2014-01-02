clear bi ed kernel basis eig eigqd phi phiqd terror terrorqd;
bi = 10;
ed = 50;
basis = 2;
ith = 8;
kernel = 2;
epsilon = 3;
terror = zeros(1,ed - bi + 1);
terrorqd = zeros(1,ed - bi + 1);
for N = bi : ed
    phi = HSeigsolve(N,kernel,basis,epsilon);
    %qdopts.RelTol = 1.e-13;
    %qdopts.AbsTol = 1.e-16;
    phiqd = HSeigsolveqd(N,kernel,basis,2,epsilon); 
    eig = phi.eigvals;
    eigqd = phiqd.eigvals;
    terror(N - bi + 1) =abs( (1/(ith^2*pi^2)-eig(ith))/(1/ith^2*pi^2));
    terrorqd(N - bi + 1) =abs((1/(ith^2*pi^2)-eigqd(ith))/(1/ith^2*pi^2));
end
hold off;
semilogy((bi : ed), terror,'sb');
hold on;
semilogy((bi : ed), terrorqd,'*r');
legend('relative error','qd relative error')
xlabel('# collocation points')
ylabel('error')
title(sprintf('eval(%d), basis(%d),epsilon(%d)',ith,basis,epsilon))
hold off
clear bi ed basis eig eigqd phi phiqd terror terrorqd;

