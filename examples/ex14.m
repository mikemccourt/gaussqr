% ex14
% This tests the speed of the fast QR factorization in 1D.

% These are the values I'm interested in testing
% Note M is stored as a percentage of N
Nvec = 50*2.^(0:7);
Mvec = .1*(1:8);

% These values are arbitrary, because I don't really care about the
% solution, just the time it takes to find it
ep = 1;
alpha = 1;

timeMat = zeros(length(Nvec),length(Mvec));

n = 1;
for N=Nvec
    b = rand(N,1);
    x = pickpoints(0,1,N);
    m = 1;
    for Mp=Mvec
        M = floor(Mp*N);
        c = zeros(M,1);
        tic
        c = computeQReig(M,x,ep,alpha,b);
        timeMat(n,m) = toc;
        fprintf('%d\t',M)
        m = m + 1;
    end
    fprintf('  |  %d\n',N)
    n = n + 1;
end