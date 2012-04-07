% ex14
% This tests the speed of the fast QR factorization in 1D.
% Note that it takes advantage of the fast phi evaluation, so as to not
% bias the cost of the QR solve
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.FAST_PHI_EVALUATION = 1;

% These are the values I'm interested in testing
% Note M is stored as a percentage of N
Nvec = 200*2.^(0:7);
Mvec = .01*2.^(0:5);

% These values are arbitrary, because I don't really care about the
% solution, just the time it takes to find it
ep = 1;
alpha = 1;

timeMat = zeros(length(Nvec),length(Mvec));
slowMat = zeros(length(Nvec),length(Mvec));

% Ignore because we're not concerned about the solution accuracy
warning off

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
        tic
        phi = rbfphi(1:M,x,ep,alpha);
        c = phi\b;
        slowMat(n,m) = toc;
        fprintf('%d\t%d\t%g\t%g\n',M,N,timeMat(n,m),slowMat(n,m))
        m = m + 1;
    end
    n = n + 1;
end

warning on

save('speedtest.mat','slowMat','timeMat')
