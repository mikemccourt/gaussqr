% ex11
% This tests the speed of the fast QR factorization in 1D.
% Note that it takes advantage of the fast phi evaluation, so as to not
% bias the cost of the QR solve
%
% These results are invalid from a stability perspective, but are relevant
% in considering the speed required to make fast QR helpful
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.FAST_PHI_EVALUATION = 1;

% These values are arbitrary, because I don't really care about the
% solution, just the time it takes to find it
ep = 1;
alpha = 1;

M = 1024;
Nvec = [3200,6400,12800,25600,51200];

warning('off','MATLAB:rankDeficientMatrix')

for N=Nvec
    b = rand(N,1);
    x = pickpoints(0,1,N);
    tic
    c = computeQReig(M,x,ep,alpha,b);
    timeMN = toc;
    tic
    phi = gqr_phi(1:M,x,ep,alpha);
    coef = phi\b;
    slowMN = toc;
    fprintf('%d\t%d\t%g\t%g\n',M,N,timeMN,slowMN)
end

fprintf('\n\n')

N = 50000;
Mvec = [4,8,16,32,64,128,256,512,1024];
b = rand(N,1);

for M=Mvec
    r = gqr_roots(M,ep,alpha);
    x = pickpoints(-r,r,N);
    tic
    c = computeQReig(M,x,ep,alpha,b);
    timeMN = toc;
    tic
    phi = gqr_phi(1:M,x,ep,alpha);
    coef = phi\b;
    slowMN = toc;
    fprintf('%d\t%d\t%g\t%g\n',M,N,timeMN,slowMN)
end

warning('on','MATLAB:rankDeficientMatrix')