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

M = 1024;
Nvec = [3200,6400,12800,25600,51200];

warning off

for N=Nvec
    b = rand(N,1);
    x = pickpoints(0,1,N);
    c = zeros(M,1);
    tic
    c = computeQReig(M,x,ep,alpha,b);
    timeMN = toc;
    tic
    phi = gqr_phi(1:M,x,ep,alpha);
    c = phi\b;
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
    c = zeros(M,1);
    tic
    c = computeQReig(M,x,ep,alpha,b);
    timeMN = toc;
    tic
    phi = gqr_phi(1:M,x,ep,alpha);
    c = phi\b;
    slowMN = toc;
%    rank(phi) % just so that I know
    fprintf('%d\t%d\t%g\t%g\n',M,N,timeMN,slowMN)
end

warning on

% save('speedtest.mat','slowMat','timeMat')
