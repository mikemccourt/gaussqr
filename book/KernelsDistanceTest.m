function KernelsDistanceTest(opt)
% This is a time test for the various distance matrix options
% available in GaussQR.
% We primarily test symmetric distance matrices, but nonsymmetric
% problems may also be relevant.
% This is defined as a function, rather than a script, to provide access to
% the various distance functions within this file, but not outside
%
% You can use opt to conduct one of the two tests that appear in the book
%    opt=1 - Fig 4.1a
%    opt=2 - Fig 4.1b

num_runs = 10;
plots_on = 'off';

if opt==0
    % Just to test that DistanceMatrixRepmat can be called here
    DistanceMatrixRepmat(rand(2^8,5));
elseif opt==1
    d = 5;
    Nvec = [100,200,500,1000,2000,5000,10000];
    time_aniso = zeros(length(Nvec),1);
    time_repmat = zeros(length(Nvec),1);
    time_bsxfun = zeros(length(Nvec),1);
    time_pdist = zeros(length(Nvec),1);
    
    pts = haltonset(d);
    fprintf('Testing Halton points in %d dimensions\n',d)
    fprintf('num_pts repmat  bsxfun  aniso   pdist\n')
    k = 1;
    for N=Nvec
        z = net(pts,N);
        
        time_repmat(k) = DistanceMatrixTimeTest(@DistanceMatrixRepmat,num_runs,z);
        time_bsxfun(k) = DistanceMatrixTimeTest(@DistanceMatrixBsxfun,num_runs,z);
        time_aniso(k) = DistanceMatrixTimeTest(@(z)AnisotropicDistanceMatrix(z,ones(size(z))),num_runs,z);
        time_pdist(k) = DistanceMatrixTimeTest(@(z)squareform(pdist(z)),num_runs,z);

        fprintf('%7d\t%3.1e\t%3.1e\t%3.1e\t%3.1e\n',N,time_repmat(k),time_bsxfun(k),time_aniso(k),time_pdist(k))
        k = k + 1;
    end

    h = figure('visible',plots_on);
    loglog(Nvec,time_pdist,'--','linewidth',2)
    hold on
    loglog(Nvec,time_repmat,'r-.','linewidth',2)
    loglog(Nvec,time_aniso,'k','linewidth',2)
    hold off
    xlabel('number of points - N')
    ylabel('time (seconds)')
    legend('pdist','repmat','bsxfun','location','northwest')
elseif opt==2
    dvec = [2.^(0:9),1000];
    N = 10000;

    time_aniso = zeros(length(dvec),1);
    time_repmat = zeros(length(dvec),1);
    time_bsxfun = zeros(length(dvec),1);
    time_pdist = zeros(length(dvec),1);
    
    fprintf('Testing %d Halton points\n',N)
    fprintf('num_dim repmat  bsxfun  aniso   pdist\n')
    k = 1;
    for d=dvec
        pts = haltonset(d);
        z = net(pts,N);
        
        time_repmat(k) = DistanceMatrixTimeTest(@DistanceMatrixRepmat,num_runs,z);
        time_bsxfun(k) = DistanceMatrixTimeTest(@DistanceMatrixBsxfun,num_runs,z);
        time_aniso(k) = DistanceMatrixTimeTest(@(z)AnisotropicDistanceMatrix(z,ones(size(z))),num_runs,z);
        time_pdist(k) = DistanceMatrixTimeTest(@(z)squareform(pdist(z)),num_runs,z);

        fprintf('%7d\t%3.1e\t%3.1e\t%3.1e\t%3.1e\n',d,time_repmat(k),time_bsxfun(k),time_aniso(k),time_pdist(k))
        k = k + 1;
    end

    h = figure('visible',plots_on);
    loglog(dvec,time_pdist,'--','linewidth',2)
    hold on
    loglog(dvec,time_repmat,'r-.','linewidth',2)
    loglog(dvec,time_aniso,'k','linewidth',2)
    hold off
    xlabel('number of dimensions - d')
    ylabel('time (seconds)')
    legend('pdist','repmat','bsxfun','location','northwest')
else
    error('You must choose opt=1 or 2')
end
    if plots_on=='off'
        saveas(h,'time_test','png')
    end
end


function DM = DistanceMatrixRepmat(z)
N = size(z,1);
sz2 = sum(z.^2,2);
DM = sqrt(repmat(sz2,1,N) + repmat(sz2',N,1) - (2*z)*z');
end

function DM = DistanceMatrixBsxfun(z)
sz2 = sum(z.^2,2);
DM = sqrt(bsxfun(@plus,sz2,sz2') - (2*z)*z');
end

% This is a simplified form of the DistanceMatrix
% function from Program 4.5, which allows for only
% symmetric distance matrices for the purposes of
% this test.
function DM = AnisotropicDistanceMatrix(z,e)
ze = bsxfun(@times,z,e);
sze2 = sum(ze.*ze,2);
DM = sqrt(bsxfun(@plus,sze2,sze2') - (2*ze)*ze');
end

function t = DistanceMatrixTimeTest(DMfunc,num_runs,z)
timevec = zeros(1,num_runs);
for m=1:num_runs
    clear DM
    tic
    DM = DMfunc(z);
    timevec(m) = toc;
end
t = mean(timevec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For reference, results obtained for the book are
% displayed below.
%%%
% Testing Halton points in 5 dimensions
% num_pts repmat  bsxfun  aniso   pdist
%     256 1.8e-03 1.1e-03 1.3e-03 1.8e-03
%     512 4.0e-03 3.7e-03 3.8e-03 6.8e-03
%    1024 1.3e-02 1.0e-02 1.0e-02 3.0e-02
%    2048 6.7e-02 3.9e-02 4.0e-02 1.1e-01
%    4096 2.6e-01 1.5e-01 1.5e-01 4.4e-01
%    8192 9.9e-01 5.9e-01 5.9e-01 1.7e+00
%   16384 3.9e+00 2.3e+00 2.3e+00 6.8e+00
%%%
% Testing 10000 Halton points
% dim repmat  bsxfun  aniso   pdist
%    1 1.5e+00 8.6e-01 8.6e-01 2.5e+00
%    2 1.5e+00 8.6e-01 8.7e-01 2.5e+00
%    4 1.5e+00 8.5e-01 8.6e-01 2.4e+00
%    8 1.5e+00 8.5e-01 9.0e-01 2.5e+00
%   16 1.5e+00 8.6e-01 8.6e-01 3.2e+00
%   32 1.5e+00 9.0e-01 9.0e-01 4.2e+00
%   64 1.6e+00 9.8e-01 9.9e-01 6.1e+00
%  128 1.8e+00 1.1e+00 1.1e+00 1.1e+01
%  256 2.2e+00 1.5e+00 1.5e+00 2.4e+01
%  512 2.9e+00 2.2e+00 2.2e+00 4.4e+01
% 1024 4.2e+00 3.6e+00 3.7e+00 8.8e+01

