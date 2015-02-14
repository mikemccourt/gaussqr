% SurrModelPiston
% This script simulates a piston's operation within a car
% Given data about the
%      W  - piston weight
%      S  - surface area
%      V0 - initial gas volume
%      k  - spring/elasticity coefficient
%      P0 - atmospheric pressure
%      Ta - ambient temperature
%      T0 - incoming gas temperature
% The circulare motion of a piston is defined with a very complicated
% function that appears in the Surrogate Modeling chapter of the text.
% It is found from [Surjanovic and Bingham (2014)]
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Define the function of interest, which is available using the pickfunc
% function available in the GaussQR library
yf = pickfunc('piston_scaled');

% If the user has the haltonset file (in Statistics and Machine Learning
% Toolbox), use that to create the Halton points; if not, use the file
% provided in GaussQR
if exist('haltonset','file')
    point_generator = haltonset(7,'Skip',1);
    haltpts7D = @(N) net(point_generator,N);
else
    haltpts7D = @(N) haltonseq(Ntot,7);
end

% All of these kernels are defined with just r, allowing for computation
% with anisotropic distance functions
rbfM0 = @(r) exp(-r);
rbfM2 = @(r) (1+r).*exp(-r);
rbfM4 = @(r) (1+r+r.^2/3).*exp(-r);
rbfIM = @(r) 1./(1+r.^2);

% Choose an anisotropic shape parameter vector
% As a rough rule of thumb, as the domain gets bigger, ep gets smaller
ep = [1.6 .8 .4 .2 .1 .05 .025];
ep = [  1  1  1  1  1   1    1]/5;
ep = 2.^(2:-1:-4);
ep = 2.^(2:-.5:-1);

% Choose points at which to study the surrogate model
% Also, choose a number of points at which to test the error
Nvec = floor(logspace(2,3.3,15));
Neval = 500;

% Choose the kernels to be analyzed
rbfvec = {rbfM0,rbfM2,rbfM4,rbfIM};

% Analyze the convergence behavior of the surrogate model
errmat = zeros(length(rbfvec),length(Nvec));
k = 1;
for N=Nvec
    % Create the data points and test points
    % All the points will be generated at first and the test points will be
    % the first Neval points
    Ntot = N+Neval;
    xtot = haltpts7D(Ntot);
    xeval = xtot(1:Neval,:);
    x = xtot(Neval+1:end,:);
    
    % Evaluate the function to create the data and test values
    y = yf(x);
    yeval = yf(xeval);
    
    % Compute the anisotropic distance matrices
    DM = DistanceMatrix(x,x,ep);
    DMeval = DistanceMatrix(xeval,x,ep);
    
    % Fit the surrogate model and test it for different kernels
    for j=1:length(rbfvec)
        K = rbfvec{j}(DM);
        Keval = rbfvec{j}(DMeval);
        seval = Keval*(K\y);
        errmat(j,k) = errcompute(seval,yeval);
    end
    
    k = k + 1
end

h = figure;
loglog(Nvec,errmat,'linewidth',2)
xlabel('number of input points')
ylabel('2-norm error')
legend('C0 Matern','C2 Matern','C4 Matern','IMQ','location','northeast')