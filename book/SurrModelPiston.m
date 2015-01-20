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

% Define the function of interest, which is available using the pickfunc
% function available in the GaussQR library
yf = pickfunc('piston');

% Define the appropriate domains for the piston function
lb = [30,.005,.002,1000,90000 ,290,340];
ub = [60,.020,.010,5000,110000,296,360];

% Choose a total number of points and a fraction of them to test on
Ntot = 2000;
test_frac = .5;
Neval = round(test_frac*Ntot);
N = Ntot-Neval;

% Define the N points that are relevant for the problem
% First create a bunch of points and then separate into samples and test
% points to determine if we have done a good job
% If the user has the haltonset file (in Statistics and Machine Learning
% Toolbox), use that to create the Halton points; if not, use the file
% provided in GaussQR
if exist('haltonset','file')
    point_generator = haltonset(7,'Skip',1);
    xhalton = net(point_generator,Ntot);
else
    xhalton = haltonseq(Ntot,7);
end
% Scale the points from [0,1]^7 to the appropriate domain
xtot = bsxfun(@plus,bsxfun(@times,ub-lb,xhalton),lb);

% Separate the locations into data sites and test sites
x = xtot(1:N,:);
xeval = xtot(N+1:end,:);

% Evaluate the function to create the data and test values
y = yf(x);
yeval = yf(xeval);

% All of these kernels are defined with just r, allowing for computation
% with anisotropic distance functions
RBFM0 = @(r) exp(-r);
RBFM2 = @(r) (1+r).*exp(-r);
RBFG = @(r) exp(-r.^2);
RBFIM = @(r) 1./(1+r.^2);

% Choose a kernel with which to perform computations
rbf = RBFM0;

% Choose an anisotropic shape parameter vector
% As a rough rule of thumb, as the domain gets bigger, ep gets smaller
ep = [1/20 8 11 1/20 1/50 1/2 1/4];

% Compute the approximation and evaluate the error
DM = DistanceMatrix(x,x,ep);
DMeval = DistanceMatrix(xeval,x,ep);
K = rbf(DM);
Keval = rbf(DMeval);
seval = Keval*(K\y);
err = errcompute(seval,yeval)