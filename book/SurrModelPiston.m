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
GAUSSQR_PARAMETERS.ERROR_STYLE = 4;
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
rbfW6 = @(r) (1+10*r+39.4*r.^2+64*r.^3).*max(1-r,0).^7;
rbfIM = @(r) 1./(1+r.^2);
rbfGA = @(r) exp(-r.^2);

% Choose an anisotropic shape parameter vector
ep = 3.^(2:-.5:-1);%ep = ep(randperm(7));

% Choose points at which to study the surrogate model
% Also, choose a number of points at which to test the error
convN = 10;
Nvec = floor(logspace(2,4,convN));
Neval = 500;

% Also monitor the predicted/test value alignment
predtestmonitor = 0;
if predtestmonitor
    h = axes;
end

% Choose the kernels to be analyzed
rbfarr = {rbfM0,rbfM2,rbfM4};
labelarr = {'C0 Matern','C2 Matern','C4 Matern'};

% Analyze the convergence behavior of the surrogate model
h_waitbar = waitbar(0,'Initializing surrogate model','Visible','on');
errmat = zeros(length(rbfarr),length(Nvec));
k = 1;
for N=Nvec
    waitbar((k-1)/convN,h_waitbar,sprintf('Preparing model, N=%d',N));
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
    
    % Fit the surrogate model for all the different kernels
    waitbar((k-1)/convN,h_waitbar,sprintf('Computing model, N=%d',N));
    errmat(:,k) = cellfun(@(rbf) ...
                      errcompute(rbf(DMeval)*(rbf(DM)\y),yeval),...
                  rbfarr);
    
    % Plot the predicted/test values on axes together
    if predtestmonitor
        plot(h,rbfarr{3}(DMeval)*(rbfarr{3}(DM)\y),yeval,'.')
        xlabel('predicted data')
        ylabel('test data')
        pause
    end
    
    k = k + 1;
end
waitbar(1,h_waitbar,'Plotting piston convergence');

h_convergence = figure;
loglog(Nvec,errmat,'linewidth',2)
xlabel('number of input points')
ylabel('relative RMS 2-norm error')
legend(labelarr,'location','northeast')

close(h_waitbar)