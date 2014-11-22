% StatFit2.m
% This is the second example of statistical data fitting involving
% topological data from Mt. Eden in New Zealand.  In this example we will
% study cross-validation schemes and the effect of splitting data into
% training, validation and test sets on the result of the fitting.
%
% The data is drawn from the main R library.  It consists of ground levels
% (heights) on an 87x61 grid.
%
% This is going to use HSSVD to perform the interpolation as needed - it
% will test first if the standard basis is acceptable and use the stable
% basis if it needs to.

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Make some plotting choices if extra output is needed
aux_plots = 1;

% Load the data into memory
%   latlong - Latitude/Longitude locations (scaled into [0,1]^2)
%   height - ground levels
load StatFit1_data.mat
x = latlong;
y = height;
N = size(x,1);

% Choose a kernel to fit to the data
rbf = @(e,r) exp(-(e*r).^2);

% Choose locations at which to make predictions for plotting
NN = 60;
xx = pick2Dpoints([-1 -1],[1 1],NN*ones(1,2));

% Choose fitting/testing fraction
% This is the proportion of data that will be used to fit the model
fitting_fraction = .7;

% Choose training fraction
% This is the proportion of the fitting data used to train the model
% This value must always be >=1/2 (leave half out) and <=1 (leave one out)
% Other fractions will be rounded to as close as possible
training_fraction = 1/5;
training_num = round(1/training_fraction);
clear x_valid y_valid x_train y_train
for k=1:training_num
    kth_valid_index = k:training_num:N;
    x_valid{k} = x(kth_valid_index,:);
    y_valid{k} = y(kth_valid_index);
    kth_train_index = setdiff(1:N,kth_valid_index);
    x_train{k} = x(kth_train_index,:);
    y_train{k} = y(kth_train_index);
end
pause

% Conduct the cross-validation training to choose an optimal kernel
epvec = logspace(-1,1,40);
if aux_plots
    h_mle = figure;
    cvvec = arrayfun(@(ep)StatFit1_CV_eval(ep,x_valid,y_valid,x_train,y_train),epvec);
    semilogx(epvec,cvvec)
end
ep_MLE = fminbnd(@(ep)StatFit1_CV_eval(ep,x_valid,y_valid,x_train,y_train),.1,10);

% Predict results from the kriging fit
DM = DistanceMatrix(x,x);
K = rbf(ep,DM);
DM_eval = DistanceMatrix(xx,x);
K_eval = rbf(ep,DM_eval);
yp = K_eval*(K\y);

% Plot the results
% We must reshape the data for a surface plot
X1 = reshape(xx(:,1),NN,NN);
X2 = reshape(xx(:,2),NN,NN);
YP = reshape(yp,NN,NN);
h = figure;
hold on
h_dots = plot3(x(:,1),x(:,2),y,'ok'); % The given data
h_surf = surf(X1,X2,YP,'edgecolor','none');
hold off
xlabel('latitude')
ylabel('longitude')
zlabel('Ferric Oxide percentage')

% Find the MLE ep
epvec = logspace(-1,1,40);
if aux_plots
    h_mle = figure;
    mlevec = arrayfun(@(ep)StatFit1_param_func(1,ep,DM,rbf,y),epvec);
    semilogx(epvec,mlevec)
end
ep_MLE = fminbnd(@(ep)StatFit1_param_func(1,ep,DM,rbf,y),.1,10);

% Find the MLE process variance, if needed 
if aux_plots
    h_cv = figure;
    pvvec = arrayfun(@(ep)StatFit1_param_func(2,ep,DM,rbf,y),epvec);
    semilogx(epvec,pvvec)
end
ep_PV = fminbnd(@(ep)StatFit1_param_func(2,ep,DM,rbf,y),.1,10);

% Find the LOOCV ep
if aux_plots
    h_cv = figure;
    cvvec = arrayfun(@(ep)StatFit1_param_func(3,ep,DM,rbf,y),epvec);
    semilogx(epvec,cvvec)
end
ep_CV = fminbnd(@(ep)StatFit1_param_func(3,ep,DM,rbf,y),.1,10);

% Find the average kriging variance at the eval points
if aux_plots
    h_cv = figure;
    kvvec = arrayfun(@(ep)StatFit1_param_func(4,ep,{DM,DM_eval},rbf,y),epvec);
    plotyy(epvec,kvvec,epvec,pvvec.*kvvec,@semilogx,@semilogx);
end
ep_KV = fminbnd(@(ep)StatFit1_param_func(4,ep,{DM,DM_eval},rbf,y),.1,10);

% Find some confidence intervals, primarily for the MLE
% First we need to compute the process variance associated with the maximum
% likelihood estimate
ep = .01
K = rbf(ep,DM);
K_eval = rbf(ep,DM_eval);
% sigma_MLE = 2*y'*(K\y)/N;
sigma_MLE = 1; % Just skip this for now
pfunc = sqrt((rbf(ep,0) - sum(K_eval'.*(K\K_eval'))));
P = reshape(pfunc,NN,NN);
h_CI = figure;
hold on
h_dots = plot3(x(:,1),x(:,2),ones(N,1),'ok'); % Data locations
h_CI = surf(X1,X2,sigma_MLE*P,'edgecolor','none');
hold off
xlabel('latitude')
ylabel('longitude')
zlabel('95% CI bound')
colorbar