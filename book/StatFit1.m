% StatFit1.m
% This is the initial example for statistical data fitting
% The data has been drawn from
%   2011 - U.S. Geological Survey Data Series 595
% It involves the 90 Animas river locations given in that report
% Some of the locations are duplicates (I'm not sure why)
% Where those duplicates occur, only the first listed data is used.
% As a result, only 81 data locations are used here.
% Plots will include confidence intervals where relevant.
%
% This example creates a variety of surfaces for the given data,
% including different shape as chosen by our parameterization
% methods (MLE and LOOCV)
% We also will use this example to demonstrate the difficulty of
% dealing with multiscale data, where some data is very closely
% clustered and other data is spread out.

% Make some plotting choices if extra output is needed
aux_plots = 1;

% Load the data into memory
%   latlong - Latitude/Longitude locations
%   FOpct - Ferric Oxide percentage in sample
load StatFit1_data.mat
N = size(latlong,1);
y = FOpct;

% We rescale the data locations to [-1,1] for simplicity
% This is not required, but helps computation
latlong_shift = min(latlong);
latlong_scale = max(latlong) - min(latlong);
x = 2*(latlong - ones(N,1)*latlong_shift)./(ones(N,1)*latlong_scale) - 1;

% Choose a kernel to fit to the data
% Is DistanceMatrix returning complex numbers?
ep = 1;
rbf = @(e,r) exp(-(e*r).^2);
rbf = @(e,r) exp(-(e*r));

% Choose locations at which to make predictions
NN = 50;
xx = pick2Dpoints([-1 -1],[1 1],NN*ones(1,2));

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
h_stuf = surf(X1,X2,YP,'edgecolor','none');
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

% Find the LOOCV ep
epvec = logspace(-1,1,40);
if aux_plots
    h_cv = figure;
    cvvec = arrayfun(@(ep)StatFit1_param_func(2,ep,DM,rbf,y),epvec);
    semilogx(epvec,cvvec)
end
ep_CV = fminbnd(@(ep)StatFit1_param_func(2,ep,DM,rbf,y),.1,10);

% Find some confidence intervals