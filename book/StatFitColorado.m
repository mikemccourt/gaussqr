% StatFitColorado
% The purpose of this example is to study the use of cross-validation in
% choosing good epsilon and mu parameters in a smoothing spline fit.
%
% We use the noisy data from the Colorado Animas River locations, which
% measures the parts per million of Lanthanum in the sample.
%
% Some of the locations are duplicates (I'm not sure why)
% Where those duplicates occur, only the first listed data is used.
% Also, any NaNs in the data are removed
% As a result, only 81 data locations are used here.
%
% Note that, although the data locations are called latlong, they are
% actually stored as [longitude,latitude] to match our normal orientation
% of x & y in mathematics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, we must grab the data, which is stored on the GaussQR server
% We thank the USGS for making this data available to us
% http://pubs.usgs.gov/ds/595/
gqr_downloaddata('AnimasRiver_data.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Once you have the data, load it into memory
% This will create data locations latlong and values Lsppm
load AnimasRiver_data

% We need to clean the data for NaN and inf
% We also choose to cut any duplicate data locations
% The data is then scaled into [-1,1] to make our lives easier
[x,y,latlong_shift,latlong_scale] = rescale_data(latlong,Lappm,1);
N = size(x,1);

% Choose the C0 Matern to fit the data
rbf = @(e,r) exp(-(e*r));

% Choose a range of epsilon and mu parameters to consider
epvec = logspace(-1,1,20);
muvec = logspace(-6,-2,15);

% Choose locations at which to make predictions
Neval = 60;
xeval = pick2Dpoints(-1,1,Neval);

% Define a function which makes the Kriging predictions
sf = @(ep,mu,x,y,xeval) rbf(ep,DistanceMatrix(xeval,x))*...
                    (rbf(ep,DistanceMatrix(x,x)+mu*eye(length(y)))\y);

% Define the LOOCV computation as a function of epsilon
cv = @(ep,mu) crossval('mse',x,y,...
                       'Predfun',@(x,y,xeval) sf(ep,mu,x,y,xeval),...
                       'leaveout',1);
                
% Evaluate the cross-validation for our epsilon values
cvepvec = arrayfun(@(ep) cv(ep,0),epvec);

% Plot the cross-validation MSE residual
h_cv = figure;
semilogx(epvec,cvepvec,'linewidth',2)
xlabel('$\varepsilon$','interpreter','latex')
ylabel('mean square LOOCV residual')
ylim([350 400])

% Find the LOOCV minimizing epsilon and evaluate with it
epopt = fminbnd(@(ep) cv(ep,0),1,2);
seval = sf(ep,0,x,y,xeval);

% Plot the results
X1 = reshape(xeval(:,1),Neval,Neval);
X2 = reshape(xeval(:,2),Neval,Neval);
S = reshape(seval,Neval,Neval);
h = figure;
h_surf = surf(X1,X2,S,'edgealpha',.6);
hold on
h_dots = plot3(x(:,1),x(:,2),y,'ok'); % The given data
hold off
view([-152 42])
xlabel('longitude')
ylabel('latitude')
zlabel('Ferric Oxide percentage')

% Simplest possible plotting of the data
h_griddata = figure;
[X,Y] = meshgrid(pickpoints(-1,1,50),pickpoints(-1,1,50));
L = griddata(x(:,1),x(:,2),y,X,Y);
surf(X,Y,L),hold on,plot3(x(:,1),x(:,2),y,'og'),hold off
view([-152 42])