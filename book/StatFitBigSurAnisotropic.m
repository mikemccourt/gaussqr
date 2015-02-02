% StatFitBigSur
% This example loads the Big Sur data set from 
%           Interpolation of track data with radial basis methods
%           R. E. Carlson and T. A. Foley
%           Comput. Math. Appl., 24:27-34, 1992.
% We can access this from the GaussQR data directory
% Based on our isotropic results, we will try to use an anisotropic kernel
% to fit this data
% We will use the anisotropic C2 Matern kernel

% Load the data in from the repository as needed
% This data comes with sizes Nx and Ny where length(x) = Nx*Ny
gqr_downloaddata('BigSur_data.mat')
load BigSur_data
[x,t] = rescale_data(locations,temperatures);
N = size(x,1);

% Define the anisotropic C2 Matern kernel
rbf = @(r) (1+r).*exp(-r);

% Choose some points at which to evaluate the power function
% Make sure they are in the convex hull
cind = convhull(x(:,1),x(:,2));
xeval = pick2Dpoints(-1,1,6,'halt');
xeval = xeval(inpolygon(xeval(:,1),xeval(:,2),x(cind,1),x(cind,2)),:);

% Consider a range of epsilon values for each dimension
% ep2dvec = [eparr1 , eparr2 , ... ]
Nep = 25;
ep2dvec = 10.^pick2Dpoints(-1,0,Nep)';

% Loop through the 2D epsilon values for each eparr
gwvec = cellfun(@(eparr)CGW_eval(x,t,rbf,eparr',xeval),num2cell(ep2dvec,1));

% Create a surface plot of the data
h_ep = figure;
E1 = reshape(ep2dvec(1,:),Nep,Nep);
E2 = reshape(ep2dvec(2,:),Nep,Nep);
GW = reshape(gwvec,Nep,Nep);
surf(E1,E2,GW)
xlabel('$\varepsilon_1$','interpreter','latex')
ylabel('$\varepsilon_2$','interpreter','latex')
zlabel('C_{GW}')
set(get(h_ep,'currentaxes'),'xscale','log')
set(get(h_ep,'currentaxes'),'yscale','log')
set(get(h_ep,'currentaxes'),'zscale','log')
zlim([8 100])

% Find the minimum of the tested values
[~,imin] = min(gwvec);
epguess = ep2dvec(:,imin)';

% Compute the minimum with fminunc
minval = fminunc(@(eparr)CGW_eval(x,t,rbf,exp(eparr),xeval),log(epguess));
epopt = exp(minval);

% Create a surface plot of the optimal result
Nplot = 55;
xplot = pick2Dpoints(-1,1,Nplot);
splot = rbf(DistanceMatrix(xplot,x,epopt))*(rbf(DistanceMatrix(x,x,epopt))\t);
X = reshape(xplot(:,1),Nplot,Nplot);
Y = reshape(xplot(:,2),Nplot,Nplot);
S = reshape(splot,Nplot,Nplot);

% Only plot points in the convex hull
% Probably unnecessary, but looks cool
S(not(inpolygon(X,Y,x(cind,1),x(cind,2)))) = NaN;

h_pred = figure;
h_surf = surf(X,Y,S,'edgealpha',.6);
hold on
plot3(x(:,1),x(:,2),t,'or','linewidth',2)
hold off
xlabel('x')
ylabel('y')
zlabel('temperature')
view([-34.5 42])