close all
clear all
warning off
addpath('common_routines')
addpath('codes_for')
M = 2;                                  % space dimension 
n = 20;                                 % number of tracks
m = 50;                                 % number of points on tracks
[N dsites yy] = TrackData2D(n,m);       % generate N = n*m track data in 2D
neval = 30;                             % parameter for evaluation points
npu = floor(((N)./(4))^(1/M));          % parameter for PU centres
xx = linspace(0,1,n);                   % subdomains in one direction
[X Y] = meshgrid(xx,yy);                % patches centred at tracks
puctrs = [X(:) Y(:)];                   % define the PU centres
rbf_aniso = @(r) 1./sqrt(1+r.^2);       % define the RBF
% The test function
f = @(x) 0.75*exp(-((9*x(:,1)-2).^2+(9*x(:,2)-2).^2)/4) ...
            + 0.75*exp(-((9*x(:,1)+1).^2/49+(9*x(:,2)+1)/10)) ...
            + 0.5*exp(-((9*x(:,1)-7).^2+(9*x(:,2)-3).^2)/4) ...
            - 0.2*exp(-((9*x(:,1)-4).^2+(9*x(:,2)-7).^2));
rhs = f(dsites);                        % function values
h = 2;                                  % upper bound for the radius
r_min = 12;                             % minimum cardinality of patches
P1 = 4;                                 % number of testing radii
ep = [3 3];                             % guess for the shape parameters

tic
[epoints Pf] = PU_for(M,dsites,neval,npu,rbf_aniso,f,rhs,...
    r_min,h,P1,ep,puctrs);
toc