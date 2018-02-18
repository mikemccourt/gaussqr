close all
warning off
addpath('common_routines')
addpath('codes_op')
M = 2;                                  % space dimension 
n = 20;                                 % number of tracks
m = 50;                                 % number of points on tracks
[N, dsites, yy] = TrackData2D(n,m);       % generate N = n*m track data in 2D
neval = 30;                             % parameter for evaluation points
npu = floor(((N)./(4))^(1/M));          % parameter for PU centres
xx = linspace(0,1,n);                   % subdomains in one direction
[X, Y] = meshgrid(xx,yy);                % patches centred at tracks
puctrs = [X(:) Y(:)];                   % define the PU centres
rbf_aniso = @(r) 1./sqrt(1+r.^2);       % define the RBF

% The test function
f = @(x) 0.75*exp(-((9*x(:,1)-2).^2+(9*x(:,2)-2).^2)/4) ...
            + 0.75*exp(-((9*x(:,1)+1).^2/49+(9*x(:,2)+1)/10)) ...
            + 0.5*exp(-((9*x(:,1)-7).^2+(9*x(:,2)-3).^2)/4) ...
            - 0.2*exp(-((9*x(:,1)-4).^2+(9*x(:,2)-7).^2));
rhs = f(dsites);                        % function values
h = 2;                                  % upper bound for the radius
param = [2./npu 1/npu 3 3];             % initial values for the parameters 
                                        % to be optimized
                                        
global time1 time2 time3
time1 = 0;
time2 = 0;
time3 = 0;
tic
[epoints, Pf, ~, ~] = PU_op(M,dsites,neval,npu,rbf_aniso,f,rhs,h,param,puctrs);
toc
fprintf('Times in Cost_ep %g %g %g\n', time1, time2, time3)