% Test comparing the PUM method to a global method
addpath('common_routines')
addpath('codes_op')

% The test function
yf = @(x) franke(x(:, 1), x(:, 2));

M = 2;                                  % space dimension 
n = 25;                                 % number of tracks
m = 80;                                 % number of points on tracks
[N, dsites, yy] = TrackData2D(n,m);       % generate N = n*m track data in 2D
neval = 30;                             % parameter for evaluation points
npu = floor(((N)/(4))^(1/M));          % parameter for PU centres
xx = linspace(0,1,n);                   % subdomains in one direction
[X, Y] = meshgrid(xx,yy);                % patches centered at tracks
puctrs = [X(:) Y(:)];                   % define the PU centres
rbf_aniso = @(r) 1./sqrt(1+r.^2);       % define the RBF
rhs = yf(dsites);                        % function values
h = 2;                                  % upper bound for the radius
param = [2/npu 1/npu 3 3];             % initial values for the parameters 
                                       % to be optimized
          
global time1 time2 time3
time1 = 0;
time2 = 0;
time3 = 0;
% Maybe the constant warning issues are slowing things down ???
warning off
outer_timer = tic;
[epoints, Pf, radii, epsilon] = PU_op(M, dsites, neval, npu, rbf_aniso, yf, rhs, h, param, puctrs);
toc(outer_timer)
warning on
fprintf('Times in Cost_ep %g %g %g\n', time1, time2, time3)

% What about if we just do interpolation with a C2 Wendland kernel
rbf = @(r) max(1-r,0).^4.*(4*r+1);
% This is a C6 Wendland kernel below, if we prefer
% k = 3;
% d = 2;
% l = floor(d / 2 + k + 1);
% rbf = @(r) (1 - r) .^ (l + 3) .* (1+(l+3)*r + (6*l^2+36*l+45)/15*r.^2 + (l^3+9*l^2+23*l+15)/15*r.^3);
% We would have to choose some kernel parameters
epvec = [10, 5];
tic
K = DistanceMatrix(dsites, dsites, epvec, rbf);
Keval = DistanceMatrix(epoints, dsites, epvec, rbf);
yeval = Keval * (K \ rhs);

exact = yf(epoints);
maxerr = norm(yeval - exact,inf);
rms_err = norm(yeval - exact)/sqrt(length(exact));
toc
fprintf('RMS error:       %e\n', rms_err);
fprintf('Maximum error:   %e\n', maxerr);
