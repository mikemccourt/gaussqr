% GSSEx3.m
% This example computes the modification of TPS
% pg. 24 of Fasshauer/Ye(2011) [Reproducing kernels of Sobolev spaces via a
% Green kernel approach with differential operators & boundary operators
% This involves the greens kernel of the biharmonic operator
% We are considering interpolation on the [-1,1]^2 domain

% The fundamental solution of the biharmonic operator
% The boundary operator term is accessed with R_func
% Notice the eps term added to prevent log(0)
G = @(x,z) DistanceMatrix(x,z).^2.*log(DistanceMatrix(x,z)+eps)/(8*pi);
Gx = @(x,z) (2*log(DistanceMatrix(x,z)+eps)+1).*(repmat(x(:,1),1,length(z))-repmat(z(:,1)',length(x),1))/(8*pi);
Gy = @(x,z) (2*log(DistanceMatrix(x,z)+eps)+1).*(repmat(x(:,2),1,length(z))-repmat(z(:,2)',length(x),1))/(8*pi);

% Choose the coefficients for the boundary kernel term
avec = [1 1 1];

% Consider data locations/kernel centers and evaluation points
N = 10;
x = pick2Dpoints([-1 -1],[1 1],N,'halton');
NN = 20;
xx = pick2Dpoints([-1 -1],[1 1],N);

% Think of some function we want to interpolate
yf = @(x) cos(x(:,1)+x(:,2));
y = yf(x);
yy = yf(xx);

% Design the MFS problem we need to solve to create the corrector
% Our centers are z_MFS and our boundary collocation points are x_MFS
N_x_MFS = 200;
N_z_MFS = 40;
N1d = ceil(N_x_MFS/4);
x1d = pickpoints(-1,1,N1d);
xb = ones(N1d,1);
x_MFS = unique([[x1d,xb];[xb,x1d];[x1d,-xb];[-xb,x1d]],'rows');
t_circle = pickpoints(0,2*pi,N_z_MFS+1);
t_circle = t_circle(1:end-1); % So that there is no duplicate center
z_MFS = 2*[cos(t_circle),sin(t_circle)];

% Identify the collocation points for each component of the problem
% x_d  = x_MFS(1:3:end,:);
% x_nx = x_MFS(2:3:end,:);
% x_ny = x_MFS(3:3:end,:);
x_d  = x_MFS;
x_nx = x_MFS(abs(x_MFS(:,1))==1,:);
x_ny = x_MFS(abs(x_MFS(:,2))==1,:);

% Create the MFS matrices
Gbar_d = G(x_d,z_MFS);
Grhs_d = G(x_d,x);
% Gbar_nx = Gx(x_nx,z_MFS);
% Grhs_nx = Gx(x_nx,x);
% Gbar_ny = Gy(x_ny,z_MFS);
% Grhs_ny = Gy(x_ny,x);
% Gbar = [Gbar_d;Gbar_nx;Gbar_ny];
% Grhs = [Grhs_d;Grhs_nx;Grhs_ny];
Gbar = [Gbar_d];
Grhs = [Grhs_d];

% Solve for the MFS coefficients
coef_MFS = Gbar\Grhs;

% Form the interpolation and evaluation matrix
K_int  = G(x,x) - G(x,z_MFS)*coef_MFS + R_func(avec,x,x);
K_eval = G(xx,x) - G(xx,z_MFS)*coef_MFS + R_func(avec,xx,x);

% Create the interpolant and check the error
yp = K_eval*(K_int\y);
errcompute(yp,yy)