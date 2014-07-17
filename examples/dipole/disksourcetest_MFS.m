% disksourcetest_MFS.m
% This example studies the convergence rate for a 2D Laplace problem when
% the boundary conditions (or interface conditions) are tainted by noise.
% The problem we want to solve is an easy one:
%     Lap(u) = 0,     interior
%          u = v,     boundary

% The Dirichlet boundary condition function
v = @(x) 1 + x(:,1) + x(:,2);
% v = @(x) x(:,1) + x(:,2) + 1./((x(:,1)-2).^2 + x(:,2).^2);

% The fundamental solution
fs = @(x,z) log(DistanceMatrix(x,z));

% Error to introduce into the boundary condition (0 for no noise)
noise = 0.01;

% Circle radius
% Fictitious boundary radius
r = 1;
R = 1.5;

% Number of collocation and source points and test points
NN = 100;
tt = linspace(0,2*pi,NN+1)';tt = tt(1:NN);
xx = [cos(tt),sin(tt)];
yy = v(xx);

% Define the problem size vector
Nvec = floor(logspace(1,3,30));

errvec = zeros(size(Nvec));
errvec_noise = zeros(size(Nvec));
coefvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Define the collocation and source points, x and z, and test points xx
    % t is the angle; the last value is dumped to avoid duplication
    t = linspace(0,2*pi,N+1)';t = t(1:N);
    x = [cos(t),sin(t)];
    z = R*x;

    % Form the linear system, with noise if requested
    K = fs(x,z);
    K_eval = fs(xx,z);
    y = v(x);
    y_noise = v(x) + noise*randn(N,1);

    % Solve the system with and without noise and compute the errors
    yp = K_eval*(K\y);
    errvec(k) = errcompute(yp,yy);
    yp_noise = K_eval*(K\y_noise);
    coefvec(k) = norm(yp_noise);
    errvec_noise(k) = errcompute(yp_noise,yy);
    k = k + 1;
end

h = figure;
loglog(Nvec,[errvec;errvec_noise],'linewidth',3)
xlabel('Collocation/Source points')
ylabel('Error')
title(sprintf('Noise = %g',noise))
legend('noiseless','noisy','location','southwest')

h = figure;
loglog(Nvec,coefvec)