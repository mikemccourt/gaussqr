% disksourcetest_MFS.m
% This example studies the convergence rate for a 2D Laplace problem when
% the boundary conditions (or interface conditions) are tainted by noise.
% The problem we want to solve is an easy one:
%     Lap(u) = 0,     interior
%          u = v,     boundary

% The Dirichlet boundary condition function
%v = @(x) 1 + x(:,1) + x(:,2);
v = @(x) x(:,1) + x(:,2) + 1./((x(:,1)-2).^2 + x(:,2).^2);

% The fundamental solution
fs = @(x,z) log(DistanceMatrix(x,z));

% Error to introduce into the boundary condition (0 for no noise)
noise = 0.0001;

% Circle radius
% Fictitious boundary radius
r = 1;
R = 2.5;

% Number of collocation and source points and test points
NN = 100;
tt = linspace(0,2*pi,NN+1)';tt = tt(1:NN);
xx = [cos(tt),sin(tt)];
yy = v(xx);

% Define the problem size vector
Nvec = floor(logspace(1,3,30));

errvec = zeros(size(Nvec));
errvec_reg = zeros(size(Nvec));
errvec_noise = zeros(size(Nvec));
errvec_noise_reg = zeros(size(Nvec));
coefvec = zeros(size(Nvec));
coefvec_reg = zeros(size(Nvec));
coefvec_noise = zeros(size(Nvec));
coefvec_noise_reg = zeros(size(Nvec));
condvec = zeros(size(Nvec));
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
    % no noise, no regularization
    yp = K_eval*(K\y);
    errvec(k) = errcompute(yp,yy);
    coefvec(k) = norm(K\y);
    % no noise, with regularization
    pinvK = pinv(K);
    c_reg = pinvK*y;
    yp_reg = K_eval*c_reg;
    errvec_reg(k) = errcompute(yp_reg,yy);
    coefvec_reg(k) = norm(c_reg);
    % noise, no regularization
    yp_noise = K_eval*(K\y_noise);
    errvec_noise(k) = errcompute(yp_noise,yy);
    coefvec_noise(k) = norm(K\y_noise);
    % noise, with regularization
    c_noise_reg = pinvK*y_noise;
    yp_noise_reg = K_eval*c_noise_reg;
    errvec_noise_reg(k) = errcompute(yp_noise_reg,yy);
    coefvec_noise_reg(k) = norm(c_noise_reg);
    % condition number (always the same)
    condvec(k) = cond(K);
    k = k + 1;
end

h1 = figure;
loglog(Nvec,[errvec;errvec_reg;errvec_noise;errvec_noise_reg],'linewidth',3)
xlabel('Collocation/Source points')
ylabel('Error')
title(sprintf('Noise = %g',noise))
legend('noiseless','noiseless pinv','noisy','noisy pinv','location','southwest')

h2 = figure;
%loglog(Nvec,[coefvec;coefvec_noise;condvec],'linewidth',3)
[AX,H1,H2] = plotyy(Nvec,[coefvec;coefvec_reg;coefvec_noise;coefvec_noise_reg],Nvec,condvec,'loglog','loglog','linewidth',3);
set(H1,'linewidth',3)
xlabel('Collocation/Source points')
ylabel('2-norm of coefficient vector')
set(H2,'linewidth',3)
ylabel(AX(2),'cond(K)')
title(sprintf('Noise = %g',noise))
legend('noiseless','noiseless pinv','noisy','noisy pinv','cond','location','northwest')
