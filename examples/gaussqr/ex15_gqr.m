% ex15.m
% This example solves a nonlinear BVP in 1D with GQRr collocation
% The problem is the linear critical gradient
%       u_t - (k(u_x)u_x)_x = f   x=[-1,1]
% The diffusivity term is nonlinear and a little confusing
% You can read about it, and all the constants, in my thesis
rbfsetup
global GAUSSQR_PARAMETERS

% True solution, also Dirichlet boundary condition
usol = @(x,t) exp(-t)*(1-x.^2);

% Choose parameters for the simulation
dt = .01;
T = .01; % Final time
ep = .01;
alpha = 1;
N = 40;
NN = 100; % Error evaluation points

% Choose the boundary conditions
%   [0 0] - Dirichlet/Dirichlet
%   [1 1] - Neumann/Neumann
%   [0 1] - Dirichlet/Neumann
%   [1 0] - Neumann/Dirichlet
BC = [0 0];

% Set up the spacial discretization
x = pickpoints(-1,1,N,'cheb');
uold = usol(x,0);
xx = pickpoints(-1,1,NN);

% First we must interpolate the initial condition for a guess of
% the coefficients for the time stepping
% This provides us an opportunity to test the choices of ep and alpha
GQRold = gqr_rsolve(x,uold,ep,alpha);
up = gqr_eval(GQRold,x);
errinit = errcompute(uold,up);
fprintf('error of initial condition interpolant : %g\n',errinit)

% Need to perform the time stepping
for t=dt:dt:T
    fu = ex15_gqr_residual(GQRold,x,uold,dt)
    c = GQRold.coef;
%     newcoef = fminsearch(@(coef) ex15_gqr_resnorm(coef,GQRold,x,uold,dt,BC,t),c);
%     newcoef = fminunc(@(coef) ex15_gqr_resnorm(coef,GQRold,x,uold,dt,BC,t),c);
    newcoef = lsqnonlin(@(coef) ex15_gqr_resnorm(coef,GQRold,x,uold,dt,BC,t),c);
    [newcoef,c]
    GQR = GQRold;
    GQR.coef = newcoef;
    up = gqr_eval(GQR,x);
    plot(x,usol(x,t),'or',x,up)
end