% Finance1.m
% Basic finance example
% Drawn from Larsson et al 2013
%
% The problem is a parabolic problem,
%     du/dt = Lu
%     u(t,boundary) = bc(boundary)
%     u(0,x) = payout(x)
%
% Our differential operator for this problem has 3 pieces:
%      Lu = L1u + L2u + L3u
% L1u = r*x'*grad(u)
% L2u = (1/2)*x'*((B*B').*hess(u))*x
% L3u = -r*u
%
% The domain for this problem is all positive points in R^d such that their
% average value is less than 4K.  This domain is like a d dimensional
% triangle with 1-norm of the points less than 4Kd.
%
% Right now, we are just going to work this in 1D, but I will consider
% spicing it up later
clear all

% K is the scaled exercise price, which we set as 1
K = 1;

% T is the exercise time, chosen to be 1
T = 1;

% d is the number of assets in the option
d = 1;

% B is the dxd volatility matrix (nonsingular)
% It has diagonal elements .3, off diagonal elements .05
B = .3*eye(d) + .05*(ones(d)-eye(d));

% r is the risk-free interest rate
r = .05;

% tol is our computational tolerance
tol = 1e-4;

% payout(x) is the contract function which describes the payout
% This is also used to generate initial conditions
payout = @(x) max(0,sum(x-K,2)/d);

% bc(x,t) defines the boundary conditions which are used at the near-field
% and far-field locations.  The near-field location is the 0 point and the
% far-field location is the plane sum(x)=4Kd.  Other boundary conditions
% are not required to provide well-posedness.
% Note: This function will return zero everywhere except along the plane 
% sum(x)==4*K*d, meaning the x=0 BC is automatically built-in
bc = @(x,t) K*(4-exp(-r*t))*(sum(x,2)==4*K*d);

% Choose some collocation points in the domain
% N is the total number of points to compute with
% x_bc and x_int are the boundary and interior points
% x_eval is a set of points to evaluate the solution on
pt_opt = 'cheb';
N = 29;
x = pickpoints(0,4*K,N,pt_opt);
x_int = x(2:end-1);
x_bc = x([1,end]);
x_all = [x_int;x_bc];
N_eval = 300;
x_eval = pickpoints(0,4*K,N_eval);

% Choose an RBF to solve with, for us it is the Gaussian
% The derivatives are needed to create the collocation matrix
rbf = @(e,r) exp(-(e*r).^2);
rbfdx = @(e,r,dx) 2*e^2*dx.*exp(-(e*r).^2);
rbfdxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
% We must choose a shape parameter ep>0
ep = 3;
% We also form distance matrices which we need for evaluation
DM_all = DistanceMatrix(x_all,x_all);
DM_bc = DistanceMatrix(x_bc,x_all);
DM_int = DistanceMatrix(x_int,x_all);
DiffM_int = DifferenceMatrix(x_int,x_all);
DM_eval = DistanceMatrix(x_eval,x_all);
% The basis is evaluated using those distance matrices
RM_all = rbf(ep,DM_all);
RM_bc = rbf(ep,DM_bc);
RM_int = rbf(ep,DM_int);
RxM_int = rbfdx(ep,DM_int,DiffM_int);
RxxM_int = rbfdxx(ep,DM_int);
RM_eval = rbf(ep,DM_eval);

% We can define our differential operator using differentiation matrices
% The differentiation matrices can be defined now and used in perpetuity
Dx = RxM_int/RM_all;
Dxx = RxxM_int/RM_all;
L = @(u,x,t) r*x.*(Dx*u) + 1/2*B^2*(x.^2).*(Dxx*u)-r*u(1:end-2);

% We need to set up a time stepping scheme
% Let's just try something simple for now
% We'll use Euler's method:
%      u_{k+1} = u_k + dt*Lu_k
% We choose a time step that is really small
dt = 1e-4;
t_vec = 0:dt:T;

% u_coef will contain the coefficients for our solution basis
% u_sol will contain the solution at the collocation points
u_sol = zeros(N,length(t_vec));
u_coef = zeros(N,length(t_vec));

% We need to record our initial condition
% We must also interpolate our initial condition
u_sol(:,1) = payout(x_all);
u_coef(:,1) = RM_all\u_sol(:,1);

% Perform the time stepping
k = 1;
for t=t_vec(2:end)
    % On the interior
    u_sol(1:end-2,k+1) = u_sol(1:end-2,k) + dt*L(u_sol(:,k),x_int,t);
    % On the boundary
    u_sol(end-1:end,k+1) = bc(x_bc,t);
    k = k + 1;
end

u_sol = u_sol(1:end-2,:);

[XX,TT] = meshgrid(x_int,t_vec);
h = surf(XX,TT,u_sol');
set(h,'edgecolor','none')