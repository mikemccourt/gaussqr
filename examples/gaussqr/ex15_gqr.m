% ex15.m
% This example solves a nonlinear BVP in 1D with GQRr collocation
% The problem is the linear critical gradient
%       u_t - (k(u_x)u_x)_x = f   x=[-1,1]
% The diffusivity term is nonlinear and a little confusing
% You can read about it, and all the constants, in my thesis
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% True solution, also Dirichlet boundary condition
% usol = @(x,t) exp(-t)*(1-x.^2);
% usol = @(x,t) exp(-t)*cos(pi*x/2);
usol = @(x,t) erf((1-exp(-t))*4*x)+1;

% Choose parameters for the simulation
dt = .1/(2^2);
T = 1; % Final time (T=dt is one time step)
ep = 3.7;
ep = .1; % What do I want to say ...
alpha = 1;
N = 25;
NN = 100; % Error evaluation points

% Choose when to save the solution
t_save = .1:.1:1;
% t_save = .001:.001:.01;
err_save = zeros(size(t_save));
save_count = 1;

% Choose physical parameters for the diffusivity
% Setting kk = 0 will convert the problem back to a linear problem
% Defaults: kk=1,z=2,C=1,k0=1
DIFF_kk = 10;
DIFF_z = 1;
DIFF_C = .5;
DIFF_k0 = 1;

% Choose the boundary conditions
%   [0 0] - Dirichlet/Dirichlet
%   [1 1] - Neumann/Neumann
%   [0 1] - Dirichlet/Neumann
%   [1 0] - Neumann/Dirichlet
BC = [0 0];

% Choose lsqnonlin solve parameters
opts = optimset('Display','off');

% Set up the spacial discretization
point_spacing = 'even';
x = pickpoints(-1,1,N,point_spacing);
uold = usol(x,0);
xx = pickpoints(-1,1,NN); % In case you want to test error elsewhere

% First we must interpolate the initial condition for a guess of
% the coefficients for the time stepping
% This provides us an opportunity to test the choices of ep and alpha
GQR = gqr_rsolve(x,uold,ep,alpha);
up = gqr_eval(GQR,x);
errinit = errcompute(up,uold);
fprintf('error of initial condition interpolant : %g\n\n',errinit)

% Store the diffusivity coefficients for use in the residual 
GQR.DIFF_kk = DIFF_kk;
GQR.DIFF_z = DIFF_z;
GQR.DIFF_C = DIFF_C;
GQR.DIFF_k0 = DIFF_k0;

% Need to perform the time stepping
for t=dt:dt:T
    % This is computed to know how good we could do
    utrue = usol(x,t);
    GQRtrue = gqr_rsolve(x,utrue,ep,alpha);
    up = gqr_eval(GQRtrue,x);
    errtrue = errcompute(up,utrue);
    intres = ex15_gqr_resBC(GQRtrue.coef,GQR,x,uold,dt,BC,t);
    fprintf('At t=%g, error of interp : %g\t residual : %g\n',t,errtrue,norm(intres))
    
    % Consider the linear version, with k(u_x) = 1
    % This provides an initial guess for the nonlinear solver
    [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
    phi = gqr_phi(Marr,x,ep,alpha);
    phixx = gqr_phi(Marr,x,ep,alpha,2);
    A = phi/dt - phixx;
    A([1,end],:) = phi([1,end],:)/dt;
    
    % Compute the source term (need to encapsulate this)
%     S_u_xx = exp(-t)*(-2);
%     S_u_t = -exp(-t)*(1-x.^2);
    S_u_xx = 256/sqrt(pi)*(exp(-t)-1)^3*x.*exp(-16*x.^2*(1-exp(-t))^2);
    S_u_t = 8/sqrt(pi)*exp(-t)*x.*exp(-16*x.^2*(1-exp(-t))^2);
    S_f = S_u_t-S_u_xx;
    
    rhs = S_f + uold/dt;
    rhs([1,end]) = utrue([1,end]); % Apply Dirichlet BC
    
    GQRlin = GQR;
    c = A\rhs;
    GQRlin.coef = c;
    
    up = gqr_eval(GQRlin,x);
    errlin = errcompute(up,utrue);
    linres = ex15_gqr_resBC(GQRlin.coef,GQRlin,x,uold,dt,BC,t);
    fprintf('\t\t\t error of linear : %g\t residual : %g\n',errlin,norm(linres))
    
    % Try the nonlinear solve, using initial guess from linear solve
%     c = GQR.coef; % Previous time step solution
%     c = GQRlin.coef; % Current linear solve solution
    c = GQRtrue.coef;
    newcoef = lsqnonlin(@(coef) ex15_gqr_resBC(coef,GQR,x,uold,dt,BC,t),c,[],[],opts);
    GQR.coef = newcoef;
    
    ur = gqr_eval(GQR,x);
    errnln = errcompute(ur,utrue);
    nlnres = ex15_gqr_resBC(newcoef,GQR,x,uold,dt,BC,t);
    fprintf('\t\t\t error of nonlin : %g\t residual : %g\n',errnln,norm(nlnres))
    
%     plot(x,abs(utrue-up),'or',x,abs(utrue-ur))
%     pause
    
    % Consider saving the error if required
    if abs(t-t_save(save_count))<1e-10
        err_save(save_count) = errnln;
        save_count = save_count + 1;
    end
    uold = ur;
end