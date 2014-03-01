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
% The exact solution is the formula for the value of a European call option
% as defined on Wikipedia:
%      C(x,t) = F_Z(d1)*u - F_Z(d2)*K*exp(-r*(T-t)) 
% where
%      F_Z(z) = P(Z<z) for Z~N(0,1)
%      d1(x) = 1/(B*sqrt(T-t))*(log(x/K)+(r+B^2/2)*(T-t))
%      d2(x) = d1(x) - B*sqrt(T-t)
% Our solution below uses t in place of T-t because we know the final value
% and not the initial value.  This means that when solving the problem we
% are really stepping backwards in time.
%
% Right now, we are just going to work this in 1D, but I will consider
% spicing it up later
% NOTE: This requires the statistics toolbox, for now
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

% The true solution is described above, and defined below
% Note that we have substituted t for T-t since we are solving an initial
% value problem and not a final value problem
d1_truesol = @(x,t) 1./(B*sqrt(t)).*(log(x/K)+(r+B^2/2)*t);
d2_truesol = @(x,t) d1_truesol(x,t) - B*sqrt(t);
C_truesol = @(x,t) normcdf(d1_truesol(x,t)).*x - K*normcdf(d2_truesol(x,t)).*exp(-r*t);

% Choose some collocation points in the domain
% N is the total number of points to compute with
% x_bc and x_int are the boundary and interior points
% x_eval is a set of points to evaluate the solution on
pt_opt = 'cheb';
N = 40;
x = pickpoints(0,4*K,N,pt_opt);
x_int = x(2:end-1);
x_bc = x([1,end]);
x_all = [x_int;x_bc];
i_int = 1:length(x_int);
i_bc = length(i_int)+1:length(i_int)+length(x_bc);
N_int = length(i_int);
N_bc = length(i_bc);
N_eval = 300;
x_eval = pickpoints(0,4*K,N_eval);

% If we want to test with finite differences
% Note that we MUST have pt_opt='even' to do this
% This also will ONLY work with one of the BDF methods
finite_diff_test = 0;
finite_diff_test = finite_diff_test*(strcmp(pt_opt,'even'));

% Choose an RBF to solve with, for us it is the Gaussian
% The derivatives are needed to create the collocation matrix
rbf = @(e,r) exp(-(e*r).^2);
rbfdx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
rbfdxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
% We must choose a shape parameter ep>0
ep = 2;
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
% D0 just picks out the elements from the interior region
D0 = RM_int/RM_all;
Dx = RxM_int/RM_all;
Dxx = RxxM_int/RM_all;
% Should we choose to instead run with finite differences
if finite_diff_test==1
    delta_x = 4*K/(N-1);
    D0 = eye(N_int,N);
    Dx = ([zeros(N_int,N_bc),eye(N_int)]-D0)/(2*delta_x);
    Dxx = (D0-2*[zeros(N_int,1),eye(N_int),zeros(N_int,1)]+[zeros(N_int,N_bc),eye(N_int)])/delta_x^2;
    I = eye(N);
    swap = I([N-1,1:N-2,N],:);
end

% We need to set up a time stepping scheme
% Choices that are available
%   1) Euler's Method: u_{k+1} = u_k + dt*Lu_k
%   2) Backward Euler: u_{k+1} = u_k + dt*Lu_{k+1}
%   3) BDF2 + 1 BE   : u_{k+2} - 4/3*u_{k+1} + 1/3*u_k = 2/3*dt*Lu_{k+2}
% Note here that t is actually measuring "time to expiry" not "time from
% right now".  As a result, we are kind of solving this problem backwards
ts_scheme = 3;
dt = 1e-3;
t_vec = 0:dt:T;
% Compute necessary time stepping components
if ts_scheme==1 || ts_scheme==4
    L = @(u,x,t) r*x.*(Dx*u) + 1/2*B^2*(x.^2).*(Dxx*u)-r*D0*u;
    if ts_scheme==4
        odefun = @(t,u) [L(u,x_int,t);u(i_bc)-bc(x_bc,t)];
        jac = [r*diag(x_int)*Dx+1/2*B^2*diag(x_int.^2)*Dxx-r*D0;...
            zeros(N_bc,N_int),eye(N_bc)];
        mass = [eye(N_int),zeros(N_int,N_bc);zeros(N_bc,N)];
        odeopt = odeset('InitialStep',dt,'Jacobian',jac,'Mass',mass,...
            'MStateDependence','none','MassSingular','yes');
    end
elseif ts_scheme==2 || ts_scheme==3
    L_mat = r*diag(x_int)*Dx+1/2*B^2*diag(x_int.^2)*Dxx-r*D0;
    if finite_diff_test==1
        L_mat = L_mat*swap;
    end
    BE_mat = eye(N) - dt*[L_mat;zeros(N_bc,N)];
    if ts_scheme==3
        BDF_mat = eye(N) - 2/3*dt*[L_mat;zeros(N_bc,N)];
    end
end

% u_coef will contain the coefficients for our solution basis
% u_sol will contain the solution at the collocation points
u_sol = zeros(N,length(t_vec));
u_coef = zeros(N,length(t_vec));

% We need to record our initial condition
% We must also interpolate our initial condition
u_sol(:,1) = payout(x_all);
u_coef(:,1) = RM_all\u_sol(:,1);

% Perform the time stepping
% If using builtin ODE solver, just call it
% Otherwise, manually perform the iterations
if ts_scheme==4
    [t_ret,u_ret] = ode15s(odefun,t_vec,u_sol(:,1),odeopt);
    u_sol = u_ret';
else
    k = 1;
    for t=t_vec(2:end)
        switch ts_scheme
            case 1
                % On the interior
                u_sol(i_int,k+1) = u_sol(i_int,k) + dt*L(u_sol(:,k),x_int,t);
                % On the boundary
                u_sol(i_bc,k+1) = bc(x_bc,t);
            case 2
                % Form the linear system [Int equations;BC equations]
                % We only need to form it once, so do that once and save it
                rhs = [u_sol(i_int,k);bc(x_bc,t)];
                u_sol(:,k+1) = BE_mat\rhs;
            case 3
                % Take one BE step to produce the necessary two old values
                if k==1
                    rhs = [u_sol(i_int,1);bc(x_bc,t)];
                    u_sol(:,2) = BE_mat\rhs;
                else
                    rhs = [4/3*u_sol(i_int,k)-1/3*u_sol(i_int,k-1);bc(x_bc,t)];
                    u_sol(:,k+1) = BDF_mat\rhs;
                end
        end
        k = k + 1;
    end
end

% Dump the boundary conditions, which are not fun
u_sol = u_sol(i_int,:);

% Plot the error in the solution
figure
[XX,TT] = meshgrid(x_int,t_vec);
h = surf(XX,TT,abs(u_sol' - C_truesol(XX,TT)));
set(h,'edgecolor','none')
title(sprintf('dt=%g,\tN=%d,\tspace=%s,\tts=%d,\tep=%g\t',dt,N,pt_opt,ts_scheme,ep))
xlabel('Spot price')
ylabel('time to expiry')
zlabel('option value error')