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
% NOTE: This is inactive right now
tol = 1e-4;

% We need to set up a time stepping scheme
% Choices that are available
%   1) Euler's Method: u_{k+1} = u_k + dt*Lu_k
%   2) Backward Euler: u_{k+1} = u_k + dt*Lu_{k+1}
%   3) BDF2 + 1 BE   : u_{k+2} - 4/3*u_{k+1} + 1/3*u_k = 2/3*dt*Lu_{k+2}
% Note here that t is actually measuring "time to expiry" not "time from
% right now".  As a result, we are kind of solving this problem backwards
ts_scheme = 4;
dt = 1e-2;

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

% We must choose a shape parameter ep>0
% This if-block will allow you to either choose epsilon here (if you
% haven't yet) or run with an existing epsilon (for batch jobs)
% If you have already run this script (and thus ep exists) the value
% provided here is not accessed at all
if not(exist('ep','var'))
    ep = 1.25;
end

% This is an option as to what spatial solver to use
%    hssvd = -2 for polynomial collocation
%          = -1 for finite differences
%          = 0  for RBF-direct
%          = 1  for HS-SVD
% NOTE: FD only allowed for evenly spaced points
% NOTE: polynomials allowed only on Chebyshev spaced points
solver = -2;

% The following are HS-SVD parameters
% The alpha value determines eigenfunction locality
% reg = 1 asks for a low rank eigenfunction expansion
reg = 1;
alpha = 3;

% If hssvd = 0, you can choose what RBF you want to run with
%    rbf_choice = 1 is Gaussian
%               = 2 is Multiquadric
rbf_choice = 1;

% Because of the discontinuity in the solution at t=0, we also have the
% option of trying to solve the problem in a coupled sense
% Choosing coupling=1 will activate this option and solve two problems: one
% on [0,K] and the other on [K,4K]
% At x=K, the solution will have matching values and first derivatives, but
% nothing else is required
% NOTE: Not active yet
coupling = 0;

% Choose some collocation points in the domain
% N is the total number of points to compute with
% pt_opt is the distribution of points in the domain
%     pt_opt = 'even' is uniformly space
%            = 'cheb' has Chebyshev spacing
%            = 'cp_even' has coupling (see below), even interiors
%            = 'cp_cheb' has coupling (see below), Cheb interiors
% NOTE: pt_opt may be overwritten below if incompatible with solver
N = 41;
pt_opt = 'cheb';

% Overwrite point selection if the solver is incompatible
if solver==-1 && ~strcmp(pt_opt(end-3:end),'even')
    pt_opt = 'even';
elseif solver==-2 && ~strcmp(pt_opt(end-3:end),'chebs')
    pt_opt = 'cheb';
end

% Select some points based on the user requests
x = pickpoints(0,4*K,N,pt_opt);
N = length(x);

% Cut up the domain into pieces to work with
% x_bc and x_int are the boundary and interior points
% x_eval is a set of points to evaluate the solution on
x_int = x(2:end-1);
x_bc = x([1,end]);
x_all = [x_int;x_bc];
i_int = 1:length(x_int);
i_bc = length(i_int)+1:length(i_int)+length(x_bc);
N_int = length(i_int);
N_bc = length(i_bc);
N_eval = 300;
x_eval = pickpoints(0,4*K,N_eval);

% Set up the solver we are using
switch solver
    case 1
    % Note the dividing factor required to account for the change in scale
    % we are imposing to normalize our spatial domain to [-1,1]
        GQR = gqr_solveprep(reg,(x_all-2)/2,ep,alpha);
        if reg==0
            Rmat = [eye(N);GQR.Rbar];
        else
            Rmat = 1;
        end
        RM_all = gqr_phi(GQR,(x_all-2)/2)*Rmat;
        RM_int = gqr_phi(GQR,(x_int-2)/2)*Rmat;
        RxM_int = gqr_phi(GQR,(x_int-2)/2,1)*Rmat/2;
        RxxM_int = gqr_phi(GQR,(x_int-2)/2,2)*Rmat/4;
    case 0
        % Create a function for the RBF
        % The derivatives are needed to create the collocation matrix
        switch rbf_choice
            case 1
                rbf = @(e,r) exp(-(e*r).^2);
                rbfdx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
                rbfdxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);
            case 2
                rbf = @(e,r) sqrt(1+(e*r).^2);
                rbfdx = @(e,r,dx) e^2*dx./sqrt(1+(e*r).^2);
                rbfdxx = @(e,r) e^2./(1+(e*r).^2).^(3/2);
        end
        % Could also consider the Multiquadrics
        % We also form distance matrices which we need for evaluation
        DM_all = DistanceMatrix(x_all,x_all);
        DM_bc = DistanceMatrix(x_bc,x_all);
        DM_int = DistanceMatrix(x_int,x_all);
        DiffM_int = DifferenceMatrix(x_int,x_all);
        % The basis is evaluated using those distance matrices
        RM_all = rbf(ep,DM_all);
        RM_int = rbf(ep,DM_int);
        RxM_int = rbfdx(ep,DM_int,DiffM_int);
        RxxM_int = rbfdxx(ep,DM_int);
    case -1
        delta_x = 4*K/(N-1);
        RM_all = 1; % This isn't needed for FD
        RM_int = eye(N_int,N);
        RxM_int = ([zeros(N_int,N_bc),eye(N_int)]-eye(N_int,N))/(2*delta_x);
        RxxM_int = (eye(N_int,N)-2*[zeros(N_int,1),eye(N_int),zeros(N_int,1)]+[zeros(N_int,N_bc),eye(N_int)])/delta_x^2;
    case -2
        % Form the Chebyshev differentiation matrices
        D_cheb = cheb(N);
        D2_cheb = D_cheb^2;
        RM_all = 1; % This isn't needed for polynomial collocation
        RM_int = eye(N_int,N);
        % Note the division by 2 or 2^2 for the change in scale
        % cheb expects [-1,1] but we're on [0,4]
        RxM_int = D_cheb(2:N-1,:)/2;
        RxxM_int = D2_cheb(2:N-1,:)/4;
end

% Form the differentiation matrices
% Note that only L_mat is used after this block
D0 = RM_int/RM_all;
Dx = RxM_int/RM_all;
Dxx = RxxM_int/RM_all;
% Form the differential operator using the differentiation matrices
L_mat = r*diag(x_int)*Dx+1/2*B^2*diag(x_int.^2)*Dxx-r*D0;
% Swap the columns, if needed for finite differences
if solver<0
    L_mat = L_mat(:,[2:N-1,1,N]);
end


% Set up the time stepping mechanism
t_vec = 0:dt:T;
% Compute necessary time stepping components
if ts_scheme==1 || ts_scheme==4
    L = @(u) L_mat*u;
    if ts_scheme==4
        odefun = @(t,u) [L(u);u(i_bc)-bc(x_bc,t)];
        jac = [L_mat;zeros(N_bc,N_int),eye(N_bc)];
        mass = [eye(N_int),zeros(N_int,N_bc);zeros(N_bc,N)];
        odeopt = odeset('InitialStep',dt,'Jacobian',jac,'Mass',mass,...
            'MStateDependence','none','MassSingular','yes');
    end
elseif ts_scheme==2 || ts_scheme==3
    BE_mat = eye(N) - dt*[L_mat;zeros(N_bc,N)];
    if ts_scheme==3
        BDF_mat = eye(N) - 2/3*dt*[L_mat;zeros(N_bc,N)];
    end
end

% u_coef will contain the coefficients for our solution basis
% u_sol will contain the solution at the collocation points
u_sol = zeros(N,length(t_vec));

% We need to record our initial condition
% We must also interpolate our initial condition
u_sol(:,1) = payout(x_all);

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
                u_sol(i_int,k+1) = u_sol(i_int,k) + dt*L(u_sol(:,k));
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

% The true solution is described above, and defined below
% Note that we have substituted t for T-t since we are solving an initial
% value problem and not a final value problem
d1_truesol = @(x,t) 1./(B*sqrt(t)).*(log(x/K)+(r+B^2/2)*t);
d2_truesol = @(x,t) d1_truesol(x,t) - B*sqrt(t);
C_truesol = @(x,t) normcdf(d1_truesol(x,t)).*x - K*normcdf(d2_truesol(x,t)).*exp(-r*t);

% Choose whether or not to plot the output
plot_sol = 1;

% Plot the error in the solution
if plot_sol
    figure
    [XX,TT] = meshgrid(x_int,t_vec);
    h = surf(XX,TT,abs(u_sol' - C_truesol(XX,TT)));
    set(h,'edgecolor','none')
    title(sprintf('dt=%g,\tN=%d,\tspace=%s,\tts=%d,solver=%d,\tep=%g\t',dt,N,pt_opt,ts_scheme,solver,ep))
    xlabel('Spot price')
    ylabel('time to expiry')
    zlabel('option value error')
end