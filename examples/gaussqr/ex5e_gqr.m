% ex5e_gqr
% This problem compares GaussQR collocation with a very refined finite
% difference solution
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .6;

clear functions % reset persistent variables

% Discretization choices
N = 35;
dt = .1;
T = 30;

% When do I want to compute the error in the solution
t_store = 0:.5:T;

% Compute finite difference solution, if needed
gmFDsols_file = 'newsols.mat'; % File to store results at
recompute_FDsols = 2; % 1 - Crank-Nicolson, 2 - ode15s
switch recompute_FDsols
    case 1
        ex5e_gqr_FD(N,dt,t_store,gmFDsols_file);
    case 2
        ex5e_gqr_FD_ode(N,dt,t_store,gmFDsols_file);
end

% Allocate storage for the solution
sol_store = cell(size(t_store));
count_store = 1;

% Solution parameters
h = 1e-7; % Jacobian-vector finite difference parameter
gmres_restart = 30;
Nit_max = 20; % Max Newton iterations
du_tol = 1e-5; % Newton step length for convergence

% Set up spatial, temporal discretization
time_steps = dt:dt:T;
x = pick2Dpoints([-1,-1],[1,1],N);

% Line search values
tls = [0 .25 .75 1];
Fls = zeros(1,4);

% User preferences
plot_initial_conditions = 0;
plot_solutions = 1;
store_solutions = 1;

% -1 - print nothing to the screen
% 0 - only print when a solution is stored
% 1 - print outputs at the Newton solve level
% 2 - print outputs at the linear solve level
keep_me_notified = 0;

% 0 - only compute the preconditioner once per time step
% 1 - direct solve at each Newton step
recompute_preconditioner = 1; 
% 0 - use the full Jacobian as the preconditioner
% 1 - use the ILUdt as the preconditioner
use_ilu_factorization = 0; % iterative solve using ILU


if plot_initial_conditions
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(N*N+1:end),N,N);
    subplot(1,2,1)
    surf(X,Y,A),title('Activator'),xlabel('x'),ylabel('y')
    subplot(1,2,2)
    surf(X,Y,H),title('Inhibitor'),xlabel('x'),ylabel('y')
    pause
end

if plot_solutions
    load(gmFDsols_file(1:end-4))
    u = sol_store{end};
    X = reshape(x(:,1),N,N);
    Y = reshape(x(:,2),N,N);
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(1+N*N:end),N,N);
    subplot(1,2,1)
    surf(X,Y,A),title('Activator'),xlabel('x'),ylabel('y')
    subplot(1,2,2)
    surf(X,Y,H),title('Inhibitor'),xlabel('x'),ylabel('y')
end
