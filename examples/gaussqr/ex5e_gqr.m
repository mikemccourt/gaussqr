% ex5e_gqr
% This problem compares GaussQR collocation with a very refined finite
% difference solution
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = 15;

clear functions % reset persistent variables

% When do I want to compute the error in the solution
T = 1;
t_store = 0:.1:T;

% Compute finite difference solution, if needed
gmFDsols_file = 'newsols.mat'; % File to store results at
recompute_FDsols = 0; % 1 - Crank-Nicolson, 2 - ode15s
if recompute_FDsols>0
    fprintf('Computing Finite Difference solution for comparison')

    % Finite difference discretization choices
    N = 25;
    switch recompute_FDsols
        case 1
            dt = .1;
            [computing_time,time_steps] = ex5e_gqr_FD(N,dt,t_store,gmFDsols_file);
        case 2
            [computing_time,time_steps] = ex5e_gqr_FD_ode(N,t_store,gmFDsols_file);
    end
end

% Compute the GaussQR solution
GAUSSQR_PARAMETERS.MY_EPSILON = 1e-4;
GAUSSQR_PARAMETERS.MY_ALPHA = 1;
N = 5;
gmRBFsols_file = 'rbfsols.mat';

[computing_time,time_steps] = ex5e_gqr_RBF_ode(N,t_store,gmRBFsols_file);