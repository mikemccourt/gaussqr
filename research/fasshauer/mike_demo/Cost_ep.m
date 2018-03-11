%-------------------------------------------------------------------------%
%
% File: Cost_ep(rbf,param,idx_ds,dsites,rhs)
%
% Goal: script that approximates the optimal parameters via LOOCV
%
% Inputs:   rbf:        radial basis function
%           param:      guess for the shape parameters 
%           idx_ds:     vector containing the indices of the data points
%                       located in a given PU subdomain
%           dsites:     NXM matrix representing a set of N data sites
%           rhs:        function values
%          strat:        which strategy to use: 'loocv' or 'mle'
%
% Outputs:  error:      estimate the error via LOOCV on the j-th subdomain
%
% Calls on: DistanceMatrixAniso
%
%-------------------------------------------------------------------------%
function error = Cost_ep(rbf, param, idx_ds, dsites, rhs, strat, tikhonov)

% Prevent ourselves from considering situations with insufficiently many points in the patch
if length(idx_ds) < 5
    error = 1e60;
    return
end

K = rbf(real(DistanceMatrixAniso(dsites(idx_ds, :), dsites(idx_ds, :), param)));

% We can add a Tikhonov regularization term if we want (though we have to
% do it elsewhere as well
K = K + tikhonov * eye(size(K));

% Check if the K matrix is even symmetric PD
% If it is not, crash out of this function and return a lousy value
% We could add an amount to the diagonal if we wanted to counter ill-conditioning?
try
    R = chol(K);
catch exception
    if strcmp(exception.identifier, 'MATLAB:posdef')
        error = 1e30;
        return
    end
    rethrow(exception)
end

% Simultaneously solve for both the inverse and the interpolant
% Use the Cholesky factorization which already exists: R'*R*X = B
%
% It's possible we're wasting juice because we're not reusing this Cholesky
% factorization later on when we actually need to do predictions (I assume
% we just invert the matrix again).  That's another possible optimization.
if strcmp(strat, 'loocv')
    rhs_and_I = [rhs(idx_ds), eye(size(K,1))];
    opts.UT = false;
    opts.LT = true;
    opts.TRANSA = true;
    temp = linsolve(R, rhs_and_I, opts);
    
    opts.UT = true;
    opts.LT = false;
    opts.TRANSA = false;
    sol_and_Kinv = linsolve(R, temp, opts);
    error = norm(sol_and_Kinv(:, 1) ./ diag(sol_and_Kinv(:, 2:end)), inf);
else
    % This is the MLE instead of the LOOCV
    % Not sure if it is better but it should be considered
    N = length(K);
    y = rhs(idx_ds);
    opts.UT = false;
    opts.LT = true;
    opts.TRANSA = true;
    temp = linsolve(R, y, opts);

    opts.UT = true;
    opts.LT = false;
    opts.TRANSA = false;
    K_inv_y = linsolve(R, temp, opts);
    error = N * log(max(y' * K_inv_y, 1e-16)) + sum(log(max(diag(R), 1e-16)));
end