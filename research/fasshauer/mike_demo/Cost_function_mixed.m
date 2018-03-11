%-------------------------------------------------------------------------%
%
% File: Cost_function(rbf,param,puradius,idx_ds,index,q,M,puctrs,...
%       dsites,rhs)
%
% Goal: script that approximates the optimal parameters via LOOCV
%
% Inputs:  rbf:          radial basis function
%          param:        parameter values (remember this is only for the isotropic case
%          puradius:     minimum radius for the subdomains
%          idx_ds:       vector containing the indices of the data points
%                        located in a given PU subdomain
%          index:        index of the block
%          q:            number of blocks in one direction
%          M:            space dimension
%          puctrs:       the centre of the subdomain
%          dsites:       NXM matrix representing a set of N data sites
%          rhs:          the function values
%
% Outputs:  error:       error approximation on the j-th subdomain
%
% Calls on: IntegerBased_MD_Neighbourhood, IntegerBased_MD_RangeSearchAniso
%
%-------------------------------------------------------------------------%
function error = Cost_function_mixed(rbf, param, puradius, idx_ds, index, q, M, puctrs, dsites, rhs, strat, tikhonov)

% Determine if we are running in the isotropic (2 param) or anisotropic (4 param) version
assert(length(param) == 2 || length(param) == 4)
if length(param) == 2
    maxparam = param(1);  % Trying to modify to allow for just isotropic problem
    domain_params = [param(1), param(1)];
    shape_params = [param(2), param(2)];
else
    maxparam = max(param(1:M));
    domain_params = param(1:M);
    shape_params = param(M+1:2*M);
end

% Find the number of blocks needed to determine the points lying on the subdomain
global time1 time2 time3
ltimer = tic;
n = 1;
while n < 99999
    if maxparam > n * max(puradius)
        n = n + 1;
    else
        t = n;
        [dxx, dx] = IntegerBased_MD_Neighbourhood(dsites, idx_ds, index, q, M, t);
        break
    end
end
time1 = time1 + toc(ltimer);

% Find data sites located in the j-th subdomain
ltimer = tic;
idx = IntegerBased_MD_RangeSearchAniso(puctrs, domain_params, dxx, dx);
time2 = time2 + toc(ltimer);

% Estimate the error via the desired strategy
ltimer = tic;
error = Cost_ep(rbf, shape_params, idx, dsites, rhs, strat, tikhonov);
time3 = time3 + toc(ltimer);