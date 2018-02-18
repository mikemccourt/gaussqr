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
function error = Cost_function_mike(rbf, param, puradius, idx_ds, index, q, M, puctrs, dsites, rhs)

% Find the number of blocks needed to determine the points lying on the subdomain
n = 1;
assert(length(param) == 2)
maxparam = param(1);  % Trying to modify to allow for just isotropic problem
global time1 time2 time3

ltimer = tic;
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
% idx = IntegerBased_MD_RangeSearchAniso(puctrs, param(1:M), dxx, dx);
idx = IntegerBased_MD_RangeSearchAniso(puctrs, [param(1), param(1)], dxx, dx);
time2 = time2 + toc(ltimer);

% Estimate the error via LOOCV
ltimer = tic;
% error = Cost_ep(rbf, param(M+1:2*M), idx, dsites, rhs);
error = Cost_ep(rbf, [param(2), param(2)], idx, dsites, rhs);
time3 = time3 + toc(ltimer);