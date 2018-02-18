%-------------------------------------------------------------------------%
%
% File: Cost_function_con(rbf,param,puradius,idx_ds,index,q,M,puctrs,...
%       dsites,rhs,h)
%
% Goal: script that imposes constraints on the optimal parameters
%       and approximates the error
%
% Inputs:  rbf:          radial basis function
%          param:        guess for parameters (radii and shape parameters)
%          puradius:     minimum radius for the subdomains
%          idx_ds:       vector containing the indices of the data points
%                        located in a given PU subdomain
%          index:        index of the block
%          q:            number of blocks in one direction
%          M:            space dimension
%          puctrs:       centre of the subdomain
%          dsites:       NXM matrix representing a set of N data sites
%          rhs:          function values
%          h:            upper bound for the radius
%
% Outputs:  error:       error approximation on the j-th subdomain
%
% Calls on: Cost_function: computes the error estimates
%
%-------------------------------------------------------------------------%
function [error] = Cost_function_con(rbf, param, puradius, idx_ds, index, q, M, puctrs, dsites, rhs, h)
% Impose constraints and evaluate the cost function
a = (param(1:M) < puradius);
b = (param(1:M) > h * puradius);
c = (param(M+1:2*M) < 0);

% !!!!!!!!!!!!!!! What does this max do??
if max(a + b + c) > 0
    fprintf('\t\t%g %g %g %g %g %g\n', a, b, c)
%     error = 1e40;
    error = Cost_function(rbf, param, puradius, idx_ds, index, q, M, puctrs, dsites, rhs);
else
    % Evaluate the error
    error = Cost_function(rbf, param, puradius, idx_ds, index, q, M, puctrs, dsites, rhs);
end

if isnan(error) || isinf(error)
    error = 1e50;
end