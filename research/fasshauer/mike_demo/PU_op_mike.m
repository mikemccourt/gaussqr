%-------------------------------------------------------------------------%
%
% File: PU_op_mike(M,dsites,neval,npu,rbf,wf,f,rhs,puctrs,isotropic,strat,tikhonov);
%
% Goal: script that performs partition of unity with variable patches
%       and shape parameters
%
% Inputs:     M:          space dimension
%             dsites:     NXM matrix representing a set of N data sites
%             neval:      number of evaluation points in one direction
%             npu:        number of PU subdomains in one direction
%             rbf:        radial basis function
%             wf:         weight function
%             f:          test function
%             rhs:        function values
%             puctrs:     the locations of the PU centers
%             isotropic:  true/false
%             strat:      what parametrization strategy to use
%             tikhonov:   a Tikhonov regularization parameter
%
% Outputs:    epoints:    evaluation points
%             Pf:         interpolant computed at the evaluation points
%
% Calls on:   IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%             IntegerBased_MD_RangeSearchAniso, MakeSDGrid,
%             IntegerBased_MD_ContainingQuery, DistanceMatrixAniso
%             PuweightAniso, Cost_function_con
%
% Remarks: 1) DistanceMatrixAniso, MakeSDGrid come from the books:
%             [G.E. Fasshauer, Meshfree Approximation Methods with Matlab,
%             World Scientific, Singapore, 2007]
%             [G.E. Fasshauer, M.J. McCourt, Kernel-based Approximation 
%             Methods using Matlab, World Scientific, Singapore, 2015]
%          2) IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%             IntegerBased_MD_ContainingQuery come from the paper:
%             [R. Cavoretto, A. De Rossi, E. Perracchione, Optimal 
%             selection of local approximants in RBF-PU interpolation,
%             J. Sci. Comput. (2017), in press. 
%             DOI: 10.1007/s10915-017-0418-7]
%
%             fminsearch is a Matlab routine used to perform uncostrained,
%             derivative-free optimization
%
%-------------------------------------------------------------------------%
function [epoints, Pf, radii, epsilon] = PU_op_mike(M, dsites, neval, npu, rbf, f, rhs, puctrs, isotropic, strat, tikhonov, verbose)

% Create neval^M equally spaced evaluation points
epoints = MakeSDGrid(M, neval);

% Initialization
puradius = ones(1, M) * (1 ./ npu);
npu_M = size(puctrs, 1);
neval_M = size(epoints, 1);

% Parameter for integer-based partitioning structure
q = ceil(1 ./ puradius(1, 1));

% Build the partitioning structure for data sites and evaluation points
idx_ds = IntegerBased_MD_Structure(dsites, q, puradius(1, 1), M);
idx_ep = IntegerBased_MD_Structure(epoints, q, puradius(1, 1), M);

% Build the partitioning structure for data sites and evaluation points
radii = zeros(npu_M, M);
rn = zeros(npu_M, 1); 
epsilon = zeros(npu_M, M);

search_only = 0;
loop_timer = tic;
for j = 1:npu_M
    % Find the block containing the j-th subdomain centre
    index = IntegerBased_MD_ContainingQuery(puctrs(j, :), q, puradius(1, 1), M);
    
    % Determine the optimal parameters (through a random search)
    search_timer = tic;
    fun = @(param2M) Cost_function_mixed( ...
        rbf, param2M, puradius, idx_ds, index, q, M, ...
        puctrs(j, :), dsites, rhs, strat, tikhonov);
    opts = struct('logspace', false, 'num_points', 125);
    bounds = [[.1, .3]; [.1, .3]; [1, 900]; [1, 900]];
    if strcmp(strat, 'no_opt')
        minval = [2 / npu, 1 / npu, 3, 3];  % Used as initial guess in earlier versions of this
    else
        if isotropic
            bounds = [bounds(1, :); bounds(3, :)];
        end
        [minval, fminval, ~, ~] = fminrnd(fun, bounds, opts);
        if isotropic
            minval = [minval(1), minval(1), minval(2), minval(2)];
        end
    end
    search_only = search_only + toc(search_timer);

    % Store the results from the optimization
    radii(j,:) = minval(1:M)'; % The optimal semi-axes
    epsilon(j,:) = minval(M + 1:2 * M)'; % The optimal shape parameters
    if verbose
        fprintf('%d %g %g %g %g %g\n', j, radii(j,1), radii(j,2), epsilon(j,1), epsilon(j,2), fminval)
    end
    
    % Find the mumber of blocks to search for the points on the patches
    n = 1;
    while n < 99999
        if max(radii) > n*max(puradius)
            n = n + 1;
        else
            t = n;
            break
        end
    end
    rn(j) = t;
    
    % Construct anisotropic weights
    % Find data sites located in the j-th subdomain
    [dxx, dx] = IntegerBased_MD_Neighbourhood(dsites, idx_ds, index, q, M, rn(j));
    locpts(j).ind = IntegerBased_MD_RangeSearchAniso(puctrs(j, :), radii(j, :), dxx, dx);
    
    [edxx, edx] = IntegerBased_MD_Neighbourhood(epoints, idx_ep, index, q, M, rn(j));
    elocpts(j).ind = IntegerBased_MD_RangeSearchAniso(puctrs(j, :), radii(j, :), edxx, edx);
end
this_time = toc(loop_timer);
fprintf('First loop time %g, Search component %g\n', this_time, search_only)

% Define the pu weight
epu = PuweightAniso(epoints, elocpts, puctrs, radii);

% Initialize the prediction vector
Pf = zeros(neval_M, 1);

loop_timer = tic;
for j = 1:npu_M
    eval_points = epoints(elocpts(j).ind, :);
    
    if (~isempty(eval_points))
        centers = dsites(locpts(j).ind, :);
        ep = epsilon(j, :);
        
        % Compute the interpolation matrix
        IM = rbf(DistanceMatrixAniso(centers, centers, ep)) + tikhonov * eye(length(centers));
        
        % Compute local evaluation matrix
        EM = rbf(DistanceMatrixAniso(eval_points, centers, ep));
        
        if verbose
            fprintf('j=%g, local N-%g\n', j, length(IM));
        end
        
        % Compute the local interpolant
        localfit = EM * (IM \ rhs(locpts(j).ind));
        
        % Accumulate the local interpolant into the global fit
        Pf(elocpts(j).ind) = Pf(elocpts(j).ind) + localfit .* epu(j).w;
    end
end
this_time = toc(loop_timer);
fprintf('Second loop time %g\n', this_time)

% Compute errors on evaluation grid
exact = f(epoints);
maxerr = norm(Pf - exact, inf);
rms_err = norm(Pf - exact) / sqrt(neval_M);
fprintf('RMS error:       %e\n', rms_err);
fprintf('Maximum error:   %e\n', maxerr);