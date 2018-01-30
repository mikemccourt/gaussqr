% Build experiments with SVD + regularization if desired
% One set of experiments should involve sampling from LOOCV chi-squared
% Experiments involving uniform distribution
% Sampling from the likelihood
% 
% Strategy: Find a point which is "representative" of the density
%   Consider a density expanded around that point
%   Randomly sample with each dimension and point
% Analysis: Consider the impact of the spread on the quality of the
%   approximation and/or the stability of the system
%
% Will need to run experiments in different dimensions
% Any differential equations examples?

% First, a 1D problem
N = 11;
yf = @(x) sin(6 * x) - x;
x = pickpoints(-1, 1, N, 'halton');
y = yf(x);
xeval = pickpoints(-1, 1, 100);
yeval = yf(xeval);

% Only working with the Gaussian RBF for now
rbf = @(r) exp(-r.^2);

% Define an epsilon center manually and some log region around it
ep = .05 * sqrt(N);
DM = DistanceMatrix(x, x);
DMeval = DistanceMatrix(xeval, x);
ypred = rbf(ep * DMeval) * random_svd_solve(rbf(ep * DM), y);
base_acc = errcompute(ypred, yeval);

num_runs = 30;
t_vals = logspace(0, 3, 40);
t_acc = zeros(size(t_vals));
j = 1;
for t=t_vals;
    ep_rand_min = ep / t;
    ep_rand_max = ep * t;
    results = zeros(1, num_runs);
    for k=1:num_runs
        epvec = exp(pickpoints(log(ep_rand_min), log(ep_rand_max), N, 'rand'))';
        ypred_rnd = rbf(bsxfun(@times, DMeval, epvec)) * random_svd_solve(rbf(bsxfun(@times, DM, epvec)), y);
        results(k) = errcompute(ypred_rnd, yeval);
    end
    t_acc(j) = median(results);
    j = j + 1;
end

loglog(t_vals, ones(size(t_vals)) * base_acc, '--k', 'linewidth', 3)
hold on
plot(t_vals, t_acc, 'linewidth', 3)
hold off
xlabel(sprintf('t, %3.2f / t < ep < %3.2f * t', ep, ep))
ylabel('MSE')
legend(sprintf('ep=%4.2f', ep), 'random', 'location', 'west')
