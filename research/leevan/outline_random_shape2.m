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

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

% First, a 1D problem
Nvec = floor(linspace(10, 500, 13));
yf = @(x) 1 ./ (1 + 25 * (x(:, 1) - .5) .^ 2);
% yf = @(x) sin(6 * x) - x;
xeval = pickpoints(-1, 1, 100);
yeval = yf(xeval);

% Only working with the Gaussian RBF for now
tvec = [2, 6, 20];
rbf = @(r) exp(-r.^2);
num_runs = 30;
t_acc = zeros(length(tvec), length(Nvec));
t_bot = zeros(length(tvec), length(Nvec));
t_top = zeros(length(tvec), length(Nvec));
base_acc = zeros(size(Nvec));

j = 1;
for N=Nvec
    x = pickpoints(-1, 1, N, 'cheb');  % Could allow this to change
    y = yf(x);
    
    % Define an epsilon center and some log region around it
    opt = struct('strat', 'loocv', 'output', 'opt', 'svd_tol', 1e-13);
    
    epbase = sqrt(N);
    bounds = [.05 * epbase, 20 * epbase];
    ep = choose_ep_rand(x, y, rbf, bounds, opt);
    fprintf('ep=%g, [%g, %g]\n', ep, bounds(1), bounds(2))
    
    DM = DistanceMatrix(x, x);
    DMeval = DistanceMatrix(xeval, x);
    ypred = rbf(ep * DMeval) * random_svd_solve(rbf(ep * DM), y);
    base_acc(j) = errcompute(ypred, yeval);
    
    tcount = 1;
    for t=tvec
        ep_rand_min = ep / t;
        ep_rand_max = ep * t;
        results = zeros(1, num_runs);
        for k=1:num_runs
            epvec = exp(pickpoints(log(ep_rand_min), log(ep_rand_max), N, 'rand'))';
            ypred_rnd = rbf(bsxfun(@times, DMeval, epvec)) * random_svd_solve(rbf(bsxfun(@times, DM, epvec)), y);
            results(k) = errcompute(ypred_rnd, yeval);
        end
        t_acc(tcount, j) = median(results);
        t_bot(tcount, j) = prctile(results, 25);
        t_top(tcount, j) = prctile(results, 75);
        tcount = tcount + 1;
    end
    j = j + 1;
end

hbase = loglog(Nvec, base_acc, '--k', 'linewidth', 3);
hold on
colors = ['b', 'r', 'g'];  % Need to generalize this somehow
handles = zeros(size(colors));
for k=1:length(colors)
    handles(k) = plot(Nvec, t_acc(k, :), colors(k), 'linewidth', 3);
    
    hh = fill_between_lines(Nvec, t_bot(k, :), t_top(k, :), colors(k));
    set(hh, 'FaceAlpha', .1);
    set(hh, 'EdgeAlpha', .1);
    hhtop = plot(Nvec, t_top(k, :), ':', 'color', colors(k), 'linewidth', 3);
    hhbot = plot(Nvec, t_bot(k, :), ':', 'color', colors(k), 'linewidth', 3);
end
hold off

xlabel(sprintf('t, %3.2f / t < ep < %3.2f * t', ep, ep))
ylabel('RMSE')
legend_vals = arrayfun(@(tt)sprintf('t=%g',tt), [1, tvec], 'uniformoutput', 0);
legend([hbase, handles], legend_vals, 'location', 'southwest', 'fontsize', 14)
