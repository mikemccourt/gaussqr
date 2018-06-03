% Trying to redo some shit I already did
% Create image for singular value distribution

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

clf reset
fontsize = 14;
N = 131;
tvec = logspace(.001, 1, 16);
x = pickpoints(-1, 1, N);
xeval = pickpoints(-1, 1, 100);
DM = DistanceMatrix(x, x);
DMeval = DistanceMatrix(xeval, x);
rbf = @(r) exp(-r.^2);
num_runs = 200;

stuff = { ...
    {'loocv', @(x) sin(6 * x) - x}, ...
    {'mle', @(x) 1 - tanh(3 * x) + exp(x)}, ...
    {'acc', @(x) 1 ./ (1 + 25 * (x - .3) .^ 2)}, ...
    };

for k=1:length(stuff)
    strat = stuff{k}{1};
    yf = stuff{k}{2};
    
    t_med = zeros(1, length(tvec));
    t_bot = zeros(1, length(tvec));
    t_top = zeros(1, length(tvec));
    j = 1;
    
    y = yf(x);
    yeval = yf(xeval);
    opt = struct('strat', strat, 'output', 'opt', 'svd_tol', 1e-13, 'xeval', xeval, 'yeval', yeval);
    bounds = [.001, 100];
    ep = choose_ep_rand(x, y, rbf, bounds, opt);
    
    ypred = rbf(ep * DMeval) * random_svd_solve(rbf(ep * DM), y);
    base_acc = errcompute(ypred, yeval);
    
    tcount = 1;
    for t=tvec
        ep_rand_min = ep / t;
        ep_rand_max = ep * t;
        results = zeros(1, num_runs);
        for m=1:num_runs
            epvec = exp(pickpoints(log(ep_rand_min), log(ep_rand_max), N, 'rand'))';
            ypred_rnd = rbf(bsxfun(@times, DMeval, epvec)) * random_svd_solve(rbf(bsxfun(@times, DM, epvec)), y);
            results(m) = errcompute(ypred_rnd, yeval);
        end
        t_bot(tcount) = prctile(results, 25);
        t_med(tcount) = prctile(results, 50);
        t_top(tcount) = prctile(results, 75);
        tcount = tcount + 1;
    end
    
    subplot(1, length(stuff), k)
    hfill = fill_between_lines(tvec, t_bot, t_top, 'r');
    hold on
    set(hfill, 'linestyle', 'none')
    set(hfill, 'facealpha', .2)
    plot(tvec, t_med, 'r', 'linewidth', 3)
    plot(tvec, t_bot, ':r', 'linewidth', 2)
    plot(tvec, t_top, ':r', 'linewidth', 2)
    hsol = plot(tvec, base_acc * ones(1, length(tvec)), '--k', 'linewidth', 3);
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlim([1, 10])
    ylim([1e-15, 1e-6])
    xlabel('$\tau$', 'fontsize', fontsize, 'interpreter', 'latex')
    ylabel('singular values', 'fontsize', fontsize)
    xticks([1, 10])
    yticks([1e-14, 1e-10, 1e-6])
    set(gca, 'fontsize', fontsize)
    legend([hsol, hfill], {sprintf('$\\varepsilon=%3.2f$', ep), 'random'}, 'location', 'best', 'interpreter', 'latex', 'fontsize', fontsize)
    hold off
    
    j = j + 1;
end

savefig('examples_1d')
saveas(gcf, 'examples_1d', 'png')