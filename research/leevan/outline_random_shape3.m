% Trying to redo some shit I already did
% Create image for singular value distribution

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

clf reset
fontsize = 14;
N = 31;
x = pickpoints(-1, 1, N);
xeval = pickpoints(-1, 1, 100);
DM = DistanceMatrix(x, x);
DMeval = DistanceMatrix(xeval, x);
rbf = @(r) exp(-r.^2);
num_runs = 100;

stuff = { ...
    {'loocv', @(x) sin(6 * x) - x}, ...
    {'mle', @(x) 1 - tanh(3 * x) + exp(x)}, ...
    {'acc', @(x) 1 ./ (1 + 25 * (x - .3) .^ 2)}, ...
    };

for k=1:length(stuff)
    strat = stuff{k}{1};
    yf = stuff{k}{2};
    y = yf(x);
    yeval = yf(xeval);
    tau = sqrt(10);
    opt = struct('strat', strat, 'output', 'opt', 'svd_tol', 1e-13, 'xeval', xeval, 'yeval', yeval);
    bounds = [.001, 100];
    ep = choose_ep_rand(x, y, rbf, bounds, opt);
    ep_rand_min = ep / tau;
    ep_rand_max = ep * tau;
    
    K = rbf(ep * DM);
    [~, S, ~] = svd(K);
    Svals = diag(S);
    
    Sresults = zeros(N, num_runs);
    results = zeros(num_runs);
    for m=1:num_runs
        epvec = exp(pickpoints(log(ep_rand_min), log(ep_rand_max), N, 'rand'))';
        Kr = rbf(bsxfun(@times, DM, epvec));
        [~, S, ~] = svd(Kr);
        Sresults(:, m) = diag(S);
    end
    
    Sbot = prctile(Sresults, 25, 2);
    Smed = prctile(Sresults, 50, 2);
    Stop = prctile(Sresults, 75, 2);
    
    subplot(1, length(stuff), k)
    hfill = fill_between_lines(0:N-1, Sbot', Stop', 'r');
    hold on
    set(hfill, 'linestyle', 'none')
    set(hfill, 'facealpha', .2)
    semilogy(0:N-1, Smed', 'r', 'linewidth', 3)
    plot(0:N-1, Sbot', ':r', 'linewidth', 2)
    plot(0:N-1, Stop', ':r', 'linewidth', 2)
    hsol = semilogy(0:N-1, Svals, '--k', 'linewidth', 3);
    set(gca, 'yscale', 'log')
    xlim([0, N-1])
    ylim([1e-16, 1e2])
    ylabel('singular values', 'fontsize', fontsize)
    xticks([0, 10, 20, 30])
    yticks([1e-15, 1e-10, 1e-5, 1e0])
    set(gca, 'fontsize', fontsize)
    legend([hsol, hfill], {sprintf('$\\varepsilon=%3.2f$', ep), 'random, $\tau=\sqrt{10}$'}, 'location', 'southwest', 'interpreter', 'latex', 'fontsize', fontsize)
    hold off
end

savefig('examples_1d_TINY_svd')
saveas(gcf, 'examples_1d_TINY_svd', 'png')