% Numerical results for different distributions on a small problem

% These are standard colors I plot with now (SigOpt brand colors)
sigopt_medium_blue = [0, 71, 187] / 255;
sigopt_aqua = [0, 134, 191] / 255;
sigopt_blue = [27, 54, 93] / 255;
sigopt_purple = [135, 24, 157] / 255;
sigopt_magenta = [187, 41, 187] / 255;
sigopt_orange = [255, 130, 0] / 255;
sigopt_yellow = [253, 218, 36] / 255;
sigopt_light_green = [151, 215, 0] / 255;
sigopt_green = [0, 177, 64] / 255;
sigopt_teal = [0, 164, 153] / 255;
sigopt_dark_gray = [83, 86, 90] / 255;

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

clf reset
fontsize = 14;
N = 131;
vvec = logspace(0, 2, 20) - 1 + .1;
x = pickpoints(-1, 1, N);
xeval = pickpoints(-1, 1, 100);
DM = DistanceMatrix(x, x);
DMeval = DistanceMatrix(xeval, x);
rbf = @(r) exp(-r.^2);
num_runs = 250;

stuff = { ...
    {'loocv', @(x) sin(6 * x) - x}, ...
    {'mle', @(x) 1 - tanh(3 * x) + exp(x)}, ...
    {'acc', @(x) 1 ./ (1 + 25 * (x - .3) .^ 2)}, ...
    };

ep_from_other = [6.01, 16.76, 10.59];
v = 10;

for k=1:length(stuff)
    strat = stuff{k}{1};
    yf = stuff{k}{2};
    
    ep = ep_from_other(k);
    IM = rbf(ep * DM);
    base_singular_values = svd(IM);
    
    results_uniform = zeros(N, num_runs);
    results_normal = zeros(N, num_runs);
    results_chi2 = zeros(N, num_runs);
    results_ga = zeros(N, num_runs);
    for m=1:num_runs
        % Log-uniform
        t = fzero(invt(v, ep), [.00001, .9999]);
        nu = 2 * ep * log(t) / (t - 1 / t);
        invt = @(v,m) (@(t) log(t) .* (t + 1./t)./(t - 1./t) - 1 - v/m^2);
        epvec = exp(pickpoints(log(nu / t), log(nu * t), N, 'rand'))';
        IM = rbf(bsxfun(@times, DM, epvec));
        results_uniform(:, m) = svd(IM);
        
        % Log-normal
        sigma = sqrt(log(v / ep^2 + 1));
        mu = log(ep^2 / sqrt(v + ep^2));
        epvec = lognrnd(mu, sigma, 1, N);
        IM = rbf(bsxfun(@times, DM, epvec));
        results_normal(:, m) = svd(IM);
        
        % Chi-squared
        epvec = chi2rnd(ep, 1, N);
        IM = rbf(bsxfun(@times, DM, epvec));
        results_chi2(:, m) = svd(IM);
        
        % Gamma
        a = ep^2 / v;
        b = v / ep;
        epvec = gamrnd(a, b, 1, N);
        IM = rbf(bsxfun(@times, DM, epvec));
        results_ga(:, m) = svd(IM);
    end
    
    t_bot_lu = prctile(results_uniform, 25, 2);
    t_med_lu = prctile(results_uniform, 50, 2);
    t_top_lu = prctile(results_uniform, 75, 2);
    t_bot_ln = prctile(results_normal, 25, 2);
    t_med_ln = prctile(results_normal, 50, 2);
    t_top_ln = prctile(results_normal, 75, 2);
    t_bot_c2 = prctile(results_chi2, 25, 2);
    t_med_c2 = prctile(results_chi2, 50, 2);
    t_top_c2 = prctile(results_chi2, 75, 2);
    t_bot_ga = prctile(results_ga, 25, 2);
    t_med_ga = prctile(results_ga, 50, 2);
    t_top_ga = prctile(results_ga, 75, 2);
    
    subplot(1, length(stuff), k)
    hbase = plot(1:N, base_singular_values, '--k', 'linewidth', 3);
    hold on
    hfill_lu = fill_between_lines(1:N, t_bot_lu', t_top_lu', sigopt_medium_blue);
    set(hfill_lu, 'linestyle', 'none')
    set(hfill_lu, 'facealpha', .2)
    hfill_ln = fill_between_lines(1:N, t_bot_ln', t_top_ln', sigopt_magenta);
    set(hfill_ln, 'linestyle', 'none')
    set(hfill_ln, 'facealpha', .2)
    hfill_c2 = fill_between_lines(1:N, t_bot_c2', t_top_c2', sigopt_orange);
    set(hfill_c2, 'linestyle', 'none')
    set(hfill_c2, 'facealpha', .2)
    hfill_ga = fill_between_lines(1:N, t_bot_ga', t_top_ga', sigopt_green);
    set(hfill_ga, 'linestyle', 'none')
    set(hfill_ga, 'facealpha', .2)
    set(gca, 'yscale', 'log')
    xlim([1, 131])
    ylim([1e-15, 1e2])
    xlabel('singular value index', 'fontsize', fontsize)
    ylabel('singular value', 'fontsize', fontsize)
    xticks([1, 65, 131])
    yticks([1e-15, 1e-10, 1e-5, 1e0])
    set(gca, 'fontsize', fontsize)
    legend([hbase, hfill_lu, hfill_ln, hfill_c2, hfill_ga], ...
        {  sprintf('$\\varepsilon=%3.2f$', ep), ...
           'log-uniform', ...
           'log-normal', ...
           'chi-square', ...
           'gamma', ...
        }, ...
        'location', 'northeast', 'interpreter', 'latex', 'fontsize', fontsize)
    hold off
    
    fig = gcf;
    ax = gca;
    bottom = .15;
    ax_height = .8;
    ax_width = .26;
    if k == 1
        left = .06;
    elseif k == 2
        left = .06 + .07 + ax_width;
        legend('location', 'southwest')
    elseif k == 3
        left = .06 + 2 * .07 + 2 * ax_width;
    end
    ax.Position = [left bottom ax_width ax_height];
    fig.Position = [560         528        1080         320];
end

savefig('examples_1d_svd_LARGE')
saveas(gcf, 'examples_1d_svd_LARGE', 'png')