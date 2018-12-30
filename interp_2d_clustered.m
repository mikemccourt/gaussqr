% Consider the impact of clustered data on the isotropic symmetric setting

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
fontsize = 18;
Nvec = logspace(1, 4, 23);
tau = sqrt(10);
rbf = @(r) exp(-r.^2);
num_splits = 15;

yf = @(x) ( ...
    exp(-3 * (x(:, 1) .^ 2 + x(:, 2) .^ 2)) + ...
    1.2 * exp(-4 * (.7 * (x(:, 1) - .5) .^ 2 + 1.3 * (x(:, 1) - .5) .^ 2))...
);
results = zeros([num_splits, length(Nvec)]);

Ncount = 1;
for N=Nvec
    xtight1 = pick2Dpoints(-.2, .2, [floor(N/3), 1], 'halton');
    xtight2 = pick2Dpoints(-.75, -.35, [floor(N/3), 1], 'halton');
    x = [
        xtight1;
        xtight2;
        pick2Dpoints(-1, 1, [N - length(xtight1) - length(xtight2), 1], 'halton');
    ];
    xeval = [
        pick2Dpoints(-.15, .15, [200, 1], 'halton');
        pick2Dpoints(-.7, -.4, [200, 1], 'halton');
    ];
    DM = DistanceMatrix(x, x);
    DMeval = DistanceMatrix(xeval, x);
    
    y = yf(x);
    yeval = yf(xeval);
    ep_base = log(N);
    epvec = exp(linspace(log(ep_base / tau), log(ep_base * tau), num_splits));
    
    fprintf('%d\t', Ncount)
    
    tcount = 1;
    for ep=epvec
        warning('off', 'MATLAB:nearlySingularMatrix')
        ypred = rbf(ep * DMeval) * (rbf(ep * DM) \ y);
        warning('on', 'MATLAB:nearlySingularMatrix')
        results(tcount, Ncount) = errcompute(ypred, yeval);
        fprintf('%d ', tcount)
        tcount = tcount + 1;
    end
    
    fprintf('\n')
    
    Ncount = Ncount + 1;
end
    
clf reset
handles = zeros([1, num_splits]);
hold on
for tcount=1:num_splits
    color = sigopt_teal * (1 - (tcount - 1) / (num_splits - 1)) + sigopt_purple * (tcount - 1) / (num_splits - 1);
    handles(tcount) = plot(Nvec, results(tcount, :), 'linewidth', 3, 'color', color);
end
hold off

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([1e1, 1e4])
ylim([1e-12, 1e0])
xlabel('N - number of points sampled', 'fontsize', fontsize, 'interpreter', 'tex')
ylabel('RMSE', 'fontsize', fontsize)
xticks([1e1, 1e2, 1e3, 1e4])
yticks([1e-10, 1e-5, 1e0])
set(gca, 'fontsize', fontsize)
legend([handles(end), handles(1)], {'largest shape parameter', 'smallest shape parameter'}, ...
    'location', 'southwest', 'fontsize', fontsize)

filename = 'examples_2d_clustered_spread';
savefig(filename);saveas(gcf, filename, 'png')