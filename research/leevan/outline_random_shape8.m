% Trying to redo some shit I already did
% Create image for singular value distribution

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

clf reset
fontsize = 14;
Nvec = floor(logspace(1, 4, 23));
Neval = 400;
tau1 = 1.1;
tau2 = sqrt(10);
rbf = @(r) exp(-r.^2);
num_runs = 30;

yf = @(x) sin(6 * (x(:, 1) .^ 2 + x(:, 2) .^ 2));
results = zeros([1, length(Nvec)]);
results1 = zeros([num_runs, length(Nvec)]);
results2 = zeros([num_runs, length(Nvec)]);

Ncount = 1;
warning('off', 'MATLAB:nearlySingularMatrix')
for N=Nvec
    x = pick2Dpoints(-1, 1, [N, 1], 'halton');
    xeval = pick2Dpoints(min(x), max(x), [400, 1], 'halton');
    DM = zeros(N);
    DMeval = zeros(Neval, N);
    
    y = yf(x);
    yeval = yf(xeval);
    ep = log(N);
    ypred = rbf(ep * DistanceMatrix(xeval, x)) * (rbf(ep * DistanceMatrix(x, x)) \ y);
    results(Ncount) = errcompute(ypred, yeval);
    
    for tcount=1:num_runs
        for col=1:N
            epvec = exp(pick2Dpoints(log(ep / tau1), log(ep * tau1), 1, 'rand'))';
            DM(:, col) = DistanceMatrix(x, x(col, :), epvec');
            DMeval(:, col) = DistanceMatrix(xeval, x(col, :), epvec');
        end
        ypred = rbf(DMeval) * (rbf(DM) \ y);
        results1(tcount, Ncount) = errcompute(ypred, yeval);
        for col=1:N
            epvec = exp(pick2Dpoints(log(ep / tau2), log(ep * tau2), 1, 'rand'))';
            DM(:, col) = DistanceMatrix(x, x(col, :), epvec');
            DMeval(:, col) = DistanceMatrix(xeval, x(col, :), epvec');
        end
        ypred = rbf(DMeval) * (rbf(DM) \ y);
        results2(tcount, Ncount) = errcompute(ypred, yeval);
    end
    
    Ncount = Ncount + 1;
end
warning('on', 'MATLAB:nearlySingularMatrix')

t1_bot = prctile(results1, 25, 1);
t1_med = prctile(results1, 50, 1);
t1_top = prctile(results1, 75, 1);
t2_bot = prctile(results2, 25, 1);
t2_med = prctile(results2, 50, 1);
t2_top = prctile(results2, 75, 1);
    
clf reset
hold on
color_blue = [31,64,125] / 255;
color_orange = [248,155,32] / 255;

hfill1 = fill_between_lines(Nvec, t1_bot, t1_top, color_orange);
set(hfill1, 'linestyle', 'none')
set(hfill1, 'facealpha', .2)
hfill1 = plot(Nvec, t1_med, 'color', color_orange, 'linewidth', 3);
plot(Nvec, t1_bot, ':', 'color', color_orange, 'linewidth', 2)
plot(Nvec, t1_top, ':', 'color', color_orange, 'linewidth', 2)

hfill2 = fill_between_lines(Nvec, t2_bot, t2_top, color_blue);
set(hfill2, 'linestyle', 'none')
set(hfill2, 'facealpha', .2)
hfill2 = plot(Nvec, t2_med, 'color', color_blue, 'linewidth', 3);
plot(Nvec, t2_bot, ':', 'color', color_blue, 'linewidth', 2)
plot(Nvec, t2_top, ':', 'color', color_blue, 'linewidth', 2)

hbase = plot(Nvec, results, '--k', 'linewidth', 3);

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([1e1, 1e4])
ylim([1e-12, 1e0])
xlabel('N - number of points sampled', 'fontsize', fontsize, 'interpreter', 'tex')
ylabel('normalized RMSE', 'fontsize', fontsize)
xticks([1e1, 1e2, 1e3, 1e4])
yticks([1e-10, 1e-5, 1e0])
set(gca, 'fontsize', fontsize)
legend([hbase, hfill1, hfill2], {'$\tau=1$ ($\varepsilon=\log(N)$)', '$\tau=1.1$', '$\tau=\sqrt{10}$'}, 'location', 'southwest', 'fontsize', fontsize, 'interpreter', 'latex')
hold off

filename = 'examples_2d_random';
savefig(filename)
saveas(gcf, filename, 'png')