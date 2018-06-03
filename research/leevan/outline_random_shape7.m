% Trying to redo some shit I already did
% Create image for singular value distribution

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C);

clf reset
fontsize = 14;
Nvec = logspace(1, 4, 23);
tau = sqrt(10);
rbf = @(r) exp(-r.^2);
num_splits = 15;
num_runs = 30;

yf = @(x) sin(6 * (x(:, 1) .^ 2 + x(:, 2) .^ 2));
results = zeros([num_splits, length(Nvec)]);

Ncount = 1;
for N=Nvec
    x = pick2Dpoints(-1, 1, [N, 1], 'halton');
    xeval = pick2Dpoints(min(x), max(x), [400, 1], 'halton');
    DM = DistanceMatrix(x, x);
    DMeval = DistanceMatrix(xeval, x);
    
    y = yf(x);
    yeval = yf(xeval);
    ep_base = log(N);
    epvec = exp(linspace(log(ep_base / tau), log(ep_base * tau), num_splits));
    
    tcount = 1;
    for ep=epvec
        warning('off', 'MATLAB:nearlySingularMatrix')
        ypred = rbf(ep * DMeval) * (rbf(ep * DM) \ y);
        warning('on', 'MATLAB:nearlySingularMatrix')
        results(tcount, Ncount) = errcompute(ypred, yeval);
        tcount = tcount + 1;
    end
    
    Ncount = Ncount + 1;
end
    
clf reset
handles = zeros([1, num_splits]);
red = [1, 0, 0];
blue = [0, 0, 1];
hold on
for tcount=1:num_splits
    color = red * (1 - (tcount - 1) / (num_splits - 1)) + blue * (tcount - 1) / (num_splits - 1);
    handles(tcount) = plot(Nvec, results(tcount, :), 'linewidth', 3, 'color', color);
end
hold off

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([1e1, 1e4])
ylim([1e-12, 1e0])
xlabel('N - number of points sampled', 'fontsize', fontsize, 'interpreter', 'tex')
ylabel('normalized RMSE', 'fontsize', fontsize)
xticks([1e1, 1e2, 1e3, 1e4])
yticks([1e-10, 1e-5, 1e0])
set(gca, 'fontsize', fontsize)
legend([handles(end), handles(1)], {'largest shape parameter', 'smallest shape parameter'}, 'location', 'southwest', 'fontsize', fontsize)

savefig('examples_2d_spread')
saveas(gcf, 'examples_2d_spread', 'png')