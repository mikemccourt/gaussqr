% Allow for using rotated kernels in interpolation

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
Nvec = floor(logspace(1, 4, 23));
Neval = 400;
varvec = [.1, 1, 10, 100];
rbf = @(r) exp(-r.^2);
rbf2 = @(r2) exp(-r2);
num_runs = 25;

yf = @(x) sin(6 * (x(:, 1) .^ 2 + x(:, 2) .^ 2));
results = zeros([1, length(Nvec)]);
resultsall = cell(1, length(varvec));
for k=1:length(varvec)
    resultsall{k} = zeros([num_runs, length(Nvec)]);
end

Ncount = 1;
warning('off', 'MATLAB:nearlySingularMatrix')
for N=Nvec
    x = pick2Dpoints(-1, 1, [N, 1], 'halton');
    xeval = pick2Dpoints(min(x), max(x), [Neval, 1], 'halton');
    DM2 = zeros(N);
    DMeval2 = zeros(Neval, N);
    
    y = yf(x);
    yeval = yf(xeval);
    ep = log(N);
    Keval = rbf(ep * DistanceMatrix(xeval, x));
    K = rbf(ep * DistanceMatrix(x, x));
    ypred = Keval * (K \ y);
    results(Ncount) = errcompute(ypred, yeval);
    
    fprintf('%d\t', Ncount)
    
    for tcount=1:num_runs
        vcount = 1;
        for v=varvec
            a = ep^2 / v;
            b = v / ep;
            epvec = gamrnd(a, b, N, 2);
            anglevec = pickpoints(-pi, pi, N, 'rand');
            
            for col=1:N
                V = [[cos(anglevec(col)), -sin(anglevec(col))];
                    [sin(anglevec(col)), cos(anglevec(col))]];
                DM2(:, col) = sum((((x - x(col, :)) * V) .* epvec(col, :)) .^ 2, 2);
                DMeval2(:, col) = sum((((xeval - x(col, :)) * V) .* epvec(col, :)) .^ 2, 2);
            end
            ypred = rbf2(DMeval2) * (rbf2(DM2) \ y);
            resultsall{vcount}(tcount, Ncount) = errcompute(ypred, yeval);
            vcount = vcount + 1;
        end
        fprintf('%d ', tcount)
    end
    
    fprintf('\n')
    Ncount = Ncount + 1;
end
warning('on', 'MATLAB:nearlySingularMatrix')
    
clf reset

colors = [sigopt_green; sigopt_medium_blue; sigopt_orange; sigopt_dark_gray];
hcount = 1;
handles = zeros(1, length(varvec) + 1);
handles(hcount) = plot(Nvec, results, '--k', 'linewidth', 3);
hold on

for vcount=1:length(varvec)
    hcount = hcount + 1;
    r = resultsall{vcount};
    color = colors(vcount, :);
    bot = prctile(r, 25, 1);
    med = prctile(r, 50, 1);
    top = prctile(r, 75, 1);
    handles(hcount) = fill_between_lines(Nvec, bot, top, color);
    set(handles(hcount), 'linestyle', 'none')
    set(handles(hcount), 'facealpha', .2)
    handles(hcount) = plot(Nvec, med, 'color', color, 'linewidth', 3);
    plot(Nvec, bot, ':', 'color', color, 'linewidth', 2)
    plot(Nvec, top, ':', 'color', color, 'linewidth', 2)
end

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim([1e1, 1e4])
ylim([1e-12, 1e0])
xlabel('N - number of points sampled', 'fontsize', fontsize, 'interpreter', 'tex')
ylabel('RMSE', 'fontsize', fontsize)
xticks([1e1, 1e2, 1e3, 1e4])
yticks([1e-10, 1e-5, 1e0])
set(gca, 'fontsize', fontsize)
labels = arrayfun(@(v) sprintf('Var=%4.1f', v), [0, varvec], 'uniformoutput', 0);
legend(handles, labels, ...
    'location', 'southwest', 'fontsize', fontsize)
hold off

filename = 'examples_2d_rotated';
savefig(filename)
saveas(gcf, filename, 'png')