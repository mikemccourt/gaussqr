function [x_sol, f_sol, x, f] = fminrnd(fun, bounds, opts)

logspace = true;
num_points = 100;
polish = false;
if nargin == 3
    if isfield(opts, 'logspace')
        logspace = opts.logspace;
    end
    if isfield(opts, 'num_points')
        num_points = opts.num_points;
    end
    if isfield(opts, 'polish')
        polish = opts.polish;
        error('Not working yet')
    end
end

if logspace
    bounds = log(bounds);
end

dim = length(bounds);
p = haltonset(dim);
unscaled_points = mod(bsxfun(@plus, net(p, num_points), rand(1, dim)), 1);

x_width = bounds(:, 2) - bounds(:, 1);
x_min = bounds(:, 1);
x = bsxfun(@plus, bsxfun(@times, unscaled_points, x_width'), x_min');

if logspace
    x = exp(x);
end

f = zeros(num_points, 1);
for k=1:num_points
    f(k) = fun(x(k, :));
end

[f_sol, ind_sol] = min(f);
x_sol = x(ind_sol, :);

end