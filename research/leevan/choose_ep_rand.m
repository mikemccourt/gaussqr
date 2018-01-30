function [ep, optvals] = choose_ep_rand(x, y, rbf, bounds, opt)
assert(strcmp(opt.strat, 'loocv') || strcmp(opt.strat, 'acc') || strcmp(opt.strat, 'mle'))
assert(strcmp(opt.output, 'opt') || strcmp(opt.output, 'plot'))

DM = DistanceMatrix(x, x);
if strcmp(opt.strat, 'loocv')
    fopt = @(ep) choose_ep_rand_loocv(ep, DM, y, rbf, opt.svd_tol);
elseif strcmp(opt.strat, 'mle')
    fopt = @(ep) choose_ep_rand_mle(ep, DM, y, rbf, opt.svd_tol);
else
    DMeval = DistanceMatrix(opt.xeval, x);
    fopt = @(ep) choose_ep_rand_acc(ep, DM, y, rbf, opt.svd_tol, DMeval, opt.yeval);
end

if strcmp(opt.output, 'opt')
    [ep, optvals] = fminbnd(fopt, bounds(1), bounds(2));
else
    ep = logspace(log10(bounds(1)), log10(bounds(2)), opt.plotnum);
    optvals = arrayfun(fopt, ep);
end
end

function loocv = choose_ep_rand_loocv(ep, DM, y, rbf, svd_tol)
N = length(DM);
I_and_y = [eye(N), y];
Kinv_and_c = random_svd_solve(rbf(ep * DM), I_and_y, svd_tol);
Kinv = Kinv_and_c(:, 1:N);
c = Kinv_and_c(:, N + 1);
loocv = sum(abs(c ./ diag(Kinv)));
end

function log_like = choose_ep_rand_mle(ep, DM, y, rbf, svd_tol)
N = length(DM);
[c, logdet] = random_svd_solve(rbf(ep * DM), y, svd_tol);
log_like = N * log(y' * c) + logdet;
end

function acc = choose_ep_rand_acc(ep, DM, y, rbf, svd_tol, DMeval, yeval)
N = length(DM);
c = random_svd_solve(rbf(ep * DM), y, svd_tol);
acc = errcompute(rbf(ep * DMeval) * c, yeval);
end