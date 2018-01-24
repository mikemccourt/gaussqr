% GP Sampling Fundamentals
% It requires the software with Fasshauer/McCourt (could be replaced):
%   - pick2Dpoints
%   - DistanceMatrix
%   - DifferenceMatrix
yf = @(x) 1 ./ (1 + 25 * (x(:, 1) - .5) .^ 2);
yfx = @(x) -50 * (x(:, 1) - .5) ./ (1 + 25 * (x(:, 1) - .5) .^ 2) .^ 2;

x = [.1, .28, .44, .61, .77, .92]';
y = yf(x);
xplt = pickpoints(0, 1, 200, 'even');
yplt = yf(xplt);

% It is assumed that some amount of noise polluted the observations
% This noise is assumed (for now) to be Gaussian for simplicity
noise_level = 1e-3;
y = y + sqrt(noise_level) * randn(length(y), 1);

% We will restrict ourselves to radial kernels for now ... can switch later
% Also we will only have a single shape parameter for simplicity
% The "process variance" is only useful in a statistical setting
process_variance = 10.0;
ep = 1.5;
rbf = @(r) process_variance * exp(-(ep * r) .^ 2);

% Consider the standard numerical analysis computations
% We have added the noise_level to the diagonal
Kint = rbf(DistanceMatrix(x, x)) + noise_level * eye(length(y));
Keval = rbf(DistanceMatrix(xplt, x));
yeval = Keval * (Kint \ y);

% This is the kriging variance (power function)
% This will be more complicated for a nonradial kernel
% We want the standard deviation so we square root it here
ksd = sqrt(max(rbf(0) - sum((Keval / Kint) .* Keval, 2), 0));

% This is a function that helps for plotting confidence intervals
fill_between_lines = @(X,Y1,Y2,C) fill( [X' fliplr(X')],  [Y1' fliplr(Y2')], C);

subplot(2, 2, 1)
h1 = plot(xplt, yplt, 'r', 'linewidth', 3);
hold on
h2 = plot(x, y, 'or', 'linewidth', 3);
h3 = plot(xplt, yeval, 'b', 'linewidth', 3);

h4 = fill_between_lines(xplt, yeval - 2 * ksd, yeval + 2 * ksd, 'b');
set(h4, 'FaceAlpha', .1);
set(h4, 'EdgeAlpha', .1);
h4top = plot(xplt, yeval + 2 * ksd, ':b', 'linewidth', 3);
h4top.Color = [0,0,1,0.2];
h4bot = plot(xplt, yeval - 2 * ksd, ':b', 'linewidth', 3);
h4bot.Color = [0,0,1,0.2];
legend([h1, h2, h3], {'true', 'data', 'pred'}, 'location', 'south')
title('2 SD confidence intervals')
hold off
ylim([0, 1])

% This is a sampling of draws from the posterior
% yeval is the posterior mean at the points xplt
% We need to compute the posterior covariance at xplt
% The wiggle room allows for ill-conditioning (could also use SVD)
wiggle_room = 1e-10;
Kplt = rbf(DistanceMatrix(xplt, xplt));
Kpostcov = Kplt - (Keval / Kint) * Keval';
Kpostcov = Kpostcov + wiggle_room * eye(size(xplt, 1));
Kpostcov_chol = chol(Kpostcov, 'lower');

n_draws = 100;
normal_draws = randn(size(xplt, 1), n_draws);
ypost = bsxfun(@plus, yeval, Kpostcov_chol * normal_draws);

subplot(2, 2, 2)
h21 = plot(xplt, ypost, 'b');
hold on
for k=1:length(h21)
    h21(k).Color = [0, 0, 1, 0.1];
end
title(sprintf('%d posterior draws', n_draws))
h22 = plot(x, y, 'or', 'linewidth', 3);
hold off
ylim([0, 1])

% If we want to consider the gradient we will need rbf gradients
% Since this is in 1D it is much easier to manage
% We will, however, need derivatives of both the x and z points
%   This will yield a system similar to Greg's symmetric collocation
% Because we are using radial kernels, we know that
%   K(x,z) = phi(0) and thus d/dx K(x,x) = 0
% For nonradial kernels this may require more care
% In higher dimensions this requires MUCH more care
rbfx = @(r, dx) process_variance * -2 * ep ^ 2 * dx .* exp(-(ep * r) .^ 2);
rbfxz = @(r, dx) process_variance * 2 * ep ^ 2 * (1 - 2 * dx .^ 2 * ep ^ 2) .* exp(-(ep * r) .^ 2);

% We consider the standard gradient from numerical analysis
% We also want to consider the variance for that prediction
yxplt = yfx(xplt);

Kevalx = rbfx(DistanceMatrix(xplt, x), DifferenceMatrix(xplt(:, 1), x(:, 1)));
yxeval = Kevalx * (Kint \ y);

ksdx = sqrt(max(rbfxz(0, 0) - sum((Kevalx / Kint) .* Kevalx, 2), 0));

subplot(2, 2, 3)
h31 = plot(xplt, yxplt, 'r', 'linewidth', 3);
hold on
h32 = plot(xplt, yxeval, 'b', 'linewidth', 3);

h4 = fill_between_lines(xplt, yxeval - 2 * ksdx, yxeval + 2 * ksdx, 'b');
set(h4, 'FaceAlpha', .1);
set(h4, 'EdgeAlpha', .1);
h4top = plot(xplt, yxeval + 2 * ksdx, ':b', 'linewidth', 3);
h4top.Color = [0,0,1,0.2];
h4bot = plot(xplt, yxeval - 2 * ksdx, ':b', 'linewidth', 3);
h4bot.Color = [0,0,1,0.2];
legend([h1, h2, h3], {'true', 'data', 'pred'}, 'location', 'south')
title('Deriv, 2 SD confidence intervals')
hold off
ylim([-4, 4])

% Same idea as before where we try to take draws from the posterior
% This time we draw from the posterior of the gradient
Kpltxx = rbfxz(DistanceMatrix(xplt, xplt), DifferenceMatrix(xplt(:, 1), xplt(:, 1)));
Kpostcovx = Kpltxx - (Kevalx / Kint) * Kevalx';
Kpostcovx = Kpostcovx + wiggle_room * eye(size(xplt, 1));
Kpostcovx_chol = chol(Kpostcovx, 'lower');

yxpost = bsxfun(@plus, yxeval, Kpostcovx_chol * normal_draws);

subplot(2, 2, 4)
h41 = plot(xplt, yxpost, 'b');
hold on
for k=1:length(h41)
    h41(k).Color = [0, 0, 1, 0.1];
end
title(sprintf('Deriv, %d posterior draws', n_draws))
hold off
ylim([-4, 4])




