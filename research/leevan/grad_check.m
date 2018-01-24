% This is a quick demo as to what it looks like to try and identify an
% approximation and then use that to estimate the gradient of the norm.
% This is adapted from PDEFiniteDiff1D.m in Fasshauer/McCourt
% It requires the software with that book (could be replaced):
%   - pick2Dpoints
%   - DistanceMatrix
%   - DifferenceMatrix

% Produce repeatable results
rng(1)

% Create and plot the data that will be present
% Need to include some noise for the article to make sense
noise_level = .01;
subplot(1, 2, 1);
N = 300;
x = pick2Dpoints(0, 1, sqrt(N), 'halton');
xplt = linspace(0, 1, 35);
yplt = linspace(0, 1, 36);

y = zeros(length(x), 1);
y(sqrt(sum((x - .4).^2, 2)) < .3) = 1;
y(sqrt(sum((x - .6).^2, 2)) < .2) = 2;
y = y + noise_level * randn(length(y), 1);

[xx1, xx2] = meshgrid(xplt, yplt);
zz = griddata(x(:,1), x(:, 2), y, xx1, xx2);
x0 = x(y < 0.5, :);
x1 = x((y > 0.5) & (y < 1.5), :);
x2 = x(y > 1.5, :);
contour(xx1, xx2, zz, 2)
hold on
h0 = plot(x0(:, 1), x0(:, 2), 'or');
h1 = plot(x1(:, 1), x1(:, 2), 'xb');
h2 = plot(x2(:, 1), x2(:, 2), 'g+');
title('Data')
hold off
legend([h0, h1, h2], {'value = 0', 'value = 1', 'value = 2'})

% Choose a kernel to study interpolant gradients
% The Tikhonov level is the quantity being added to the diagonal
tikhonov_level = .001;
subplot(1, 2, 2);
K = @(e, r) (1 + (e * r) + (e * r).^2 / 3) .* exp(-e * r);
Kd = @(e, r, dx) -e^2 / 3 * dx .* (1 + e * r) .* exp(-e * r);
Neval = 25;
xeval = pick2Dpoints(0, 1, Neval, 'even');
ep = 15;
Kint = K(ep, DistanceMatrix(x, x)) + tikhonov_level * eye(length(x));
Keval = K(ep, DistanceMatrix(xeval, x));
Kevalx = Kd(ep, DistanceMatrix(xeval, x), DifferenceMatrix(xeval(:, 1), x(:, 1)));
Kevaly = Kd(ep, DistanceMatrix(xeval, x), DifferenceMatrix(xeval(:, 2), x(:, 2)));
c = Kint \ y;
ypred = Keval * c;
ypredx = Kevalx * c;
ypredy = Kevaly * c;
ypredgrad = [ypredx, ypredy];
ypredgradnorm = sqrt(sum(ypredgrad.^2, 2));

X1 = reshape(xeval(:, 1), Neval, Neval);
X2 = reshape(xeval(:, 2), Neval, Neval);
YG = reshape(ypredgradnorm, Neval, Neval);
contour(X1, X2, YG, [10, 20, 30]);
title('Gradient Norm')

% Now we can consider taking draws from the posterior 
% We consider the joint distribution of [Y_x1, Y_x2, Y]
%   where the derivative terms are at xeval and Y is at x.
% The covariance is then known to be
%   [[C_11_ee, C_12_ee', C_1_ex'],
%    [C_12_ee, C_22_ee,  C_2_ex'],
%    [C_1_xe,  C_2_xe,   C_xx]   ]
% e.g., C_11 is the derivative wrt the first dimension of both the
%   evaluation point and the center, both evaluated and centered at xeval
% Also, C_xx is the normal covariance evaluated and centered at x
% Note -- I am doing this with the Gaussian because I have already found
%   these gradients.  It could be done with something else as well.
tikhonov_level = .001;
ep = 6;
a = 1.0; % Process variance ... just in case we want to change it later
diagonal_safety = a * 1e-11;
K = @(e, r) a * exp(-(e * r) .^ 2);
Kd = @(e, r, dx) a * -2 * e ^ 2 * dx .* exp(-(e * r) .^ 2);
Kd11 = @(e, r, dx) a * 2 * ep ^ 2 * (1 - 2 * dx .^ 2 * ep ^ 2) .* exp(-(ep * r) .^ 2);
Kd12 = @(e, r, d1, d2) a * -4 * ep ^ 4 * d1 .* d2 .* exp(-(ep * r) .^ 2);

Kint = K(ep, DistanceMatrix(x, x)) + tikhonov_level * eye(length(x));
Kevalx1 = Kd(ep, DistanceMatrix(xeval, x), DifferenceMatrix(xeval(:, 1), x(:, 1)));
Kevalx2 = Kd(ep, DistanceMatrix(xeval, x), DifferenceMatrix(xeval(:, 2), x(:, 2)));
Kevalx1x1 = Kd11(ep, DistanceMatrix(xeval, xeval), DifferenceMatrix(xeval(:, 1), xeval(:, 1)));
Kevalx2x2 = Kd11(ep, DistanceMatrix(xeval, xeval), DifferenceMatrix(xeval(:, 2), xeval(:, 2)));
Kevalx1x2 = Kd12(ep, DistanceMatrix(xeval, xeval), DifferenceMatrix(xeval(:, 1), xeval(:, 1)), DifferenceMatrix(xeval(:, 2), xeval(:, 2)));

Kpriorcovgrad = [Kevalx1x1, Kevalx1x2'; Kevalx1x2, Kevalx2x2];
Kposteriorcovgrad = Kpriorcovgrad - ([Kevalx1; Kevalx2] / Kint) * [Kevalx1; Kevalx2]' + diagonal_safety * eye(size(Kpriorcovgrad));
Kposteriorcovgrad_chol = chol(Kposteriorcovgrad, 'lower');

ypredx1 = Kevalx1 * (Kint \ y);
ypredx2 = Kevalx2 * (Kint \ y);

n_draws = 100;
normal_draws = randn(length(Kposteriorcovgrad_chol), n_draws);
ygradvar_draws = Kposteriorcovgrad_chol * normal_draws;
ygradx1posterior_draws = bsxfun(@plus, ypredx1, ygradvar_draws(1:length(ypredx), :));
ygradx2posterior_draws = bsxfun(@plus, ypredx2, ygradvar_draws(length(ypredy)+1:end, :));
ygradposteriornorm_draws = sqrt(ygradx1posterior_draws .^ 2 + ygradx2posterior_draws .^ 2);

figure
axis_count = 1;
for ycheck=[5, 10, 15, 20]
    ycount = sum(ygradposteriornorm_draws > ycheck, 2) / n_draws;
    YC = reshape(ycount, Neval, Neval);
    subplot(2, 2, axis_count)
    contour(X1, X2, YC, 10)
    title(sprintf('Prob norm > %g', ycheck))
    caxis([0, 1])
    colorbar
    axis_count = axis_count + 1;
end

% Uncomment this to see each of the individual outcomes
% figure
% fprintf('Press enter to step through the random draws, ^c to exit\n')
% tt = linspace(0, 2 * pi, 200);
% xt1 = .4 + .3 * cos(tt);yt1 = .4 + .3 * sin(tt);
% xt2 = .6 + .2 * cos(tt);yt2 = .6 + .2 * sin(tt);
% for k=1:n_draws;
%     YGX1D = reshape(ygradx1posterior_draws(:, k), Neval, Neval);
%     YGX2D = reshape(ygradx2posterior_draws(:, k), Neval, Neval);
%     YGD = reshape(ygradposteriornorm_draws(:, k), Neval, Neval);
%     subplot(1, 3, 1);
%     contour(X1, X2, YGX1D, 10);
%     hold on
%     plot(xt1, yt1, 'k', 'linewidth', 2)
%     plot(xt2, yt2, 'k', 'linewidth', 2)
%     hold off
%     subplot(1, 3, 2);
%     contour(X1, X2, YGX2D, 10);
%     hold on
%     plot(xt1, yt1, 'k', 'linewidth', 2)
%     plot(xt2, yt2, 'k', 'linewidth', 2)
%     hold off
%     subplot(1, 3, 3);
%     contour(X1, X2, YGD, [10, 20, 30]);
%     hold on
%     plot(xt1, yt1, 'k', 'linewidth', 2)
%     plot(xt2, yt2, 'k', 'linewidth', 2)
%     hold off
%     title(sprintf('Image %d/%d', k, n_draws))
%     colorbar
%     pause
% end