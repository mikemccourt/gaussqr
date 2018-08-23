% This is a workplace for the Laguerre kernels (name pending)
x = linspace(0, 10, 200);

% This is a plot showing what the impact of alpha is
avec = [0, 1, 2, 3];
ac = 1;
for a=avec
    subplot(2, 2, ac)
    cla reset
    hold on
    for k=0:5
        plot(x, LaguerrePoly(k, x, a), 'linewidth', 3);
    end
    title(sprintf('alpha = %g', a));
    ylim([-40, 60])
    hold off
    ac = ac + 1;
end

% This should check for normality of the weight function
a = .78;
d = .17;
n = 3;
rho = @(x,a,d) x.^a .* exp((2*d - 1)*x) * (1 - 2 * d)^(a+1) / gamma(a + 1);
check = quadgk(@(x) rho(x,a,d), 0, Inf);
if abs(check - 1) > 1e-10
    fprintf('normality of weight function failed')
end

% This defines the normalization term for the eigenfunctions
g = @(n,a,d) sqrt(gamma(n+1)/gamma(n+a+1)*gamma(a+1)/(1-2*d)^(a+1));

% And the eigenfunctions themselves
phi = @(n,x,a,d) g(n,a,d)*exp(-d*x).*LaguerrePoly(n,x,a);

% With their orthonormality
check = quadgk(@(x)phi(n,x,a,d).*phi(n+1,x,a,d).*rho(x,a,d),0,Inf);
if abs(check) > 1e-10
    fprintf('orthogonality of eigenfunctions failed')
end
check = quadgk(@(x)phi(n,x,a,d).*phi(n,x,a,d).*rho(x,a,d),0,Inf);
if abs(check - 1) > 1e-10
    fprintf('normality of eigenfunctions failed')
end

% Now I also want to be able to consider the impact of d
x = linspace(0, 100, 200);
a = .1;
figure
dvec = [.1, .15, .2, .24];
dc = 1;
for d=dvec
    subplot(2, 2, dc)
    cla reset
    hold on
    for k=0:5
        plot(x, phi(k, x, a, d), 'linewidth', 3);
    end
    title(sprintf('alpha=%g, delta = %g', a, d));
    hold off
    dc = dc + 1;
end

% Now for some kernels ... first we check that the closed form works
x = linspace(0, 16, 200);
a = .3;
d = .22;
t = .48;
nvec = 0:20;
figure
clf reset
hold on
zvec = [0, 1, 2, 3, 4] * 2;
h = zeros(size(zvec));
hl = cell(size(zvec));
zc = 1;
for z=zvec
    Kvals = zeros(size(x));
    Kc = 1;
    for xx=x
        phix = arrayfun(@(n)phi(n,xx,a,d), nvec);
        phiz = arrayfun(@(n)phi(n,z,a,d), nvec);
        lamvec = (1 - t) * t .^ nvec;
        Kvals(Kc) = phix * (lamvec .* phiz)';
        Kc = Kc + 1;
    end
    K_closed = laguerre_kernel(x, z, a, d, t);
    h(zc) = plot(x, Kvals, 'linewidth', 3);
    hl{zc} = sprintf('z=%g', z);
    hc = plot(x, K_closed, '--k', 'linewidth', 3);
    hc.Color(4) = .3;
    zc = zc + 1;
end
hold off
title(sprintf('a=%g, d=%g, t=%g', a, d, t))
legend(h, hl, 'location', 'northeast')

% So, I think that we can conclude that the closed form is okay
% Maybe some issues in the context of really large values??
% Will need to think about that ...
% Let's try a basic interpolation problem and see what we get
% Also, maybe should use a better name than t for the free parameter
rbf = @(x, z, e) exp(-(e * DistanceMatrix(x, z)).^2);

% Just randomly creating points now -- likely to be structured in time
N = 30;
yf = @(x, t) cos(5 * x) .* exp(-t / 3);
xx = rand(N, 1);
xt = 5 * rand(N, 1);
y = yf(xx, xt);

% Maybe try a little log-likelihood while we're trying stuff ??
opts.logspace = false;opts.num_points = 1000;
bounds = [[-.9, .9]; [.01, .24]; [.01, .99]; [.8, 5]];
kernel_eval = @(a, d, t, e) rbf(xx, xx, e) .* laguerre_kernel(xt, xt, a, d, t);
opt_adte = fminrnd(@(v) kernel_mle(kernel_eval(v(1), v(2), v(3), v(4)), y), bounds, opts);
a = opt_adte(1);d = opt_adte(2);t = opt_adte(3);e = opt_adte(4);

x1d = linspace(0, 1, 60);
t1d = linspace(0, 5, 61);
[XX, TT] = meshgrid(x1d, t1d);
YY = yf(XX, TT);
xxeval = XX(:);
xteval = TT(:);

Kint_x = rbf(xx, xx, e);
Kint_t = laguerre_kernel(xt, xt, a, d, t);
K = Kint_x .* Kint_t;
Keval_x = rbf(xxeval, xx, e);
Keval_t = laguerre_kernel(xteval, xt, a, d, t);
Keval = Keval_x .* Keval_t;
ypred = Keval * (K \ y);
YP = reshape(ypred, length(t1d), length(x1d));

figure
subplot(1, 2, 1)
contour(XX, TT, YY, 'linewidth', 2)
c = colorbar;
hold on
plot(xx, xt, 'ok', 'markersize', 6, 'MarkerFaceColor', 'k')
hold off
title('Original data')

subplot(1, 2, 2)
contour(XX, TT, log10(abs(YP - YY)), 'linewidth', 2)
c = colorbar;
c.Label.String = 'log_{10}-error';
hold on
plot(xx, xt, 'ok', 'markersize', 6, 'MarkerFaceColor', 'k')
hold off
title(sprintf('a=%g, d=%g, t=%g; e=%g', a, d, t, e))

% That error plot is kinda cool, right?  Like the data isn't talking in a
% purely isotropic fashion.  The error "ridges", so to speak, have a
% structure in them which is maybe the result of the kernel choice?

