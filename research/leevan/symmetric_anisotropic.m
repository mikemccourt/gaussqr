% This is some experimentation involving the kernels you are working on
% In this setting the goal is to work with nonradial kernels
% By keeping them as symmetric, though, there is still a chance to
% leverage the standard theory
% This is section 4 associated with 14-Dec-2017

% All kernels of the variety (4.1) start from an intermediate function
% which we will call J(x,z).  This can be any symmetric function of two
% variables but there are probably good choices available.
% Using this function, we can define a new kernel
%   K(x,y) = sum_i=1^N J(x,z_i)J(y,z_i)
% Let's look at some examples of what this looks like

plot_index = 1;

z = linspace(.2, .9, 10);
y = linspace(0.05, .95, 6);
xplt = linspace(0, 1, 222);

% J(x,z) = min(x,z) - xz
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
Jx = bsxfun(@min, xplt', z) - xplt' * z;
Jy = bsxfun(@min, y', z) - y' * z;
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = min(x, z) - xz')
leg_names = arrayfun(@(x)sprintf('y=%g',x), y, 'UniformOutput', 0);
legend(leg_names, 'Location', 'northwest')

%{
% Special kernels only defined on (0, 1)
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
ep = @(z) -2 ./ (.1 + abs(z - .88) .^ 2);
Jx = bsxfun(@rdivide,...
     sinh(bsxfun(@times,ep(z),bsxfun(@min, xplt', z))) .* ...
     sinh(bsxfun(@times,ep(z),1-bsxfun(@max, xplt', z))), ...
     ep(z) .* sinh(ep(z)));
Jy = bsxfun(@rdivide,...
     sinh(bsxfun(@times,ep(z),bsxfun(@min, y', z))) .* ...
     sinh(bsxfun(@times,ep(z),1-bsxfun(@max, y', z))), ...
     ep(z) .* sinh(ep(z)));
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = min(x, z)')
leg_names = arrayfun(@(x)sprintf('y=%g',x), y, 'UniformOutput', 0);
legend(leg_names, 'Location', 'northwest')
%}

% J(x,z) = exp(-3 * |x - z|)
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
Jx = exp(-3 * abs(bsxfun(@minus, xplt', z)));
Jy = exp(-3 * abs(bsxfun(@minus, y', z)));
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = exp(-3 * |x - z|)')
% leg_names = arrayfun(@(x)sprintf('y=%g',x), y, 'UniformOutput', 0);
% legend(leg_names, 'Location', 'south')

%{
% J(x,z) = exp(-10 * |x - z|^2)
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
Jx = exp(-10 * bsxfun(@minus, xplt', z) .^ 2);
Jy = exp(-10 * bsxfun(@minus, y', z) .^ 2);
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = exp(-10 * |x - z|^2)')
% leg_names = arrayfun(@(x)sprintf('y=%g',x), y, 'UniformOutput', 0);
% legend(leg_names, 'Location', 'northwest')
%}

% J(x,z) = exp(-5 / (.1 + |z - .88|^2) * |x - z|^2)
z = linspace(.1, .9, 10);
y = linspace(0.05, .95, 6);
xplt = linspace(0, 1, 222);
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
Jx = exp(bsxfun(@times, -5 ./ (.1 + (z - .88) .^ 2), bsxfun(@minus, xplt', z) .^ 2));
Jy = exp(bsxfun(@times, -5 ./ (.1 + (z - .88) .^ 2), bsxfun(@minus, y', z) .^ 2));
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = exp(-1 / (1 + |z - .85|) * |x - z|^2)')

% Special kernels only defined on [-1, 1] from Fasshauer/McCourt
subplot(2, 2, plot_index);
plot_index = plot_index + 1;
z = linspace(-.8, .9, 10)';
y = linspace(-1, 1, 6)';
xplt = linspace(-1, 1, 222)';
ep = @(x, z) repmat(.95./(1 + abs(z-.7)'),length(x), 1);
J = @(b,x,z) .4 + (1-b).* ...
    (b.*(1-b.^2) - 2*b.*bsxfun(@plus,x.^2,z.^2') + (1+3*b.^2).*(x*z'))./ ...
    ((1-b.^2).^2 + 4*b.*(b.*bsxfun(@plus,x.^2,z.^2')-(1+b.^2).*(x*z')));
Jx = J(ep(xplt,z),xplt,z);
Jy = J(ep(y,z),y,z);
Kxy = Jx * Jy' / length(z);
plot(xplt, Kxy, 'linewidth', 3)
title('J(x,z) = Cinf Chebyshev kernel')
% leg_names = arrayfun(@(x)sprintf('y=%g',x), y, 'UniformOutput', 0);
% legend(leg_names, 'Location', 'south')

%%%
% This is an interpolation example
figure;
f = @(x) sin(4 * x)';
y = linspace(0.05, .95, 6);
xplt = linspace(0, 1, 222);
fint = f(y);
feval = f(xplt);

hold on
plot(xplt, feval, '--', 'linewidth', 2);
plot(y, fint, 'or');
leg_names = {'true', 'data'};
for z_num=[10, 30, 90];
    z = linspace(.1, .9, z_num);
    J = @(x, z) exp(bsxfun(@times, -5 ./ (.1 + (z - .88) .^ 2), bsxfun(@minus, x', z) .^ 2));
    K = @(x, y, z) J(x, z) * J(y, z)' / length(z);
    Kint = K(y, y, z);
    Keval = K(xplt, y, z);
    fpred = Keval * (Kint \ fint);
    plot(xplt, fpred, 'linewidth', 2);
    leg_names{length(leg_names) + 1} = sprintf('z_{num} = %d', z_num);
end
legend(leg_names)
hold off