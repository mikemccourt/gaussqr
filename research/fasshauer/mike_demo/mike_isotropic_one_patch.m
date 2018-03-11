% The goal of this is to study just an isotropic PUM to understand how well
% the LOOCV or MLE criteria works in that situation
%
% I'm using the same code that we have used for the anisotropic setting
% because I have already been working on it, but it should be 100%
% computing isotropically.

% The test function
yf = @(x) franke(x(:, 1), x(:, 2));

M = 2;                                  % space dimension 
n = 25;                                 % number of tracks
m = 80;                                 % number of points on tracks
[N, dsites, yy] = TrackData2D(n, m);       % generate N = n*m track data in 2D
rhs = yf(dsites);                        % function values

rbf = @(r) 1 ./ sqrt(1 + r .^ 2);       % define the RBF
neval = 30;                             % parameter for evaluation points

npu = floor((N / 4) ^ (1 / M));          % parameter for PU centres
xx = linspace(0, 1, n);                   % subdomains in one direction
[X, Y] = meshgrid(xx, yy);                % patches centered at tracks
puctrs = [X(:) Y(:)];                   % define the PU centres

strat = 'loocv';                        % how to find an optimal choice
tikhonov_parameter = 0;                 % a regularization parameter

% Here I choose a PU center for which to consider the parametrization surface
point_to_test = puctrs(23, :);
index = IntegerBased_MD_ContainingQuery(point_to_test, npu, puradius(1, 1), M);

% This initializes some of the indexing tools
puradius = ones(1, M) * (1 ./ npu);
idx_ds = IntegerBased_MD_Structure(dsites, npu, puradius(1, 1), M);

% This is the choices for the graph that we are going to create
radii_vec = linspace(.07, .35, 25);
ep_vec = logspace(0, 2, 26);
[R, E] = meshgrid(radii_vec, ep_vec);
L = zeros(size(E));

for row=1:size(R, 1)
    for col=1:size(E, 2)
        L(row, col) = Cost_function_mixed( ...
            rbf, [R(row, col), E(row, col)], puradius, idx_ds, index, ...
            npu, M, puctrs, dsites, rhs, strat, tikhonov_parameter);
    end
end

L(L == 1e60) = nan;
L(L == 1e30) = nan;
if strcmp(strat, 'loocv')
    contourf(R, log10(E), log10(L), 20);
    color_label = 'log10 loocv residual';
else
    contourf(R, log10(E), L, 20);
    color_label = 'MLE';
end
c = colorbar;
xlabel('radius size')
ylabel('log10 shape parameter')
c.Label.String = color_label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Notes from Mike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Admittedly, I do not understand everything that is going on in the
% integer-based search, but the point of this demo is to look at the
% surface of the parametrization strategy.  I want to understand what it is
% that we are optimizing when we iterate over each patch.
%
% When I look at the surface for a given point_to_test (as chosen above) it
% seems like the LOOCV surface is pointing rather directly away from what
% someone might consider the "good" region (remember that the goal is to
% minimize the LOOCV residual).  It seems that the LOOCV resuidual can be
% minimized by letting the shape parameter become incredibly large.  This
% basically has the effect of ignoring the points that one is trying to
% fit, I think.
%
% I am mentioning this because I think this matches the structure of what
% the original research on this topic consisted of.  I may not have fully
% understood it, but in looking at these experiments it seems, at first
% glance, as though there is little correlation between the LOOCV value and
% the quality of the resulting PUM-interpolation.  After all, if the shape
% parameter becomes TOO small there's no way the resulting interpolant will
% have no predictive capacity.
%
% The reason I ran this example is because I was trying to understand what
% I was seeing in the anisotropic case.  I was seeing behavior somewhat
% similar to this when running tests on anisotropic problems and I could
% not understand exactly what was happening.  I thought that maybe if I
% looked at this simpler version of the problem (with only 2 free
% parameters instead of 4) I would be able to track the behavior better.