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

rbf_aniso = @(r) 1 ./ sqrt(1 + r .^ 2);       % define the RBF
neval = 30;                             % parameter for evaluation points

npu = floor((N / 4) ^ (1 / M));          % parameter for PU centres
xx = linspace(0, 1, n);                   % subdomains in one direction
[X, Y] = meshgrid(xx, yy);                % patches centered at tracks
puctrs = [X(:) Y(:)];                   % define the PU centres


% This is the choices for the graph that we are going to create
radii_vec = linspace(.07, .3, 25);
ep_vec = logspace(0, 3, 26);
[R, E] = meshgrid(radii_vec, ep_vec);
L = zeros(size(E));

puradius = ones(1, M) * (1 ./ npu);
idx_ds = IntegerBased_MD_Structure(dsites, npu, puradius(1, 1), M);
point_to_test = puctrs(106, :);
index = IntegerBased_MD_ContainingQuery(point_to_test, npu, puradius(1, 1), M);

for row=1:size(R, 1)
    for col=1:size(E, 2)
        L(row, col) = Cost_function_mike(rbf, [R(row, col), E(row, col)], puradius, idx_ds, index, npu, M, puctrs, dsites, rhs);
    end
end

L(L == 1e60) = nan;
L(L == 1e30) = nan;
contourf(R, log10(E), log10(L), 20);
colorbar





h = 2;                                  % upper bound for the radius
param = [.3, .3, 3, 3];             % initial values for the parameters 
                                       % to be optimized


% Alternate patches ?????
[X, Y] = meshgrid(linspace(0, 1, 7), linspace(.02, .98, 6));

% What about if we just do interpolation with a C2 Wendland kernel
rbf = @(r) max(1-r,0).^4.*(4*r+1);
% This is a C6 Wendland kernel below, if we prefer
% k = 3;
% d = 2;
% l = floor(d / 2 + k + 1);
% rbf = @(r) (1 - r) .^ (l + 3) .* (1+(l+3)*r + (6*l^2+36*l+45)/15*r.^2 + (l^3+9*l^2+23*l+15)/15*r.^3);
% We would have to choose some kernel parameters
epvec = [10, 5];
tic
K = DistanceMatrix(dsites, dsites, epvec, rbf);
Keval = DistanceMatrix(epoints, dsites, epvec, rbf);
yeval = Keval * (K \ rhs);

exact = yf(epoints);
maxerr = norm(yeval - exact,inf);
rms_err = norm(yeval - exact)/sqrt(length(exact));
toc
fprintf('RMS error:       %e\n', rms_err);
fprintf('Maximum error:   %e\n', maxerr);
