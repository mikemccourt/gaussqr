clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Potential interpolation problem for error bound
%
% Last modified: 2012/09/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data
%--------------------------------------------------------------------------

% Medium data
R = 0.1;                      % Sphere radius [m]
sig = 0.2;                    % Electric conductivity [S/m]
% Sources data
dipmom = 2.7e-12.*[1, 0, 0];  % Dipole moment [Am]
srcpnts = [0, 0, 0.6*R];      % Dipole position [m]

% Parameters for numerical computation
radbasfun = 'imq';            % Radial basis function
Npnts_surf = 120;             % Number of desired points on sphere's surface
ep = 5;                       % RBFs shape parameter


% RBF definition and derivatives
%--------------------------------------------------------------------------
% rbf     = Radial basis function;
rbf = pickRBF(radbasfun);


% Collocation matrix and known-terms vector
%--------------------------------------------------------------------------

% Collocation points
[bdydata, intdata, dist] = SphereRegGoldPoints(Npnts_surf, R);
data = [bdydata; intdata];
% Centers
ctrs = data;

% Evaluation points
evalpnts = SphereRegUnifPoints(dist/3, R); % dist/3 controls the distance 
                                           % between evaluation points
neval = size(evalpnts,1);

% Compute exact solution at data sites (function values -> rhs)
rhs = HomSpherePotential(R, sig, srcpnts, dipmom, data);

% Compute exact solution at evaluation points 
exact = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

% Compute distance matrices
DM_eval = DistanceMatrix(evalpnts, ctrs);
DM_data = DistanceMatrix(data, ctrs);

% Compute evaluation matrix EM
EM = rbf(ep, DM_eval);

% Compute interpolation matrix IM
IM = rbf(ep, DM_data);

% Compute interpolant
c = IM\rhs;
Pf = EM * c;%(IM\rhs);

% Comparison and maximum errors
%--------------------------------------------------------------------------
max_err_phi = norm(Pf-exact,inf);
normdiff = norm(Pf-exact);
rel_err_phi = normdiff/norm(exact);
rms_err_phi = norm(Pf-exact)/neval;
COND = cond(IM);