clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MEG meshfree forward solver for a single sphere model
%                            - Kansa's method -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calls on:
%   SphereRegGoldPoints.m
%   SphereSurfGoldPoints.m
%   SphereRegUnifPoints.m
%   DistanceMatrix.m
%   DifferenceMatrix.m
%   InfMediumPotential.m
%   HomSpherePotential.m
%   HomSphereInduction.m
%
% To do:
%   -- Computation of J vector as -sig*grad(phi).
%   -- Computation of B (magnitude?normal component?) with Ampere-Laplace
%      formula.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: singularities of analytic potential formula for dipole located
% at origin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Magnetic induction field observation points
% obspnts = ; ...
% Parameters for numerical computation
radbasfun = 'gaussian';            % Radial basis function
Npnts_surf = 300;             % Number of desired points on sphere's surface
ep = 50;                    % RBFs shape parameter


% RBF definition and derivatives
%--------------------------------------------------------------------------
% rbf     = Radial basis function;
% dxrbf   = component along x of the gradient of the RBF
% dyrbf   = component along y of the gradient of the RBF
% dzrbf   = component along z of the gradient of the RBF
% Lrbf    = Laplacian of the RBF in 3D
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);


% Collocation matrix and known-terms vector
%--------------------------------------------------------------------------

% Collocation points
[bdydata, intdata, dist] = SphereRegGoldPoints(Npnts_surf, R);
% Boundary centers outside the domain
bdyctrs = SphereSurfGoldPoints(Npnts_surf, R+0*dist);
% Centers
ctrs = [intdata; bdyctrs];
% Evaluation points
%evalpnts = SphereRegUnifPoints(dist/3, R); % Dist/3 controls the distance 
                                           % between evaluation points
[evalpnts, tri_evalpnts] = SphereSurfUniformPoints(R,3);
neval = size(evalpnts,1);

% Compute evaluation matrix EM
DM_eval = DistanceMatrix(evalpnts, ctrs);
EM = rbf(ep, DM_eval);

% Compute blocks for collocation matrix
% Interior points
DM_intdata = DistanceMatrix(intdata,ctrs);
LCM = Lrbf(ep,DM_intdata);
% Boundary points
DM_bdydata = DistanceMatrix(bdydata,ctrs);
dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));

NV = bdydata/R;   % Unit vectors normal to sphere surface

A = bsxfun(@times,NV(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
B = bsxfun(@times,NV(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
C = bsxfun(@times,NV(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));

BCM = bsxfun(@plus,A,B);
BCM = bsxfun(@plus,BCM,C);

% Collocation matrix
CM = [LCM; BCM];

% Compute known-terms vector (a.k.a. righthand side vector)
gradphi_F = gradphiF_dip(bdydata, srcpnts, dipmom, sig);% Gradient of the 
                                                    % potential at boundary
                                                    % in the unbound case
rhs = [ zeros(size(intdata,1),1); -sum(NV.*gradphi_F,2) ];


% Numerical solution for the potential
%--------------------------------------------------------------------------

% Potential at evalpnts in the source free case
phi0 = EM * (CM\rhs);
% Potential at evalpnts in the unbound domain case
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);
% Potential at evalpnts (superposition of effects)
phi = phi0 + phi_F;


% Numerical solution for the magnetic induction field
%--------------------------------------------------------------------------


%  Analytic solution for the potential
%--------------------------------------------------------------------------
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);


% Analytic solution for the magnetic induction field
%--------------------------------------------------------------------------
% B_exact = HomSphereInduction(srcpnts, dipmom, obspnts);


% Comparison and maximum errors
%--------------------------------------------------------------------------
% Potential
max_err_phi = norm(phi-phi_an,inf);
normdiff_phi = norm(phi-phi_an);
rel_err_phi = normdiff_phi/norm(phi_an)
rms_err_phi = norm(phi-phi_an)/neval;
COND = cond(CM)
% Magnetic induction field
% max_err_B = norm(B-B_an,inf);
% normdiff_B = norm(B-B_an);
% rel_err_B = normdiff_B/norm(B_an);
% rms_err_B = norm(B-B_an)/size(obspnts,1);

% Plots
%--------------------------------------------------------------------------

