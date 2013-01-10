clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MEG meshfree forward solver for a single sphere model
%                            - Kansa's method -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calls on:
%   BallGeometry.m
%   SphereSurfGoldPoints.m
%   DistanceMatrix.m
%   DifferenceMatrix.m
%   gradphiF_dip.m
%   phiF_dip.m
%   HomSpherePotential.m
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: singularities of analytic potential formula for dipole located
% at origin
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
radbasfun = 'gaussian';     % Radial basis function
Npnts = 700;                % Number of desired interior points
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
[POINTS, NORMALS] = BallGeometry(R, Npnts, 'kansa');
intdata = POINTS.int1;
bdydata = POINTS.bdy11;

% Boundary centers (inside the domain)
bdyctrs = bdydata;

% Centers
ctrs = [intdata; bdydata];

% Evaluation points
evalpnts = SphereSurfGoldPoints(1000, R);
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

A = bsxfun(@times,NORMALS.n11(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
B = bsxfun(@times,NORMALS.n11(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
C = bsxfun(@times,NORMALS.n11(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));

BCM = bsxfun(@plus,A,B);
BCM = bsxfun(@plus,BCM,C);

% Collocation matrix
CM = [LCM; BCM];

% Compute known-terms vector (a.k.a. righthand side vector)
gradphi_F = gradphiF_dip(bdydata, srcpnts, dipmom, sig);% Gradient of the 
                                                    % potential at boundary
                                                    % in the unbound case
rhs = [ zeros(size(intdata,1),1); -sum(NORMALS.n11.*gradphi_F,2) ];


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
% Delaunay triangulation of the sphere
d_evalpnts = delaunayn(evalpnts);
tr = TriRep(d_evalpnts, evalpnts);
tr_surf = freeBoundary(tr);

subplot(1,2,1)
SurfacePlot_dip(evalpnts, tr_surf, phi_an, 'Analytic potential [V]')
subplot(1,2,2)
SurfacePlot_dip(evalpnts, tr_surf, phi, 'Computed potential [V]')