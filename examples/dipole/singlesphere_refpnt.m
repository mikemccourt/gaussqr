% singlesphere_refpnt.m
%
%  In this experiment we are interested in studying the effect of the 
%  choice of the reference point on the error on the solution (potential 
%  difference).
%  The script generates a plot showing how the error on the solution 
%  behaves when we move the reference point over the sphere.
%
%  We consider the solution to the Laplace equation on a
%  sphere with Neumann boundary conditions.
%
%  The problem has several physical parameters relating to the
%  underlying EEG/MEG physical system.  These parameters are:
%    R - Sphere radius [dm] <default = 1>
%    sig - Electric conductivity [S/dm] <default = .02>
%    dipmom - Dipole moment [x10^-12 Am] <default = 2.7e*[1,0,0]>
%    srcpnts - Dipole position [dm] <default = [0,0,0.6*R]>
%
%  This script allows you to test the error with respect to the reference
%  point for the potential.
%
%  The solution parameters to be considered are
%     sol_type - How you want to solve the system <default = 'kansa'>
%                'kansa' : Nonsymmetric collocation
%                          rbf_choice and ep must also be chosen
%                'mfs' : Method of fundamental solutions
%                        MFS_frac and ctr_sphere must also be chosen
%     rbf_choice - RBF for collocation <default = 'imq'>
%     ep - RBF shape parameter <default = 10>
%     mfs_frac - How many centers for MFS, in [0.0,1.0]*N <default = 1.0>
%     mfs_sphere - Fraction beyond R (eg, 1.3R) for centers <default = 1.3>
%
%  Other parameters to be set:
%     Npnts - Number of desired collocation points <default = 1000>
%     dip_cushion - How much space should be given around the dipole where
%               no RBF centers are allowed <default = .005>
%     N_eval - # evaluation points <default = 1001>
%
%
%  The results of these experiments are stored in
%     errvec - Errors computed considering each evaluation point over the 
%              sphere as a reference point
%
%  May need to consider the results of this experiment on the choice of
%  reference point for future experiment.  It seems to have quite an
%  effect, as should be expected.


R = 1;
sig = 0.02;
dipmom = 2.7*[1, 0, 0];
srcpnts = [0, 0, 0.6*R];

sol_type = 'kansa';
radbasfun = 'imq';
ep = 10;
mfs_frac = 1.0;
mfs_sphere = 1.3;
eval_diff = 1;

Npnts = 1000;
dip_cushion = .005; % Not yet implemented
N_eval = 1001;


%%%%%%%%%%%%%%%%%%%%%
% Basic setup stuff independent of this problem

% Consider the standard GQR parameters for the errcompute function
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Set random number generator to constant
% This is used in choosing which BC points are Dirichlet in the mixed case
rng(0);


%%%%%%%%%%%%%%%%%%%%%
% This is the start of the solver

% RBF definition and derivatives
if strcmp(sol_type,'kansa')
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);
else
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');
end

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval, R);

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);

% Analytic solution for the potential
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

% Loop through the requested N values

% Determine collocation points
[POINTS, NORMALS] = BallGeometry(R,Npnts,sol_type,[],srcpnts);
intdata = POINTS.int1;
bdydata = POINTS.bdy11;
N_int = size(intdata,1);
N_bdy = size(bdydata,1);
N_tot = N_int + N_bdy;

% Compose a vector of all the RBF centers
% In the MFS setting, these are chosen in a sphere around the ball
if strcmp(sol_type,'mfs')
    ctrs = SphereSurfGoldPoints(floor(mfs_frac*Npnts), mfs_sphere*R);
else % For kansa, the centers and collocation points coincide
    ctrs = [intdata; bdydata];
end


% Compute the collocation block for the interior
DM_intdata = DistanceMatrix(intdata,ctrs);
LCM = Lrbf(ep,DM_intdata);
rhs_int = zeros(N_int,1);

% Compute the evaluation matrix
DM_eval = DistanceMatrix(evalpnts, ctrs);
EM = rbf(ep, DM_eval);

%   Notice the use of zeros(0,3), not []
%   To allow for bdydata_neu(:,1) calls later

bdydata_neu = bdydata;
normvecs = NORMALS.n11;

% Compute the collocation block for the boundary conditions
% This also computes the RHS for the problem
% First we consider the Neumann BC component
DM_bdydata_neu = DistanceMatrix(bdydata_neu,ctrs);

% Find all the necessary difference matrices
dx_bdydata_neu = DifferenceMatrix(bdydata_neu(:,1),ctrs(:,1));
dy_bdydata_neu = DifferenceMatrix(bdydata_neu(:,2),ctrs(:,2));
dz_bdydata_neu = DifferenceMatrix(bdydata_neu(:,3),ctrs(:,3));

% Compute normal derivative collocation matrix for boundary
A = repmat(normvecs(:,1),1,N_tot).*dxrbf(ep,DM_bdydata_neu,dx_bdydata_neu);
B = repmat(normvecs(:,2),1,N_tot).*dyrbf(ep,DM_bdydata_neu,dy_bdydata_neu);
C = repmat(normvecs(:,3),1,N_tot).*dzrbf(ep,DM_bdydata_neu,dz_bdydata_neu);
BCM_neu = A + B + C;

% Compute known-terms vector (a.k.a. righthand side vector)
% This requires the gradient of the unbounded potential at boundary
gradphi_F_neu = gradphiF_dip(bdydata_neu, srcpnts, dipmom, sig);
rhs_bdy_neu = -sum(normvecs.*gradphi_F_neu,2);


% Create the full linear system from the blocksand solve it
% Compose rhs
rhs = [rhs_int;rhs_bdy_neu];
% Compose collocation matrix in same order as rhs
CM = [LCM;BCM_neu];
% Coefficients for evaluation
coefs = linsolve(CM,rhs);

% Potential at evalpnts in the source free case
phi0 = EM * coefs;
% Potential at evalpnts (superposition of effects)
phi = phi0 + phi_F;

% If requested, compute the difference of the solution with a\
% reference point, arbitrarily chosen as evalpnts(1)

errvec = zeros(N_eval,1);
for i = 1:N_eval
    phi_comp = phi - phi(i);
    phi_true = phi_an - phi_an(i);
    errvec(i) = errcompute(phi_comp,phi_true);
end

figure
SurfacePlot_dip(evalpnts, errvec)
title('Error on the solution (potential difference) vs. reference point')