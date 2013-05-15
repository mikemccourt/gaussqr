% singlesphere_N.m
%
%  For this problem, we consider the solution to the Laplace equation on a
%  sphere with Neumann boundary conditions.
%
%  The problem has several physical parameters relating to the
%  underlying EEG/MEG physical system.  These parameters are:
%    R - Sphere radius [dm] <default = 1>
%    sig - Electric conductivity [S/dm] <default = .02>
%    dipmom - Dipole moment [x10^-12 Am] <default = 2.7*[1,0,0]>
%    srcpnts - Dipole position [dm] <default = [0,0,0.6*R]>
%
%  This script allows you to test the convergence rate (with respect to N
%  of different RBFs and different epsilon values.
%
%  The solution parameters to be considered are
%     sol_type - How you want to solve the system <default = 'kansa'>
%                'kansa' : Nonsymmetric collocation
%                          rbf_choice and ep must also be chosen
%                'mfs' : Method of fundamental solutions
%                        MFS_frac and ctr_sphere must also be chosen
%     rbf_choice - RBF for collocation <default = 'imq'>
%     int_point_dist - How the (collocation) points are spread in the ball
%                  'halton' : Halton cube, restricted to ball <default>
%                  'even' : Gridden points, restricted to ball
%                  'random' : Uniform random
%                  'cheb' : Chebyshev spaced along radii of the ball
%     bdy_point_dist - How the (collocation) points are spread on the surface
%                  'spiral' : evenly (golden spiral) <default> 
%                  'halton' : quasi-randomly (Halton)
%     ep - RBF shape parameter <default = 10>
%     mfs_frac - How many centers for MFS, in [0.0,1.0]*N <default = 0.4>
%     mfs_sphere - Fraction beyond R (eg, 1.3R) for centers <default = 1.3>
%     BC_choice - How to choose the boundary conditions <default = 1>
%                 1 : Neumann
%                 2 : Dirichlet
%                 3 : Mixed (6 random Dirichlet, the rest Neumann)
%                 4 : Neumann and one Dirichlet at a reference point (zero 
%                     potential point)
%     refpnt - Reference point (zero potential); this is used only if
%              BC_choice = 4 <default = [0,0,R]>
%     eval_diff - Consider the solution as the difference between all
%                 values and a reference point <default = 1>
%
%  The value we are interested in studying is the effect of increasing N,
%  so you must specify a vector of N values that you want to study
%     Npnts - Number of collocation points <default = 1000>
%     Nvec_samples - Row vector of sample sizes for integration 
%                    <default = 100:2000:31000;>
%     BC_frac - The fraction of the total points to be used to enforce
%               boundary conditions <default = .3>
%     dip_cushion - How much space should be given around the dipole where
%               no RBF centers are allowed <default = .005>
%     N_eval_phi - # points to evaluate error on the potential <default = 1001>
%     N_eval_B - # points to evaluate error on the magnetic field <default = 400>
%
%  Some outputs are available if you would like them
%     iter_out - Print output during the solves <default = 1>
%     plot_err - log-log plot of error vs. sample size <default = 1>
%     errcolor - Color for error line in log-log plot <default = 'b'>
%
%  The results of these experiments are stored in
%     Nvec_true - Actual number of collocation points, because the point
%                 distribution in a sphere is tricky
%     errvec_phi - Error on potential computed at Nvec_true
%     errvec_B - Errors on magnetic field computed at Nvec_samples
%
%  Note that if MFS with fewer centers than collocation points is chosen,
%  condition doesn't make sense (rectangular, not square, system)

R = 1;
sig = 0.02;
dipmom = 2.7.*[1, 0, 0];
srcpnts = [0, 0, 0.6*R];

sol_type = 'mfs';
radbasfun = 'imq';
int_point_dist = 'halton';
bdy_point_dist = 'spiral';
ep = 1;
mfs_frac = 0.4;
mfs_sphere = 1.3;
BC_choice = 1;
refpnt = [0, 0, R];
eval_diff = 1;

Npnts = 1000;
Nvec_samples = 100:2000:31000;
BC_frac = .3; % Not yet implemented
dip_cushion = .05;
N_eval_phi = 1001;
N_eval_B = 300;

iter_out = 1;
plot_err = 1;
errcolor = 'b';


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

% Determine the evaluation points for phi and B
evalpnts_phi = SphereSurfGoldPoints(N_eval_phi, R);
evalpnts_B = SphereSurfGoldPoints(N_eval_B, 1.05*R);

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts_phi,srcpnts,dipmom,sig);

% Analytic solution for the potential
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts_phi);

% If requested, compute the difference of the solution with a reference
% point, arbitrarily chosen as evalpnts(1)
% Similarly, if Neumann BC with a Dirichlet condition at the reference 
% point are used, compute the difference of the solution with the reference
% point 
if eval_diff || BC_choice == 4
    if BC_choice == 4 
        phi_ref = HomSpherePotential(R, sig, srcpnts, dipmom, refpnt);
        phi_true = phi_an - phi_ref;
    else
        phi_true = phi_an - phi_an(1);
    end
else
    phi_true = phi_an;
end

% Determine collocation points
[POINTS, NORMALS] = BallGeometry(R,Npnts,sol_type,int_point_dist,bdy_point_dist,srcpnts,dip_cushion);
intdata = POINTS.int1;
bdydata = POINTS.bdy11;
N_int = size(intdata,1);
N_bdy = size(bdydata,1);
N_tot = N_int + N_bdy;  % Total number of collocation points

% Compose a vector of all the RBF centers
% In the MFS setting, these are chosen in a sphere around the ball
if strcmp(sol_type,'mfs')
    N_ctrs = floor(mfs_frac*Npnts);
    ctrs = SphereSurfGoldPoints(N_ctrs, mfs_sphere*R);
else % For kansa, the centers and collocation points coincide
    if BC_choice == 4 % In this case we need an "extra" center
        % (reference point)
        ctrs = [intdata; bdydata; refpnt];
        N_tot = N_tot + 1;
    else
        ctrs = [intdata; bdydata];
    end
    N_ctrs = N_tot;
end

% Compute the collocation block for the interior
DM_intdata = DistanceMatrix(intdata,ctrs);
LCM = Lrbf(ep,DM_intdata);
rhs_int = zeros(N_int,1);

% Compute the evaluation matrix
DM_eval = DistanceMatrix(evalpnts_phi, ctrs);
EM = rbf(ep, DM_eval);

% Determine which points are Neumann and which are Dirichlet
%   Notice the use of zeros(0,3), not []
%   To allow for bdydata_neu(:,1) calls later
if BC_choice==1 % Do the standard Neumann BC
    bdydata_neu = bdydata;
    normvecs = NORMALS.n11;
    bdydata_dir = zeros(0,3);
elseif BC_choice==2 % Run a test with Dirichlet BC
    bdydata_neu = zeros(0,3);
    normvecs = zeros(0,3);
    bdydata_dir = bdydata;
elseif BC_choice==3 % Run a test with Mixed BC
    % Right now, fixed at 6 Dirichlet BC points
    % Could be variable, but not important
    N_dir = min(6,N_bdy);
    i_dir = randperm(N_bdy,N_dir);
    i_neu = setdiff(1:N_bdy,i_dir);
    
    bdydata_neu = bdydata(i_neu,:);
    normvecs = NORMALS.n11(i_neu,:);
    bdydata_dir = bdydata(i_dir,:);
else % Run a test with Neumann BC + Dirichlet BC at the reference point
    % (zero potential point)
    bdydata_neu = bdydata;
    normvecs = NORMALS.n11;
    bdydata_dir = refpnt;
end

% Compute the collocation block for the boundary conditions
% This also computes the RHS for the problem
% First we consider the Neumann BC component
DM_bdydata_neu = DistanceMatrix(bdydata_neu,ctrs);

% Find all the necessary difference matrices
dx_bdydata_neu = DifferenceMatrix(bdydata_neu(:,1),ctrs(:,1));
dy_bdydata_neu = DifferenceMatrix(bdydata_neu(:,2),ctrs(:,2));
dz_bdydata_neu = DifferenceMatrix(bdydata_neu(:,3),ctrs(:,3));

% Compute normal derivative collocation matrix for boundary
A = repmat(normvecs(:,1),1,N_ctrs).*dxrbf(ep,DM_bdydata_neu,dx_bdydata_neu);
B = repmat(normvecs(:,2),1,N_ctrs).*dyrbf(ep,DM_bdydata_neu,dy_bdydata_neu);
C = repmat(normvecs(:,3),1,N_ctrs).*dzrbf(ep,DM_bdydata_neu,dz_bdydata_neu);
BCM_neu = A + B + C;

% Compute known-terms vector (a.k.a. righthand side vector)
% This requires the gradient of the unbounded potential at boundary
gradphi_F_neu = gradphiF_dip(bdydata_neu, srcpnts, dipmom, sig);
rhs_bdy_neu = -sum(normvecs.*gradphi_F_neu,2);


% Now we consider the Dirichlet BC component
DM_bdydata_dir = DistanceMatrix(bdydata_dir,ctrs);
BCM_dir = rbf(ep,DM_bdydata_dir);

if BC_choice == 4
    % Dirichlet condition at the reference point (zero potential)
    rhs_bdy_dir = 0;
else
    % Compute the true solution to be used as Dirichlet BC
    phi_F_bdy_dir = phiF_dip(bdydata_dir,srcpnts,dipmom,sig);
    phi_bdy_dir = HomSpherePotential(R, sig, srcpnts, dipmom, bdydata_dir);
    rhs_bdy_dir = phi_bdy_dir - phi_F_bdy_dir;
end

% Create the full linear system from the blocksand solve it
% Compose rhs
rhs = [rhs_int;rhs_bdy_neu;rhs_bdy_dir];
% Compose collocation matrix in same order as rhs
CM = [LCM;BCM_neu;BCM_dir];
% Coefficients for evaluation
[coefs,recip_cond] = linsolve(CM,rhs);

% Potential at evalpnts in the source free case
phi0 = EM * coefs;
% Potential at evalpnts (superposition of effects)
phi = phi0 + phi_F;

% If requested, compute the difference of the solution with a
% reference point, arbitrarily chosen as evalpnts(1)
if eval_diff && BC_choice ~= 4
    phi_comp = phi - phi(1);
else
    phi_comp = phi;
end

% Compute the total errors
errvec = errcompute(phi_comp,phi_true);
% Store the condition of the system
% For a low-rank system, instead store the rank
% This may happen for some MFS problems
if floor(recip_cond)==recip_cond
    condvec = recip_cond;
else
    condvec = 1/recip_cond;
end
Nvec_true = N_tot;

fprintf('err_phi = %g\ncond = %g\nN = %d\n\n',errvec,condvec,Nvec_true);

% Compute analytic magnetic (induction) field
B_an = HomSphereInduction(srcpnts,dipmom,evalpnts_B);

% Loop through the requested N values
errvec_B = [];
k=1;
for N = Nvec_samples
    if iter_out
        fprintf('k=%d\n',k)
    end
    B = Induction(evalpnts_B,srcpnts,dipmom,sol_type,radbasfun,ep,coefs,ctrs,sig,R,N );
    errvec_B(k) = norm (B - B_an)/norm(B_an);
    if iter_out
        fprintf('\terr_B = %g\n\tN = %d\n',errvec_B(k),Nvec_samples(k));
    end
    k=k+1;
end

if plot_err
    H1 = loglog(Nvec_samples,errvec_B);
    xlabel('No. of samples')
    ylabel('Relative error')
    set(H1,'LineWidth',3,'Color',errcolor)
    title('Magnetic field convergence test')
end