% multispheresourcetest_MFS.m
% This example studies the convergence rate for the potential problem in a
% three-layered sphere setting when the centers for the kernels are all
% placed outside the outermost sphere.

% Physical parameters
R = [0.087, 0.092, 0.1];
sig = [0.33, 0.0125, 0.33];
dipmom = [1, 0, 0];
srcpnts = [0, 0, 0.6*R(1)];
reference = [0,0,-R(end)];

% Number of evaluation points on the outer sphere
N_eval = 100;

% Number of terms for the series expansion of the semi-analytical solution
sol_acc = 200;

% Both interface conditions are imposed at each collocation point (1) or
% two differt sets of collocation points are considered (0)
match_couple = 0;

% Number of collocation points
Nvec = logspace(1,3,15);

% Ratio between the number of centers and the number of collocation points
mfs_frac = 1.0;

% Coefficient for the radius of the sphere the centers are located on
mfs_sphere = 1.2;

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
% rng(0);
rand('state',0);

% Find the sigma value where the dipole is located
% This is needed to evaluate phi_F and gradphi_F
src_loc = DistanceMatrix(srcpnts,[0,0,0]);
sig_dip = sig(find(src_loc<R,1,'first'));

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval-1, R(end));
evalpnts = [reference;evalpnts];

% Analytic solution for the potential
phi_an = MultiSpherePotential(R, sig, srcpnts, dipmom, evalpnts, sol_acc);

% Compute the difference of the solution with a reference point, which was 
% attached at the top of evalpnts earlier
phi_true = phi_an - phi_an(1);

% RBF definition and derivatives
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');


errvec = zeros(size(Nvec));
condvec = zeros(size(Nvec));
coefvec = zeros(size(Nvec));
Nvec_true = zeros(size(Nvec));
oldtest = zeros(N_eval,1);
newtest = zeros(N_eval,1);
errdvec = zeros(size(Nvec));
for k=1:length(Nvec)
    fprintf('k=%d\n',k)
    
    % Choose all points in the geometry
    [POINTS, NORMALS] = BallGeometry(R,round( [0.8 1 1.2]/3*Nvec(k) ),'mfs',[],'spiral');
    
    % Cut up the domain into the appropriate sections
    % Section A is the inner ball, B is the layer in the middle, C in the
    % outer layer
    A_cpl_out = POINTS.bdy11;
    B_cpl_in = POINTS.bdy12;
    B_cpl_out = POINTS.bdy22;
    C_cpl_in = POINTS.bdy23;
    C_bdy = POINTS.bdy33;
    
    %%% These are the normal vectors to the interface and boundary
    % Note that the interface normals both point toward the boundary, even
    % though that is technically the 'negative' normal for the outer ball,
    % since it is pointing in to rather than out of the outer domain
    A_cpl_out_nv = NORMALS.n11;
    B_cpl_in_nv = NORMALS.n12;
    B_cpl_out_nv = NORMALS.n22;
    C_cpl_in_nv = NORMALS.n23;
    C_bdy_nv = NORMALS.n33;
    
    % If the user wants, we can require the values and derivatives to be
    % matched at the same points rather than different points.  I don't
    % think it matters, but it's nice to have this option.
    if match_couple
        B_cpl_in = A_cpl_out;
        C_cpl_in = B_cpl_out;
        B_cpl_in_nv = A_cpl_out_nv;
        C_cpl_in_nv = B_cpl_out_nv;
    end
    
    % Classify all the sizes of each domain/boundary
    N_A_cpl_out = size(A_cpl_out,1);
    N_B_cpl_in = size(B_cpl_in,1);
    N_B_cpl_out = size(B_cpl_out,1);
    N_C_cpl_in = size(C_cpl_in,1);
    N_C_bdy = size(C_bdy,1);
    
    % Compose a vector of all the RBF centers
    N_ctrs = floor(mfs_frac*(N_A_cpl_out+N_B_cpl_in+N_B_cpl_out+N_C_cpl_in+N_C_bdy));
    all_ctrs = SphereSurfGoldPoints(N_ctrs, mfs_sphere*R(end));
    
    ctrs_perm = randperm(N_ctrs);
    A_ctrs_ind = ctrs_perm(1:floor(mfs_frac*N_A_cpl_out));
    B_ctrs_ind = ctrs_perm(length(A_ctrs_ind)+1:length(A_ctrs_ind)+floor(mfs_frac*(N_B_cpl_in+N_B_cpl_out)));
    C_ctrs_ind = ctrs_perm(length(B_ctrs_ind)+1:end);
    
%     A_ctrs_ind = randperm(N_ctrs,floor(mfs_frac*N_A_cpl_out));
%     
%     BC_ctrs_ind = setdiff(1:N_ctrs,A_ctrs_ind);
%     B_ctrs_ind = randperm(length(BC_ctrs_ind),floor(mfs_frac*(N_B_cpl_in+N_B_cpl_out)));
%     B_ctrs_ind = BC_ctrs_ind(B_ctrs_ind);
%     
%     C_ctrs_ind = setdiff(1:N_ctrs,[A_ctrs_ind B_ctrs_ind]);
    
    A_ctrs = all_ctrs(A_ctrs_ind,:);
    B_ctrs = all_ctrs(B_ctrs_ind,:);
    C_ctrs = all_ctrs(C_ctrs_ind,:);
    
    N_A = size(A_ctrs,1);
    N_B = size(B_ctrs,1);
    N_C = size(C_ctrs,1);
    
    % Compute the evaluation matrix
    DM_eval = DistanceMatrix(evalpnts, C_ctrs);
    EM = rbf([], DM_eval);
    
    % Compute the collocation block for the boundary conditions
    % This also computes the RHS for the problem
    DM_C_bdy = DistanceMatrix(C_bdy,C_ctrs);
    
    % Find all the necessary difference matrices
    C_bdy_dx = DifferenceMatrix(C_bdy(:,1),C_ctrs(:,1));
    C_bdy_dy = DifferenceMatrix(C_bdy(:,2),C_ctrs(:,2));
    C_bdy_dz = DifferenceMatrix(C_bdy(:,3),C_ctrs(:,3));
    
    % Compute normal derivative collocation matrix for boundary
    D1 = repmat(C_bdy_nv(:,1),1,N_C).*dxrbf([],DM_C_bdy,C_bdy_dx);
    D2 = repmat(C_bdy_nv(:,2),1,N_C).*dyrbf([],DM_C_bdy,C_bdy_dy);
    D3 = repmat(C_bdy_nv(:,3),1,N_C).*dzrbf([],DM_C_bdy,C_bdy_dz);
    BCM_C_bdy = D1 + D2 + D3;
    
    % Compute known-terms vector (a.k.a. righthand side vector)
    % This requires the gradient of the unbounded potential at boundary
    rhs_C_bdy = zeros(N_C_bdy,1);
    
    % Deal with the coupling region
    % First we consider the Dirichlet coupling condition
    % This requires values from B to be mapped to points on A
    DM_A_cpl_out_dir = DistanceMatrix(A_cpl_out,A_ctrs);
    CCM_A_cpl_out_dir = rbf([],DM_A_cpl_out_dir);
    
    DM_B_cpl_in_dir = DistanceMatrix(A_cpl_out,B_ctrs);
    CCM_B_cpl_in_dir = rbf([],DM_B_cpl_in_dir);
    
    DM_B_cpl_out_dir = DistanceMatrix(B_cpl_out,B_ctrs);
    CCM_B_cpl_out_dir = rbf([],DM_B_cpl_out_dir);
    
    DM_C_cpl_in_dir = DistanceMatrix(B_cpl_out,C_ctrs);
    CCM_C_cpl_in_dir = rbf([],DM_C_cpl_in_dir);
    
    rhs_A_cpl_out_dir = -phiF_dip(A_cpl_out,srcpnts,dipmom,sig_dip);
    rhs_B_cpl_out_dir = zeros(N_B_cpl_out,1);
    
    % Now we consider the Neumann coupling condition
    % Here, normal derivatives need to be mapped from A to points on B
    DM_A_cpl_out_neu = DistanceMatrix(B_cpl_in,A_ctrs);
    A_cpl_out_dx_neu = DifferenceMatrix(B_cpl_in(:,1),A_ctrs(:,1));
    A_cpl_out_dy_neu = DifferenceMatrix(B_cpl_in(:,2),A_ctrs(:,2));
    A_cpl_out_dz_neu = DifferenceMatrix(B_cpl_in(:,3),A_ctrs(:,3));
    D1 = repmat(B_cpl_in_nv(:,1),1,N_A).*dxrbf([],DM_A_cpl_out_neu,A_cpl_out_dx_neu);
    D2 = repmat(B_cpl_in_nv(:,2),1,N_A).*dyrbf([],DM_A_cpl_out_neu,A_cpl_out_dy_neu);
    D3 = repmat(B_cpl_in_nv(:,3),1,N_A).*dzrbf([],DM_A_cpl_out_neu,A_cpl_out_dz_neu);
    CCM_A_cpl_out_neu = sig(1)*(D1 + D2 + D3);
    
    DM_B_cpl_in_neu = DistanceMatrix(B_cpl_in,B_ctrs);
    B_cpl_in_dx_neu = DifferenceMatrix(B_cpl_in(:,1),B_ctrs(:,1));
    B_cpl_in_dy_neu = DifferenceMatrix(B_cpl_in(:,2),B_ctrs(:,2));
    B_cpl_in_dz_neu = DifferenceMatrix(B_cpl_in(:,3),B_ctrs(:,3));
    D1 = repmat(B_cpl_in_nv(:,1),1,N_B).*dxrbf([],DM_B_cpl_in_neu,B_cpl_in_dx_neu);
    D2 = repmat(B_cpl_in_nv(:,2),1,N_B).*dyrbf([],DM_B_cpl_in_neu,B_cpl_in_dy_neu);
    D3 = repmat(B_cpl_in_nv(:,3),1,N_B).*dzrbf([],DM_B_cpl_in_neu,B_cpl_in_dz_neu);
    CCM_B_cpl_in_neu = sig(2)*(D1 + D2 + D3);
    
    DM_B_cpl_out_neu = DistanceMatrix(C_cpl_in,B_ctrs);
    B_cpl_out_dx_neu = DifferenceMatrix(C_cpl_in(:,1),B_ctrs(:,1));
    B_cpl_out_dy_neu = DifferenceMatrix(C_cpl_in(:,2),B_ctrs(:,2));
    B_cpl_out_dz_neu = DifferenceMatrix(C_cpl_in(:,3),B_ctrs(:,3));
    D1 = repmat(C_cpl_in_nv(:,1),1,N_B).*dxrbf([],DM_B_cpl_out_neu,B_cpl_out_dx_neu);
    D2 = repmat(C_cpl_in_nv(:,2),1,N_B).*dyrbf([],DM_B_cpl_out_neu,B_cpl_out_dy_neu);
    D3 = repmat(C_cpl_in_nv(:,3),1,N_B).*dzrbf([],DM_B_cpl_out_neu,B_cpl_out_dz_neu);
    CCM_B_cpl_out_neu = sig(2)*(D1 + D2 + D3);
    
    DM_C_cpl_in_neu = DistanceMatrix(C_cpl_in,C_ctrs);
    C_cpl_in_dx_neu = DifferenceMatrix(C_cpl_in(:,1),C_ctrs(:,1));
    C_cpl_in_dy_neu = DifferenceMatrix(C_cpl_in(:,2),C_ctrs(:,2));
    C_cpl_in_dz_neu = DifferenceMatrix(C_cpl_in(:,3),C_ctrs(:,3));
    D1 = repmat(C_cpl_in_nv(:,1),1,N_C).*dxrbf([],DM_C_cpl_in_neu,C_cpl_in_dx_neu);
    D2 = repmat(C_cpl_in_nv(:,2),1,N_C).*dyrbf([],DM_C_cpl_in_neu,C_cpl_in_dy_neu);
    D3 = repmat(C_cpl_in_nv(:,3),1,N_C).*dzrbf([],DM_C_cpl_in_neu,C_cpl_in_dz_neu);
    CCM_C_cpl_in_neu = sig(3)*(D1 + D2 + D3);
    
    gradphiF_B_cpl_in = gradphiF_dip(B_cpl_in,srcpnts,dipmom,sig_dip);
    rhs_B_cpl_in_neu = -sig(1)*sum(B_cpl_in_nv.*gradphiF_B_cpl_in, 2);
    
    rhs_C_cpl_in_neu = zeros(N_C_cpl_in,1);
    
    % Create the full linear system from the blocks
    % Compose rhs
    rhs = [rhs_A_cpl_out_dir; ...
           rhs_B_cpl_in_neu; ...
           rhs_B_cpl_out_dir; ...
           rhs_C_cpl_in_neu; ...
           rhs_C_bdy];
    % Compose collocation matrix in same order as rhs
    % Notice the coupling conditions are of the form A-B
    % That's why the (-) appears before the B coupling components
    CM = [CCM_A_cpl_out_dir,-CCM_B_cpl_in_dir,zeros(N_A_cpl_out,N_C); ...
          CCM_A_cpl_out_neu,-CCM_B_cpl_in_neu,zeros(N_B_cpl_in,N_C); ...
          zeros(N_B_cpl_out,N_A),CCM_B_cpl_out_dir,-CCM_C_cpl_in_dir; ...
          zeros(N_C_cpl_in,N_A),CCM_B_cpl_out_neu,-CCM_C_cpl_in_neu; ...
          zeros(N_C_bdy,N_A+N_B),BCM_C_bdy];
    
    % Solve the BVP linear system
    coefs = linsolve(CM,rhs);
    
    % Extract the oefficients for evaluation
    C_coefs = coefs(N_A+N_B+1:end);
    % Potential at evalpnts in the source free case
    phi = EM * C_coefs;
    phi_comp(:,k) = phi - phi(1);
    
    % Compute the total error
    errvec(k) = errcompute(phi_comp(:,k),phi_true);
    
    condvec(k) = cond(CM);
    
    coefvec(k) = norm(C_coefs);
    
    Nvec_true(k) = N_A_cpl_out + N_B_cpl_in + ...
                   N_B_cpl_out + N_C_cpl_in + N_C_bdy;
    newtest = phi_comp(:,k);
    if k>1
        errdvec(k) = errcompute(oldtest,newtest);
        if k==2
            errdvec(1) = errdvec(2);
        end
    end
    oldtest = newtest;
end

h = figure;
loglog(Nvec_true,[errvec;condvec;coefvec;errdvec],'linewidth',3)
xlabel('Collocation points')
legend('Error','Cond','Coef Norm','Test Diff')
title('Potential Problem - Dipole in a Three-layered Sphere')