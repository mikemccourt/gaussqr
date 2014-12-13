% realgeom_test.m is a test script for the MFS solution of the EEG forward
% problem with realistic head geometries. 
%
% The accuracy of the solution is evaluated with respect to a BEM reference 
% solution obtained with a fine mesh by means of the 2-norm relative 
% difference.
% The potential is evaluated at a set of points almost uniformly
% distributed on the scalp.
%
% For details on the meaning of the parameters, take the script
% 'threesphere_paper.m' as a reference since this code is mostly based on
% that file.
%
% WARNING: this script needs a set of MAT files in the path to run. See the
% following for further details.

%%%%% Input data load and preliminary operations
% Reset the random numbers generator
rng(0);

% Define the vector containing the number of collocation points
% You need a MAT file containing the geometry data for each element in this
% vector. The file name must be in the form 'geometryMFS_xxx.mat', where
% xxx is the number of collocation points.
Nvec = [500 600 700 800 900 1000 1200 1500 2000 3000 6000 8000 9000];
N_Nvec = length(Nvec);

% Physical parameters
% If this parameters are changed, the reference solution must be 
% re-evaluated accordingly!
sig = [0.33, 0.0125, 0.33];
sig_dip  = sig(1);
srcpnt = [0 0 0.055];
dipmom = [1 0 0];

% Load the BEM reference solution from the appropriate MAT file
load reference_BEM_sol.mat

% MFS solver parameters
% mfs_frac = 0.8; Defined in loop below
coeff = [1.3 0.85 1.3 0.85 1.3];
match_couple = 1;

% Load the functions the MFS solver needs
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');

% Set to 1 if you want to check what goes on during the simulation
iter_out = 1;

% Load the evaluation points on the scalp from the appropriate MAT file
load eval_points.mat; N_eval = length(evalpnts);

%%%%% MFS solution
for mfs_frac = [.4 .6 .8]
phi_MFS = zeros(N_eval,N_Nvec);
N_rows = zeros(N_Nvec,1); N_ctrs = zeros(N_Nvec,1);
elapsed_t_MFS = zeros(N_Nvec,1); errvec = zeros(N_Nvec,1);
for k = 1:N_Nvec
    if iter_out
        fprintf('k=%d\n',k)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the geometry from the appropriate file
    filename = strcat('geometryMFS_',num2str(Nvec(k)),'.mat');
    load(filename);
    
    % Cut up the domain into the appropriate sections
    % Section A is the brain, B is the skull in the middle, C is the
    % scalp
    % If the user wants, we can require the values and derivatives to be
    % matched at the same points rather than different points.
    C_bdy = GEOMETRY_MFS.scalp.Points;
    C_bdy_nv = GEOMETRY_MFS.scalp.Normals;
    if match_couple
        A_cpl_out = GEOMETRY_MFS.inner_skull.Points;
        B_cpl_in = GEOMETRY_MFS.inner_skull.Points;
        B_cpl_out = GEOMETRY_MFS.outer_skull.Points;
        C_cpl_in = GEOMETRY_MFS.outer_skull.Points;
        
        A_cpl_out_nv = GEOMETRY_MFS.inner_skull.Normals;
        B_cpl_in_nv = GEOMETRY_MFS.inner_skull.Normals;
        B_cpl_out_nv = GEOMETRY_MFS.outer_skull.Normals;
        C_cpl_in_nv = GEOMETRY_MFS.outer_skull.Normals;
        
        % Classify all the sizes of each domain/boundary
        N_A_cpl_out = size(A_cpl_out,1);
        N_B_cpl_in = size(B_cpl_in,1);
        N_B_cpl_out = size(B_cpl_out,1);
        N_C_cpl_in = size(C_cpl_in,1);
        N_C_bdy = size(C_bdy,1);
        
        N_innsk = N_A_cpl_out; % # Collocation points on the inner skull
        N_outsk = N_B_cpl_out;  % # Collocation points on the outer skull
        N_scalp = N_C_bdy; % # Collocation points on the scalp
    else
        % Half for the continuity of the potential, half for the
        % continuity of the normal component of the current density
        ind1 = randperm(VpS,round(VpS/2));
        ind2 = setdiff(1:VpS,ind1);
        
        A_cpl_out = GEOMETRY_MFS.inner_skull.Points(ind1,:);
        B_cpl_in = GEOMETRY_MFS.inner_skull.Points(ind2,:);
        B_cpl_out = GEOMETRY_MFS.outer_skull.Points(ind1,:);
        C_cpl_in = GEOMETRY_MFS.outer_skull.Points(ind2,:);
        
        A_cpl_out_nv = GEOMETRY_MFS.inner_skull.Normals(ind1,:);
        B_cpl_in_nv = GEOMETRY_MFS.inner_skull.Normals(ind2,:);
        B_cpl_out_nv = GEOMETRY_MFS.outer_skull.Normals(ind1,:);
        C_cpl_in_nv = GEOMETRY_MFS.outer_skull.Normals(ind2,:);
        
        % Classify all the sizes of each domain/boundary
        N_A_cpl_out = size(A_cpl_out,1);
        N_B_cpl_in = size(B_cpl_in,1);
        N_B_cpl_out = size(B_cpl_out,1);
        N_C_cpl_in = size(C_cpl_in,1);
        N_C_bdy = size(C_bdy,1);
        
        N_innsk = N_A_cpl_out + N_B_cpl_in; % # Collocation points on the inner skull
        N_outsk = N_B_cpl_out + N_C_cpl_in;  % # Collocation points on the outer skull
        N_scalp = N_C_bdy; % # Collocation points on the scalp
    end
    
    % Define the RBF centers
    % Inflate/deflate collocation point randomly reduced clouds
    ind1 = randperm(N_innsk,round(mfs_frac*N_innsk));
    A_ctrs = centers(GEOMETRY_MFS.inner_skull.Points(ind1,:), coeff(1));
    
    ind1 = randperm(N_innsk,round(mfs_frac*N_innsk));
    B_ctrs_in = centers(GEOMETRY_MFS.inner_skull.Points(ind1,:), coeff(2));
    
    ind1 = randperm(N_outsk,round(mfs_frac*N_outsk));
    B_ctrs_out = centers(GEOMETRY_MFS.outer_skull.Points(ind1,:), coeff(3));
    B_ctrs = [B_ctrs_in; B_ctrs_out];
    
    ind1 = randperm(N_outsk,round(mfs_frac*N_outsk));
    C_ctrs_in = centers(GEOMETRY_MFS.outer_skull.Points(ind1,:), coeff(4));
    ind1 = randperm(N_scalp,round(mfs_frac*N_scalp));
    C_ctrs_out = centers(GEOMETRY_MFS.scalp.Points(ind1,:), coeff(5));
    C_ctrs = [C_ctrs_in; C_ctrs_out];
    
    N_A = size(A_ctrs,1);
    N_B = size(B_ctrs,1);
    N_C = size(C_ctrs,1);
    
    tic;
    
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
    D1 = repmat(C_bdy_nv(:,1),1,N_C).*dxrbf(1,DM_C_bdy,C_bdy_dx);
    D2 = repmat(C_bdy_nv(:,2),1,N_C).*dyrbf(1,DM_C_bdy,C_bdy_dy);
    D3 = repmat(C_bdy_nv(:,3),1,N_C).*dzrbf(1,DM_C_bdy,C_bdy_dz);
    BCM_C_bdy = D1 + D2 + D3;
    
    % Compute known-terms vector (a.k.a. righthand side vector)
    % This requires the gradient of the unbounded potential at boundary
    % gradphi_F_bdy = gradphiF_dip(C_bdy, srcpnt, dipmom, sig_dip);
    rhs_C_bdy = zeros(N_C_bdy,1);
    
    % Deal with the coupling region
    % First we consider the Dirichlet coupling condition
    % This requires values from B to be mapped to points on A
    DM_A_cpl_out_dir = DistanceMatrix(A_cpl_out,A_ctrs);
    CCM_A_cpl_out_dir = rbf(1,DM_A_cpl_out_dir);
    
    DM_B_cpl_in_dir = DistanceMatrix(A_cpl_out,B_ctrs);
    CCM_B_cpl_in_dir = rbf(1,DM_B_cpl_in_dir);
    
    DM_B_cpl_out_dir = DistanceMatrix(B_cpl_out,B_ctrs);
    CCM_B_cpl_out_dir = rbf(1,DM_B_cpl_out_dir);
    
    DM_C_cpl_in_dir = DistanceMatrix(B_cpl_out,C_ctrs);
    CCM_C_cpl_in_dir = rbf(1,DM_C_cpl_in_dir);
    
    rhs_A_cpl_out_dir = -phiF_dip(A_cpl_out,srcpnt,dipmom,sig_dip);
    rhs_B_cpl_out_dir = zeros(N_B_cpl_out,1);
    
    % Now we consider the Neumann coupling condition
    % Here, normal derivatives need to be mapped from A to points on B
    DM_A_cpl_out_neu = DistanceMatrix(B_cpl_in,A_ctrs);
    A_cpl_out_dx_neu = DifferenceMatrix(B_cpl_in(:,1),A_ctrs(:,1));
    A_cpl_out_dy_neu = DifferenceMatrix(B_cpl_in(:,2),A_ctrs(:,2));
    A_cpl_out_dz_neu = DifferenceMatrix(B_cpl_in(:,3),A_ctrs(:,3));
    D1 = repmat(B_cpl_in_nv(:,1),1,N_A).*dxrbf(1,DM_A_cpl_out_neu,A_cpl_out_dx_neu);
    D2 = repmat(B_cpl_in_nv(:,2),1,N_A).*dyrbf(1,DM_A_cpl_out_neu,A_cpl_out_dy_neu);
    D3 = repmat(B_cpl_in_nv(:,3),1,N_A).*dzrbf(1,DM_A_cpl_out_neu,A_cpl_out_dz_neu);
    CCM_A_cpl_out_neu = sig(1)*(D1 + D2 + D3);
    
    DM_B_cpl_in_neu = DistanceMatrix(B_cpl_in,B_ctrs);
    B_cpl_in_dx_neu = DifferenceMatrix(B_cpl_in(:,1),B_ctrs(:,1));
    B_cpl_in_dy_neu = DifferenceMatrix(B_cpl_in(:,2),B_ctrs(:,2));
    B_cpl_in_dz_neu = DifferenceMatrix(B_cpl_in(:,3),B_ctrs(:,3));
    D1 = repmat(B_cpl_in_nv(:,1),1,N_B).*dxrbf(1,DM_B_cpl_in_neu,B_cpl_in_dx_neu);
    D2 = repmat(B_cpl_in_nv(:,2),1,N_B).*dyrbf(1,DM_B_cpl_in_neu,B_cpl_in_dy_neu);
    D3 = repmat(B_cpl_in_nv(:,3),1,N_B).*dzrbf(1,DM_B_cpl_in_neu,B_cpl_in_dz_neu);
    CCM_B_cpl_in_neu = sig(2)*(D1 + D2 + D3);
    
    DM_B_cpl_out_neu = DistanceMatrix(C_cpl_in,B_ctrs);
    B_cpl_out_dx_neu = DifferenceMatrix(C_cpl_in(:,1),B_ctrs(:,1));
    B_cpl_out_dy_neu = DifferenceMatrix(C_cpl_in(:,2),B_ctrs(:,2));
    B_cpl_out_dz_neu = DifferenceMatrix(C_cpl_in(:,3),B_ctrs(:,3));
    D1 = repmat(C_cpl_in_nv(:,1),1,N_B).*dxrbf(1,DM_B_cpl_out_neu,B_cpl_out_dx_neu);
    D2 = repmat(C_cpl_in_nv(:,2),1,N_B).*dyrbf(1,DM_B_cpl_out_neu,B_cpl_out_dy_neu);
    D3 = repmat(C_cpl_in_nv(:,3),1,N_B).*dzrbf(1,DM_B_cpl_out_neu,B_cpl_out_dz_neu);
    CCM_B_cpl_out_neu = sig(2)*(D1 + D2 + D3);
    
    DM_C_cpl_in_neu = DistanceMatrix(C_cpl_in,C_ctrs);
    C_cpl_in_dx_neu = DifferenceMatrix(C_cpl_in(:,1),C_ctrs(:,1));
    C_cpl_in_dy_neu = DifferenceMatrix(C_cpl_in(:,2),C_ctrs(:,2));
    C_cpl_in_dz_neu = DifferenceMatrix(C_cpl_in(:,3),C_ctrs(:,3));
    D1 = repmat(C_cpl_in_nv(:,1),1,N_C).*dxrbf(1,DM_C_cpl_in_neu,C_cpl_in_dx_neu);
    D2 = repmat(C_cpl_in_nv(:,2),1,N_C).*dyrbf(1,DM_C_cpl_in_neu,C_cpl_in_dy_neu);
    D3 = repmat(C_cpl_in_nv(:,3),1,N_C).*dzrbf(1,DM_C_cpl_in_neu,C_cpl_in_dz_neu);
    CCM_C_cpl_in_neu = sig(3)*(D1 + D2 + D3);
    
    % Gotta check to make sure these are correct for different sig values
    % in the different domains.  No issues if they are all equal though.
    gradphiF_B_cpl_in = gradphiF_dip(B_cpl_in,srcpnt,dipmom,sig_dip);
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
    
    [N_rows(k), N_ctrs(k)] = size(CM);
    
    % Solve the BVP linear system
    coefs = linsolve(CM,rhs);
    
    % Extract the oefficients for evaluation
    C_coefs = coefs(N_A+N_B+1:end);
    
    % Potential at evalpnts in the source free case
    phi_MFS(:,k) = EM * C_coefs;
    % Refer the potential to its mean value on the scalp
    phi_MFS(:,k) = phi_MFS(:,k) - sum(phi_MFS(:,k))/N_eval;
    
    elapsed_t_MFS(k) = toc;
    
    errvec(k) = norm(phi_MFS(:,k)-phi_BEM(:,end))/norm(phi_BEM(:,end));
    
    if iter_out
        fprintf('Elapsed time=%d s\n2-norm relative error=%d\n\n',...
                elapsed_t_MFS(k),errvec(k))
    end
end
save(sprintf('dipole_time_tests%s.mat',int2str(floor(10*mfs_frac))),'elapsed_t_MFS','errvec')
end