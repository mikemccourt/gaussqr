%% threesphere_paper.m
% This script runs the simulation we need for the paper (homogeneous
% sphere case): it is based on multispheresphere_N.m, so look at the help
% of that file for further information.
% Here we investigate convergence, conditioning and cost of our MFS solver
% (with different mfs_frac) and make comparisons with BEM.
%

% Choose the solvers you want to use
BEM = 1;
MFS = 1;
kansa = 0;

R = [0.087, 0.092, 0.1];
sig = [0.33, 0.0125, 0.33];
dipmom = [1, 0, 0];
srcpnts = [0, 0, 0.6*R(1)];

radbasfun = 'imq';
ep = 10;

int_point_dist = 'halton';
bdy_point_dist = 'spiral';

mfs_frac = [1.0 0.8 0.4 0.2];
mfs_sphere = [1.5 0.8 1.5 0.8 1.5];

reference = [0,0,-R(end)];

match_couple = 0;
sol_acc = 200;

Nvec = 350:500:5350;
N_eval = 1000;

iter_out = 1;
plot_sol = 1;
sol_err_style = 1;

% N_plotsol = length(Nvec);
N_plotsol = 3;
mfs_frac_plotsol = 1;

resultsfilename = 'three_sphere_tests';


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
% Eventually, I need to come up with a way to call this regardless of what
% version of Matlab you are running.  Think about it ...
rng(0);

% Find the sigma value where the dipole is located
% This is needed to evaluate phi_F and gradphi_F
src_loc = DistanceMatrix(srcpnts,[0,0,0]);
sig_dip = sig(find(src_loc<R,1,'first'));


%%%%%%%%%%%%%%%%%%%%%
% This is the start of the solver

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval-1, R(end));
evalpnts = [reference;evalpnts];

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
tic;
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig_dip);
tt = toc;

% Analytic solution for the potential
phi_an = MultiSpherePotential(R, sig, srcpnts, dipmom, evalpnts, sol_acc);
% G = MultiSpherePotentialLF(srcpnts,evalpnts,[0 0 0],R,sig); 
% phi_an = G*dipmom';

% Compute the difference of the solution with a reference point, which was 
% attached at the top of evalpnts earlier
phi_true = phi_an - phi_an(1);

%% BEM solutions
N_elements = zeros(length(Nvec),3);
if BEM
    phi_comp_BEM = zeros(N_eval,length(Nvec));
    for k=1:length(Nvec)
        if iter_out
            fprintf('k=%d\n',k)
        end
        [vol, N_elements(k,:)] = prepare_multisphere_mesh(Nvec(k),R);
        [phi_comp_BEM(:,k),~,~,elapsed_t_BEM(k)] = BEM_potential(vol,sig,evalpnts,srcpnts,dipmom);
        errvec_BEM(k) = errcompute(phi_comp_BEM(:,k),phi_true);
        Nvec_true_BEM(k) = sum(N_elements(k,:));
        if iter_out
            fprintf('\terr = %g\n\tN = %d\n',errvec_BEM(k),Nvec_true_BEM(k));
        end
    end
else
    for k=1:length(Nvec)
        [~, N_elements(k,:)] = prepare_multisphere_mesh(Nvec(k),R);
    end
end

%% MFS solutions
if MFS
    % RBF definition and derivatives
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');
    
    phi_comp_MFS = zeros(N_eval,length(Nvec),length(mfs_frac));
    for i=1:length(mfs_frac)
        
        for k=1:size(N_elements,1)
            
            if iter_out
                fprintf('k=%d\n',k)
            end
            
            % Choose all points in the geometry
            [POINTS, NORMALS] = BallGeometry(R,N_elements(k,:),'mfs',int_point_dist,bdy_point_dist);
            
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
            A_ctrs = SphereSurfGoldPoints(floor(mfs_frac(i)*N_A_cpl_out), mfs_sphere(1)*R(1) );
            B_ctrs_in = SphereSurfGoldPoints(floor(mfs_frac(i)*N_B_cpl_in), mfs_sphere(2)*R(1) );
            B_ctrs_out = SphereSurfGoldPoints(floor(mfs_frac(i)*N_B_cpl_out), mfs_sphere(3)*R(2) );
            B_ctrs = [B_ctrs_in; B_ctrs_out];
            C_ctrs_in = SphereSurfGoldPoints(floor(mfs_frac(i)*N_C_cpl_in), mfs_sphere(4)*R(2) );
            C_ctrs_out = SphereSurfGoldPoints(floor(mfs_frac(i)*N_C_bdy), mfs_sphere(5)*R(3) );
            C_ctrs = [C_ctrs_in; C_ctrs_out];
            N_A = size(A_ctrs,1);
            N_B = size(B_ctrs,1);
            N_C = size(C_ctrs,1);
            
            tic;
            
            % Compute the evaluation matrix
            DM_eval = DistanceMatrix(evalpnts, C_ctrs);
            EM = rbf(ep, DM_eval);
            
            % Compute the collocation block for the boundary conditions
            % This also computes the RHS for the problem
            DM_C_bdy = DistanceMatrix(C_bdy,C_ctrs);
            
            % Find all the necessary difference matrices
            C_bdy_dx = DifferenceMatrix(C_bdy(:,1),C_ctrs(:,1));
            C_bdy_dy = DifferenceMatrix(C_bdy(:,2),C_ctrs(:,2));
            C_bdy_dz = DifferenceMatrix(C_bdy(:,3),C_ctrs(:,3));
            
            % Compute normal derivative collocation matrix for boundary
            D1 = repmat(C_bdy_nv(:,1),1,N_C).*dxrbf(ep,DM_C_bdy,C_bdy_dx);
            D2 = repmat(C_bdy_nv(:,2),1,N_C).*dyrbf(ep,DM_C_bdy,C_bdy_dy);
            D3 = repmat(C_bdy_nv(:,3),1,N_C).*dzrbf(ep,DM_C_bdy,C_bdy_dz);
            BCM_C_bdy = D1 + D2 + D3;
            
            % Compute known-terms vector (a.k.a. righthand side vector)
            % This requires the gradient of the unbounded potential at boundary
            gradphi_F_bdy = gradphiF_dip(C_bdy, srcpnts, dipmom, sig_dip);
            rhs_C_bdy = -sum(C_bdy_nv.*gradphi_F_bdy,2);
            
            % Deal with the coupling region
            % First we consider the Dirichlet coupling condition
            % This requires values from B to be mapped to points on A
            DM_A_cpl_out_dir = DistanceMatrix(A_cpl_out,A_ctrs);
            CCM_A_cpl_out_dir = rbf(ep,DM_A_cpl_out_dir);
            
            DM_B_cpl_in_dir = DistanceMatrix(A_cpl_out,B_ctrs);
            CCM_B_cpl_in_dir = rbf(ep,DM_B_cpl_in_dir);
            
            DM_B_cpl_out_dir = DistanceMatrix(B_cpl_out,B_ctrs);
            CCM_B_cpl_out_dir = rbf(ep,DM_B_cpl_out_dir);
            
            DM_C_cpl_in_dir = DistanceMatrix(B_cpl_out,C_ctrs);
            CCM_C_cpl_in_dir = rbf(ep,DM_C_cpl_in_dir);
            
            rhs_A_cpl_out_dir = zeros(N_A_cpl_out,1);
            rhs_B_cpl_out_dir = zeros(N_B_cpl_out,1);
            
            % Now we consider the Neumann coupling condition
            % Here, normal derivatives need to be mapped from A to points on B
            % The RHS is confusing, need to ask Sal for comments
            DM_A_cpl_out_neu = DistanceMatrix(B_cpl_in,A_ctrs);
            A_cpl_out_dx_neu = DifferenceMatrix(B_cpl_in(:,1),A_ctrs(:,1));
            A_cpl_out_dy_neu = DifferenceMatrix(B_cpl_in(:,2),A_ctrs(:,2));
            A_cpl_out_dz_neu = DifferenceMatrix(B_cpl_in(:,3),A_ctrs(:,3));
            D1 = repmat(B_cpl_in_nv(:,1),1,N_A).*dxrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dx_neu);
            D2 = repmat(B_cpl_in_nv(:,2),1,N_A).*dyrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dy_neu);
            D3 = repmat(B_cpl_in_nv(:,3),1,N_A).*dzrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dz_neu);
            CCM_A_cpl_out_neu = sig(1)*(D1 + D2 + D3);
            
            DM_B_cpl_in_neu = DistanceMatrix(B_cpl_in,B_ctrs);
            B_cpl_in_dx_neu = DifferenceMatrix(B_cpl_in(:,1),B_ctrs(:,1));
            B_cpl_in_dy_neu = DifferenceMatrix(B_cpl_in(:,2),B_ctrs(:,2));
            B_cpl_in_dz_neu = DifferenceMatrix(B_cpl_in(:,3),B_ctrs(:,3));
            D1 = repmat(B_cpl_in_nv(:,1),1,N_B).*dxrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dx_neu);
            D2 = repmat(B_cpl_in_nv(:,2),1,N_B).*dyrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dy_neu);
            D3 = repmat(B_cpl_in_nv(:,3),1,N_B).*dzrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dz_neu);
            CCM_B_cpl_in_neu = sig(2)*(D1 + D2 + D3);
            
            DM_B_cpl_out_neu = DistanceMatrix(C_cpl_in,B_ctrs);
            B_cpl_out_dx_neu = DifferenceMatrix(C_cpl_in(:,1),B_ctrs(:,1));
            B_cpl_out_dy_neu = DifferenceMatrix(C_cpl_in(:,2),B_ctrs(:,2));
            B_cpl_out_dz_neu = DifferenceMatrix(C_cpl_in(:,3),B_ctrs(:,3));
            D1 = repmat(C_cpl_in_nv(:,1),1,N_B).*dxrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dx_neu);
            D2 = repmat(C_cpl_in_nv(:,2),1,N_B).*dyrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dy_neu);
            D3 = repmat(C_cpl_in_nv(:,3),1,N_B).*dzrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dz_neu);
            CCM_B_cpl_out_neu = sig(2)*(D1 + D2 + D3);
            
            DM_C_cpl_in_neu = DistanceMatrix(C_cpl_in,C_ctrs);
            C_cpl_in_dx_neu = DifferenceMatrix(C_cpl_in(:,1),C_ctrs(:,1));
            C_cpl_in_dy_neu = DifferenceMatrix(C_cpl_in(:,2),C_ctrs(:,2));
            C_cpl_in_dz_neu = DifferenceMatrix(C_cpl_in(:,3),C_ctrs(:,3));
            D1 = repmat(C_cpl_in_nv(:,1),1,N_C).*dxrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dx_neu);
            D2 = repmat(C_cpl_in_nv(:,2),1,N_C).*dyrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dy_neu);
            D3 = repmat(C_cpl_in_nv(:,3),1,N_C).*dzrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dz_neu);
            CCM_C_cpl_in_neu = sig(3)*(D1 + D2 + D3);
            
            % Gotta check to make sure these are correct for different sig values
            % in the different domains.  No issues if they are all equal though.
            gradphiF_B_cpl_in = gradphiF_dip(B_cpl_in,srcpnts,dipmom,sig_dip);
            rhs_B_cpl_in_neu = -(sig(1)-sig(2))*sum(B_cpl_in_nv.*gradphiF_B_cpl_in, 2);
            
            gradphiF_C_cpl_in = gradphiF_dip(C_cpl_in,srcpnts,dipmom,sig_dip);
            rhs_C_cpl_in_neu = -(sig(2)-sig(3))*sum(C_cpl_in_nv.*gradphiF_C_cpl_in, 2);
            
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
            % coefs = linsolve(CM,rhs);
            [U,S,V] = svd(CM,0);
            sing_val = diag(S);
            inv_sing_val_TSVD = 1./sing_val;
            
            indices = find(sing_val < 1e-5*max(sing_val));
            inv_sing_val_TSVD(indices) = 0;
            % inv_sing_val_TSVD(floor(0.2*length(inv_sing_val)):end) = 0;
            cbeta = U'*rhs;
            coefs = V*(inv_sing_val_TSVD.*cbeta(1:length(sing_val)));
            
            % Extract the oefficients for evaluation
            C_coefs = coefs(N_A+N_B+1:end);
            % Potential at evalpnts in the source free case
            phi0 = EM * C_coefs;
            % Potential at evalpnts (superposition of effects)
            phi = phi0 + phi_F;
            phi_comp_MFS(:,k,i) = phi - phi(1);
            
            elapsed_t_MFS(k,i) = tt + toc;
            
            % Compute the total errors
            errvec_MFS(k,i) = errcompute(phi_comp_MFS(:,k,i),phi_true);
            
            % Store the condition of the system
            % [U,S,~] = svd(CM);
            % sing_val = diag(S);
            % cbeta = U'*rhs;
            condvec_MFS(k,i) = norm(rhs)/min(sing_val)/...
                norm(cbeta(1:length(sing_val))./sing_val);
            condstr = 'Effective condition number';
            
            Nvec_true_MFS(k) = N_A_cpl_out + N_B_cpl_in + ...
                N_B_cpl_out + N_C_cpl_in + N_C_bdy;
            
            if iter_out
                fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec_MFS(k,i),condvec_MFS(k,i),Nvec_true_MFS(k));
            end
            
        end
    end
end
%% Kansa solutions
if kansa
    % RBF definition and derivatives
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);
    phi_comp_kansa = zeros(N_eval,length(Nvec));
    for k=1:length(Nvec)
        
        if iter_out
            fprintf('k=%d\n',k)
        end
        
        % Choose all points in the geometry
        [POINTS, NORMALS] = BallGeometry(R,sum(N_elements(k,:)),'kansa',int_point_dist,bdy_point_dist);
        
        % Cut up the domain into the appropriate sections
        % Section A is the inner ball, B is the layer in the middle, C in the
        % outer layer
        A_int = POINTS.int1;
        A_cpl_out = POINTS.bdy11;
        B_cpl_in = POINTS.bdy12;
        B_int = POINTS.int2;
        B_cpl_out = POINTS.bdy22;
        C_cpl_in = POINTS.bdy23;
        C_int = POINTS.int3;
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
        N_A_int = size(A_int,1);
        N_A_cpl_out = size(A_cpl_out,1);
        N_B_cpl_in = size(B_cpl_in,1);
        N_B_int = size(B_int,1);
        N_B_cpl_out = size(B_cpl_out,1);
        N_C_cpl_in = size(C_cpl_in,1);
        N_C_int = size(C_int,1);
        N_C_bdy = size(C_bdy,1);
        
        % Compose a vector of all the RBF centers
        A_ctrs = [A_int;A_cpl_out];
        B_ctrs = [B_cpl_in;B_int;B_cpl_out];
        C_ctrs = [C_cpl_in;C_int;C_bdy];
        
        N_A = size(A_ctrs,1);
        N_B = size(B_ctrs,1);
        N_C = size(C_ctrs,1);
        
        tic;
        
        % Compute the collocation block for the interior portions
        DM_A_int = DistanceMatrix(A_int,A_ctrs);
        LCM_A_int = Lrbf(ep,DM_A_int);
        rhs_A_int = zeros(N_A_int,1);
        
        DM_B_int = DistanceMatrix(B_int,B_ctrs);
        LCM_B_int = Lrbf(ep,DM_B_int);
        rhs_B_int = zeros(N_B_int,1);
        
        DM_C_int = DistanceMatrix(C_int,C_ctrs);
        LCM_C_int = Lrbf(ep,DM_C_int);
        rhs_C_int = zeros(N_C_int,1);
        
        % Compute the evaluation matrix
        DM_eval = DistanceMatrix(evalpnts, C_ctrs);
        EM = rbf(ep, DM_eval);
        
        % Compute the collocation block for the boundary conditions
        % This also computes the RHS for the problem
        DM_C_bdy = DistanceMatrix(C_bdy,C_ctrs);
        
        % Find all the necessary difference matrices
        C_bdy_dx = DifferenceMatrix(C_bdy(:,1),C_ctrs(:,1));
        C_bdy_dy = DifferenceMatrix(C_bdy(:,2),C_ctrs(:,2));
        C_bdy_dz = DifferenceMatrix(C_bdy(:,3),C_ctrs(:,3));
        
        % Compute normal derivative collocation matrix for boundary
        D1 = repmat(C_bdy_nv(:,1),1,N_C).*dxrbf(ep,DM_C_bdy,C_bdy_dx);
        D2 = repmat(C_bdy_nv(:,2),1,N_C).*dyrbf(ep,DM_C_bdy,C_bdy_dy);
        D3 = repmat(C_bdy_nv(:,3),1,N_C).*dzrbf(ep,DM_C_bdy,C_bdy_dz);
        BCM_C_bdy = D1 + D2 + D3;
        
        % Compute known-terms vector (a.k.a. righthand side vector)
        % This requires the gradient of the unbounded potential at boundary
        gradphi_F_bdy = gradphiF_dip(C_bdy, srcpnts, dipmom, sig_dip);
        rhs_C_bdy = -sum(C_bdy_nv.*gradphi_F_bdy,2);
        
        % Deal with the coupling region
        % First we consider the Dirichlet coupling condition
        % This requires values from B to be mapped to points on A
        DM_A_cpl_out_dir = DistanceMatrix(A_cpl_out,A_ctrs);
        CCM_A_cpl_out_dir = rbf(ep,DM_A_cpl_out_dir);
        
        DM_B_cpl_in_dir = DistanceMatrix(A_cpl_out,B_ctrs);
        CCM_B_cpl_in_dir = rbf(ep,DM_B_cpl_in_dir);
        
        DM_B_cpl_out_dir = DistanceMatrix(B_cpl_out,B_ctrs);
        CCM_B_cpl_out_dir = rbf(ep,DM_B_cpl_out_dir);
        
        DM_C_cpl_in_dir = DistanceMatrix(B_cpl_out,C_ctrs);
        CCM_C_cpl_in_dir = rbf(ep,DM_C_cpl_in_dir);
        
        rhs_A_cpl_out_dir = zeros(N_A_cpl_out,1);
        rhs_B_cpl_out_dir = zeros(N_B_cpl_out,1);
        
        % Now we consider the Neumann coupling condition
        % Here, normal derivatives need to be mapped from A to points on B
        % The RHS is confusing, need to ask Sal for comments
        DM_A_cpl_out_neu = DistanceMatrix(B_cpl_in,A_ctrs);
        A_cpl_out_dx_neu = DifferenceMatrix(B_cpl_in(:,1),A_ctrs(:,1));
        A_cpl_out_dy_neu = DifferenceMatrix(B_cpl_in(:,2),A_ctrs(:,2));
        A_cpl_out_dz_neu = DifferenceMatrix(B_cpl_in(:,3),A_ctrs(:,3));
        D1 = repmat(B_cpl_in_nv(:,1),1,N_A).*dxrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dx_neu);
        D2 = repmat(B_cpl_in_nv(:,2),1,N_A).*dyrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dy_neu);
        D3 = repmat(B_cpl_in_nv(:,3),1,N_A).*dzrbf(ep,DM_A_cpl_out_neu,A_cpl_out_dz_neu);
        CCM_A_cpl_out_neu = sig(1)*(D1 + D2 + D3);
        
        DM_B_cpl_in_neu = DistanceMatrix(B_cpl_in,B_ctrs);
        B_cpl_in_dx_neu = DifferenceMatrix(B_cpl_in(:,1),B_ctrs(:,1));
        B_cpl_in_dy_neu = DifferenceMatrix(B_cpl_in(:,2),B_ctrs(:,2));
        B_cpl_in_dz_neu = DifferenceMatrix(B_cpl_in(:,3),B_ctrs(:,3));
        D1 = repmat(B_cpl_in_nv(:,1),1,N_B).*dxrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dx_neu);
        D2 = repmat(B_cpl_in_nv(:,2),1,N_B).*dyrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dy_neu);
        D3 = repmat(B_cpl_in_nv(:,3),1,N_B).*dzrbf(ep,DM_B_cpl_in_neu,B_cpl_in_dz_neu);
        CCM_B_cpl_in_neu = sig(2)*(D1 + D2 + D3);
        
        DM_B_cpl_out_neu = DistanceMatrix(C_cpl_in,B_ctrs);
        B_cpl_out_dx_neu = DifferenceMatrix(C_cpl_in(:,1),B_ctrs(:,1));
        B_cpl_out_dy_neu = DifferenceMatrix(C_cpl_in(:,2),B_ctrs(:,2));
        B_cpl_out_dz_neu = DifferenceMatrix(C_cpl_in(:,3),B_ctrs(:,3));
        D1 = repmat(C_cpl_in_nv(:,1),1,N_B).*dxrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dx_neu);
        D2 = repmat(C_cpl_in_nv(:,2),1,N_B).*dyrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dy_neu);
        D3 = repmat(C_cpl_in_nv(:,3),1,N_B).*dzrbf(ep,DM_B_cpl_out_neu,B_cpl_out_dz_neu);
        CCM_B_cpl_out_neu = sig(2)*(D1 + D2 + D3);
        
        DM_C_cpl_in_neu = DistanceMatrix(C_cpl_in,C_ctrs);
        C_cpl_in_dx_neu = DifferenceMatrix(C_cpl_in(:,1),C_ctrs(:,1));
        C_cpl_in_dy_neu = DifferenceMatrix(C_cpl_in(:,2),C_ctrs(:,2));
        C_cpl_in_dz_neu = DifferenceMatrix(C_cpl_in(:,3),C_ctrs(:,3));
        D1 = repmat(C_cpl_in_nv(:,1),1,N_C).*dxrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dx_neu);
        D2 = repmat(C_cpl_in_nv(:,2),1,N_C).*dyrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dy_neu);
        D3 = repmat(C_cpl_in_nv(:,3),1,N_C).*dzrbf(ep,DM_C_cpl_in_neu,C_cpl_in_dz_neu);
        CCM_C_cpl_in_neu = sig(3)*(D1 + D2 + D3);
        
        % Gotta check to make sure these are correct for different sig values
        % in the different domains.  No issues if they are all equal though.
        gradphiF_B_cpl_in = gradphiF_dip(B_cpl_in,srcpnts,dipmom,sig_dip);
        rhs_B_cpl_in_neu = -(sig(1)-sig(2))*sum(B_cpl_in_nv.*gradphiF_B_cpl_in, 2);
        
        gradphiF_C_cpl_in = gradphiF_dip(C_cpl_in,srcpnts,dipmom,sig_dip);
        rhs_C_cpl_in_neu = -(sig(2)-sig(3))*sum(C_cpl_in_nv.*gradphiF_C_cpl_in, 2);
        
        % Create the full linear system from the blocks
        % Compose rhs
        rhs = [rhs_A_int; ...
            rhs_A_cpl_out_dir; ...
            rhs_B_cpl_in_neu; ...
            rhs_B_int; ...
            rhs_B_cpl_out_dir; ...
            rhs_C_cpl_in_neu; ...
            rhs_C_int; ...
            rhs_C_bdy];
        % Compose collocation matrix in same order as rhs
        % Notice the coupling conditions are of the form A-B
        % That's why the (-) appears before the B coupling components
        CM = [LCM_A_int,zeros(N_A_int,N_B+N_C); ...
            CCM_A_cpl_out_dir,-CCM_B_cpl_in_dir,zeros(N_A_cpl_out,N_C); ...
            CCM_A_cpl_out_neu,-CCM_B_cpl_in_neu,zeros(N_B_cpl_in,N_C); ...
            zeros(N_B_int,N_A),LCM_B_int,zeros(N_B_int,N_C); ...
            zeros(N_B_cpl_out,N_A),CCM_B_cpl_out_dir,-CCM_C_cpl_in_dir; ...
            zeros(N_C_cpl_in,N_A),CCM_B_cpl_out_neu,-CCM_C_cpl_in_neu; ...
            zeros(N_C_int,N_A+N_B),LCM_C_int; ...
            zeros(N_C_bdy,N_A+N_B),BCM_C_bdy];
        
        % Solve the BVP linear system
        %     coefs = linsolve(CM,rhs);
        [U,S,V] = svd(CM,0);
        sing_val = diag(S);
        inv_sing_val_TSVD = 1./sing_val;
        
        indices = find(sing_val < 1e-10);
        inv_sing_val_TSVD(indices) = 0;
        % inv_sing_val_TSVD(floor(0.2*length(inv_sing_val)):end) = 0;
        cbeta = U'*rhs;
        coefs = V*(inv_sing_val_TSVD.*cbeta(1:length(sing_val)));
        
        % Extract the oefficients for evaluation
        C_coefs = coefs(N_A+N_B+1:end);
        % Potential at evalpnts in the source free case
        phi0 = EM * C_coefs;
        % Potential at evalpnts (superposition of effects)
        phi = phi0 + phi_F;
        phi_comp_kansa(:,k) = phi - phi(1);
        
        elapsed_t_kansa(k) = tt + toc;
        
        % Compute the total errors
        errvec_kansa(k) = errcompute(phi_comp_kansa(:,k),phi_true);
        
        % Store the condition of the system
        [U,S,~] = svd(CM);
        sing_val = diag(S);
        beta = U'*rhs;
        condvec_kansa(k) = norm(rhs)/min(sing_val)/...
            norm(beta(1:length(sing_val))./sing_val);
        condstr = 'Effective condition number';
        
        Nvec_true_kansa(k) = N_A_cpl_out + N_B_cpl_in + ...
            N_B_cpl_out + N_C_cpl_in + N_C_bdy;
        
        if iter_out
            fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec_kansa(k),condvec_kansa(k),Nvec_true_kansa(k));
        end
    end
end

%% Save results in a mat file
save(resultsfilename,'R','sig','evalpnts','dipmom','srcpnts',...
     'phi_comp_*','phi_true');
if BEM
    save(resultsfilename,'*_BEM','-append');
end
if MFS
    save(resultsfilename,'mfs_sphere','mfs_frac','*_MFS','-append');
end
if kansa
    save(resultsfilename,'ep','*_kansa','-append');
end

%% Plot simulation results
% Error vs No. of collocation points or elements
h1 = figure;
axes1 = axes('Parent',h1,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if BEM
    convtest_BEM = semilogy(Nvec_true_BEM, errvec_BEM,...
                            'Parent',axes1,...
                            'MarkerSize',8,...
                            'LineWidth',3);
    set(convtest_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    convtest_kansa = semilogy(Nvec_true_kansa, errvec_kansa,...
                              'Parent',axes1,...
                              'MarkerSize',8,...
                              'LineWidth',3);
    set(convtest_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    convtest_MFS = semilogy(Nvec_true_MFS,errvec_MFS,...
                            'Parent',axes1,...
                            'MarkerSize',8,...
                            'LineWidth',3);
    set(convtest_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(convtest_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(convtest_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(convtest_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Three Layered Sphere - Convergence test',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points or Elements',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('Relative Error on Surface Potential',...
       'FontWeight','bold','FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthWest',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% CPU time vs No. of collocation points or elements
h2 = figure;
axes1 = axes('Parent',h2,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if BEM
    time_BEM = semilogy(Nvec_true_BEM, elapsed_t_BEM,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);...
    set(time_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    time_kansa = semilogy(Nvec_true_kansa, elapsed_t_kansa,...
                          'Parent',axes1,...
                          'MarkerSize',8,...
                          'LineWidth',3);
    set(time_kansa,...
                          'Marker','diamond',...
                          'Color',[0 0 1],...
                          'DisplayName','Kansa (IMQ)');
end
if MFS
    time_MFS = semilogy(Nvec_true_MFS,elapsed_t_MFS,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);
    set(time_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(time_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(time_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(time_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Three Layered Sphere - CPU time',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points or Elements',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('CPU Time [s]',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthEast',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% Effective condition number vs No. of collocation points or elements
h3 = figure;
axes1 = axes('Parent',h3,'FontSize',16,'YScale','log','XLim',[100 6000]);
hold(axes1,'all');
if kansa
    cond_kansa = semilogy(Nvec_true_kansa, condvec_kansa,...
                          'Parent',axes1,...
                          'MarkerSize',8,...
                          'LineWidth',3);
    set(cond_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    cond_MFS = semilogy(Nvec_true_MFS,condvec_MFS,...
                        'Parent',axes1,...
                        'MarkerSize',8,...
                        'LineWidth',3);
    set(cond_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(cond_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(cond_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(cond_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Three Layered Sphere - System Condition',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Total Collocation Points',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('Effective Condition Number',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthEast',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% CPU time vs Error
h4 = figure;
axes1 = axes('Parent',h4,'FontSize',16,'XScale','log','YScale','log');
hold(axes1,'all');
if BEM
    timeerr_BEM = semilogy(errvec_BEM, elapsed_t_BEM,...
                           'Parent',axes1,...
                           'MarkerSize',8,...
                           'LineWidth',3);
    set(timeerr_BEM,...
        'Marker','hexagram',...
        'Color',[0 0 0],...
        'DisplayName','Symmetric BEM');
end
if kansa
    timeerr_kansa = semilogy(errvec_kansa, elapsed_t_kansa,...
                             'Parent',axes1,...
                             'MarkerSize',8,...
                             'LineWidth',3);
    set(timeerr_kansa,...
        'Marker','diamond',...
        'Color',[0 0 1],...
        'DisplayName','Kansa (IMQ)');
end
if MFS
    timeerr_MFS = semilogy(errvec_MFS,elapsed_t_MFS,...
                           'Parent',axes1,...
                           'MarkerSize',8,...
                           'LineWidth',3);
    set(timeerr_MFS(1),...
        'Marker','square',...
        'Color',[1 0 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(1))));
    set(timeerr_MFS(2),...
        'Marker','o',...
        'Color',[0 1 1],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(2))));
    set(timeerr_MFS(3),...
        'Marker','*',...
        'Color',[1 0 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(3))));
    set(timeerr_MFS(4),...
        'Marker','^',...
        'Color',[0 1 0],...
        'DisplayName',strcat('MFS, R_{MFS} = ',num2str(mfs_frac(4))));
end
title('Current Dipole in Three Layered Sphere - Cost per Accuracy',...
      'FontWeight','bold',...
      'FontSize',20)
xlabel('Relative Error on Surface Potential',...
       'FontWeight','bold',...
       'FontSize',16);
ylabel('CPU Time [s]',...
       'FontWeight','bold',...
       'FontSize',16);
legend1 = legend(axes1,'show');
set(legend1,...
    'Location','SouthWest',...
    'FontWeight','bold',...
    'FontSize',14);
grid on

% I requested, plot also the analitycal solution, the computed solution and 
% the error distribution for a specific No. of collocation points or 
% elements
if plot_sol
    switch sol_err_style
        case 1
            if BEM
                sol_err_BEM = abs(phi_true - phi_comp_BEM(:,N_plotsol));
                plotstr_BEM = 'BEM Absolute Error [V]';
            end
            if MFS
                sol_err_MFS = abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol));
                plotstr_MFS = 'MFS Absolute Error [V]';
            end
            if kansa
                sol_err_kansa = abs(phi_true - phi_comp_kansa(:,N_plotsol));
                plotstr_kansa = 'KM Absolute Error [V]';
            end
        case 2
            if BEM
                sol_err_BEM = log10(abs(phi_true - phi_comp_BEM(:,N_plotsol)));
                plotstr_BEM = 'BEM Log10 of Absolute Error';
            end
            if MFS
                sol_err_MFS = log10(abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol)));
                plotstr_MFS = 'MFS Log10 of Absolute Error';
            end
            if kansa
                sol_err_kansa = log10(abs(phi_true - phi_comp_kansa(:,N_plotsol)));
                plotstr_kansa = 'KM Log10 of Absolute Error';
            end
        case 3
            if BEM
                sol_err_BEM = log10(abs(phi_true - phi_comp_BEM(:,N_plotsol))./...
                                    (abs(phi_true)+eps)+eps);
                plotstr_BEM = 'BEM Log10 of Pointwise Relative Error';
            end
            if MFS
                sol_err_MFS = log10(abs(phi_true - phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol))./...
                                    (abs(phi_true)+eps)+eps);
                plotstr_MFS = 'MFS Log10 of Pointwise Relative Error';
            end
            if kansa
                sol_err_kansa = log10(abs(phi_true - phi_comp_kansa(:,N_plotsol))./...
                                      (abs(phi_true)+eps)+eps);
                plotstr_kansa = 'KM Log10 of Pointwise Relative Error';
            end
        otherwise
            error('Unknown 3D plot error style %g',sol_err_style)
    end
    
    h5 = figure; % Analytic
    SurfacePlot_dip(evalpnts, phi_true)
    title('Analytic Potential [V]','FontWeight','bold','FontSize',16)
    if BEM
        h6 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_BEM(:,N_plotsol));
        title('BEM Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_BEM);
        title(plotstr_BEM,'FontWeight','bold','FontSize',16)
    end
    if MFS
        h7 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_MFS(:,N_plotsol,mfs_frac_plotsol))
        title('MFS Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_MFS);
        title(plotstr_MFS,'FontWeight','bold','FontSize',16)
    end
    if kansa
        h8 = figure;
        subplot(1,2,1)
        SurfacePlot_dip(evalpnts, phi_comp_kansa(:,N_plotsol));
        title('Kansa Computed Potential [V]','FontWeight','bold','FontSize',16)
        subplot(1,2,2)
        SurfacePlot_dip(evalpnts, sol_err_kansa)
        title(plotstr_kansa,'FontWeight','bold','FontSize',16)
    end
end