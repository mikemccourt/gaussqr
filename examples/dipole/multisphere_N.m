% multisphere_N.m
%
%  For this problem, we consider the solution to the Laplace equation on a
%  sphere with Neumann boundary conditions.  We will consider 2 spheres
%  (eventually 3), one inside the other, with 2nd order (values and fluxes)
%  coupling between the domains.
%
%  The coupling will be carried out with mismatched points (different
%  coupling locations for the inner and outer spheres) in general.  We will
%  include an option to fix them to be the same points, perhaps.
%
%  The problem has several physical parameters relating to the
%  underlying EEG/MEG physical system.  These parameters are:
%    R - Sphere radius [dm] <default = 1>
%    sig - Electric conductivity [S/dm] <default = [.02,.02]>
%    dipmom - Dipole moment [x10^-12 Am] <default = 2.7*[1,0,0]>
%    srcpnts - Dipole position [dm] <default = [0,0,0.6*R]>
%
%  This script allows you to test the convergence rate with respect to N
%  of different RBFs and different epsilon values.
%
%  The solution parameters to be considered are
%     sol_type - How you want to solve the system <default = 'kansa'>
%                'kansa' : Nonsymmetric collocation
%                          rbf_choice and ep must also be chosen
%                'mfs' : Method of fundamental solutions
%                        MFS_frac and ctr_sphere must also be chosen
%     rbf_choice - RBF for collocation <default = 'imq'>
%     ep - RBF shape parameter <default = 10>
%     mfs_frac - How many centers for MFS, in (0.0->1.0)*N <default = 1.0>
%     mfs_sphere - Fraction beyond R (eg, 1.3R) for centers <default = 1.4>
%     BC_choice - How to choose the boundary conditions <default = 1>
%                 1 : Neumann
%                 2 : Dirichlet
%                 3 : Mixed (6 random Dirichlet, the rest Neumann)
%                 4 : Neumann + 1 Dirichlet at reference, set to 0
%     eval_diff - Consider the solution as the difference between all
%                 values and a reference point <default = 1>
%     reference - The point chosen to evaluate the difference at, or set to
%                 0 for BC_choice=4 <default = [R(end),0,0]>
%                 Use reference=[] for no additional reference point
%
%  These are parameters directly related to the coupling between the
%  concentric balls
%     match_couple - Use the same points on the interface to couple the
%                    values and the fluxes <default = 0>
%     sol_acc - How long a series should be used in computing the true
%               solution on the boundary <default = 100>
%
%  The value we are interested in studying is the effect of increasing N,
%  so you must specify a vector of N values that you want to study
%     Nvec - Row vector of N values <default = 100:50:500>
%     BC_frac - The fraction of the total points to be used to enforce
%               boundary conditions <default = .3>
%     dip_cushion - How much space should be given around the dipole where
%               no RBF centers are allowed <default = .005>
%     N_eval - # points to evaluate error <default = 1001>
%
%  Some outputs are available if you would like them
%     iter_out - Print output during the solves <default = 0>
%     plot_sol - 3D surface plot of boundary solution <default = 0>
%     sol_err_style - How do you want the 3D solution error displayed
%                     0 : No error computed, just the solution
%                     1 : Absolute error <default>
%                     2 : Log absolute error
%                     3 : Log pointwise relative error
%     plot_err - log-log plot of error vs. N <default = 1>
%     errcolor - Color for error line in log-log plot <default = 'b'>
%     condcolor - Color for condition line in log-log plot <default = 'r'>
%
%  The results of these experiments are stored in
%     Nvec_true - Actual number of collocation points, because the point
%                 distribution in a sphere is tricky
%     errvec - Errors computed at Nvec_true
%     condvec - Collocation matrix condition numbers at Nvec_true
%
%  Note that if MFS with fewer centers than collocation points is chosen,
%  the value in condvec is the rank of the matrix, not the condition

R = [0.7, 1];
sig = [0.02, 0.02];
dipmom = 2.7*[1, 0, 0];
srcpnts = [0, 0, 0.6*R(1)];

sol_type = 'mfs';
radbasfun = 'imq';
ep = 1;
mfs_frac = .4;
mfs_sphere = 1.4;
BC_choice = 1;
eval_diff = 1;
reference = [R(end),0,0];
%reference = [-0.5955 0.0699 0.8003]; % Minimum potential

match_couple = 1;
sol_acc = 100;

Nvec = 100:100:2000;
BC_frac = .3; % Not yet implemented
dip_cushion = .01;
N_eval = 1001;

iter_out = 1;
plot_sol = 1;
sol_err_style = 1;
plot_err = 1;
errcolor = 'b';
condcolor = 'r';


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

% RBF definition and derivatives
if strcmp(sol_type,'kansa')
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);
else
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF('fundamental_3d');
end

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval, R(end));
evalpnts = [reference;evalpnts];

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig_dip);

% Analytic solution for the potential
phi_an = MultiSpherePotential(R, sig, srcpnts, dipmom, evalpnts, sol_acc);

% If requested, compute the difference of the solution with a reference
% point, which was attached at the top of evalpnts earlier
if eval_diff
    phi_true = phi_an - phi_an(1);
else
    phi_true = phi_an;
end

% Loop through the requested N values
errvec = [];
condvec = [];
Nvec_true = [];
k = 1;
for Npnts = Nvec
    if iter_out
        fprintf('k=%d\n',k)
    end
    
    % Choose all points in the geometry
    [POINTS, NORMALS] = BallGeometry(R,Npnts,sol_type,[],srcpnts,dip_cushion);
    
    % Cut up the domain into the appropriate sections
    % Section A is the inner ball, B is the outer ball
    A_int = POINTS.int1;
    A_cpl_out = POINTS.bdy11;
    B_cpl_in = POINTS.bdy12;
    B_int = POINTS.int2;
    B_bdy = POINTS.bdy22;
    
    % These are the normal vectors to the interface and boundary
    % Note that the interface normals both point toward the boundary, even
    % though that is technically the 'negative' normal for the outer ball,
    % since it is pointing in to rather than out of the outer domain
    A_cpl_out_nv = NORMALS.n11;
    B_cpl_in_nv = NORMALS.n12;
    B_bdy_nv = NORMALS.n22;
    
    % If the user wants, we can require the values and derivatives to be
    % matched at the same points rather than different points.  I don't
    % think it matters, but it's nice to have this option.
    if match_couple
        B_cpl_in = A_cpl_out;
        B_cpl_in_nv = A_cpl_out_nv;
    end
    
    % Classify all the sizes of each domain
    N_A_int = size(A_int,1);
    N_A_cpl_out = size(A_cpl_out,1);
    N_B_cpl_in = size(B_cpl_in,1);
    N_B_int = size(B_int,1);
    N_B_bdy = size(B_bdy,1);
    
    % Compose a vector of all the RBF centers
    % In the MFS setting, these are chosen in a sphere around the ball
    % The MFS center choices may not be the best, but it will work for now
    if strcmp(sol_type,'mfs')
        N_ctrs = floor(mfs_frac*(N_A_cpl_out+N_B_cpl_in+N_B_bdy));
        all_ctrs = SphereSurfGoldPoints(N_ctrs, mfs_sphere*R(end));
        A_ctrs_ind = randperm(N_ctrs,floor(mfs_frac*N_A_cpl_out));
        B_ctrs_ind = setdiff(1:N_ctrs,A_ctrs_ind);
        A_ctrs = all_ctrs(A_ctrs_ind,:);
        B_ctrs = all_ctrs(B_ctrs_ind,:);
    else % For kansa, the centers and collocation points coincide
        A_ctrs = [A_int;A_cpl_out];
        if BC_choice~=4
            B_ctrs = [B_cpl_in;B_int;B_bdy];
        else
            B_ctrs = [B_cpl_in;B_int;B_bdy;reference];
        end
    end
    N_A = size(A_ctrs,1);
    N_B = size(B_ctrs,1);
    
    
    % Compute the collocation block for the interior portions
    DM_A_int = DistanceMatrix(A_int,A_ctrs);
    LCM_A_int = Lrbf(ep,DM_A_int);
    rhs_A_int = zeros(N_A_int,1);
    
    DM_B_int = DistanceMatrix(B_int,B_ctrs);
    LCM_B_int = Lrbf(ep,DM_B_int);
    rhs_B_int = zeros(N_B_int,1);
    
    
    % Compute the evaluation matrix
    DM_eval = DistanceMatrix(evalpnts, B_ctrs);
    EM = rbf(ep, DM_eval);
    
    
    % Determine which points are Neumann and which are Dirichlet
    % This only applies to the boundary, not the interface
    %   Notice the use of zeros(0,3), not []
    %   To allow for bdydata_neu(:,1) calls later
    switch BC_choice
        case 1 % Do the standard Neumann BC
            B_bdy_neu = B_bdy;
            B_bdy_neu_nv = B_bdy_nv;
            B_bdy_dir = zeros(0,3);
        case 2 % Run a test with Dirichlet BC
            B_bdy_neu = zeros(0,3);
            B_bdy_neu_nv = zeros(0,3);
            B_bdy_dir = B_bdy;
        case 3 % Run a test with Mixed BC
            % Right now, fixed at 6 Dirichlet BC points
            % Could be variable, but not important
            N_B_dir = min(6,N_B_bdy);
            i_dir = randperm(N_B_bdy,N_B_dir);
            i_neu = setdiff(1:N_B_bdy,i_dir);
            
            B_bdy_neu = B_bdy(i_neu,:);
            B_bdy_neu_nv = B_bdy_nv(i_neu,:);
            B_bdy_dir = B_bdy(i_dir,:);
        case 4
            B_bdy_neu = B_bdy;
            B_bdy_neu_nv = B_bdy_nv;
            B_bdy_dir = reference;
            N_B_bdy = N_B_bdy + 1;
        otherwise
            error('What was this BC_choice %d that you passed?',BC_choice)
    end
    
    
    % Compute the collocation block for the boundary conditions
    % This also computes the RHS for the problem
    % First we consider the Neumann BC component
    DM_B_bdy_neu = DistanceMatrix(B_bdy_neu,B_ctrs);

    % Find all the necessary difference matrices
    B_bdy_dx_neu = DifferenceMatrix(B_bdy_neu(:,1),B_ctrs(:,1));
    B_bdy_dy_neu = DifferenceMatrix(B_bdy_neu(:,2),B_ctrs(:,2));
    B_bdy_dz_neu = DifferenceMatrix(B_bdy_neu(:,3),B_ctrs(:,3));
    
    % Compute normal derivative collocation matrix for boundary
    D1 = repmat(B_bdy_neu_nv(:,1),1,N_B).*dxrbf(ep,DM_B_bdy_neu,B_bdy_dx_neu);
    D2 = repmat(B_bdy_neu_nv(:,2),1,N_B).*dyrbf(ep,DM_B_bdy_neu,B_bdy_dy_neu);
    D3 = repmat(B_bdy_neu_nv(:,3),1,N_B).*dzrbf(ep,DM_B_bdy_neu,B_bdy_dz_neu);
    BCM_B_bdy_neu = D1 + D2 + D3;
    
    % Compute known-terms vector (a.k.a. righthand side vector)
    % This requires the gradient of the unbounded potential at boundary
    gradphi_F_bdy_neu = gradphiF_dip(B_bdy_neu, srcpnts, dipmom, sig_dip);
    rhs_B_bdy_neu = -sum(B_bdy_neu_nv.*gradphi_F_bdy_neu,2);
    
    % Now we consider the Dirichlet BC component
    DM_B_bdy_dir = DistanceMatrix(B_bdy_dir,B_ctrs);
    BCM_B_bdy_dir = rbf(ep,DM_B_bdy_dir);
    
    % Compute the true solution to be used as Dirichlet BC
    if BC_choice==4
        rhs_B_bdy_dir = 0;
    else
        phi_F_bdy_dir = phiF_dip(B_bdy_dir,srcpnts,dipmom,sig_dip);
        phi_bdy_dir = MultiSpherePotential(R, sig, srcpnts, dipmom, B_bdy_dir, sol_acc);
        rhs_B_bdy_dir = phi_bdy_dir - phi_F_bdy_dir;
    end
    
    % Combine the BC components
    rhs_B_bdy = [rhs_B_bdy_neu;rhs_B_bdy_dir];
    BCM_B_bdy = [BCM_B_bdy_neu;BCM_B_bdy_dir];
    
    
    % Deal with the coupling region
    % First we consider the Dirichlet coupling condition
    % This requires values from B to be mapped to points on A
    DM_A_cpl_out_dir = DistanceMatrix(A_cpl_out,A_ctrs);
    CCM_A_cpl_out_dir = rbf(ep,DM_A_cpl_out_dir);
    
    DM_B_cpl_in_dir = DistanceMatrix(A_cpl_out,B_ctrs);
    CCM_B_cpl_in_dir = rbf(ep,DM_B_cpl_in_dir);
    
    rhs_A_cpl_out_dir = zeros(N_A_cpl_out,1);
    
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
    
    % Gotta check to make sure these are correct for different sig values
    % in the different domains.  No issues if they are all equal though.
    gradphiF_B_cpl_in = gradphiF_dip(B_cpl_in,srcpnts,dipmom,sig_dip);
    rhs_B_cpl_in_neu = -(sig(1)-sig(2))*sum(B_cpl_in_nv.*gradphiF_B_cpl_in, 2);
    
    % Create the full linear system from the blocks
    % Compose rhs
    rhs = [rhs_A_int; ...
           rhs_A_cpl_out_dir; ...
           rhs_B_cpl_in_neu; ...
           rhs_B_int; ...
           rhs_B_bdy];
    % Compose collocation matrix in same order as rhs
    % Notice the coupling conditions are of the form A-B
    % That's why the (-) appears before the B coupling components
    CM = [LCM_A_int,zeros(N_A_int,N_B); ...
          CCM_A_cpl_out_dir,-CCM_B_cpl_in_dir; ...
          CCM_A_cpl_out_neu,-CCM_B_cpl_in_neu; ...
          zeros(N_B_int,N_A),LCM_B_int; ...
          zeros(N_B_bdy,N_A),BCM_B_bdy];
      
    % Solve the BVP linear system
    [coefs,recip_cond] = linsolve(CM,rhs);
    % Extract the oefficients for evaluation
    %%% There might be a smarter way to extract this
    B_coefs = coefs(N_A+1:end);
    % Potential at evalpnts in the source free case
    phi0 = EM * B_coefs;
    % Potential at evalpnts (superposition of effects)
    phi = phi0 + phi_F;

    % If requested, compute the difference of the solution with a
    % reference point, which was put at the top earlier
    if eval_diff
        phi_comp = phi - phi(1);
    else
        phi_comp = phi;
    end
    
    % Compute the total errors
    errvec(k) = errcompute(phi_comp,phi_true);
    % Store the condition of the system
    % For a low-rank system, instead store the rank
    % This may happen for some MFS problems
    if floor(recip_cond)==recip_cond
        condvec(k) = recip_cond;
    else
        condvec(k) = 1/recip_cond;
    end
    Nvec_true(k) = N_A + N_B;
    
    if iter_out
        fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec(k),condvec(k),Nvec_true(k));
    end
    
    k = k + 1;
end

if plot_err
    clf reset
    
    switch GAUSSQR_PARAMETERS.ERROR_STYLE
        case 1
            errstr = 'Pointwise Rel Err';
        case 2
            errstr = 'Absolute Error';
        case 3
            errstr = 'Relative Error';
        case 4
            errstr = 'RMS Relative Error';
    end
    switch BC_choice
        case 1
            bcstr = 'Neumann BC';
        case 2
            bcstr = 'Dirichlet BC';
        case 3
            bcstr = 'Mixed BC';
        case 4
            bcstr = 'Zero Reference BC';
    end
    epstr = sprintf(', \\epsilon=%g',ep);
    
    [AX,H1,H2] = plotyy(Nvec_true,errvec,Nvec_true,condvec,@loglog);
    xlabel('Total collocation points')
    set(AX(1),'Xlim',[Nvec_true(1),Nvec_true(end)])
    set(get(AX(1),'Ylabel'),'String',errstr)
    set(AX(1),'Ycolor',errcolor)
    set(AX(2),'Xlim',[Nvec_true(1),Nvec_true(end)])
    set(get(AX(2),'Ylabel'),'String','Matrix Condition','Color',condcolor)
    set(AX(2),'Ycolor',condcolor)
    set(H1,'LineWidth',3,'Color',errcolor)
    set(H2,'LineWidth',3,'Color',condcolor)
    title(strcat(bcstr,epstr))
end

if plot_sol
    figure
    switch sol_err_style
        case 0
            sol_err = phi_comp;
            plotstr = 'Computed Solution';
        case 1
            sol_err = abs(phi_true - phi_comp);
            plotstr = 'Absolute Error';
        case 2
            sol_err = log10(abs(phi_true - phi_comp));
            plotstr = 'Log10 of Absolute Error';
        case 3
            sol_err = log10(abs(phi_true - phi_comp)./(abs(phi_true)+eps)+eps);
            plotstr = 'Log10 of Pointwise Relative Error';
        otherwise
            error('Unknown 3D plot error style %g',sol_err_style)
    end
    
    subplot(1,2,1)
    SurfacePlot_dip(evalpnts, phi_true)
    title('Analytic potential','FontWeight','bold','FontSize',12)
    subplot(1,2,2)
    SurfacePlot_dip(evalpnts, sol_err);
    title(plotstr,'FontWeight','bold','FontSize',12)
end