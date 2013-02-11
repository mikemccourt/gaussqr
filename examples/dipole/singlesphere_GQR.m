% singlesphere_GQR.m
%
%  For this problem, we consider the solution to the Laplace equation on a
%  sphere with Neumann boundary conditions.
%
%  The problem has several physical parameters relating to the
%  underlying EEG/MEG physical system.  These parameters are:
%    R - Sphere radius [dm] <default = 1>
%    sig - Electric conductivity [S/dm] <default = .02>
%    dipmom - Dipole moment [x10^-12 Am] <default = 2.7e*[1,0,0]>
%    srcpnts - Dipole position [dm] <default = [0,0,0.6*R]>
%
%  This script allows you to test the convergence rate (with respect to N
%  of different RBFs and different epsilon values.
%
%  The solution parameters to be considered are
%     ep - Gaussian shape parameter <default = 1e-5>
%     alpha - GaussQR scale parameter <default = 1>
%     lowrank - Use the low rank GaussQR method <default = 1>
%     BC_choice - How to choose the boundary conditions <default = 1>
%                 1 : Neumann
%                 2 : Dirichlet
%                 3 : Mixed (6 random Dirichlet, the rest Neumann)
%     eval_diff - Consider the solution as the difference between all
%                 values and a reference point <default = 1>
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
%  condition doesn't make sense (rectangular, not square, system)

R = 1;
sig = 0.02;
dipmom = 2.7*[1, 0, 0];
srcpnts = [0, 0, 0.6*R];

ep = 1e-5;
alpha = 1;
lowrank = 1;
BC_choice = 1;
eval_diff = 1;

Nvec = 100:100:1500;
BC_frac = .3; % Not yet implemented
dip_cushion = .01;
N_eval = 1001;

iter_out = 1;
plot_sol = 1;
plot_err = 1;
errcolor = 'b';
condcolor = 'r';


%%%%%%%%%%%%%%%%%%%%%
% Basic setup stuff for GaussQR

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = -200;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .4;

% Set random number generator to constant
% This is used in choosing which BC points are Dirichlet in the mixed case
rng(0);


%%%%%%%%%%%%%%%%%%%%%
% This is the start of the solver

% Determine the evaluation points (all on the boundary)
evalpnts = SphereSurfGoldPoints(N_eval, R);

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);

% Analytic solution for the potential
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

% If requested, compute the difference of the solution with a reference
% point, arbitrarily chosen as evalpnts(1)
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
for N_requested = Nvec
    if iter_out
        fprintf('k=%d\n',k)
    end
    
    % Determine collocation points
    [POINTS, NORMALS] = BallGeometry(R,N_requested,'kansa',[],srcpnts,dip_cushion);
    intdata = POINTS.int1;
    bdydata = POINTS.bdy11;
    N_int = size(intdata,1);
    N_bdy = size(bdydata,1);
    N = N_int + N_bdy;
    
    % Compose a vector of all the kernel centers
    ctrs = [intdata; bdydata];
    
    % Set up the GaussQR solver
    GQR = gqr_solveprep(lowrank,ctrs,ep,alpha);
    M = size(GQR.Marr,2);
    
    
    % Compute the Laplacian block and RHS for the interior
    Phi_int = gqr_phi(GQR,intdata,[2,0,0]) + ...
              gqr_phi(GQR,intdata,[0,2,0]) + ...
              gqr_phi(GQR,intdata,[0,0,2]);
    rhs_int = zeros(N_int,1);
    
    % Compute the evaluation matrix
    Phi_eval = gqr_phi(GQR,evalpnts);
    
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
    else % Run a test with Mixed BC
        % Right now, fixed at 6 Dirichlet BC points
        % Could be variable, but not important
        N_dir = min(6,N_bdy);
        i_dir = randperm(N_bdy,N_dir);
        i_neu = setdiff(1:N_bdy,i_dir);
        
        bdydata_neu = bdydata(i_neu,:);
        normvecs = NORMALS.n11(i_neu,:);
        bdydata_dir = bdydata(i_dir,:);
    end
    
    % First we consider the Neumann BC component
    % We need one derivative in each direction
    Phi_neu_dx = gqr_phi(GQR,bdydata_neu,[1,0,0]);
    Phi_neu_dy = gqr_phi(GQR,bdydata_neu,[0,1,0]);
    Phi_neu_dz = gqr_phi(GQR,bdydata_neu,[0,0,1]);
    
    % Compute normal derivative collocation matrix for boundary
    A = repmat(normvecs(:,1),1,M).*Phi_neu_dx;
    B = repmat(normvecs(:,2),1,M).*Phi_neu_dy;
    C = repmat(normvecs(:,3),1,M).*Phi_neu_dz;
    Phi_neu = A + B + C;
    
    % Compute known-terms vector (a.k.a. righthand side vector)
    % This requires the gradient of the unbounded potential at boundary
    gradphi_F_neu = gradphiF_dip(bdydata_neu, srcpnts, dipmom, sig);
    rhs_bdy_neu = -sum(normvecs.*gradphi_F_neu,2);
    
    
    % Now we consider the Dirichlet BC component
    Phi_dir = gqr_phi(GQR,bdydata_dir);
    
    % Compute the true solution to be used as Dirichlet BC
    phi_F_bdy_dir = phiF_dip(bdydata_dir,srcpnts,dipmom,sig);
    phi_bdy_dir = HomSpherePotential(R, sig, srcpnts, dipmom, bdydata_dir);
    rhs_bdy_dir = phi_bdy_dir - phi_F_bdy_dir;
    
    
    % Create the full linear system from the blocksand solve it
    % Compose rhs
    rhs = [rhs_int;rhs_bdy_neu;rhs_bdy_dir];
    % Compose collocation matrix in same order as rhs
    Phi = [Phi_int;Phi_neu;Phi_dir];
    % Choose to solve with low rank or with full QR
    if lowrank
        CM = Phi;
        EM = Phi_eval;
    else
        CM = Phi(:,1:N) + Phi(:,N+1:end)*GQR.Rbar;
        EM = Phi_eval(:,1:N) + Phi_eval(:,N+1:end)*GQR.Rbar;
    end
    
    % Coefficients for evaluation
    [coefs,recip_cond] = linsolve(CM,rhs);
    
    % Potential at evalpnts in the source free case
    phi0 = EM * coefs;
    % Potential at evalpnts (superposition of effects)
    phi = phi0 + phi_F;

    % If requested, compute the difference of the solution with a\
    % reference point, arbitrarily chosen as evalpnts(1)
    if eval_diff
        phi_comp = phi - phi(1);
    else
        phi_comp = phi;
    end
    
    % Compute the total errors
    errvec(k) = errcompute(phi_comp,phi_true);
    % Store the condition of the system
    % For a low-rank system, instead store the rank
    if floor(recip_cond)==recip_cond
        condvec(k) = recip_cond;
    else
        condvec(k) = 1/recip_cond;
    end
    Nvec_true(k) = N;
    
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
            bcstr = sprintf('Mixed BC');
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
    
    subplot(1,2,1)
    SurfacePlot_dip(evalpnts, phi_true)
    title('Analytic potential')
    subplot(1,2,2)
    SurfacePlot_dip(evalpnts, abs(phi_true-phi_comp))
    title('Absolute error')
end