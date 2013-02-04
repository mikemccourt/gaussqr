% singlesphere_ep.m
%
%  For this problem, we consider the solution to the Laplace equation on a
%  sphere with Neumann boundary conditions.
%  This script allows you to test the error with respect to epsilon of
%  different RBFs.
%
%  The problem has several physical parameters relating to the
%  underlying EEG/MEG physical system.  These parameters are:
%    R - Sphere radius [dm] <default = 1>
%    sig - Electric conductivity [S/dm] <default = .02>
%    dipmom - Dipole moment [x10^-12 Am] <default = 2.7e*[1,0,0]>
%    srcpnts - Dipole position [dm] <default = [0,0,0.6*R]>
%
%  The solution parameters to be considered are
%     rbfset - set of RBFs for collocation to be tested
%               <default = { 'Gaussian' 'IMQ' 'MQ' ...
%                            'Wendland_C4' 'Wendland_C6' }>
%     BC_choice - How to choose the boundary conditions <default = 1>
%                 1 : Neumann
%                 2 : Dirichlet
%                 3 : Mixed (6 random Dirichlet, the rest Neumann)
%     eval_diff - Consider the solution as the difference between all
%                 values and a reference point <default = 1>
%
%  The value we are interested in studying is the effect of varying epsilon
%  so you must specify a vector of epsilon values that you want to study
%     epvec - Row vector of epsilon values <default = logspace(0,2,100)>
%
%     Npnts - desired number of interior points <default = 500>
%     BC_frac - The fraction of the total points to be used to enforce
%               boundary conditions <default = .3>
%     dip_cushion - How much space should be given around the dipole where
%               no RBF centers are allowed <default = .005>
%     N_eval - # points to evaluate error <default = 1001>
%
%  Some outputs are available if you would like them
%     iter_out - Print output during the solves <default = 0>
%     plot_err - log-log plot of error vs. epsilon <default = 1>
%     errcolor - Color for error line in log-log plot <default = 'b'>
%     condcolor - Color for condition line in log-log plot <default = 'r'>
%     save_fig - Save figures (.fig) (for plot_err = 1) <default = 0>
%     save_data - Save epvec, errvec, condvec <default = 0>
%
%  The results of these experiments are stored in
%     errvec - Errors computed at Nvec_true
%     condvec - Collocation matrix condition numbers at Nvec_true


R = 1;
sig = 0.02;
dipmom = 2.7*[1, 0, 0];
srcpnts = [0, 0, 0.6*R];
 
rbfset = { 'Gaussian' 'IMQ' 'MQ' 'Wendland_C4' 'Wendland_C6' };
epvec = logspace(0,2,100);
BC_choice = 1;
eval_diff = 1;

Npnts = 500;
BC_frac = .3; % Not yet implemented
dip_cushion = .01;
N_eval = 1001;

iter_out = 1;
plot_err = 1;
errcolor = 'b';
condcolor = 'r';
save_fig = 0;
save_data = 0;


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

% Determine collocation points
[POINTS, NORMALS] = BallGeometry(R,Npnts,'kansa',[],srcpnts,dip_cushion);
intdata = POINTS.int1;
bdydata = POINTS.bdy11;
N_int = size(intdata,1);
N_bdy = size(bdydata,1);

% Compose a vector of all the RBF centers
% For kansa, the centers and collocation points coincide
ctrs = [intdata; bdydata];

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

% Find all the necessary distance matrices
% For interior block
DM_intdata = DistanceMatrix(intdata,ctrs);
% For Dirichlet BC block
DM_bdydata_dir = DistanceMatrix(bdydata_dir,ctrs);
% For Neumann BC block
DM_bdydata_neu = DistanceMatrix(bdydata_neu,ctrs);
% For evaluation matrix
DM_eval = DistanceMatrix(evalpnts, ctrs);

% Find all the necessary difference matrices
dx_bdydata_neu = DifferenceMatrix(bdydata_neu(:,1),ctrs(:,1));
dy_bdydata_neu = DifferenceMatrix(bdydata_neu(:,2),ctrs(:,2));
dz_bdydata_neu = DifferenceMatrix(bdydata_neu(:,3),ctrs(:,3));

% Compute known-terms vector (a.k.a. righthand side vector)
% This requires the gradient of the unbounded potential at boundary
rhs_int = zeros(N_int,1);
gradphi_F_neu = gradphiF_dip(bdydata_neu, srcpnts, dipmom, sig);
rhs_bdy_neu = -sum(normvecs.*gradphi_F_neu,2);

% Compute the true solution to be used as Dirichlet BC
phi_F_bdy_dir = phiF_dip(bdydata_dir,srcpnts,dipmom,sig);
phi_bdy_dir = HomSpherePotential(R, sig, srcpnts, dipmom, bdydata_dir);
rhs_bdy_dir = phi_bdy_dir - phi_F_bdy_dir;

% Compose rhs
rhs = [rhs_int;rhs_bdy_neu;rhs_bdy_dir];


for j = 1:length(rbfset)
    [rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(rbfset{j});
    if iter_out
        fprintf('RBF: %s\n',rbfset{j})
    end
    % Loop through the requested epsilon values
    errvec = [];
    condvec = [];
    k = 1;
    for ep=epvec
        
        % Compute the collocation block for the interior
        LCM = Lrbf(ep,DM_intdata);
        
        % Compute normal derivative collocation matrix for boundary
        A = bsxfun(@times,normvecs(:,1),dxrbf(ep,DM_bdydata_neu,dx_bdydata_neu));
        B = bsxfun(@times,normvecs(:,2),dyrbf(ep,DM_bdydata_neu,dy_bdydata_neu));
        C = bsxfun(@times,normvecs(:,3),dzrbf(ep,DM_bdydata_neu,dz_bdydata_neu));
        BCM_neu = A + B + C;
        
        % Now we consider the Dirichlet BC component
        BCM_dir = rbf(ep,DM_bdydata_dir);
        
        % Create the full linear system from the blocksand solve it
        % Compose collocation matrix in same order as rhs
        CM = [LCM;BCM_neu;BCM_dir];
        % Coefficients for evaluation
        [coefs,recip_cond] = linsolve(CM,rhs);
        
        % Compute the evaluation matrix
        EM = rbf(ep, DM_eval);
        
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
        condvec(k) = 1/recip_cond;
        
        if iter_out
            fprintf('\tep = %g\n\t\terr = %g\n\t\tcond = %g\n\n',epvec(k),errvec(k),condvec(k));
        end
        k = k+1;
    end
    
    if plot_err
%         clf reset
        
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
                
        if save_fig
            h = figure('Name', strcat(rbfset{j}, ' RBF, ', ...
                num2str(N_int), ' interior points, ', ...
                num2str(N_bdy), ' boundary points'));
        end
        
        [AX,H1,H2] = plotyy(epvec,errvec,epvec,condvec,@loglog);
        xlabel('\epsilon')
        set(AX(1),'Xlim',[epvec(1),epvec(end)])
        set(get(AX(1),'Ylabel'),'String',errstr)
        set(AX(1),'Ycolor',errcolor)
        set(AX(2),'Xlim',[epvec(1),epvec(end)])
        set(get(AX(2),'Ylabel'),'String','Matrix Condition','Color',condcolor)
        set(AX(2),'Ycolor',condcolor)
        set(H1,'LineWidth',3,'Color',errcolor)
        set(H2,'LineWidth',3,'Color',condcolor)
        title(strcat(bcstr,', ',rbfset{j}))
        
        if save_fig
            filename = strcat(rbfset{j},'_',...
                              num2str(N_int),'int_',...
                              num2str(N_bdy),'bdy');
            saveas(h, filename,'fig')
        end
        
    end
    
    if save_data
        filename = strcat(rbfset{j},'_',...
            num2str(N_int),'int_',...
            num2str(N_bdy),'bdy');
        save(filename, 'errvec', 'condvec', 'epvec')
    end
    
end