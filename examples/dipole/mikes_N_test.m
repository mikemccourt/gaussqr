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
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
GAUSSQR_PARAMETERS.ERROR_STYLE = 3;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

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
radbasfun = 'imq';     % Radial basis function
Nvec = 100:50:1000;                % Number of desired interior points
ep = 30;                % RBF shape parameter
BC_choice = 1;          % 1 - neumann, 2 - dirichlet

% RBF definition and derivatives
%--------------------------------------------------------------------------
% rbf     = Radial basis function;
% dxrbf   = component along x of the gradient of the RBF
% dyrbf   = component along y of the gradient of the RBF
% dzrbf   = component along z of the gradient of the RBF
% Lrbf    = Laplacian of the RBF in 3D
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);

% Evaluation points
evalpnts = SphereSurfGoldPoints(1000, R);
neval = size(evalpnts,1);

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);

%  Analytic solution for the potential
%--------------------------------------------------------------------------
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);


% Loop through the requested N values
errvec = [];
condvec = [];
k = 1;
for Npnts = Nvec
    fprintf('k=%d\n',k)
    
    % Determine collocation points
    [POINTS, NORMALS] = BallGeometry(R, Npnts, 'kansa');
    intdata = POINTS.int1;
    bdydata = POINTS.bdy11;
    N_int = size(intdata,1);
    N_bdy = size(bdydata,1);
    
    % Compose a vector of all the RBF centers
    ctrs = [intdata; bdydata];
    
    
    % Compute the collocation block for the interior
    DM_intdata = DistanceMatrix(intdata,ctrs);
    LCM = Lrbf(ep,DM_intdata);
    rhs_int = zeros(N_int,1);
    
    % Compute the evaluation matrix
    DM_eval = DistanceMatrix(evalpnts, ctrs);
    EM = rbf(ep, DM_eval);
    
    % Compute the collocation block for the boundary conditions
    % This also computes the RHS for the problem
    % Find the distance matrix for the boundary conditions
    DM_bdydata = DistanceMatrix(bdydata,ctrs);
    if BC_choice==1 % Do the standard Neumann BC
        % Find all the necessary difference matrices
        dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
        dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
        dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));
        
        % Compute normal derivative collocation matrix for boundary
        A = bsxfun(@times,NORMALS.n11(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
        B = bsxfun(@times,NORMALS.n11(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
        C = bsxfun(@times,NORMALS.n11(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));
%         A = diag(NORMALS.n11(:,1))*dxrbf(ep,DM_bdydata,dx_bdydata);
%         B = diag(NORMALS.n11(:,2))*dyrbf(ep,DM_bdydata,dy_bdydata);
%         C = diag(NORMALS.n11(:,3))*dzrbf(ep,DM_bdydata,dz_bdydata);
        BCM = A + B + C;
    
        % Compute known-terms vector (a.k.a. righthand side vector)
        % This requires the gradient of the unbounded potential at boundary
        gradphi_F = gradphiF_dip(bdydata, srcpnts, dipmom, sig);
        rhs_bdy = -sum(NORMALS.n11.*gradphi_F,2);
    else % Run a test with Dirichlet BC
        BCM = rbf(ep,DM_bdydata);
        
        % Compute the true solution to be used as BC
        phi_F_bdy = phiF_dip(bdydata,srcpnts,dipmom,sig);
        phi_bdy = HomSpherePotential(R, sig, srcpnts, dipmom, bdydata);
        rhs_bdy = phi_bdy - phi_F_bdy;
    end
    
    
    % Compose rhs
    rhs = [rhs_int;rhs_bdy];
    % Compose collocation matrix in same order as rhs
    CM = [LCM;BCM];
    % Coefficients for evaluation
    [coefs,recip_cond] = linsolve(CM,rhs);
    
    % Potential at evalpnts in the source free case
    phi0 = EM * coefs;
    % Potential at evalpnts (superposition of effects)
    phi = phi0 + phi_F;
    
    % Comparison and maximum errors
    %--------------------------------------------------------------------------
    % Potential
    errvec(k) = errcompute(phi,phi_an);
    condvec(k) = 1/recip_cond;
    Nvec(k) = N_int + N_bdy;
    fprintf('\terr = %g\n\tcond = %g\n\tN = %d\n',errvec(k),condvec(k),Nvec(k));
    k = k + 1;
end

clf reset

[AX,H1,H2] = plotyy(Nvec,errvec,Nvec,condvec,@semilogy);
xlabel('Total collocation points')