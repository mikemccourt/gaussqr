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
Npnts = 300;                % Number of desired interior points
epvec = logspace(0,2,20);  % Vector of epsilon values for study


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

% Find all the necessary distance matrices
DM_eval = DistanceMatrix(evalpnts, ctrs);
DM_intdata = DistanceMatrix(intdata,ctrs);
DM_bdydata = DistanceMatrix(bdydata,ctrs);

% Find all the necessary difference matrices
dx_bdydata = DifferenceMatrix(bdydata(:,1),ctrs(:,1));
dy_bdydata = DifferenceMatrix(bdydata(:,2),ctrs(:,2));
dz_bdydata = DifferenceMatrix(bdydata(:,3),ctrs(:,3));

% Compute known-terms vector (a.k.a. righthand side vector)
gradphi_F = gradphiF_dip(bdydata, srcpnts, dipmom, sig);% Gradient of the 
                                                    % potential at boundary
                                                    % in the unbound case
rhs = [ zeros(size(intdata,1),1); -sum(NORMALS.n11.*gradphi_F,2) ];

% Potential at evalpnts in the unbound domain case
% This is the analytic component of the computed solution
phi_F = phiF_dip(evalpnts,srcpnts,dipmom,sig);

%  Analytic solution for the potential
%--------------------------------------------------------------------------
phi_an = HomSpherePotential(R, sig, srcpnts, dipmom, evalpnts);

errvec = [];
condvec = [];
k = 1;
for ep=epvec
    fprintf('ep=%g\n',ep)
    
    % Compute blocks for collocation matrix using this ep
    % Interior points
    LCM = Lrbf(ep,DM_intdata);
    % Boundary points
    A = bsxfun(@times,NORMALS.n11(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
    B = bsxfun(@times,NORMALS.n11(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
    C = bsxfun(@times,NORMALS.n11(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));
    BCM = bsxfun(@plus,A,B);
    BCM = bsxfun(@plus,BCM,C);
    % Evaluation matrix
    EM = rbf(ep, DM_eval);
    
    % Collocation matrix
    CM = [LCM; BCM];
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
    k = k + 1;
end

clf reset

[AX,H1,H2] = plotyy(epvec,errvec,epvec,condvec,@loglog);
[min_err,min_ep_ind] = min(errvec);
fprintf('\tbest_err = %g\n\tbest_ep = %g\n\tN = %d\n',min_err,epvec(min_ep_ind),length(ctrs));