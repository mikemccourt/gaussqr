clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MEG meshfree forward solver for a multisphere model
%                           - Kansa's method -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Requires:
%   SphereRegGoldPoints2.m
%   SphereSurfGoldPoints.m
%   SphereRegUnifPoints.m
%   DistanceMatrix.m
%   DifferenceMatrix.m
%   gradphiF_dip.m
%   phiF_dip.m
%   MultiSpherePotential.m
%
% To do:
%   -- Computation of J vector as -sig*grad(phi).
%   -- Computation of B (magnitude?normal component?) with Ampere-Laplace
%      formula.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: singularity of analytic potential formula for dipole at origin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified: 2012/10/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Input data
%--------------------------------------------------------------------------

% Medium data
R = [0.07, 0.1];              % Sphere radii [m]
sig = [0.2, 0.2];           % Electric conductivities [S/m]
% Sources data
dipmom = 2.7e-12.*[1, 0, 0];  % Dipole moment [Am]
srcpnts = [0, 0, 0.6*R(1)];   % Dipole position [m]
% Magnetic induction field observation points
% obspnts = ; ...
% Parameters for numerical computation
radbasfun = 'gaussian';            % Radial basis function
% Npnts_surf = [800, 1000];             % Number of desired points on sphere's surface
Npnts_surf = 300;    
ep = 50;                    % RBFs shape parameter

% RBF definition and derivatives
%--------------------------------------------------------------------------
% rbf     = Radial basis function;
% dxrbf   = component along x of the gradient of the RBF
% dyrbf   = component along y of the gradient of the RBF
% dzrbf   = component along z of the gradient of the RBF
% Lrbf    = Laplacian of the RBF in 3D
[rbf, dxrbf, dyrbf, dzrbf, Lrbf] = pickRBF(radbasfun);


% Collocation points
[dsites, Npnts, NV, dist] = SphereRegGoldPoints2(Npnts_surf, R, 'kansa');
% [dsites, Npnts, NV] = SphereRegGoldPoints3(R, Npnts_surf, 'kansa');

Na_int = Npnts(1);      % No. of interior points of the innermost layer (A)
Na_bdy = Npnts(2);      % No. of points on internal interface considered as
                        % belonging to the innermost layer
Nb_bdy = Npnts(3);      % No. of points on internal interface considered as
                        % belonging to the outermost layer
Nb_int = Npnts(4);      % No. of interior points of the outermost layer (B)
N0     = Npnts(5);      % No. of points on the external boundary
Npnts_tot = sum(Npnts,1); % Total no. of points

a_interior   = 1 : Na_int;
a_interface  = Na_int + 1 : Na_int + Na_bdy;
a_all        = 1 : Na_int + Na_bdy;
b_interior   = Na_int + Na_bdy + Nb_bdy + 1 : Na_int + Na_bdy + Nb_bdy + Nb_int;
b_interface  = Na_int + Na_bdy + 1 : Na_int + Na_bdy + Nb_bdy;
out_boundary = Na_int + Na_bdy + Nb_bdy + Nb_int + 1 : Npnts_tot;
b_all        = Na_int + Na_bdy + 1 : Npnts_tot;

% Centers
ctrs = dsites;

% Evaluation points
% evalpnts = SphereSurfGoldPoints(Npnts_surf/2, R(size(R,2))); 
[evalpnts, tri_evalpnts] = SphereSurfUniformPoints(R(2),3);
neval = size(evalpnts,1);

% Compute evaluation matrix EM
EM = DistanceMatrix(evalpnts, ctrs);
EM = rbf(ep, EM);

% Collocation matrix and known-terms vector
%--------------------------------------------------------------------------

L_AA = DistanceMatrix( dsites(a_all,:), ctrs(a_all,:) );
L_AA = rbf(ep,L_AA);
% A = inv(L_AA);
[L_AA,U_AA] = lu(L_AA);

L_BB = DistanceMatrix( dsites(b_all,:), ctrs(b_all,:) );
L_BB = rbf(ep,L_BB);
% B = inv(L_BB);
[L_BB,U_BB] = lu(L_BB);


% Interior points innermost layer: J1
J1 = DistanceMatrix(dsites(a_interior,:), ctrs(a_all,:));
J1 = Lrbf(ep,J1);
J1 = (J1/U_AA)/L_AA;


%   Continuity of the potential on the internal boundary: J2 and J4
% Interpolation of the potential given (Nb_bdy+Nb_int+N0) points in the 
% outermost layer (interpolation matrix B) and evaluation of the potential
% at Na_bdy points on the internal interface considered as belonging to the
% innermost layer (evaluation matrix A).

J2 = [zeros(Na_bdy,Na_int), eye(Na_bdy, Na_bdy)];
% J2 = (J2/U_AA)/L_AA;

J4 = DistanceMatrix( dsites(a_interface,:), ctrs(b_all,:) );
J4 = rbf(ep,J4);
J4 = -(J4/U_BB)/L_BB;



%   Continuity of the current density on the internal boundary: J3 and J5
% Interpolation of the potential given (Na_int+Na_bdy) points in the 
% innermost layer (interpolation matrix AA) and evaluation of the current 
% density at Nb_bdy points on the internal interface considered as 
% belonging to the outermost layer (evaluation matrix BB involving normal
% derivatives of the potential).

% J3
J3 = DistanceMatrix( dsites(b_interface,:), ctrs(a_all,:) );
dx = DifferenceMatrix( dsites(b_interface,1), ctrs(a_all,1) );
dy = DifferenceMatrix( dsites(b_interface,2), ctrs(a_all,2) );
dz = DifferenceMatrix( dsites(b_interface,3), ctrs(a_all,3) );
A = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,1),dxrbf(ep,J3,dx));
B = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,2),dyrbf(ep,J3,dy));
C = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,3),dzrbf(ep,J3,dz));
B = bsxfun(@plus,A,B);
B = bsxfun(@plus,B,C);
J3 = sig(1).*(B/U_AA)/L_AA;

%J5
J5 = DistanceMatrix( dsites(b_interface,:), ctrs(b_all,:) );
dx = DifferenceMatrix( dsites(b_interface,1), ctrs(b_all,1) );
dy = DifferenceMatrix( dsites(b_interface,2), ctrs(b_all,2) );
dz = DifferenceMatrix( dsites(b_interface,3), ctrs(b_all,3) );
A = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,1),dxrbf(ep,J5,dx));
B = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,2),dyrbf(ep,J5,dy));
C = bsxfun(@times,NV(Na_bdy+1:Na_bdy+Nb_bdy,3),dzrbf(ep,J5,dz));
B = bsxfun(@plus,A,B);
B = bsxfun(@plus,B,C);
J5 = -sig(2).*(B/U_BB)/L_BB;

% Interior points outermost layer: J6
J6 = DistanceMatrix( dsites(b_interior,:), ctrs(b_all,:) );
J6 = Lrbf(ep,J6);
J6 = (J6/U_BB)/L_BB;


% Points on the outermost surface: J7
J7 = DistanceMatrix( dsites(out_boundary,:), ctrs(b_all, :) );
dx = DifferenceMatrix( dsites(out_boundary,1), ctrs(b_all,1) );
dy = DifferenceMatrix( dsites(out_boundary,2), ctrs(b_all,2) );
dz = DifferenceMatrix( dsites(out_boundary,3), ctrs(b_all,3) );
A = bsxfun(@times,NV(Na_bdy+Nb_bdy+1:Na_bdy+Nb_bdy+N0,1),dxrbf(ep,J7,dx));
B = bsxfun(@times,NV(Na_bdy+Nb_bdy+1:Na_bdy+Nb_bdy+N0,2),dyrbf(ep,J7,dy));
C = bsxfun(@times,NV(Na_bdy+Nb_bdy+1:Na_bdy+Nb_bdy+N0,3),dzrbf(ep,J7,dz));
B = bsxfun(@plus,A,B);
B = bsxfun(@plus,B,C);
J7 = (B/U_BB)/L_BB;

% Collocation matrix
CM = [ J1                          , zeros(Na_int,Nb_bdy+Nb_int+N0);
       J2                          , J4 ;
       J3                          , J5 ;
       zeros(Nb_int,Na_int+Na_bdy) , J6 ;
       zeros(N0,Na_int+Na_bdy)     , J7                             ];
   

% Compute known-terms vector (a.k.a. righthand side vector)
rhs = [ zeros(Na_int+Na_bdy,1);
    
        -(sig(1)-sig(2)).*sum(NV(Na_bdy+1:Na_bdy+Nb_bdy, :).*...
        gradphiF_dip(dsites(b_interface,:),srcpnts,dipmom,sig(1)), 2);
        
        zeros(Nb_int,1);
        
        -sum(NV(Na_bdy+Nb_bdy+1:Na_bdy+Nb_bdy+N0, :).*...
        gradphiF_dip(dsites(out_boundary,:),srcpnts,dipmom,sig(1)), 2) ];


% Numerical solution for the potential
%--------------------------------------------------------------------------
% Potential at evalpnts in the source free case
phi0 = EM * (CM\rhs);
% Potential at evalpnts in the unbound domain case
phi_F = phiF_dip(evalpnts, srcpnts, dipmom, sig(1));
% Potential at evalpnts (superposition of effects)
phi = phi0 + phi_F;


% Numerical solution for the magnetic induction field
%--------------------------------------------------------------------------
% TO DO

%  Analytic solution for the potential
%--------------------------------------------------------------------------
phi_an = MultiSpherePotential(R, sig, srcpnts, dipmom, evalpnts, 30);

% Comparison and maximum errors
%--------------------------------------------------------------------------
% Potential
err = phi-phi_an;
max_err_phi = norm(err,inf);
normdiff_phi = norm(err);
rel_err_phi = normdiff_phi/norm(phi_an)
rms_err_phi = normdiff_phi/neval;
% COND = cond(CM)


% Plots
%--------------------------------------------------------------------------

% figure('Color',[1 1 1]);
% subplot(2,2,1)
% SurfacePlot_dip(evalpnts, tri_evalpnts, phi_an, 'Analytic potential')
% subplot(2,2,2)
% SurfacePlot_dip(evalpnts, tri_evalpnts, phi, 'Computed potential \Phi = \phi_0 + \phi_F')
% subplot(2,2,3)
 SurfacePlot_dip(evalpnts, tri_evalpnts, phi0, '\phi_0')
% subplot(2,2,4)
% SurfacePlot_dip(evalpnts, tri_evalpnts, phi_F, '\phi_F')

% figure('Color',[1 1 1]);
% SurfacePlot_dip(evalpnts, tri_evalpnts, abs(err./phi_an), 'Relative error')


% figure('Color',[1 1 1]);
% plot3(dsites(a_interior,1),...
%       dsites(a_interior,2),...
%       dsites(a_interior,3),'o'); % Interior A points
% hold on;
% plot3(dsites(a_interface,1),...
%       dsites(a_interface,2),...
%       dsites(a_interface,3),'ro'); % Interface A points
% plot3(dsites(b_interface,1),...
%       dsites(b_interface,2),...
%       dsites(b_interface,3),'go'); % Interface B points
% plot3(dsites(b_interior,1),...
%       dsites(b_interior,2),...
%       dsites(b_interior,3),'bo'); % Interior B points
% plot3(dsites(out_boundary,1),...
%       dsites(out_boundary,2),...
%       dsites(out_boundary,3),'co'); % External boundary points
% axis equal