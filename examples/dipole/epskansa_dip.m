function [rms_err_phi, rel_err_phi, max_err_phi, COND ] = epskansa(ep, NV, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf, neval)

% Compute evaluation matrix EM
EM = rbf(ep, DM_eval);

% Compute blocks for collocation matrix
% Interior points
LCM = Lrbf(ep,DM_intdata);
% Boundary points

A = bsxfun(@times,NV(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
B = bsxfun(@times,NV(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
C = bsxfun(@times,NV(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));

BCM = bsxfun(@plus,A,B);
BCM = bsxfun(@plus,BCM,C);

% Collocation matrix
CM = [LCM; BCM];


% Numerical solution for the potential
%--------------------------------------------------------------------------

% Potential at evalpnts in the source free case
phi0 = EM * (CM\rhs);

% Potential at evalpnts (superposition of effects)
phi = phi0 + phi_F;




% Analytic solution for the magnetic induction field
%--------------------------------------------------------------------------
% B_exact = HomSphereInduction(srcpnts, dipmom, obspnts);


% Comparison and maximum errors
%--------------------------------------------------------------------------
% Potential
max_err_phi = norm(phi-phi_an,inf);
normdiff_phi = norm(phi-phi_an);
rel_err_phi = normdiff_phi/norm(phi_an);
rms_err_phi = norm(phi-phi_an)/neval;
COND = cond(CM);
% Magnetic induction field
% max_err_B = norm(B-B_an,inf);
% normdiff_B = norm(B-B_an);
% rel_err_B = normdiff_B/norm(B_an);
% rms_err_B = norm(B-B_an)/size(obspnts,1);