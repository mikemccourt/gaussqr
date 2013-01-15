function [rel_err_phi, COND] = epskansa(ep, NORMALS, rhs, DM_eval, DM_intdata, DM_bdydata, dx_bdydata, dy_bdydata, dz_bdydata, phi_F, phi_an, rbf, dxrbf, dyrbf, dzrbf, Lrbf)
% Compute evaluation matrix EM
EM = rbf(ep, DM_eval);

% Compute blocks for collocation matrix
% Interior points
LCM = Lrbf(ep,DM_intdata);

% Boundary points
A = bsxfun(@times,NORMALS.n11(:,1),dxrbf(ep,DM_bdydata,dx_bdydata));
B = bsxfun(@times,NORMALS.n11(:,2),dyrbf(ep,DM_bdydata,dy_bdydata));
C = bsxfun(@times,NORMALS.n11(:,3),dzrbf(ep,DM_bdydata,dz_bdydata));
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

% Comparison and maximum errors
%--------------------------------------------------------------------------
% Potential
normdiff_phi = norm(phi-phi_an);
rel_err_phi = normdiff_phi/norm(phi_an);
COND = cond(CM);
end