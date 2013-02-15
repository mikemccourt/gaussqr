function G = gradphi0_dip(ep, c, obspnts, ctrs, rbf_type)
% gradphi0 calculates the gradient of a RBF expansion (N terms) at M 
% observation points.
%
% Input arguments:
% ep       =  RBFs shape parameter.
% c        =  N vector of expansion coefficients.
% obspnts  =  Mx3 matrix of observation points coordinates.
% ctrs     =  RBFs centers.
% rbf_type =  type of RBF.
%
% Outputs:
% G =  Mx3 matrix of gradient at observation points.
%

[~, dxrbf, dyrbf, dzrbf] = pickRBF(rbf_type);

DM_eval = DistanceMatrix(obspnts,ctrs);
dx = DifferenceMatrix(obspnts(:,1),ctrs(:,1));
dy = DifferenceMatrix(obspnts(:,2),ctrs(:,2));
dz = DifferenceMatrix(obspnts(:,3),ctrs(:,3));

Gx = dxrbf(ep, DM_eval, dx);
Gy = dyrbf(ep, DM_eval, dy);
Gz = dzrbf(ep, DM_eval, dz);

Gx = Gx * c;
Gy = Gy * c;
Gz = Gz * c;

G = [Gx, Gy, Gz];

end