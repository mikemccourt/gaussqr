% This tests to determine how to approximate the true solution in the 1
% sphere case.
% We may need this to invoke the Neumann BC.
% If possible, we will compute the interpolant to this and then evaluate
% the derivative of the interpolant.
% I'm going to try to use the stable basis, we'll see if that works out.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = -200;

N = 700;
NN = 1000;
ep = .01;
alpha = 65; % Arrived at this value (65) after some experimenting

R = 0.1;                      % Sphere radius [m]
sig = 0.2;                    % Electric conductivity [S/m]
dipmom = 2.7e-12.*[1, 0, 0];  % Dipole moment [Am]
srcpnts = [0, 0, 0.6*R];      % Dipole position [m]

[POINTS,NORMALS] = BallGeometry(R,N,'kansa');
x = [POINTS.int1;POINTS.bdy11];
u = HomSpherePotential( R, sig, srcpnts, dipmom, x );

xx = SphereSurfGoldPoints(NN, R);
uu = HomSpherePotential( R, sig, srcpnts, dipmom, xx );

% epvec = logspace(-3,0,30);
% errvec = [];
% k = 1;
% for ep=epvec
%     GQR = gqr_solve(x,u,ep,alpha);
%     us = gqr_eval(GQR,xx);
%     errvec(k) = errcompute(us,uu);
%     k = k + 1
% end
% loglog(epvec,errvec)
% pause

% Delaunay triangulation of the sphere
d_evalpnts = delaunayn(xx);
tr = TriRep(d_evalpnts, xx);
tr_surf = freeBoundary(tr);

us = gqr_eval(gqr_solve(x,u,ep,alpha),xx);
clf reset
subplot(1,3,1)
SurfacePlot_dip(evalpnts, tr_surf, uu, 'True solution')
subplot(1,3,2)
SurfacePlot_dip(evalpnts, tr_surf, us, 'Computed solution')
subplot(1,3,3)
SurfacePlot_dip(evalpnts, tr_surf, uu-us, 'Error')