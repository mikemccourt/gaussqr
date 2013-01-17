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

% Choose which function you would like to test with
%   1 - Use the standard potential function
%   2 - Use the constant value 1e-10
%   3 - Use the function cos(sqrt(x^2+y^2+z^2)*pi/R)
test_choice = 3;

switch test_choice
    case 1
        uf = @(x) HomSpherePotential(R,sig,srcpnts,dipmom,x);
    case 2
        uf = @(x) 1e-10*ones(size(x,1),1);
    case 3
        uf = @(x) cos(sqrt(sum(x.^2,2))*pi/R);
end

[POINTS,NORMALS] = BallGeometry(R,N,'kansa');
x = [POINTS.int1;POINTS.bdy11];
u = uf(x);
% u = HomSpherePotential( R, sig, srcpnts, dipmom, x );
% u = 1e-10*ones(size(x,1),1);

xx = SphereSurfGoldPoints(NN, R);
uu = uf(xx);
% uu = HomSpherePotential( R, sig, srcpnts, dipmom, xx );
% uu = 1e-10*ones(size(xx,1),1);

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

GQR = gqr_solve(x,u,ep,alpha);
us = gqr_eval(GQR,xx);
us_x = gqr_eval(GQR,xx,[1,0,0]);
us_y = gqr_eval(GQR,xx,[0,1,0]);
us_z = gqr_eval(GQR,xx,[0,0,1]);

normvecs = xx/R;
us_n = sum(normvecs.*[us_x,us_y,us_z],2);

clf reset
subplot(2,2,1)
SurfacePlot_dip(evalpnts, tr_surf, uu, 'True solution')
subplot(2,2,2)
SurfacePlot_dip(evalpnts, tr_surf, us, 'Computed solution')
subplot(2,2,3)
SurfacePlot_dip(evalpnts, tr_surf, uu-us, 'Error')
subplot(2,2,4)
SurfacePlot_dip(evalpnts, tr_surf, us_n, 'Normal Derivative')