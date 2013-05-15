function B = Induction( obspnts,  srcpnts, dipmom, ...
                        sol_type, radbasfun, ep, ...
                        coefs, ctrs, ...
                        sig, radius,  N )
% Induction calculates the magnetic (induction) field at M observation 
% points outside the head due to N current dipoles. 
% The relative magnetic permeability is assumed to be 1.
%
% The magnetic field is given by the sum of two terms:
% - one takes into account the effect of primary currents only. This term
%   is evaluated by means of an analytical formula.
% - one takes into account the effect of secondary (volume) currents only.
%   This term is evaluated by means of numerical integration.
%
% Concerning the numerical term, we make use of a corollary of the 
% Divergence Theorem to integrate the potential function over the 
% interfaces, instead of the current density (i.e. the gradient of the 
% potential function) over the entire domain.
% For details see, e.g., the review on MEG by Hamalainen et al.:
% Hämäläinen M., et al. "Magnetoencephalography—theory, instrumentation, 
% and applications to noninvasive studies of the working human brain." 
% Reviews of Modern Physics 65.2 (1993):413.
%
% Note that in MEG literature it is well established that the magnetic
% field is basically independent of currents flowing in the skull and in
% the scalp.
% This fact was observed for the first time in:
% Hämäläinen M., and Sarvas J. "Realistic conductivity geometry model of 
% the human head for interpretation of neuromagnetic data." Biomedical 
% Engineering, IEEE Transactions on 36.2 (1989): 165-171.
% and widely verified/accepted later on.
% Therefore, even though we take all the layers into account to evaluate
% the potential function on the scalp (EEG problem), we can consider only 
% the innermost interface to evaluate the magnetic field outside the head 
% (MEG problem) without an appreciable loss in accuracy.
% This is the reason why this function works with only one layer. We could
% generalize it if needed.
%
% We use (right now) a very simple quasi Monte Carlo integration strategy 
% with sampling Halton points on the interface. 
% We need to know surface area in order to evaluate the integral properly.
%
% Right now we deal with simple spherical surfaces, so the sampling is
% straightforward and the interface has known analytical surface area. 
% For real geometries we have to think about:
% - the sampling (map points from a unit square to real surface?);
% - a way of estimating surface areas.
% It will depend on how we will handle the geometry. Dealing with point
% clouds only (RBF surfaces interpolation) would be ideal. But also very 
% rough meshes could be ok since they are easy and fast to build (however 
% someone could say our method wouldn't be fully meshfree/messfree...)
%
% *****NOTE*****
% This function requires the sample size as input. Of course this "static" 
% approach is not optimal in general: we should prefer an adaptive strategy 
% with a certain tolerance provided as input. 
% It is very easy to modify the code to reach this goal taking advantage of
% the nested sampling.
%
% *****WARNING***** 
% bsxfun is used to speed up the code. There might be compatibility issues 
% with old versions of Matlab.
%
% Input:
%   obspnts    = Mx3 matrix of observation points coordinates [m].
%   srcpnts    = Nx3 matrix of source dipoles coordinates [m].
%   dipmom     = Nx3 matrix of dipoles moments components [A*m].
%
%   sol_type   = Solver type
%                   'kansa' : Nonsymmetric collocation.
%                   'mfs' : Method of fundamental solutions.
%   radbasfun  = RBF for collocation (it matters only if sol_type =
%                'kansa')
%   ep         = RBF shape parameter (it matters only if sol_type =
%                'kansa')
%
%   coefs      = coefficients of RBF expansion.
%   ctrs       = centers for RBF expansion.
%
%   sig        = electrical conductivity of the medium [S/m].
%   radius     = sphere radius [m].
%   N          = No. of points over the sphere for integration. 
%
% Output:
%   B          = Mx3 matrix of magnetic field at observation points.
%
% Calls on: 
%   pickRBF.m
%   SphereSurfHaltonPoints.m
%   DistanceMatrix.m
%   BF.m
%

% Check input data
P = size(srcpnts); Q = size(dipmom);
if not(isequal(P,Q))
    error('srcpnts and dipmom dimensions must match!')
end

M = size(obspnts,1); % No. of field observation points outside the head
integral = zeros(M,3); % Inizialize the matrix that will contain the
                       % integral (numerical) part of the magnetic field

% RBF definition
% We could think about a better way to call it...
if strcmp(sol_type,'kansa')
    rbf = pickRBF(radbasfun);
else
    rbf = pickRBF('fundamental_3d');
end

% Sphere surface area
surf = 4*pi*radius^2; 

% Source points to be selected on the surface
r_src = SphereSurfHaltonPoints(1:N,radius);
normals = r_src./radius; % Normal unit vectors

% Evaluate potentials at source points
phi_F = phiF_dip(r_src,srcpnts,dipmom,sig);
DM_eval = DistanceMatrix(r_src, ctrs);
EM = rbf(ep, DM_eval);
phi0 = EM * coefs;
phi = phi0 + phi_F;
phi = phi - phi(1); % To evaluate magnetic field the potential difference 
                    % shouldn't be necessary.

for i = 1:M
    % Distances between source points and field observation point
    R = bsxfun(@minus,obspnts(i,:),r_src);
    R_norm = R(:,1).*R(:,1)+R(:,2).*R(:,2)+R(:,3).*R(:,3);
    R_norm = R_norm.^(3/2);
    
    % Evaluate the integrand function
    F = [ normals(:,2).*R(:,3)-normals(:,3).*R(:,2), ...
          normals(:,3).*R(:,1)-normals(:,1).*R(:,3), ...
          normals(:,1).*R(:,2)-normals(:,2).*R(:,1)];
    F = bsxfun(@rdivide,F,R_norm);
    F = bsxfun(@times,F,phi);
    
    integral(i,:) = integral(i,:) + sum(F,1);
end
B = - 1e-07 * sig * surf * integral / N;

% Adding the analitycal part
B = B + BF(obspnts,srcpnts,dipmom);
end