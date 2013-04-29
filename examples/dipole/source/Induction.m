function B = Induction( r_obs, radius, sig, srcpnts, dipmom, ep, coefs, ctrs, radbasfun, N, sol_type )

if strcmp(sol_type,'kansa')
    rbf = pickRBF(radbasfun);
else
    rbf = pickRBF('fundamental_3d');
end


surf = 4*pi*radius^2; % Sphere surface
hs = haltonset(2); % 2D Halton set of points (in the unit square)
M = size(r_obs,1);
integral = zeros(M,3);

% Source points to be selected on the surface
A = hs(1:N,:);
theta = 2*pi*A(:,1);
phi = acos(2*A(:,2)-1);
r_src = [ radius.*sin(phi).*cos(theta), ...
    radius.*sin(phi).*sin(theta), ...
    radius.*cos(phi)];
normals = r_src./radius; % Normal unit vectors

% Evaluate potentials at source points
phi_F = phiF_dip(r_src,srcpnts,dipmom,sig);
DM_eval = DistanceMatrix(r_src, ctrs);
EM = rbf(ep, DM_eval);
phi0 = EM * coefs;
phi = phi0 + phi_F;
phi = phi - phi(1);

for i = 1:M
    % Distances between source points and field observation point
    R = bsxfun(@minus,r_obs(i,:),r_src);
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
% plot3(r_src(:,1),r_src(:,2),r_src(:,3),'o'); axis square
B = - 1e-07 * sig * surf * integral / N;

B = B + BF(r_obs,srcpnts,dipmom);
end


% 
% 
% function [ B, res ] = Induction( r_obs, radius, sig, srcpnts, dipmom, ep, coefs, ctrs, radbasfun, reltol, sol_type )
% 
% deltaN = 100; % Points to be added at each iteration
% 
% bdim = 2*radius; % Bounding cube dimension
% bbvol = bdim^3;
% hs = haltonset(3); % 3D Halton set of points (in the unit cube)
% 
% N1 = 1;
% N2 = 500;
% 
% M = size(r_obs,1);
% integral_old = zeros(M,3);
% integral = integral_old;
% B = integral_old;
% res = reltol+eps;
% % while res >= reltol
%     % Source points to be selected in a bounding cube containing the domain
%     r_src = hs(N1:N2,:) * bdim - radius;
%     
%     % Select source points inside the domain
%     r_src = r_src(r_src(:,1).^2 + r_src(:,2).^2 + r_src(:,3).^2 < radius^2, :); 
% 
%     % Evaluate current density at source points
%     J = -sig*gradphiF_dip(r_src, srcpnts, dipmom, sig);
%     J = J-sig*gradphi0_dip(ep, coefs, r_src, ctrs, radbasfun);
% 
%     for i = 1:M
%         % Distances between source points and field observation point
%         R = bsxfun(@minus,r_obs(i,:),r_src);
%         R_norm = R(:,1).*R(:,1)+R(:,2).*R(:,2)+R(:,3).*R(:,3);
%         R_norm = R_norm.^(3/2);
% 
%         % Evaluate the integrand function
%         F = [ J(:,2).*R(:,3)-J(:,3).*R(:,2), ...
%               J(:,3).*R(:,1)-J(:,1).*R(:,3), ...
%               J(:,1).*R(:,2)-J(:,2).*R(:,1)];
%         F = bsxfun(@rdivide,F,R_norm);
% 
%         integral(i,:) = integral(i,:) + sum(F,1);
%     end
%     
%     B_old = B;
%     B = 1e-07 * bbvol * integral / N2;
%     
%     N1 = N2 + 1;
%     N2 = N2 + deltaN;
%     
%     res = norm(B-B_old);
% % end
% B = B + BF(r_obs,srcpnts,dipmom);
% end