function V = HomSpherePotential( r, sigma, srcpnts, dipmom, obspnts )
% HomSpherePotential calculates the electric potential in M points 
% insiede a homogeneous sphere due to N current dipoles, making use of the 
% analytical formula given in:
% Yao, "Electric Potential Produced by a Dipole in a Homogeneous Conducting
% Sphere", IEEE Transactions on Biomedical Engineering, 47(7) (2000), 
% pp. 964-966.
%
% Input arguments:
% r        =  sphere radius [m].
% sigma    =  sphere electric conductivity [S/m].
% srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% dipmom   =  Nx3 matrix of dipole moments components [A*m].
% obspnts  =  Mx3 matrix of observation points coordinates [m].
%
% Output:
% V        =  Mx1 vector of electric potential at observation points [V].
%
%--------------------------------------------------------------------------
% !!WARNING!!: the formula is singular for dipoles located at sphere origin
%--------------------------------------------------------------------------
%
M = size(obspnts,1); N = size(srcpnts,1);
V = zeros(M,1);
for i = 1:M
        lungr2 = obspnts(i,1)*obspnts(i,1) + ...
                 obspnts(i,2)*obspnts(i,2) + ...
                 obspnts(i,3)*obspnts(i,3);
        lungr = sqrt(lungr2);
    for j = 1:N
        lungr0 = srcpnts(j,1)*srcpnts(j,1) + ...
                 srcpnts(j,2)*srcpnts(j,2) + ...
                 srcpnts(j,3)*srcpnts(j,3);
        lungr0 = sqrt(lungr0);
        B = lungr*lungr0;
        A = B/r/r;
        cosfi = (obspnts(i,:)*srcpnts(j,:)') / B;
        rp = lungr2 + lungr0*lungr0 - 2*B*cosfi;
        rp = sqrt(rp);
        Acosfi = A*cosfi;
        rpi = 1 + A*A - 2*Acosfi;
        rpi = sqrt(rpi);
        C = (obspnts(i,:)-srcpnts(j,:))/rp^3 + ...
            (obspnts(i,:)-srcpnts(j,:)*lungr2/r/r)/(rpi^3*r^3) + ...
            1/(r^3*rpi)*(obspnts(i,:) + ...
            (obspnts(i,:)*Acosfi-lungr2/r/r*srcpnts(j,:))/(rpi+1-Acosfi));
        V(i) = V(i) + dipmom(j,:)*C';
    end
end
a = 1/(4*pi*sigma);
V = a*V;
end