function phi = phiF( obspnts, srcpnts, dipmom, sigma )
% phiF calculates the electric potential at M observation points in a 
% homogeneus and infinite medium by superposing the effects of N current 
% dipoles.
%
% Input arguments:
% obspnts  =  Mx3 matrix of observation points coordinates [m].
% srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% dipmom   =  Nx3 matrix of dipoles moments components [A*m].
% sigma    =  electric conductivity of the medium [S/m].
%
% Outputs:
% phi      =  Mx1 vector of electric potential at observation points [V].
%
M = size(obspnts,1); N = size(srcpnts,1);
phi = zeros(M,1);
for i = 1:M
    for j = 1:N
        R = obspnts(i,:) - srcpnts(j,:);
        RP = R * dipmom(j,:)';
        a = R(1)*R(1) + R(2)*R(2) + R(3)*R(3); 
        a = a^(3/2);
        phi(i) = phi(i) + RP / a;
    end
end
a = 1/(4*pi*sigma);
phi = a*phi;
end