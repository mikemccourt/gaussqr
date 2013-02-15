function grad_phi = gradphiF_dip( obspnts, srcpnts, dipmom, sigma )
% gradphiF calculates the gradient of electric potential at M observation 
% points in a homogeneus and infinite medium by superposing the effects of 
% N current dipoles.
%
% Input arguments:
% obspnts  =  Mx3 matrix of observation points coordinates [m].
% srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% dipmom   =  Nx3 matrix of dipoles moments components [A*m].
% sigma    =  electric conductivity of the medium [S/m].
%
% Outputs:
% grad_phi =  Mx3 matrix of gradient at observation points [V/m].
%
M = size(obspnts,1); N = size(srcpnts,1);
grad_phi = zeros(M,3);
for i = 1:M
    for j = 1:N
        R = obspnts(i,:) - srcpnts(j,:);
        RP = R * dipmom(j,:)';
        a = R(1)*R(1) + R(2)*R(2) + R(3)*R(3);
        b = a^(5/2);
        grad_phi(i,:) = grad_phi(i,:) + (dipmom(j,:)*a - 3*RP.*R)/b;
    end
end
a = 1/(4*pi*sigma);
grad_phi = a*grad_phi;
end

% Formula implemented in a little bit different form... almost the same.
%
% function grad_phi = gradphiF( obspnts, srcpnts, dipmom, sigma )
% % gradphiF calculates the gradiet of electric potential at M observation 
% % points in a homogeneus and infinite medium by superposing the effects of 
% % N current dipoles.
% %
% % Input arguments:
% % obspnts  =  Mx3 matrix of observation points coordinates [m].
% % srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% % dipmom   =  Nx3 matrix of dipoles moments components [A*m].
% % sigma    =  electric conductivity of the medium [S/m].
% %
% % Outputs:
% % grad_phi =  Mx3 matrix of gradient at observation points [V/m].
% %
% M = size(obspnts,1); N = size(srcpnts,1);
% grad_phi = zeros(M,3);
% for i = 1:M
%     for j = 1:N
%         R = obspnts(i,:) - srcpnts(j,:);
%         RP = R * dipmom(j,:)';
%         a = R(1)*R(1) + R(2)*R(2) + R(3)*R(3);
%         b = a^(5/2);
%         a = a^(3/2);
%         grad_phi(i,:) = grad_phi(i,:) + dipmom(j,:)/a - 3*RP.*R/b;
%     end
% end
% a = 1/(4*pi*sigma);
% grad_phi = a*grad_phi;
% end