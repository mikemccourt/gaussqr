function BF = BF( obspnts,srcpnts,dipmom )
% BF calculates the magnetic (induction) field at M observation points in a
% homogeneus and infinite medium by superposing the effects of N current 
% dipoles. 
% The relative magnetic permeability of the medium is assumed to be 1.
%
% Input arguments:
% obspnts  =  Mx3 matrix of observation points coordinates [m].
% srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% dipmom   =  Nx3 matrix of dipoles moments components [A*m].
%
% Output:
% B        =  Mx3 matrix of magnetic (induction) field at observation 
%             points [T].
%

% Check input data
N = size(srcpnts); P = size(dipmom);
if not(isequal(N,P))
    error('srcpnts and dipmom dimensions must match')
end

M = size(obspnts,1);
BF = zeros(M,3);
for i = 1:M
    R = bsxfun(@minus,obspnts(i,:),srcpnts);
    R_norm = R(:,1).*R(:,1)+R(:,2).*R(:,2)+R(:,3).*R(:,3);
    R_norm = R_norm.^(3/2);
    d = [dipmom(:,2).*R(:,3) - dipmom(:,3).*R(:,2), ...
         dipmom(:,3).*R(:,1) - dipmom(:,1).*R(:,3), ...
         dipmom(:,1).*R(:,2) - dipmom(:,2).*R(:,1)];
    d = bsxfun(@rdivide,d,R_norm);
    BF(i,:) = sum(d,1);
end
BF = 1e-07*BF;
end