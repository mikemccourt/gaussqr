function B = HomSphereInduction ( srcpnts, dipmom, obspnts )
% HomSphereInduction calculates the magnetic induction at M points in
% the space surrounding a homogeneous sphere due to N current dipoles 
% making use of the analytical formula given in:
% Sarvas, "Basic Mathematical and Electromagnetic Concepts of the 
% Biomagnetic Inverse Problem", Phys. Med. Biol., 32(1) (1987), pp. 11-22.
%
% Input arguments:
% srcpnts  =  Nx3 matrix of source dipoles coordinates [m].
% dipmom   =  Nx3 matrix of dipole moments components [A*m].
% obspnts  =  Mx3 matrix of observation points coordinates [m].
%
% Output:
% B        =  Mx3 matrix of magnetic induction components [T].
%
M = size(obspnts,1); N = size(srcpnts,1);
B = zeros(M,3);
for i = 1:M
	lungr2 = obspnts(i,1)*obspnts(i,1) + ...
             obspnts(i,2)*obspnts(i,2) + ...
             obspnts(i,3)*obspnts(i,3);
    lungr = sqrt(lungr2);
    for j = 1:N
        a = obspnts(i,:) - srcpnts(j,:);
        lunga = a(1)*a(1) + a(2)*a(2) + a(3)*a(3);
        lunga = sqrt(lunga);
        b = a*obspnts(i,:)';
        c = srcpnts(j,:)*obspnts(i,:)';
        F = lunga*(lungr*lunga+lungr2-c);
        F2 = F*F;
        DF = (lunga^2/lungr+b/lunga+2*lunga+2*lungr)*obspnts(i,:) - ...
             (lunga+2*lungr+b/lunga)*srcpnts(j,:);
        d = [srcpnts(j,3)*dipmom(j,2) - srcpnts(j,2)*dipmom(j,3), ...
             srcpnts(j,1)*dipmom(j,3) - srcpnts(j,3)*dipmom(j,1), ...
             srcpnts(j,2)*dipmom(j,1) - srcpnts(j,1)*dipmom(j,2)];
        e = d*obspnts(i,:)';
        B(i) = B(i) + (F*d-e*DF)/F2;
    end
end
B = 1e-07*B;
end