%-------------------------------------------------------------------------%
%
% File: WeightAniso(epoints,elocpts,puctrs,radii)
%
% Goal: script that evaluates Wendland's functions
%
% Inputs:   epoints:     evaluation points
%           elocpts:     evaluation points on each subdomain
%           puctrs:      centre of the subdomain
%           radii:       semi-axes of the ellipses
%
% Outputs:  phi:         Wendland's functions
%
%-------------------------------------------------------------------------%
function [phi] = WeightAniso(epoints,elocpts,puctrs,radii)
[Np,M] = size(puctrs);
for i = 1:Np
    r = 0;
    for k = 1:M
        d(:,k) = epoints(elocpts(i).ind,k)-puctrs(i,k);
        r = r + (radii(i,k)).^2.*d(:,k).^2;
    end
    phi(elocpts(i).ind,i) = (4*sqrt(r)+1).*(max(0,1-sqrt(r))).^4;
    clear d
end










