%-------------------------------------------------------------------------%
%
% File: PuweightAniso(epoints,elocpts,puctrs,radii)
%
% Goal: script that constructs the PU weights
%
% Inputs:   epoints:     evaluation points
%           elocpts:     evaluation points on each subdomain
%           puctrs:      centre of the subdomain
%           radii:       semi-axes of the ellipses
%
% Outputs:  pu:          PU weights
%
% Calls on: WeightAniso
%
%-------------------------------------------------------------------------%
function pu = PuweightAniso(epoints,elocpts,puctrs,radii)
Np = size(puctrs,1);
[phi] = WeightAniso(epoints,elocpts,puctrs,radii);
% Compute the sums 
s = sum(phi,2);
for i = 1:Np
    loc = elocpts(i).ind;
    if isempty(loc)
        s(loc) = [];
    end
    pu(i).w = phi(loc,i)./s(loc);
end