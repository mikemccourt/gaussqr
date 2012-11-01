function [ surf_pnts, int_pnts, Dist ] = SphereRegGoldPoints( Npnts_surf, R )
% SphereRegGoldPoints creates a (almost) regular distribution of points inside 
% a sphere of given radius R starting from the desired number Npnts_bnd of 
% points equally distributed on the circumference (boundary).
%
% Criterion:
% The distances between every interior point and the closest others are
% almost equal to the distance Dist between two near points on the surface.
%
% Inputs:
% Npnts_surf    =    desired number of points on the surface;
% R             =    sphere radius.
%
% Outputs:
% surf_pnts     =    matrix of surface points' coordinates;
% int_pnts      =    matrix of interior points' coordinates.
% Dist          =    distance between points (approx)
%
% Calls on: 
%   SphereSurfGoldPoints.m
%   DistanceMatrix.m

surf_pnts = SphereSurfGoldPoints(Npnts_surf, R);

DM = DistanceMatrix(surf_pnts, surf_pnts); 
Dist = mean(min(DM+eye(size(DM))));

NL = round( R/Dist ); % number of internal spheres ("layers")
r = R/NL;
a = R^2/Npnts_surf; % number of points per unit surface / 4*pi

Npnts_i = zeros(1,NL-1); R_i = zeros(1,NL-1);
for i = 1:(NL-1)
   R_i(i) = r * i; % radius of the i-th inner sphere
   Npnts_i(i) = round( R_i(i)^2/a );
end

int_pnts = zeros(sum(Npnts_i)+1,3); % +1 to include the origin

% Evaluation of internal points' coordinates
N = Npnts_i(1) + 1;
int_pnts(2:N,:) = SphereSurfGoldPoints(Npnts_i(1), R_i(1));
for i = 2:(NL-1)
    int_pnts(N+1:N+Npnts_i(i),:) = SphereSurfGoldPoints(Npnts_i(i), R_i(i));
    N = N + Npnts_i(i);
end

end