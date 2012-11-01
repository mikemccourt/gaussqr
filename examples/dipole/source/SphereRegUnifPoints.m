function [ pnts ] = SphereRegUnifPoints( Dist, R )
%SPHEREREGUNIFPOINTS creates a uniform distribution of points inside a
% a sphere of given radius R starting from the desired distance between 
% points Dist.
%
% Inputs:
% Dist          =    desired distance between points;
% R             =    sphere radius.
%
% Outputs:
% pnts          =    matrix of points' coordinates.

[x, y, z] = meshgrid(-R:Dist:R, -R:Dist:R, -R:Dist:R);
ind = find(x(:).^2+y(:).^2+z(:).^2 < R^2);
n=length(ind);
pnts = zeros(n,3);
for i=1:n
  pnts(i,:) = [x(ind(i)), y(ind(i)), z(ind(i))];
end
end

