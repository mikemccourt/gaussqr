function p = SphereSurfGoldPoints(n, R)
%SphereSurfGoldPoints generates a set of evenly distributed points on a 
% sphere of given radius using the method of the Golden Section Spiral.
% Based on Python code by P. Boucher (http://www.xsi-blog.com/archives/115)
%
% Inputs:
% n        =    desired number of points on the surface;
% R        =    sphere radius.
%
% Outputs:
% pnts     =    nx3 matrix of points' coordinates.
%
inc = pi * (3-sqrt(5));
off = 2/n;
k = 0:n-1;
y = k * off - 1 + (off/2);
r = sqrt(1 - (y.^2));
phi = k * inc;
p = R*[cos(phi).*r; y ; sin(phi).*r]';
end