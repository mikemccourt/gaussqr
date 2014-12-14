function Y = sphHarm(l,m,c1,c2,c3)
%SPHHARM Computes the normalized, real-valued, spherical harmonic of degree
%   L order M at a given set of locations on the sphere.
%
% NOTE: This computes using the harmonics on the scale
%       [1,2,3,4,...,2*l+1], not [-l,-l+1,...,-1,0,1,...,l-1,l]
%       This is contradictory to some references, and the original form of
%       this code, so beware
%
%   Y = sphHarm(L,M,LAM,TH) returns the degree L order M normalized
%   spherical harmonic at the points (LAM,TH) on the sphere expressed in
%   longitude-latitude coordinates (or azimuthal-elevation).  Here
%            -pi <= lam <= pi   is the longitude (azimuthal) coordinate and
%           -pi/2 <= th <= pi/2 is the latitude (elevation) coordinate.
%   
%   Y = sphHarm(L,M,X,Y,Z) returns the degree L order M normalized
%   spherical harmonic at the points (x,y,z) on the sphere expressed in
%   Cartesian coordinates.
%
%   Y = sphHarm(L,M,X) This is the equivalent of calling
%   sphHarm(L,M,X(:,1),X(:,2),X(:,3)), but this only makes sense if you are
%   computing sphHarm of a vector, and X, Y, Z are not matrices
%
%   Example 1: Spherical coordinates
%       [lam,th] = meshgrid(linspace(-pi,pi,81),linspace(-pi/2,pi/2,41));
%       f = sphHarm(6,0,lam,th) + sqrt(14/11)*sphHarm(6,5,lam,th);
%       surf(lam,th,f), shading interp;
%
%   Example 2: Cartesain coordinates
%       [x,y,z] = sphere(101);
%       f = sphHarm(6,0,x,y,z) + sqrt(14/11)*sphHarm(6,5,x,y,z);
%       surf(x,y,z,f), shading interp, axis equal;
%
% Author: Grady Wright, 2014
%
% This function was written by Grady Wright and is available as part of his
% rbfsphere package which can be downloaded at
%         http://math.boisestate.edu/~wright/montestigliano/
% PLEASE consider using that package if you are interested in computing
% with RBFs on manifolds.  GaussQR has no serious capability in that
% regard.
%
% Programmer's Note: This could be adapted to accept vector l and m without
% too much difficulty and evaluate the harmonics in a loop

if ~(numel(l)==1 && numel(m)==1)
    error('Only a single harmonic can be evaluated at once, for now.  numel(l)=%d, numel(m)=%d',l,m)
end

if ~(floor(abs(l))==l && floor(abs(m))==m)
    error('Unacceptable indices: l=%g, m=%g',l,m)
end

if l<0
    error('The lowest order spherical harmonic is 0, but l=%d',l)
end

if m<1 || m>2*l+1
    error('The degree m must be satisfy 1<=m<=2*l+1, but m=%d and l=%d',m,l)
end

% Shift the input degree from [1,2*l+1] to [-l,l] because that is how the
% code below expects to recieve the index
m = m-1-l;

switch nargin
    case 3 % Cartesian coordinates in a single vector used.
        [lam,th] = cart2sph(c1(:,1),c1(:,2),c1(:,3));
        clear c1
    case 4 % Spherical coordinates used.
        lam = c1;
        th = c2;
        clear c1 c2
    case 5 % Cartesian coordinates used.
        [lam,th] = cart2sph(c1,c2,c3);
        clear c1 c2 c3
    otherwise
        error('Unacceptable inputs: nargin=%d',nargin)
end

% Flatten and transpose th and lam so they work with the legendre function
sz = size(th); th = th(:)'; lam = lam(:)';

% Normalization
a = sqrt((2*l+1)/2/pi*factorial(l-abs(m))/factorial(l+abs(m))*(2-double(m==0)));
Y = legendre(l,sin(th));
% Get the right associated legendre function
Y = squeeze(Y(abs(m)+1,:,:));
% Determine if the cos or sin term should be added.
pos = abs(max(0,sign(m+1)));
% Compute the spherical harmonic
Y = (pos*cos(m*lam) + (1-pos)*sin(m*lam)).*(a*Y);
% Reshape so it is the same size as the th and lam that were passed in.
Y = reshape(Y,sz);

end