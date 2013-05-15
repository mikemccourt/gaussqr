function pnts = SphereSurfHaltonPoints( npnts, R )
% SphereSurfHaltonPoints generates a set of points randomly distributed 
% over a sphere surface making use of a 2D Halton sequence.
%
% Input arguments:
% r      =  sphere radius.
% npnts  =  number of points (if scalar) or 
%           indexing to take the desired points of the sequence (if vector).
%
% Output:
% pnts   =  npnts x 3 matrix of points coordinates.
% 
% Required function:
%   (haltonseq.m) as an alternative to the MATLAB built-in haltonset
%

if not(isscalar(R))
    error('R must be a scalar')
end

hs = haltonset(2); % Load Halton set

if isscalar(npnts)
    % Take the first npnts points of the sequence
    A = net(hs,npnts);
elseif isvector(npnts)
    % Take points of the sequence specified in the vector npnts
    A = hs(npnts,:); 
else
    error('npnts must be a scalar or a vector')
end

% A = haltonseq(npnts,2); % This is an alternative to MATLAB haltonset
theta = 2*pi*A(:,1);
phi = acos(2*A(:,2)-1);
pnts = [ R .* sin(phi) .* cos(theta), ...
         R .* sin(phi) .* sin(theta), ...
         R .* cos(phi)];
end