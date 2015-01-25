function DiffM = DifferenceMatrix(xcoord,zcoord)
% function DiffM = DifferenceMatrix(datacoord,centercoord)
% Borrowed from Fasshauer's "Meshfree Approximation Methods in Matlab"
%
% Forms the difference matrix of two sets of points in R
% (some fixed coordinate of point in R^s), i.e.,
% DM(j,k) = datacoord_j - centercoord_k .
%
% You cannot call this with d-dimensional points - you can only call this
% to compute the difference between a single dimension

[N,xd] = size(xcoord);
[M,zd] = size(zcoord);

% First returns an empty matrix, for degenerate data
% Then confirms that the user passed a column of data, not all dimensions
% If possible, computes the difference matrix
if (isempty(xcoord) || isempty(zcoord))
    DiffM = zeros(N,M);
elseif xd~=1 || zd~=1
    error('Only 1 dimension at a time can be differenced: x dim=%d, z dim=%d',xd,zd)
else
    DiffM = bsxfun(@minus,xcoord,zcoord');
end