%-------------------------------------------------------------------------%
%
% File: DistanceMatrixAniso(x,z,epvec)
%
% Goal: script that compute anisostopic distances between x and z in R^M
%
% Inputs:   x:        NXM matrix representing a set of N points in R^M
%                     (i.e., each row contains one M-dimensional point)
%           z:        PXM matrix representing a set of P points in R^M
%                     (one point per row)
%           epvec:    row vector with M shape parameters 
%                     (one value per each space dimension)
%
% Outputs:  DM:       NXP matrix of anisotropic distances in Euclidean norm
%           
% Remark:   this function comes from the book:
%           [G.E. Fasshauer, M.J. McCourt, Kernel-based Approximation 
%           Methods using Matlab, World Scientific, Singapore, 2015]
%
%-------------------------------------------------------------------------%
function DM = DistanceMatrixAniso(x,z,epvec)
if not(exist('epvec','var'))
    epvec = ones(1,size(x,2));
end
xe = bsxfun(@times,x,epvec);
ze = bsxfun(@times,z,epvec);
sxe2 = sum(xe.^2,2);
sze2 = sum(ze.^2,2);
DM = sqrt(bsxfun(@plus,sxe2,sze2')-2*xe*ze');