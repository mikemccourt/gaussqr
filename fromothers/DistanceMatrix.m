% function DM = DistanceMatrix(dsites,ctrs,epvec)
% Forms the distance matrix of two sets of points in R^d,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
% Input
%   dsites: Nxd matrix representing a set of N data sites in R^d
%              (i.e., each row contains one s-dimensional point)
%   ctrs:   Mxd matrix representing a set of M centers in R^d
%              (one center per row)
%   epvec:  Optional weighting row vector for anisotropic kernels
%              (one ep per column of dsites)
%              <default = vector of all ones>
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (a - b)^2 = a.^2 + b.^2 - 2*a*b;
function DM = DistanceMatrix(dsites,ctrs,epvec)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end

[N,dN] = size(dsites); [M,dM] = size(ctrs);

if dM~=dN
    error('data points in %dD but centers in %dD',dN,dM)
else
    d = dM;
end

if nargin==3
    if any(abs(real(epvec))~=epvec)
        error('epvec must only have positive, real values')
    end
else
    epvec = ones(1,d);
end

% Allow user to pass [] and get back []
% In case no points fit some condition which is considered
if (M==0 || N==0)
    DM = zeros(N,M);
else % Otherwise, compute distance matrix
    try
        DM = DistanceMatrix_COMPUTE(dsites,ctrs,epvec);
    catch err
        if strcmp(err.identifier,'MATLAB:nomem')
            error('You probably do not have enough memory for this distance matrix')
        else
            rethrow(err);
        end
    end
end

DM = sqrt(DM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fastest way to compute distance matrix
%%% Allows for anisotropy
%%% NOTE: Can produce negative values on the order of
%%%       sqrt((machine precision)*(dimension))
%%%       I may need to check for that eventually
%%% NOTE: In Matlab2015, this is causing a crash
%%%       with DM assignment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM = DistanceMatrix_COMPUTE(dsites,ctrs,epvec)
xe = bsxfun(@times,dsites,epvec);
ze = bsxfun(@times,ctrs,epvec);
sxe2 = sum(xe.^2,2);
sze2 = sum(ze.^2,2);

DM = bsxfun(@plus,sxe2,sze2') - 2*xe*ze';

end