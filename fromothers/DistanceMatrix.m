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
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (a - b)^2 = a.^2 + b.^2 - 2*a*b;
%
% If epvec is passed, instead compute
%         ||x-z|| = ep1^2*(x1-z1)^2 + ... + epd^2*(xd-zd)^2
function DM = DistanceMatrix(dsites,ctrs,epvec)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

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

% Initialize so that we can check if the data is unstable
negative_values = 0;

% Allow user to pass [] and get back []
% In case no points fit some condition which is considered
if (M==0 || N==0)
    DM = zeros(N,M);
elseif all(epvec==ones(1,d)) % Otherwise, compute distance matrix
    try
        [DM,negative_values] = DistanceMatrix_BASIC(dsites,ctrs);
    catch err
        if strcmp(err.identifier,'MATLAB:nomem')
            warning('Memory allocation error, attempting less demanding strategy')
            try
                DM = DistanceMatrix_ANISOTROPIC(dsites,ctrs,epvec);
            catch err_safer
                if strcmp(err_safer.identifier,'MATLAB:nomem')
                    error('Nope ... that didn''t work either.  You have not enough memory for this matrix.')
                else
                    rethrow(err_safer)
                end
            end
        else
            rethrow(err);
        end
    end
end

if alertuser && negative_values
    warning('Data can yield negative distances because of cancelation')
end

if negative_values || any(epvec~=ones(1,d))
    DM = DistanceMatrix_ANISOTROPIC(dsites,ctrs,epvec);
end

DM = sqrt(DM);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing the distance matrix in the fastest way
%%% Checks to see if numerical cancelation is a problem
%%% Note that this return r.^2, so the sqrt is considered externally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DM,negative_values] = DistanceMatrix_BASIC(dsites,ctrs)
N = size(dsites,1);
M = size(ctrs,1);

T1 = sum(dsites.*dsites,2);
T2 = -2*dsites*ctrs';
T3 = (sum(ctrs.*ctrs,2))';
DM = T1(:,ones(M,1)) + T2 + T3(ones(N,1),:);

negative_values = any(DM(:)<0);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Allows for epvec to be something other than all ones
%%% Also safe if numerical cancelation was an issue
%%% Note that this return r.^2, so the sqrt is considered externally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM = DistanceMatrix_ANISOTROPIC(dsites,ctrs,epvec)
[N,d] = size(dsites);
M = size(ctrs,1);
DM = zeros(N,M);

% Accumulate sum of squares of coordinate differences
% Each coordinate is weighted by the epvec
% Need to convert to bsxfun if it is available
for k=1:d
    DM = DM + epvec(k)^2*bsxfun(@minus,dsites(:,k),ctrs(:,k)').^2;
end

end