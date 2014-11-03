% function DM = DistanceMatrix(dsites,ctrs)
% Forms the distance matrix of two sets of points in R^s,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
% Input
%   dsites: Mxs matrix representing a set of M data sites in R^s
%              (i.e., each row contains one s-dimensional point)
%   ctrs:   Nxs matrix representing a set of N centers in R^s
%              (one center per row)
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (x1 - x2)^2 = x1.^2 + x2.^2 - 2*x1*x2;
function DM = DistanceMatrix(dsites,ctrs)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

[M,dM] = size(dsites); [N,dN] = size(ctrs);

if dM~=dN
    error('data points in %dD but centers in %dD',dM,dN)
else
    d = dM;
end

% Allow user to pass [] and get back []
% In case no points fit some condition which is considered
if (M==0 || N==0)
    DM = zeros(M,N);
else % Otherwise, compute distance matrix
    T1 = sum(dsites.*dsites,2);
    T2 = -2*dsites*ctrs';
    T3 = (sum(ctrs.*ctrs,2))';
    DM = sqrt(T1(:,ones(N,1)) + T2 + T3(ones(M,1),:));

    % It is possible that the method above can yield negative distances
    %    which is obviously not good
    % If that happens, it is corrected below
    if any(imag(DM(:)))
        if alertuser
            warning('Data can yield negative distances because of cancelation')
        end
        DM = zeros(M,N);
        % Accumulate sum of squares of coordinate differences
        for l=1:d
            DM = DM + (repmat(dsites(:,l),1,N)-repmat(ctrs(:,l)',M,1)).^2;
        end
        DM = sqrt(DM);
    end
end