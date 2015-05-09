% function [DM,DM_time] = DistanceMatrix(dsites,ctrs,epvec,sp_opt)
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
%                         pass [] to trigger the default
%              If a single value is passed, it is used for all dimensions
%   sp_opt: Pass 1 to create a sparse distance matrix
%              <default = 0>
%           You may pass an rbf function as this argument, and it will
%           return the kernel matrix evaluated with rbf at the nonzero
%           locations rather than the distance matrix
%           Ex: sp_opt = @(r) 1-r;
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
%   time:   Time required for computing distance matrix
%
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (a - b)^2 = a.^2 + b.^2 - 2*a*b;
function [DM,DM_time] = DistanceMatrix(dsites,ctrs,epvec,sp_opt)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
statsavail = GAUSSQR_PARAMETERS.STATISTICS_TOOLBOX_AVAILABLE;

% Confirm the user passed something
if not(exist('dsites','var'))
    error('You passed nothing to this function')
else
    [N,d] = size(dsites);
end

% Allow for the user to pass one variable and match ctrs and dsites
% If the user passed something, confirm the size is acceptable
if not(exist('ctrs','var'))
    ctrs = dsites;
end
[M,dM] = size(ctrs);
if dM~=d
    error('data points in %dD but centers in %dD',d,dM)
end

% Set the default for epvec, if needed
% If it is passed, make sure it is something valid
% Allow the user to pass a single number for all dimensions
if not(exist('epvec','var')) || isempty(epvec)
    epvec = ones(1,d);
else
    if any(abs(real(epvec))~=epvec)
        error('epvec must only have positive, real values')
    end
    if length(epvec)==1
        epvec = epvec*ones(1,d);
    elseif any(size(epvec)~=[1,d])
        error('size(epvec)=[%d %d], but should be [1 %d]',size(epvec),d)
    end
end

% Set the default sparsity to none if not specified
% Also check to see if the user passed a radial kernel to be evaluated
% before forming the sparse distance matrix
if not(exist('sp_opt','var'))
    sp_opt = 0;
else
    if isa(sp_opt,'function_handle')
        sp_func = sp_opt;
        sp_opt = 1;
    else
        sp_func = @(x) x;
    end
    if sp_opt && ~statsavail
        error('Sparse distance matrix computation requires statistics toolbox')
    end
end


% Allow user to pass [] and get back []
% In case no points fit some condition which is considered
if (M==0 || N==0)
    DM = zeros(N,M);
else % Otherwise, compute distance matrix
    try
        % Scale the data here so epvec need not be passed
        tic
        xe = bsxfun(@times,dsites,epvec);
        ze = bsxfun(@times,ctrs,epvec);
        if sp_opt
            DM = DistanceMatrix_SPARSE(xe,ze,sp_func);
        else
            DM = DistanceMatrix_COMPUTE(xe,ze);
        end
        DM_time = toc;
    catch err
        if strcmp(err.identifier,'MATLAB:nomem')
            error('You probably do not have enough memory for this distance matrix')
        else
            rethrow(err);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fastest way to compute distance matrix
%%% Allows for anisotropy
%%% NOTE: Can produce negative values on the order of
%%%       sqrt((machine precision)*(dimension))
%%%       This is why the max(x,0) command is used here,
%%%       to avoid for complex square roots
%%% NOTE: In Matlab2015 and Windows 8.1, this is causing a crash
%%%       with bsxfun sometimes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM = DistanceMatrix_COMPUTE(xe,ze)
sxe2 = sum(xe.^2,2);
sze2 = sum(ze.^2,2);

DM = sqrt(max(bsxfun(@plus,sxe2,sze2') - 2*xe*ze',0));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute a sparse distance matrix
%%% This requires the statistics toolbox
%%% This makes the assumption that your CSRBF is of
%%%   the form max(1-r,0)
%%% If this is not the case, I don't know what to say
%%% sp_func is the rbf that the user wants to evaluate
%%% The default is no rbf evaluation at all
%%%   sp_func = @(x) x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM = DistanceMatrix_SPARSE(xe,ze,sp_func)
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

N = size(xe,1);
M = size(ze,1);

% Perform the distance computation and record how many points are good
[Cidx,Cdist] = rangesearch(xe,ze,1);
nnz = sum(cellfun(@length,Cidx));

% Create vectors of the rows, columns and distances
ivec = cell2mat(Cidx');
jvec = cell2mat(cellfun(@(x,k)k*ones(1,length(x)),Cidx',num2cell(1:M),'UniformOutput',false));
svec = cell2mat(Cdist');

% Add something nonzero to avoid the squeeze Matlab would otherwise apply
svec = svec + eps;

% Form the sparse matrix after evaluating any rbf that may have been passed
DM = sparse(ivec,jvec,sp_func(svec),N,M,nnz);

% Alert user if there is the matrix is not so dense
if alertuser && nnz/N^2>.333
    warning('Sparse matrix is %g dense; may be better to store as dense',nnz/N^2)
end
end