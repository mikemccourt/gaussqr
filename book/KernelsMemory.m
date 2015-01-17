% KernelsMemory
% This program demonstrates the use of try/catch statements to handle
% memory allocation issues, and how summing the system in pieces can allow
% for a reduced memory burden

% Define the points under consideration
N = 100;
x = pickpoints(-1,1,N);

% Define a necessary length for the series; this is problem dependent
M = 1000000;

% Define the phi functions and lambda vector
phi = @(n,x) bsxfun(@times,sqrt(2-(n==0)),cos(bsxfun(@times,n,acos(x))));
lam = @(n) (n==0)*.6 + (n>0).*.4./(n.^2*pi^2/6 + eps);

% Try to allocate enough space for the kernel matrix
% If this fails, there is no hope
try
    Kmat = zeros(N);
catch err
    if strcmp(err.identifier,'MATLAB:nomem')
        error('Insufficient memory for interpolation matrix')
    else
        rethrow(err);
    end
end

% Set up a loop to try to allocate memory, potentially fail, and then
% allocate less memory and build the kernel in pieces
more_computing_needed = 1;
allocation_size = M;
current_narr = 1:M;
bound_by_M = @(vecstart,vecrange) unique(min(vecstart-1+(1:vecrange),M));
while more_computing_needed
    try
        Phi = phi(current_narr,x);
        lamvec = lam(current_narr);
        Philam = bsxfun(@times,lamvec,Phi);
    catch err
        if any(strcmp(err.identifier,...
              {'MATLAB:nomem','MATLAB:array:SizeLimitExceeded'}))
            allocation_size = floor(allocation_size/5);
            current_narr = bound_by_M(current_narr(1),allocation_size);
            continue
        else
            rethrow(err);
        end
    end
    Kmat = Kmat + Philam*Phi';
    
    if any(current_narr~=M)
        current_narr = bound_by_M(current_narr(end)+1,allocation_size);
    else
        more_computing_needed = 0;
    end
end