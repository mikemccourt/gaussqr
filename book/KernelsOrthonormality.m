% KernelsOrthonormality
% This script demonstrates the shortcomings of bsxfun and the need to
% compute with loops from time to time

% Define the Chebyshev eigenfunctions and weight function
phi = @(n,x) bsxfun(@times,sqrt(2-(n==0)),cos(bsxfun(@times,n,acos(x))));
rho = @(x) 1./(pi*sqrt(1-x.^2));

% Define the orthonormality testing function
ortheval = @(n1,n2) integral(@(x)...
            bsxfun(@times,phi(n1,x).*phi(n2,x),rho(x)),-1,1);

% Define a range of eigenfunctions indices to test for orthonormality
ntest = [0,2,3:5,8];

% Causes an error
try
    orthmat = bsxfun(@(n1,n2)ortheval(n1,n2),ntest',ntest);
catch err
    if strcmp(err.identifier,'MATLAB:dimagree')
        fprintf('bsxfun calls ortheval with different sized n1 and n2\n')
    end
end

% Less elegant, but without errors
% Create a cell array with all the indices required
% Then compute orthogonality for each element of that cell array
[N1,N2] = meshgrid(ntest,ntest);
ncell = num2cell([N1(:),N2(:)],2);
orthvec = cellfun(@(nvals)ortheval(nvals(1),nvals(2)),ncell);
disp(reshape(orthvec,length(ntest),length(ntest)))