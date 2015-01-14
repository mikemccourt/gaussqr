% KernelsCSRBF.m
% This is a brief example of computing with sparse distance matrics and the
% compactly supported kernels that exist as a result.

% These functions are used to manipulate the sparse matrices without
% causing them to be dense matrices
% sppc - sparse plus constant (although the arguments are in reverse order
sppc = @(const,Mat) spfun(@(x)const + x,Mat);
spd  = @(Mat1,Mat2) Mat1.*spfun(@(x)1./x,Mat2);

% Definition of the compactly supported kernel
% Note that the shape parameter must be handled within the DistanceMatrix
% computation, thus the radial kernel is just a single argument
% Also note that, because the distance matrix will be sparse, the spfun
% computation in needed to avoid log(0) terms from occurring and to avoid
% 1+sqrt(1+r.^2) from introducing nonzero values throughout the matrix
rbf = @(r) sppc(1,2*r.^2).*sqrt(sppc(1,-r.^2)) + ...
           3*r.^2.*spfun(@log,r./sppc(1,sppc(1+eps,-r.^2)));