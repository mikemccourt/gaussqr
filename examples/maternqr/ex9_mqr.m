% First test for HS-SVD for tensor product
% I will try to create a tensor kernel with 2 compact Materns that use
% (potentially) different shape parameters

% We need functions that define the eigenfunctions and eigenvalues of
% each of the kernels that make up the tensor kernel
%
% The standard format for the eigenfunctions will be
%        (parameters,indices,locations)
% The eigenfunctions should take 1D points as a column vector x and indices
% as a row vector n.  It should return a matrix where each row is for a
% specific x value and each column is a specific n value.
%
% For the eigenvalue evaluation, the standard format will be
%        (parameters,indices)
% The input vector n should be a row vector and it should return a row
% vector of eigenvalues.
%
% The p object can be anything that contains the parameters for your
% Hilbert-Schmidt series.  Obviously here there aren't many choices but
% this could be something very general.
ef1 = @(p,n,x) sqrt(2)*sin(x*pi*n);
ev1 = @(p,n) (n.^2*pi^2 + p(1)^2).^(-p(2));
ef2 = @(p,n,x) sqrt(2)*sin(x*pi*n);
ev2 = @(p,n) (n.^2*pi^2 + p(1)^2).^(-p(2));

% Define the compact Matern parameters
p1 = [20,9];
p2 = [50,8];

% Choose some points at which we want to evaluate our HS-SVD
x = pick2Dpoints(.05,.95,5,'halton');
N = size(x,1);

% Define the total number of extra eigenfunctions to add after N
% Normally we would take this from GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNCS
% but this is just a demonstration so I'm not worried about doing
% everything the correct way.
M = N + 100;

% I don't immediately have a smart way to adaptively choose the
% eigenvalues, so I'm going to just sample a bunch and dump most of them
n1_test = 1:floor(M/5);
n2_test = 1:floor(M/5);

% Create the eigenvalue Kronecker product
% Theoretically, we could actually use gqr_formMarr to create the test
% options ... as I think about that it might actually be better, unless
% there is a significant imbalance between dimensions.  Like if one
% dimension decays a lot faster than others then we will end up looking at
% a lot of eigenvalues we do not need.   More to think on.
%
% Create a selection of Marr possibilities corresponding to the Kronecker
% product eigenvalues above.  Then choose the Marr that matches what we
% have up there.  Doing this in higher dimensions is tougher ... so we'll
% need to think on that.
Marr_test = [kron(ones(size(n2_test)),n1_test);kron(n2_test,ones(size(n1_test)))];
lamvec1_test = ev1(p1,Marr_test(1,:));
lamvec2_test = ev2(p2,Marr_test(2,:));
lamvec_test = lamvec1_test.*lamvec2_test;
[lamsort,i] = sort(lamvec_test,'descend');

% Isolate the necessary indices
Marr = Marr_test(:,i(1:M));
lamvec = lamsort(1:M);
lamvec1 = lamvec(1:N);
lamvec2 = lamvec(N+1:M);
Lambda1 = diag(lamvec1);

% Now that we have the Marr we want, we could just proceed with the
% standard computational strategy, as I do below.  We could also save some
% computational cost in evaluating the eigenfunctions by taking advantage
% of the Kronecker product structure of the eigenvalues.  This probably
% doesn't matter very much, since the cost of forming Psi is still much
% greater than forming Phi, but it's a viable optimization.
Phi_K1 = ef1(p1,Marr(1,:),x(:,1));
Phi_K2 = ef2(p2,Marr(2,:),x(:,2));
Phi = Phi_K1.*Phi_K2;
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:end);

% Form the CbarT matrix
CbarT = (Phi2'/Phi1').*bsxfun(@rdivide,lamvec2',lamvec1);

% Choose some points at which to evaluate the kernel
xeval = pick2Dpoints(0,1,6,'rand');
Neval = size(xeval,1);
Phieval_K1 = ef1(p1,Marr(1,:),xeval(:,1));
Phieval_K2 = ef2(p2,Marr(2,:),xeval(:,2));
Phieval = Phieval_K1.*Phieval_K2;
Psieval = Phieval*[eye(N);CbarT];

% Confirm that this worked
K1_dir = ibb(xeval(:,1),x(:,1),p1(1),p1(2));
K2_dir = ibb(xeval(:,2),x(:,2),p2(1),p2(2));
K_dir = K1_dir.*K2_dir;

K_hssvd = Psieval*Lambda1*Phi1';

fprintf('Error in HS-SVD is %g\n',norm(K_dir-K_hssvd)/norm(K_dir))