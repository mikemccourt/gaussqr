% ex16
% This example demonstrates a strategy for pivoting good eigenfunctions
% into your HS-SVD basis.
%
% This is necessary in a setting where we are encountering ill-conditioning
% on the way to a singular Phi_1 matrix because the ep->0 limit is a
% polynomial limit and polynomials do not always produce a unique
% interpolant on some sets of points.
%
% I don't know exactly what polynomial limit you will reach with this
% strategy - although I could probably explain it constructively using the
% lexicographic ordering that we use for the eigenfunctions.
%
% This example uses 6 points in 2D evenly spaced on the unit circle; this
% will produce a singular matrix when interpolating with a degree 2
% polynomial because it is a conic section (I think).
%
% To run this example you need the GaussQR library.  I can't immediately
% think of a way to demonstrate this without the library.
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.NORM_TYPE = inf;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

% First create the data locations
t = pickpoints(-pi,pi,7);t = t(1:end-1);
x = [cos(t),sin(t)];

% Try to do the polynomial interpolation
polyint_func = @(x) [ones(size(x,1),1),x(:,1),x(:,2),x(:,1).^2,x(:,1).*x(:,2),x(:,2).^2];
polyint_mat = polyint_func(x);
if rank(polyint_mat)<6
    fprintf('Polynomial interpolation not unique\n')
end

% Create the potential GaussQR interpolation object, GQR
ep = .1;
gqr_alpha = 1;
GQR = gqr_solveprep(-1,x,ep,gqr_alpha);
Phi = gqr_phi(GQR,x);

% Note that this object has the same difficulties
%
% This happens to be true for all epsilon because these points all have the
% same delta^2 scaling term (since they all lie on a circle)
%
% For points off a circle but still on a conic section this would only be a
% problem for small epsilon, but it would still become a problem
if rank(Phi(:,1:6))<6
    fprintf('First 6 Gaussian eigenfunctions unsuitable, cond(Phi1)=%g\n',cond(Phi(:,1:6)))
end

% Try and pivot out eigenfunctions that contribute nothing
%
% Initialize with the first eigenfunction and keep adding more until we
% reach the full rank of the matrix
%
% This creates a list of indices to remove and then eliminates them from
% the set of total indices available (as determined by the user in rbfsetup
% as well as by the ep and gqr_alpha values)
curr_rank = 1;
ind_to_remove = [];
Phi1 = [Phi(:,1),zeros(6,5)];
next_to_fill = 2;
next_to_test = 2;
while curr_rank<6
    Phi1(:,next_to_fill) = Phi(:,next_to_test);
    new_rank = rank(Phi1);
    if new_rank==curr_rank
        fprintf('Eigenfunction [%d,%d] not helping\n',GQR.Marr(:,next_to_test))
        ind_to_remove = [ind_to_remove,next_to_test];
    else
        next_to_fill = next_to_fill + 1;
    end
    curr_rank = new_rank;
    next_to_test = next_to_test + 1;
end
% Eliminate from consideration the indices we didn't like
total_ind = 1:size(GQR.Marr,2);
good_ind = setdiff(total_ind,ind_to_remove);
% Reset the GQR object with only the eigenfunctions we like
GQR.Marr = GQR.Marr(:,good_ind);

% Check to make sure Phi1 is invertible now
Phi = gqr_phi(GQR,x);
Phi1 = Phi(:,1:6);
fprintf('After pivoting, cond(Phi1)=%g\n\n',cond(Phi1))

% Create the HS-SVD basis to confirm equality with the Gaussian interpolant
% in the standard basis, which we know exists even with the points on a
% conic section
%
% Note that, in general, bsxfun should be used for applying lambda; here it
% doesn't matter though so the more "Linear Algebra"-obvious approach is
% used for clarity
Phi2 = Phi(:,7:end);
lamvec = GQR.eig(GQR.Marr);
Lambda1 = diag(lamvec(1:6));
Lambda2 = diag(lamvec(7:end));
CorrectorT = Lambda2*(Phi2'/Phi1')/Lambda1;

% Create some data to interpolate to check against RBF-Direct
yf = @(x) cos(x(:,1)+x(:,2));
y = yf(x);
xeval = pick2Dpoints(-1,1,15);
yeval = yf(xeval);

% Evaluate the interpolant with RBF-Direct
rbf = @(e,r) exp(-(e*r).^2);
K = rbf(ep,DistanceMatrix(x,x));
Keval = rbf(ep,DistanceMatrix(xeval,x));
yeval_direct = Keval*(K\y);

% Evaluate the interpolant with the HS-SVD basis
Psi = Phi1 + Phi2*CorrectorT;
Psieval = gqr_phi(GQR,xeval)*[eye(6);CorrectorT];
yeval_hssvd = Psieval*(Psi\y);

% Check that both the interpolants are good and that they match each other
fprintf('For Gaussians with ep=%g\n',ep)
fprintf('Error in RBF-Direct: %g\n',errcompute(yeval_direct,yeval))
fprintf('Error in HSSVD basis: %g\n',errcompute(yeval_hssvd,yeval))
fprintf('Max difference between solutions: %g\n',norm(yeval_direct-yeval_hssvd,'inf'))
fprintf('Max norm of function: %g\n',norm(yeval,'inf'))
fprintf('Smaller ep gives closer matching between solutions\n')