% HSSVDInterpExample
% This example deals with the Iterated Brownian Bridge kernels
% Because they do not have a closed form for large beta, but they may
% suffer from ill-conditioning, we use the HS-SVD basis
% We will plot the stable result, along with the unstable result
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = 2;

% Pick a function that satisfies the IBB boundary conditions
yf = @(x) .25^(-28)*max(x-.25,0).^14.*max(.75-x,0).^14;

% Set up some data points at which to sample
% Remember to dump the boundaries for the IBB kernel
N = 19;
x = pickpoints(0,1,N+2);  x = x(2:end-1);
y = yf(x);

% Set up some evaluation points
Neval = 300;
xeval = pickpoints(0,1,Neval);
yeval = yf(xeval);

% Choose IBB parameters for computation
ep = 1;
beta = 7;

% Define the eigenfunctions and eigenvalues
% lamfunc will return a vector if n is a vector
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lamfunc = @(b,e,n) ((pi*n).^2+e^2).^(-b);

% Determine how many eigenfunctions will be needed for accuracy
% This formula comes from the book
M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));
narr = 1:M;

% Evaluate phi for the necessary eigenfunctions
Phi = phifunc(narr,x);
Phi1 = Phi(:,1:N);
Phi2 = Phi(:,N+1:end);
Phieval1 = phifunc(narr(1:N),xeval);
Phieval2 = phifunc(narr(N+1:end),xeval);

% Evaluate lamvec for the necessary eigenvalues
lamvec = lamfunc(beta,ep,narr);
lamvec1 = lamvec(1:N);
lamvec2 = lamvec(N+1:end);

% Form the Rbar corrector matrix for the HS-SVD basis
Rbar = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');

% Compute the interpolant in the HS-SVD basis
Psi = Phi*[eye(N);Rbar];
Psieval = Phieval1 + Phieval2*Rbar;
seval = Psieval*(Psi\y);

% Create the standard basis kernel matrices
% Evaluate the interpolant using them
K = bsxfun(@times,lamvec,Phi)*Phi';
Keval = bsxfun(@times,lamvec1,Phieval1)*Phi1' + ...
        bsxfun(@times,lamvec2,Phieval2)*Phi2';
sevalstandard = Keval*(K\y);

% Display the ranks and errors of the approximations
standarderr = errcompute(sevalstandard,yeval);
hssvderr = errcompute(seval,yeval);
fprintf('rank(K) = %d, standard basis error = %2.1e\n',rank(K),standarderr)
fprintf('rank(Psi) = %d, HS-SVD basis error = %2.1e\n',rank(Psi),hssvderr)

% Plot the results on the same plot
plot(xeval,sevalstandard,'--','linewidth',2)
hold on
plot(xeval,seval,'linewidth',2)
hold off
legend('standard basis','hssvd basis','location','northeast')