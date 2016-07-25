% ex21.m
% This is an experiment involving Gaussian quadrature for eigenfunctions

% First, confirm the formula for standard Gauss-Hermite quadrature
yf = @(x) (cos(3 * x - .3)) .* exp(-x.^2);

% Known zeros of Hermite polynomials
% http://nvlpubs.nist.gov/nistpubs/jres/048/2/V48.N02.A04.pdf
Hzeros = {
    [-sqrt(2), sqrt(2)], ... # Roots of the n=2 Hermite polynomial
    [-1.224744871391589, 0, 1.224744871391589], ...
    [-1.650680123885785, -.524647623275290, .524647623275290, 1.650680123885785], ...
    [-2.020182870456086, -0.958572464613819, 0, 0.958572464613819, 2.020182870456086], ...
    [-2.350604973674492, -1.335849074013697, -.436077411927617, .436077411927617, 1.335849074013697, 2.350604973674492], ...
    [-2.651961350835233, -1.673561628767471, -.816287882858965, 0, .816287882858965, 1.673561628767471, 2.651961350835233], ...
    [-2.930637420257244, -1.981656756695843, -1.157193712446780, -.381186990207322, .381186990207322, 1.157193712446780, 1.981656756695843, 2.930637420257244], ...
    [-3.190993201781528, -2.266580584531843, -1.468553289216668, -.723651018752838, 0, .723651018752838, 1.468553289216668, 2.266580584531843, 3.190993201781528], ...
    [-3.436159118837738, -2.532731674232790, -1.756683849299882, -1.0366108297895114, -.342901327223705, .342901327223705, 1.0366108297895114, 1.756683849299882, 2.532731674232790, 3.436159118837738]
};

fprintf('True integral is %e\n',integral(yf,-Inf,Inf))

for k=1:length(Hzeros)
    xi = Hzeros{k}';
    n = k + 1;
    wi = 2^(n-1) * factorial(n) * sqrt(pi) ./ (n^2 * HermiteProd(n-1, xi).^2);
    appx_int = wi' * (yf(xi) .* exp(xi.^2));
    fprintf('%d\t%e\n',n,appx_int)
end

% Thought: using the zeros of the Hermite polynomial, which we should be
% able to compute to arbitrary accuracy using the eigenvalues of the J
% matrix, can we write the eigenfunctions as products of factors which
% would provide a more stable way to compute them.  It also could provide a
% faster methodology for stably computing things like the determinant -
% using just the highest values of the polynomial because it can be
% computed in log space.

% What is going on with the eigenvalue recurrence matrix?
% First thing to consider is that if we set ep=0 & alpha=1 we recover the
% Hermite polynomials (I think ... I need to consider the weight function)
ep = 0;
alpha = 1;
beta = @(alpha,ep) (1 + (2 * ep / alpha)^2)^.25;
Jf = @(alpha,ep,n) sqrt(diag(1:n, -1) + diag(1:n, 1)) / (alpha*beta(alpha,ep)*sqrt(2));

% Using this matrix, we can determine the roots of the Hermite polynomial
J7 = Jf(alpha,ep,6);
format long
these_should_match = [eig(J7), Hzeros{6}']
format short

% Part of what we want to accomplish here is to find a matrix W so that,
% given a Phi matrix evaluated at the roots, we have Phi'*W*Phi = I.
% I have done this before, so I know it is possible, I just can't remember
% how I did it last time.  I wrote down the answer, but not exactly how I
% got to it.
N = 6;
x = Hzeros{N-1}';
GQR = gqr_solveprep(-1,x,ep+eps,alpha,N);
hermiteweights = 2^(N-1) * factorial(N) * sqrt(pi) ./ (N^2 * HermiteProd(N-1, x).^2);
W = diag(hermiteweights / sqrt(pi));
Phi = gqr_phi(GQR, x);
fprintf('Error in inverse for ep=%g, alpha=%g is %e\n',ep,alpha,norm(Phi'*W*Phi - eye(N), 'fro'))

% Can I replicate these results with alpha ~= 1?
N = 6;
ep = 0;
alpha = 3.1;
GQR = gqr_solveprep(-1,x,ep+eps,alpha,N);
phiNroots = eig(Jf(alpha,ep,N-1)); % Definitely functions of alpha
% Can define these weights in terms of the eigenfunctions
hermiteweights = 2^(N-1) * factorial(N) * sqrt(pi) ./ (N^2 * HermiteProd(N-1, phiNroots * alpha).^2);
W = diag(hermiteweights / sqrt(pi));
Phi = gqr_phi(GQR, phiNroots);
fprintf('Error in inverse for ep=%g, alpha=%g is %e\n',ep,alpha,norm(Phi'*W*Phi - eye(N), 'fro'))

% Okay, now can I use ep > 0?
N = 6;
ep = .5;
alpha = 3.1;
GQR = gqr_solveprep(-1,x,ep+eps,alpha,N);
phiNroots = eig(Jf(alpha,ep,N-1)); % Also functions of ep
hermiteweights = 2^(N-1) * factorial(N) * sqrt(pi) ./ (N^2 * (gqr_phi(N,phiNroots,ep,alpha) / gammaN(N,alpha,ep)).^2);
W = diag(hermiteweights / sqrt(pi) / beta(alpha,ep));
W = diag(1./(N*gqr_phi(N,phiNroots,ep,alpha).^2)); % This is the correct way for the problem now
Phi = gqr_phi(GQR, phiNroots);
fprintf('Error in inverse for ep=%g, alpha=%g is %e\n',ep,alpha,norm(Phi'*W*Phi - eye(N), 'fro'))

% Can we use this to do eigenfunction approximation on a grid?
N = 20;
ep = .5;
alpha = 3.1;
x = eig(Jf(alpha,ep,N-1));
yf = @(x) cos(pi*x);
y = yf(x);
GQR = gqr_solveprep(-1,x,ep+eps,alpha,N);
Phi = gqr_phi(GQR,x);
weights = 1./(N*gqr_phi(N,x,ep,alpha).^2);
c = Phi\y;
ci = Phi'*(weights.*y);
fprintf('Applying the inverse has %g difference\n',norm(c-ci))

% What about for a 2D tensor grid?
N1 = 30;
N2 = 29;
ep1 = .5;
ep2 = .6;
alpha1 = 3.1;
alpha2 = 1.9;
x1_1d = eig(Jf(alpha1,ep1,N1-1));
x2_1d = eig(Jf(alpha2,ep2,N2-1));
[X1,X2] = meshgrid(x1_1d,x2_1d);
yf = @(x1,x2) cos((x1.^2+x2.^2)) + .3*x1 - .2*x2;
Y = yf(X1,X2);
x1 = X1(:);
x2 = X2(:);
y = Y(:);

GQR1 = gqr_solveprep(-1,x1_1d,ep1+eps,alpha1,N1);
GQR2 = gqr_solveprep(-1,x2_1d,ep2+eps,alpha2,N2);
Phi1 = gqr_phi(GQR1,x1_1d);
Phi2 = gqr_phi(GQR2,x2_1d);
weights1 = 1./(N1*gqr_phi(N1,x1_1d,ep1+eps,alpha1).^2);
weights2 = 1./(N2*gqr_phi(N2,x2_1d,ep2+eps,alpha2).^2);
Phi_big = kron(Phi1,Phi2);
c = Phi_big\y;
Phi_big_inv = kron(bsxfun(@times,Phi1',weights1'),bsxfun(@times,Phi2',weights2'));
ci = Phi_big_inv*y;
fprintf('Applying a 2D tensor inverse has %g difference\n',norm(c-ci))

surf(X1,X2,Y)
xp1_1d = pickpoints(min(x1_1d),max(x1_1d),50);
xp2_1d = pickpoints(min(x2_1d),max(x2_1d),48);
Phieval1 = gqr_phi(GQR1,xp1_1d);
Phieval2 = gqr_phi(GQR2,xp2_1d);
yp = kron(Phieval1,Phieval2)*ci;
[XP1,XP2] = meshgrid(xp1_1d,xp2_1d);
YP = reshape(yp,[48,50]);
surf(XP1,XP2,YP)