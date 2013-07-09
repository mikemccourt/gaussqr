% ex18_gqr.m
% This should compute the likelihood function for a set of given data
% We are interested in looking at the relationship between this likelihood
% value and the error
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.STORED_PHI_FOR_EVALUATION = 1;

epvec = logspace(-2,1,31);

N = 15; %number of points
NN = 200;
x = pickpoints(-1,1,N,'cheb'); %vector containing cheb. pts.
%yf = @(x) x+1./(1+x.^2); 
% fstring = 'y(x) = x + 1/(1+x^2)';
yf = @(x) x.^3-3*x.^2+2*x+1; %function is a polynomial
fstring = 'y(x) = x^3-3x^2+2x+1';

y = yf(x);
xx = pickpoints(-1,1,NN);
yy = yf(xx);
alpha = 1;
lamratio = 1e-12;
pinvtol = 1e-11;

errvec = [];
detvec = [];
mvec   = [];
lvec   = [];
derrvec = [];
ddetvec = [];
dmvec   = [];
dlvec   = [];

rbf = @(e,r) (exp(-e*r).^2);
DM = DistanceMatrix(x,x);
EM = DistanceMatrix(xx,x);

% Note that yPhi and yPsi are computed with a truncated SVD of Phi1 and Psi
% respectively.  The tolerance for this is pinvtol and can be set above.
k = 1;
for ep=epvec
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);



