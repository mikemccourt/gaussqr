% % ParameterPower14_1
% This example shows that the HS-SVD can be used to compute the power
% function stably, but that the power function computation is still subject
% to numerical cancelation (from subtracting nearby terms).  This is dealt
% with in the ParameterPower.m script.

% If this is running too slowly for you, reduce the number of epsilon
% values that are tested -- 123 is an arbitrary choice
epvec = [logspace(-2, -.8, 10), logspace(-.79, .6, 123), logspace(.6, 1, 10)];
alpha = 1;  % An arbitrary decision, we could search for a good choice
K = @(x, z, ep) exp(-(ep * DistanceMatrix(x, z)).^2);  % the Gaussian

% Consider 10 uniform samples of cosine
figure

% Create the data and choose some test points at which to estimate the
% interpolation accuracy and the norm of the power function
n = 10;
x = pickpoints(-1, 1, n);
xx = pickpoints(-1, 1, 100);
yf = @(x) cos(pi * x);
y = yf(x);
yy = yf(xx);

powervec = zeros(size(epvec));
errorvec = zeros(size(epvec));
k = 1;
for ep=epvec 
    Keval = K(xx, x, ep);

    % We use the GQR object to hide some of the computations
    GQR = gqr_solveprep(0, x, ep, alpha);
    Psi = GQR.stored_phi1 + GQR.stored_phi2 * GQR.CbarT;
    
    % We evaluate the psi function at the test points
    Phieval = gqr_phi(GQR.Marr, xx, ep, alpha);
    Phieval1 = Phieval(:, 1:n);
    Phieval2 = Phieval(:, n+1:end);
    Psieval = Phieval1 + Phieval2 * GQR.CbarT;
    
    % Compute the power at the test locations
    % Here, we hard-code the 1 as the Gaussian with K(x, x) == 1
    % We also vectorize the inner products at all the locations
    % Power(x) = K(x, x) - k(x)^T   * inv(K)   * k(x)
    %          = K(x, x) - psi(x)^T * inv(Psi) * k(x)
    cardinal_functions = Psieval / Psi;
    power = sqrt(1 - sum(cardinal_functions .* Keval, 2));
    powervec(k) = norm(power);
    
    % Compute the interpolation accuracy at test locations
    yeval = cardinal_functions * y;
    error = yy - yeval;
    errorvec(k) = norm(error);
    k = k + 1;
end

yyaxis left
loglog(epvec, powervec)
yyaxis right
loglog(epvec, errorvec, '--r')
legend({'Power function', 'Interpolation error'}, 'Location', 'northwest')

% Consider 15 chebyshev samples of the bessel function
figure

n = 15;
x = pickpoints(-1, 1, n, 'cheb');
xx = pickpoints(-1, 1, 100);
yf = @(x) besselj(0, abs(4 * x));
y = yf(x);
yy = yf(xx);

powervec = zeros(size(epvec));
errorvec = zeros(size(epvec));
k = 1;
for ep=epvec 
    Keval = K(xx, x, ep);

    GQR = gqr_solveprep(0, x, ep, alpha);
    Psi = GQR.stored_phi1 + GQR.stored_phi2 * GQR.CbarT;
    
    Phieval = gqr_phi(GQR.Marr, xx, ep, alpha);
    Phieval1 = Phieval(:, 1:n);
    Phieval2 = Phieval(:, n+1:end);
    Psieval = Phieval1 + Phieval2 * GQR.CbarT;
    
    cardinal_functions = Psieval / Psi;
    power = sqrt(1 - sum(cardinal_functions .* Keval, 2));
    powervec(k) = norm(power);
    
    yeval = cardinal_functions * y;
    error = yy - yeval;
    errorvec(k) = norm(error);
    k = k + 1;
end

yyaxis left
loglog(epvec, powervec)
yyaxis right
loglog(epvec, errorvec, '--r')
legend({'Power function', 'Interpolation error'}, 'Location', 'northwest')