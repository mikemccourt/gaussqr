% ex21b.m
% We are now considering the time benefits of using the "analytic" form of
% the inverse of the eigenfunction interpolation matrix.
%
% The results from the original test showed that things worked as they were
% supposed to, although there may be some question as to how well things
% scale, and obviously some limitation regarding the potential to scale to
% really large numbers.  There's also some question as to how well this
% could deal with sampling increasing density within a fixed domain.  For
% right now, though, we'll just stick with it in the simple way.

% We'll just work in 2D at first
Nvec = 10 * 2.^(0:3);  % Careful with the required timings
tvec = zeros(size(Nvec));
tivec = zeros(size(Nvec));

ep = .51;
alpha = 2.1;
yf = @(x1,x2) cos((x1.^2+x2.^2)) + .3*x1 - .2*x2;

k = 1;
for N=Nvec
    x1d = eig(Jf(alpha,ep,N-1));
    [X1,X2] = meshgrid(x1d,x1d);
    Y = yf(X1,X2);
    x1 = X1(:); x2 = X2(:); y = Y(:);
    
    GQR = gqr_solveprep(-1,x1d,ep+eps,alpha,N);
    Phi = gqr_phi(GQR,x1d);
    clear Phi_big c
    tic
    Phi_big = kron(Phi,Phi);
    c = Phi_big\y;
    tvec(k) = toc;
    clear Phi_big_inv ci
    tic
    weights_inv = N*gqr_phi(N,x1d,ep+eps,alpha).^2;
    Phi_inv = bsxfun(@rdivide,Phi',weights_inv');
    Phi_big_inv = kron(Phi_inv,Phi_inv);
    ci = Phi_big_inv*y;
    tivec(k) = toc;
    k = k + 1;
end

loglog(Nvec,tvec,'b',Nvec,tivec,'r')
xlabel('points per dimension')
ylabel('time to solve system')
legend('direct','tensor')