N = 9;
M = 5;
yf = @(x) exp(x);

x = pickpoints(-1,1,N);
y = yf(x);

ep = 1;
alpha = fzero(@(alpha)1-fzero(@(x)HermiteProd(M,(1+(2*ep/alpha)^2)^.25*alpha*x'),3),1);
beta = (1+(2*ep/alpha)^2)^.25;
delta2 = .5*alpha^2*(beta^2-1);

phi = rbfphi(1:M,x,ep,alpha);
X = zeros(N,M);
H = zeros(N,M);
L = zeros(M);
for k=1:M
    X(:,k) = (beta*alpha*x).^(k-1);
    H(:,k) = HermiteProd(k-1,beta*alpha*x);
    L(k,:) = [fliplr(HermitePoly(k-1)),zeros(1,M-k)];
end
D = diag(exp(-delta2*x.^2));
G = diag(sqrt(beta./(2.^(0:M-1).*gamma(1:M))));

% Sanity check 1
norm(phi - D*X*L'*G)


U = chol(phi'*phi);
S = diag(diag(U));
U = S\U;

Ut = chol(X'*D^2*X);
St = diag(diag(Ut));
Ut = St\Ut;

% Sanity check 2
norm(U'*S^2*U - G*L*Ut'*St^2*Ut*L'*G)


A = diag(diag(L'*G));

% Important result
norm(A*U - Ut*L'*G)