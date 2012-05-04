N = 9;
M = 5;
yf = @(x) exp(x);

x = pickpoints(0,1,N,'cheb');
y = yf(x);

ep = 1;
alpha = fzero(@(alpha)1-fzero(@(x)HermiteProd(M,(1+(2*ep/alpha)^2)^.25*alpha*x'),3),1);
beta = (1+(2*ep/alpha)^2)^.25;
ba = beta*alpha;
delta2 = .5*alpha^2*(beta^2-1);

phi = gqr_phi(1:M,x,ep,alpha);
X = zeros(N,M);
H = zeros(N,M);
L = zeros(M);
for k=1:M
    X(:,k) = (beta*alpha*x).^(k-1);
    H(:,k) = HermiteProd(k-1,ba*x);
    L(k,:) = [fliplr(HermitePoly(k-1)),zeros(1,M-k)];
end
D = diag(exp(-delta2*x.^2));
G = diag(sqrt(beta./(2.^(0:M-1).*gamma(1:M))));
Z = [zeros(1,M-1),0;eye(M-1),zeros(M-1,1)];

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
A = diag(sqrt(beta*2.^(0:M-1)./gamma(1:M)));

% Important result
norm(A*U - Ut*L'*G)

T = diag(sqrt(1:M-1),1) + diag(sqrt(1:M-1),-1);
Ts = T/(sqrt(2)*ba);

% Intermediate result
norm( (phi'*diag(x)-Ts*phi') - [zeros(M-1,N);gqr_phi(M+1,x,ep,alpha)'*sqrt(M/2)/(beta*alpha)] )

Hp = fliplr(HermitePoly(M));

% Next important result
% This deals with the recurrence, and the Frobenius matrix
norm( X'*D*diag(x) - ( (L\(G\Ts)*G*L)*X'*D + [zeros(M-1,N);gqr_phi(M+1,x,ep,alpha)'*sqrt(gamma(M+1)/(beta*2^M))/(beta*alpha)] ) )
norm( L\(G\Ts)*G*L - [[zeros(M-1,1),eye(M-1)];-Hp(1:M)/Hp(M+1)]/(ba) )


Z = [zeros(1,M-1),0;eye(M-1),zeros(M-1,1)];
weM = G\(L'\Hp(1:M)'/Hp(M+1))*flipud(eye(M,1))'*L'*G;

% Confirm a substitution from the paper I am using analogously
norm(G\(L'\Z)*L'*G - ba*Ts' - weM)