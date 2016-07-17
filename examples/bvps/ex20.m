% ex20 - the NETNA example (for real this time)
% This poses the appropriate PDE in the sphere and then solves the needed
% interpolation problem on the boundary.

% Choose the collocation points on the surface of the sphere
gqr_downloaddata('sphereMDpts_data.mat')
load sphereMDpts_data
x = sphereMDpts{1};
M = size(x,1);

% Create some fake data
yf = @(x) -3*cos(7*x(:,1)-2*x(:,3)) + 2*exp(-2*x(:,2).^2);
y = yf(x);

% Choose the kernel centers (probably a better strategy)
N = 25;
t = .1;  % Inflation parameter
z = (1 + t)*x(1:N,:);

% Choose the kernel to be the Green's function
rbf = @(r) 1./(4*pi*r);

% Go through and compute the the actual residual from the full
% cross-validation
cv_total = 0;
for k=1:N
    % Leave out the appropriate values
    z_train = z(setdiff(1:N,k),:);
    x_train = x(setdiff(1:M,k),:);
    y_train = y(setdiff(1:M,k));
    x_valid = x(k,:);
    y_valid = y(k);
    
    K_appx = rbf(DistanceMatrix(x_train,z_train));
    K_pred = rbf(DistanceMatrix(x_valid,z_train));
    y_pred = K_pred*(K_appx\y_train);
    
    this_cv = abs(y_valid - y_pred);
    fprintf('%d\t%f\n',k,this_cv)
    cv_total = cv_total + this_cv;
end
fprintf('True total CV: %f\n',cv_total)

% Now, we take a look at the actual pseudoinverse
K = rbf(DistanceMatrix(x,z));
A = pinv(K);
fprintf('If this is small the matrix is full rank %e\n',norm(A*K - eye(N)))

Kswap = K;
Kswap([3,6],:) = Kswap([6,3],:);
Kswap(:,[12,20]) = Kswap(:,[20,12]);
Aswap = A;
Aswap(:,[3,6]) = Aswap(:,[6,3]);
Aswap([12,20],:) = Aswap([20,12],:);
fprintf('If this is small then I understand row swaps %e\n',norm(Aswap*Kswap - eye(N)))

c = K\y;
fprintf('If this is small then I understand linear algebra %e\n',errcompute(A*y,c))

% Go through and compute the the faster residual using the full pinv
% I'm going to do this by manually pivoting the matrix rows/columns
% In reality the necessary crap can just be extracted
cv_total = 0;
for k=1:N
    Kswap = K;
    Kswap([k,M],:) = Kswap([M,k],:);
    Kswap(:,[k,N]) = Kswap(:,[N,k]);
    Aswap = A;
    Aswap([k,N],:) = Aswap([N,k],:);
    Aswap(:,[k,M]) = Aswap(:,[M,k]);
    yswap = y;
    yswap([k,M]) = yswap([M,k]);
    cswap = c;
    cswap([k,N]) = cswap([N,k]);
    % These identities seem to work:
    % Aswap*Kswap == eye(N)
    % Atv*Ktv + Avv*Kvv == 1
    % Att*Ktv + Avt*Kvv == 0
    % Atv*Ktt + Avv*Kvt == 0
    % Atv*yt + Avv*yv == cv
    % Att*yt + Avt*yv == ct
    % Also, yv - Kvt*(Ktt\yt) matches y_valid - y_pred from above
    
    Ktt = Kswap(1:M-1,1:N-1);
    Kvt = Kswap(M,1:N-1);
    Ktv = Kswap(1:M-1,N);
    Kvv = Kswap(M,N);
    Avv = Aswap(N,M);
    Atv = Aswap(N,1:M-1);
    Avt = Aswap(1:N-1,M);
    Att = Aswap(1:N-1,1:M-1);
    yt = yswap(1:M-1);
    yv = yswap(M);
    ct = cswap(1:N-1);
    cv = cswap(N);
    
    true_cv = yv - Kvt*(Ktt\yt);
    test_cv = yv + (1/Avv)*Atv*Ktt*(Ktt\yt);
    less_cv = yv + (1/Avv)*Atv*yt;
    base_cv = cv/Avv;
    
%     fprintf('\t%f\t%f\t%f\t%f\n',true_cv,test_cv,less_cv,base_cv);
    
    this_cv = abs(yv - Kvt*(Ktt\yt));
    fprintf('%d\t%f\n',k,this_cv)
    cv_total = cv_total + this_cv;
end
fprintf('Submatrix total CV: %f\n',cv_total)