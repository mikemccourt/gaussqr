% ex20b - the NETNA example (for real this time)
% Trying to study the role of the t parameter on both the actual
% cross-validation residual and the cheaper residual

% Choose the collocation points on the surface of the sphere
gqr_downloaddata('sphereMDpts_data.mat')
load sphereMDpts_data
x = sphereMDpts{3};
M = size(x,1);

% Create some fake data
yf = @(x) -3*cos(7*x(:,1)-2*x(:,3)) + 2*exp(-2*x(:,2).^2);
y = yf(x);

% Prepare the kernel centers (probably a better strategy)
N = 250;
tvec = linspace(.001,7,20);  % Possible inflation values to consider
zbase = x(1:N,:);

% Choose the kernel to be the Green's function
rbf = @(r) 1./(4*pi*r);

% Consider all possible inflation parameters
j = 1;
cvvec_true = zeros(size(tvec));
cvvec_cheap = zeros(size(tvec));
h_waitbar = waitbar(0,'Initializing');
for t=tvec
    z = (1 + t)*zbase;
    K = rbf(DistanceMatrix(x,z));
    A = pinv(K);
    c = A*y;
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
        
        Ktt = Kswap(1:M-1,1:N-1);
        Kvt = Kswap(M,1:N-1);
        Avv = Aswap(N,M);
        yt = yswap(1:M-1);
        yv = yswap(M);
        cv = cswap(N);
        
        warning('off','MATLAB:rankDeficientMatrix')
        true_cv = yv - Kvt*(Ktt\yt);
        warning('on','MATLAB:rankDeficientMatrix')
        cheap_cv = cv/Avv;
        
%         fprintf('%d\t%f\t%f\n',k,true_cv,cheap_cv)
        cvvec_true(j) = cvvec_true(j) + abs(true_cv);
        cvvec_cheap(j) = cvvec_cheap(j) + abs(cheap_cv);
    end
    progress = floor(100*j/length(tvec))/100;
    waitbar(progress,h_waitbar,sprintf('Computing, t=%g',t))
    j = j + 1;
end

waitbar(100,h_waitbar,sprintf('Plotting'))
% Plot the results
h_joint = figure;

yyaxis right
h_true = semilogy(tvec,cvvec_true,'linewidth',2);
xlabel('t','fontsize',14)
ylabel('cv','fontsize',14)

yyaxis left
h_cheap = semilogy(tvec,cvvec_cheap,'linewidth',2);
ylabel('metric value','fontsize',14)

legend([h_true,h_cheap],{'True CV','Cheap CV'}, 'location', 'north', 'fontsize', 14)
close(h_waitbar)