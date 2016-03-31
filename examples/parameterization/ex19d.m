% ex19b
% This is a 2D example for stable parametrization using both the MLE and
% the joint metric.
% We use the Gaussian in 2D and try to optimize for a single ep value.

d = 4;
yf = @(x) cos(sqrt(sum(bsxfun(@times,x,2.^(-(0:d-1)+1)).^2,2)));

% Set up some data points at which to sample
N = 25;
p = haltonset(d,'Skip',1);
x = net(p,N);
y = yf(x);

% Choose the locations at which to test the power function
Npf = 200;
p = haltonset(d,'Skip',N+1);
xpf = net(p,NN);

% Set up some evaluation points
NN = 150;
p = haltonset(d,'Skip',N+Npf+1);
xx = net(p,NN);
yy = yf(xx);

% Parameters for the individual dimensions
alpha = 1;
epvec = [logspace(-2,-1,20),logspace(-1,0,40),logspace(0,1,20)];

% The closed form of the radial Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);

dirvec = zeros(length(epvec),1);
dlivec = zeros(length(epvec),1);
djvec = zeros(length(epvec),1);

k = 1;
for ep=epvec    
    K = rbf(ep, DistanceMatrix(x, x));
    warning('off','MATLAB:nearlySingularMatrix')
    c = K\y;
    warning('on','MATLAB:nearlySingularMatrix')
    Keval = rbf(ep, DistanceMatrix(xx, x));
    dirvec(k) = errcompute(Keval*c,yy);
    
    logdetK = sum(log(svd(K)));
    log_mahaldist = log(abs(c'*y));
    dlivec(k) = N*log_mahaldist + logdetK;
    
    logdetKtilde = zeros(length(xpf),1);
    for m=1:length(xpf)
        xp = [x;xpf(m,:)];
        Kp = rbf(ep, DistanceMatrix(xp, xp));
        logdetKtilde(m) = sum(log(svd(Kp)));
    end
    log_PF = logdetKtilde - logdetK;
    KV = log_mahaldist + max(log_PF);
    djvec(k) = dlivec(k) + KV;
    
    fprintf('%d, condition of K %e\n',k,cond(K))
    k = k + 1;
end

% Plot the results
h_joint = figure;

yyaxis right
h_err = loglog(epvec,dirvec,'linewidth',2);
xlabel('\gamma','fontsize',14)
ylabel('max pointwise error','fontsize',14)
ylim([1e-2,1e1])
ax = gca;
ax.YTick = [1e-1, 1, 1e1];

yyaxis left
h_mle = semilogx(epvec,dlivec,'linewidth',3);
hold on
h_j = semilogx(epvec,djvec,'--','linewidth',3);
hold off
ylabel('metric value','fontsize',14)
ax = gca;
ax.YTick = [-50, 200];
ax.FontSize = 14;

legend([h_err,h_mle,h_j],{'Error','MLE Metric','DET Metric'}, 'location', 'north', 'fontsize', 14)