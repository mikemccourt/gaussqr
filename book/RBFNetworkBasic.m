% RBFNetworkBasic
% This first example considers the effect of the mu parameter on data that
% has been fixed and a fixed epsilon

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Testing data
N = 50;
M = 15;
yf = @(x) (1-4*x+32*x.^2).*exp(-16*x.^2);
rbf = @(e,r) exp(-(e*r).^2);

% Pick points to evaluate the function at
% Add some error to the data
x = pickpoints(-1,1,N-2,'rand');
x = [x;-1;1];
noise = .2;
y = yf(x) + noise*randn(N,1);

% Choose points at which to center the basis functions, as needed
z = pickpoints(-1,1,M,'halton');

% Pick a range of Tikhonov Regularization parameters
muvec = [logspace(-10,-5,5),logspace(-4.9,5,50)];

% For plotting purposes
xeval = pickpoints(-1,1,300);
yeval = yf(xeval);

% Pick a shape parameter and form the kernel matrix
% The relevant computations are done with the SVD
% This is only for stability and can be changed for larger sizes
ep = 8;
K_fit = rbf(ep,DistanceMatrix(x,z));
K_predict = rbf(ep,DistanceMatrix(xeval,z));
[U,S,V] = svd(K_fit,0);Sv = diag(S);

% Conduct the loop over the regularization values
errvec = zeros(size(muvec));
loovec = zeros(size(muvec));
gcvvec = zeros(size(muvec));
k = 1;
for mu=muvec
    % Solve for the newtork weights, here with the standard basis
    w = V*((U'*y)./(Sv+mu./Sv));

    % Evaluate predictions on the test/plotting points
    % Check the error
    ypred = K_predict*w;
    errvec(k) = errcompute(ypred,yeval);

    % Find the projection matrix, used to evaluated the residual
    P = eye(N) - bsxfun(@rdivide,U,(1+mu./Sv.^2)')*U';
    Py = P*y;

    % Evaluate the parameterization schemes
    % The projection matrix is needed for this as well
    loovec(k) = Py'*(Py./diag(P).^2)/N;
    gcvvec(k) = N*(Py'*Py)/trace(P)^2;
    k = k + 1;
end

[~,ie] = min(errvec);
[~,ig] = min(gcvvec);
[~,il] = min(loovec);
handles = [];
figure
handles(1) = loglog(muvec,errvec,'k','linewidth',3);
hold on
handles(2) = loglog(muvec,gcvvec,'--','linewidth',3);
handles(3) = loglog(muvec,loovec,'-.','color',[.7 .5 0],'linewidth',3);
loglog(muvec(ie),errvec(ie),'k+','linewidth',3,'markersize',16)
loglog(muvec(ig),gcvvec(ig),'x','linewidth',3,'markersize',16)
loglog(muvec(il),loovec(il),'x','color',[.7 .5 0],'linewidth',3,'markersize',16)
title(sprintf('ep=%g,N=%d,M=%d',ep,N,M))
xlabel('\mu')
legend(handles,'Error','GCV','LOOCV','location','northwest')
hold off

% Store the error (not the min) for the best GCV
muopt_gcv = muvec(ig);
ypred_gcv = K_predict*(V*((U'*y)./(Sv+muopt_gcv./Sv)));
gcv_min_err = errcompute(ypred_gcv,yeval);

% Store the error (not the min) for the best error
muopt_err = muvec(ie);
ypred_err = K_predict*(V*((U'*y)./(Sv+muopt_err./Sv)));
min_err = errcompute(ypred_err,yeval);

% Plot the results
figure
plot(x,y,'or')
hold on
plot(xeval,yeval,'r','linewidth',2)
plot(xeval,ypred_err,'k','linewidth',2)
plot(xeval,ypred_gcv,'--','linewidth',2)
hold off
ylim([-.5,2])
title(sprintf('ep=%g,mu=%g,N=%d,M=%d',ep,muopt_gcv,N,M))
legend('Data','True',...
    sprintf('Min Error err=%2.2g',min_err),...
    sprintf('Min GCV err=%2.2g',gcv_min_err),...
    'location','northeast')