function val = StatFit1_param_func(pstrat,ep,DM,rbf,y)
% function val = StatFit1_param_func(pstrat,ep,DM,rbf,y)
%
%     StatFit1_param_func(1,ep,DM,rbf,y)
% This function evaluates a version of the % likelihood of the data having
% been observed for the given shape parameter.
% This is used in a minimization algorithm to find the maximum likelihood
% estimator, e.g., fminbnd(@(e)StatFit1_MLE_func(1,e,DM,rbf,y).1,10);
% This is just the basic version of the function, not using any
% Hilbert-Schmidt strategy
%
%     StatFit1_param_func(2,ep,DM,rbf,y)
% This function evaluates the LOOCV residual.
% This is used in a minimization algorithm to find the maximum likelihood
% estimator, e.g., fminbnd(@(e)StatFit1_MLE_func(2,e,DM,rbf,y).1,10);
% This is just the basic version of the function, not using any
% Hilbert-Schmidt strategy

K = rbf(ep,DM);

if pstrat==1
    N = length(y);
    [U,S,V] = svd(K);
    Mdist = y'*(V*diag(0*(diag(S)<1e-14)+(1./diag(S)).*(diag(S)>1e-14))*U'*y);
    logdetK = sum(log(diag(S)+eps));
    val = N*log(Mdist) + logdetK;
elseif pstrat==2
    invK = pinv(K);
    val = norm((invK*y)./diag(invK));
end