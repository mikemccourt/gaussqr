% IBBKernelCardFuncs
% This function creates plots of some cardinal functions for the IBB and
% Gaussian kernel on uniformly spaced and Chebyshev points
% #texttt{Choose kernel data size and evaluation points}
N = 20;  xeval = pickpoints(0,1,300);
% #texttt{Define IBB kernel}
ep = 50;  beta = 20;  I = eye(N);
% #texttt{Choose the cardinal functions to plot}
cfuneval = [1 5 10 15 20];
% #texttt{Loop through different designs}
pts = {'even','cheb'};
for k=1:2
    % #texttt{Create the desired points}
    x = pickpoints(0,1,N+2,pts{k});  x = x(2:end-1);
    cardfuncs = HSSVD_IBBSolve(ep,beta,x,I,xeval);
    % To evaluate one column at a time ...
%     cardfuncs = cell2mat(cellfun(@(y)HSSVD_IBBSolve(ep,beta,x,y,xeval),num2cell(I,1),'UniformOutput',0));
    subplot(2,2,k), plot(xeval,cardfuncs(:,cfuneval),'linewidth',3)
end
%% Gaussian kernel
rbf = @(e,r) exp(-(e*r).^2);
ep = 5.75;                                  % to match scaling of K_{20,50}
cfuneval = cfuneval+1;          % small shift since Gaussians nonzero at boundary
for k=1:2
    x = pickpoints(0,1,N+2,pts{k});
    DM = DistanceMatrix(x,x);           % distance matrix for interpolation
    DMeval = DistanceMatrix(xeval,x);         %       and for evaluation
    K = rbf(ep,DM);                     % interpolation matrix
    Keval = rbf(ep,DMeval);             %   and for evaluation
    % Evaluate cardinal functions and plot them
    cardfuncs = Keval/K;
    subplot(2,2,k+2), plot(xeval,cardfuncs(:,cfuneval),'linewidth',3)
end