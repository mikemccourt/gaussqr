% IBBKernelConvergence
% This studies the convergence behavior of the IBB kernel for small values
% of epsilon where ill-conditioning can occur.  For larger values of beta,
% the smoothness parameter, the HSSVD basis is employed.

% Pick a function to be studied
f = @(x) 5.0887^6.*(max(0,x-.0567)).^6.*(max(0,(1-.0567)-x)).^6;

% Define the study to be run: parameters and data
ep = 1;  betavec = 1:5;  Nvec = 2.^(3:8);

% Choose test data
Neval = 397;
xeval = linspace(0,1,Neval)';  yeval = f(xeval);

% Loop through beta value
errvec = zeros(length(betavec),length(Nvec));
errorForBetas = zeros(length(betavec),Neval);
k = 1; j = 1;
for N=Nvec
    % Create the data for the desired N value
    x = pickpoints(0,1,N+2);  x = x(2:end-1);  y = f(x);
    
    j = 1;
    for beta=betavec
        % Use the closed form if stable enough
        if beta==1 || beta==2
            K = cell2mat(arrayfun(@(z) ...
                  IBBkernel(x,z,1,ep,beta),x','UniformOutput',0))';
            Keval = cell2mat(arrayfun(@(z) ...
                  IBBkernel(x,z,1,ep,beta),xeval','UniformOutput',0))';
            c = K\y;
            seval = Keval*c;
        else
            seval = HSSVD_IBBSolve(ep,beta,x,y,xeval);
        end
        errvec(j,k) = errcompute(seval,yeval);
        errorForBetas(beta,:) = yeval-seval;
        j = j + 1;
    end
    k = k + 1;
end

% Plot RMS error as beta and N vary
h_conv = figure;  loglog(Nvec,errvec,'linewidth',2);
% Plot individual error
h_errs = figure;  title('Error as \beta varies');
for beta = betavec;
    errorFig = subplot(ceil(length(betavec)/2),2,beta);
    semilogy(xeval,abs(errorForBetas(beta,:)));
end
ylim([1e-18 1e-3]);  set(gca,'YTick',[1e-15 1e-10 1e-5])