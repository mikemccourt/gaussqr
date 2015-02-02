% StatFitVolcano
% This example loads the Volcano data set from the main R example database
% https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/volcano.html
% We can access this from the GaussQR data directory
% Because we have so much data, we are going to use a C2 Wendland kernel
% We will allow it to be anisotropic

h_waitbar = waitbar(0,'Initializing Volcano data');

% Load the data in from the repository as needed
% This data comes with sizes Nx and Ny where length(x) = Nx*Ny
gqr_downloaddata('volcano_data.mat')
load volcano_data
x = locations;
h = heights;
N = Nx*Ny;

% If we wanted to plot just the data, we could execute the following
% X = reshape(x(:,1),Ny,Nx);
% Y = reshape(x(:,2),Ny,Nx);
% H = reshape(h,Ny,Nx);
% surf(X,Y,H),view([-211 28])

% Define xeval, some locations at which to test the power function
xeval = pick2Dpoints(0,1,10);

% Define the PD anisotropic C2 Wendland kernel in 2D
% Define it in its dense form, for passing to DistanceMatrix
rbf = @(r) max(1-r,0).^4.*(4*r+1);

% Consider a range of epsilon values to compute the kriging variance
epvec = [.01,.025,.05,.1,.25,.5,1:20];
% epvec = linspace(.5,1.5,24);

% Choose a cutoff point, beyond which we work with dense matrices
dense_cutoff = .2;

% Loop through the epsilon values
gwvec = zeros(size(epvec));
mplevec = zeros(size(epvec));
densevec = zeros(size(epvec));
DMtimevec = zeros(size(epvec));
solvetimevec = zeros(size(epvec));
k = 1;
for ep=epvec
    waitbar(k/length(epvec),h_waitbar,sprintf('\\epsilon=%g',ep));
    % Compute the kernel matrix
    [K,DMtimevec(k)] = DistanceMatrix(x,x,ep,rbf);
    densevec(k) = length(find(K))/N^2;
%     K = rbf(DistanceMatrix(x,x,ep));
%     densevec(k) = 1;
    if densevec(k)>dense_cutoff
        K = full(K);
        p = 1:N;
    else
        p = symamd(K);
    end
    
    tic
    % Factor the kernel matrix: Lp*Lp' = Kp
    Kp = K(p,p);
    Lp = chol(Kp,'lower');
    
    % Compute the solution to the system
    hp = h(p);
    ctemp = Lp'\(Lp\hp);
    c = zeros(N,1);c(p) = ctemp;
    
    % Compute the Hilbert space norm
    hsnorm = sqrt(c'*h);
    
    % Compute the log(det(K))
    logdetK = 2*sum(log(full(diag(Lp))));
    
    % Compute the power function
    if densevec(k)>dense_cutoff
        tic
        Keval = rbf(DistanceMatrix(xeval,x,ep));
        evaltime = toc;
    else
        [Keval,evaltime] = DistanceMatrix(xeval,x,ep,rbf);
    end
    pow = norm(sqrt(1 - sum(full(Keval(:,p)/Lp').^2,2)),'inf');
    
    solvetimevec(k) = toc - evaltime;
    DMtimevec(k) = DMtimevec(k) + evaltime;
    gwvec(k) = hsnorm*pow;
    mplevec(k) = 2*N*log(hsnorm) + logdetK;
    k = k + 1;
end

waitbar(1,h_waitbar,'Plotting')

h_yy = figure;
[AX,h1,h2] = plotyy(epvec,gwvec,epvec,mplevec,'loglog','semilogx');
set(get(AX(1),'xlabel'),'string','$\varepsilon$','interpreter','latex')
set(get(AX(1),'ylabel'),'string','$C_{GW}(\varepsilon,\infty)$','interpreter','latex')
set(get(AX(2),'ylabel'),'string','$C_{MPLE}(\varepsilon)$','interpreter','latex')
set(h1,'linewidth',2,'linestyle','--')
set(h2,'linewidth',2)
legend([h1 h2],'GW criterion','MPLE criterion','location','northwest')
set(AX,'xlim',[min(epvec),max(epvec)])
% Makes the axes black, for the book
% set(AX,{'ycolor'},{'k';'k'})

h_dense = figure;
semilogx(epvec,densevec,'linewidth',2)
xlabel('$\varepsilon$','interpreter','latex')
ylabel('density')
xlim([min(epvec),max(epvec)])

h_time = figure;
loglog(densevec,DMtimevec,'linewidth',2)
hold on
loglog(densevec,solvetimevec,'--','linewidth',2)
loglog([.2 .2],[.1 100],':','linewidth',2)
hold off
xlim([.01,1])
xlabel('density')
ylabel('time')
legend('distance matrix','C_{GW},C_{LW} evaluation','sparse-dense transition','location','northwest')

close(h_waitbar)