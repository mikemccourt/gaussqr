% This runs time tests on the SVM solver with the full rank and low rank
% representations.  It also compares the cost of the sparse SVM evaluations
% to ignorantly performing a full summation, which would be the associated
% cost of using an RBF network.  There is also an option to study the cost
% of different parameterizations for solving the quadratic program.

% To allow for the low-rank expansion parameter to be set
global GAUSSQR_PARAMETERS

% Choose a range of parameters to test over, or fixed parameters if testing
% over something else
epvec = logspace(-2,2,31);
bcvec = logspace(-2,4,30);
ep = .01;
bc = 1;

% Choose whether or not to plot the output, which you may not want to do if
% running the file remotely or within a batch script
% Set to 'off' to prevent plotting, 'on' to plot
% You can see a non-visible figure with handle h using figure(h)
% If this is 'off', instead of displaying the plot it automatically saves
% the figure using gqr_savefig
plots_on = 'on';

% Set the low rank parameter if desired
low_rank = 0;
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .1;

% Choose which test you want to run
% 1 - Fixed parameterization, increasing size, low-rank vs. standard
% 2 - Range of SV amounts and associated summation cost
% 3 - Cost of solving optimization problem with various ep and bc values
test_opt = 3;

switch test_opt
    case 1
        % Create random training and test data
        train_N_vec = round(logspace(1,4,30));
        
        lowvec = zeros(length(train_N_vec),1);
        fullvec = zeros(length(train_N_vec),1);
        k = 1;
        h_waitbar = waitbar(0,'Initializing','Visible',plots_on);
        for train_N=train_N_vec
            [train_data,train_class] = SVM_setup(1,train_N,10);
            tic
            SVMlow = gqr_fitsvm(train_data,train_class,ep,bc,1);
            lowvec(k) = toc;
            tic
            SVMfull = gqr_fitsvm(train_data,train_class,ep,bc,0);
            fullvec(k) = toc;
            
            progress = floor(100*k/length(train_N_vec))/100;
            waitbar(progress,h_waitbar,sprintf('num points=%d, low rank time=%5.2f, full rank time=%5.2f',train_N,lowvec(k),fullvec(k)))
            if strcmp(plots_on,'off')
                fprintf('%5.2f%% complete, num points=%d, low rank time=%5.2f, full rank time=%5.2f\n',progress,train_N,lowvec(k),fullvec(k))
            end
            k = k + 1;
        end
        
        waitbar(1,h_waitbar,'Plotting')
        h = figure('Visible',plots_on);
        loglog(train_N_vec,lowvec,'linewidth',2)
        hold on
        loglog(train_N_vec,fullvec,'--','linewidth',2)
        xlabel('number of training points')
        ylabel('training time')
        legend('Low rank','Full rank','location','southeast')
        hold off
        if strcmp(plots_on,'off')
            gqr_savefig(h,sprintf('SVMtimetests%d',test_opt));
        end
        
        close(h_waitbar)
    case 2
    case 3
        % Create random training and test data
        [train_data,train_class] = SVM_setup(1,100,10);

        timemat = zeros(length(epvec),length(bcvec));
        k_ep = 1;
        h_waitbar = waitbar(0,'Initializing','Visible',plots_on);
        for ep=epvec
            % Run here to set up SVM persistent varaibles
            SVM = gqr_fitsvm(train_data,train_class,ep,bc,low_rank);
            k_bc = 1;
            for bc=bcvec
                tic
                SVM = gqr_fitsvm(train_data,train_class,ep,bc,low_rank);
                timemat(k_ep,k_bc) = toc;
                
                progress = floor(100*((k_ep-1)*length(bcvec)+k_bc)/(length(epvec)*length(bcvec)))/100;
                waitbar(progress,h_waitbar,sprintf('compute time=%5.2f, \\epsilon=%5.2f C=%5.2f',timemat(k_ep,k_bc),ep,bc))
                if strcmp(plots_on,'off')
                    fprintf('%5.2f%% complete, ep=%5.2f C=%5.2f, time=%5.2f\n',progress,ep,bc,timemat(k_ep,k_bc))
                end
                k_bc = k_bc + 1;
            end
            k_ep = k_ep + 1;
        end
        
        waitbar(1,h_waitbar,'Plotting')
        [E,B] = meshgrid(epvec,bcvec);
        
        h = figure('Visible',plots_on);
        h_ev = surf(E,B,timemat');
        set(h_ev,'edgecolor','none')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca,'ytick',[1e-2,1e1,1e4])
        xlabel('\epsilon')
        ylabel('C')
        zlabel('SVM training time')
        shading interp
        grid off
        view([-.7,1,1])
        if strcmp(plots_on,'off')
            gqr_savefig(h,sprintf('SVMtimetests%d',test_opt));
        end
        
        close(h_waitbar)
    otherwise
        error('Unacceptable test case test_opt=%d',test_opt)
end