% SurrModelEDF1D
% This example considers data drawn from a Generalized Pareto distribution
% and computes the associated EDF.  Then it fits a CDF to that data and
% computes the PDF.


% Standardize the random results
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define the RBF we use for this problem
rbf = @(r) (1+r).*exp(-r);
rbfx = @(ep,r,dx) -ep^2*exp(-r).*dx;

% This function allows you to evaluate the EDF
% Here, xe are the evaluation points, x are the observed locations
Fhat = @(xe,x) reshape(sum(all(repmat(x,[1,1,size(xe,1)])<=repmat(reshape(xe',[1,size(xe,2),size(xe,1)]),[size(x,1),1,1]),2),1),size(xe,1),1)/size(x,1);

% Defined the Generalized Pareto distribution of interest
% With these parameters, the distribution only exists on [0,2]
gp_k = -1/2;
gp_sigma = 1;
gp_theta = 0;

% Create the data locations for this problem
N = 700;Ncuts = 5; % Make sure N/Ncuts is an integer
xall = icdf('gp',rand(N,1),gp_k,gp_sigma,gp_theta);
xcell = arrayfun(@(n)xall((n-1)*N/Ncuts+1:n*N/Ncuts,:),... 
                     (1:Ncuts)','UniformOutput',0);

% Evaluate the EDF at the given data locations
ycell = cellfun(@(x)Fhat(x,x),xcell,'UniformOutput',0);

% Organize into one vector
[x,isort] = sort(cell2mat(xcell));
yunsorted = cell2mat(ycell);
y = yunsorted(isort);

% Conduct the KS test for normality on the difference
ydiff = y - cdf('gp',x,gp_k,gp_sigma,gp_theta);
[h,p] = kstest((ydiff-mean(ydiff))/std(ydiff))

% Plot the true CDF for comparison purposes
Nplot = 500;
xplot = pickpoints(gp_theta,gp_theta-gp_sigma/gp_k,Nplot);
h_cdf = figure;
plot(x,y,'.b');
hold on
plot(xplot,cdf('gp',xplot,gp_k,gp_sigma,gp_theta),'r','linewidth',2);
hold off

% Create the surrogate model
ep = .4;
mu = 3e-3;
K_cdf = rbf(DistanceMatrix(x,x,ep));
cdf_coef = (K_cdf+mu*eye(N))\y;
cdf_eval = @(xeval) rbf(DistanceMatrix(xeval,x,ep))*cdf_coef;
pdf_eval = @(xeval) rbfx(ep,DistanceMatrix(xeval,x,ep),DifferenceMatrix(xeval,x))*cdf_coef;

% Evaluate and plot the surrogate CDF
cplot = cdf_eval(xplot);

h_edf = figure;
plot(x,y,'.','linewidth',2)
hold on
plot(xplot,cplot,'r','linewidth',2)
hold off
legend('EDF data','Model','location','southeast')

% Evaluate and plot the surrogate PDF
h_pdf = figure;
pplot = pdf_eval(xplot);
plot(xplot,pdf('gp',xplot,gp_k,gp_sigma,gp_theta),'r','linewidth',2)
hold on
plot(xplot,pplot,'--','linewidth',2)
legend('PDF','Model')
hold off