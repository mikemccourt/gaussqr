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
rbfdx = @(r,dx,ep) -ep^2*exp(-r).*dx;

% This function allows you to evaluate the EDF
% Here, xe are the evaluation points, x are the observed locations
Fhat = @(xe,x) reshape(sum(all(repmat(x,[1,1,size(xe,1)])<=repmat(reshape(xe',[1,size(xe,2),size(xe,1)]),[size(x,1),1,1]),2),1),size(xe,1),1)/size(x,1);

% Defined the Generalized Pareto distribution of interest
% With these parameters, the distribution only exists on [0,2]
gp_k = -1/2;
gp_sigma = 1;
gp_theta = 0;

% Create the data locations for this problem
N = 700;
x = sort(icdf('gp',rand(N,1),gp_k,gp_sigma,gp_theta));

% Evaluate the EDF at the given data locations
y = Fhat(x,x);

% Plot the true CDF for comparison purposes
Nplot = 500;
xplot = pickpoints(gp_theta,gp_theta-gp_sigma/gp_k,Nplot);
h_cdf = figure;
plot(xplot,cdf('gp',xplot,gp_k,gp_sigma,gp_theta),'r','linewidth',2);

% Choose an RBF to work with
rbf = rbfM4;
rbfdx = rbfM4dx;

% Create the surrogate model
ep = .3;
mu = 1e-2;
K_cdf = rbf(DistanceMatrix(x,x,ep));
cdf_coef = (K_cdf+mu*eye(N))\y;
cdf_eval = @(xeval) rbf(DistanceMatrix(xeval,x,ep))*cdf_coef;
pdf_eval = @(xeval) rbfdx(DistanceMatrix(xeval,x,ep),DifferenceMatrix(xeval,x),ep)*cdf_coef;

% Evaluate and plot the surrogate CDF
cplot = cdf_eval(xplot);

h_edf = figure;
plot(x,y,'linewidth',2)
hold on
plot(xplot,cplot,'--','linewidth',2);
hold off
legend('EDF','Model','location','southeast')

% Evaluate and plot the surrogate PDF
h_pdf = figure;
pplot = pdf_eval(xplot);
plot(xplot,pdf('gp',xplot,gp_k,gp_sigma,gp_theta),'r','linewidth',2);
hold on
plot(xplot,pplot,'--','linewidth',2);
title('Surrogate PDF')
legend('PDF','Model')
hold off