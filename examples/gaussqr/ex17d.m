% ex17d
% This example approximates the posterior distribution for the epsilon
% given an approximation problem
%     f(x) = sin(x/2)-2*cos(x)+4*sin(pi*x)
%          N=30 points at the Chebyshev nodes
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 120;

rbf = @(e,r) exp(-(e*r).^2);
N = 40;
NN = 100;

epvec = logspace(-1,1,75);

alpha_base = 2;
yf = @(x) tanh(x);

x = pickpoints(-3,3,N,'cheb');
y = yf(x);
xx = pickpoints(-3,3,NN);
yy = yf(xx);

DM = DistanceMatrix(x,x);
EM = DistanceMatrix(xx,x);

% This is first used to create the error graph for comparison
errvec = [];
edvec = [];
k = 1;
for ep=epvec
    if ep<1
        alpha = alpha_base;
    else
        alpha = alpha_base/(6^log10(ep));
    end
    GQR = gqr_solve(x,y,ep,alpha);
    yp = gqr_eval(GQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    A = rbf(ep,DM);
    warning off
    yp = rbf(ep,EM)*(A\y);
    warning on
    edvec(k) = errcompute(yp,yy);
    
    fprintf('ep=%g\tM=%d\terr=%g\n',ep,length(GQR.Marr),log10(errvec(k)))
    k = k + 1;
end

loglog(epvec,edvec,'linewidth',2),hold on
loglog(epvec,errvec,'r','linewidth',3),hold off
xlabel('\epsilon')
ylabel('error')
legend('Direct method','Stable method','location','southeast')
ylim([1e-12,1e-1])

pause

alpha = 2;

% Define the prior distribution
galpha = 1.5;
gbeta = 2.5;
gammadist = @(e,a,b) 1/(b^a*gamma(a))*e.^(a-1).*exp(-e/b);

% Define the Bayesian parameters
walk_sd = sqrt(.002);
Nsteps = 1000;
Nburnin = 100;
rand('state',0);
randn('state',0);

% Initialize the problem
ep = 1;
Lhood = gqr_likelihood(x,y,ep,alpha);
Lprior = log(gammadist(ep,galpha,gbeta));
Lpost = Lhood + Lprior;
num_accepted = 0;

ephist = zeros(Nsteps,1);
for step=1:Nsteps
    ep_new = ep + walk_sd*randn;
    Lhood_new = gqr_likelihood(x,y,ep_new,alpha);
    Lprior_new = log(gammadist(ep_new,galpha,gbeta));
    Lpost_new = Lhood_new + Lprior_new;
    
    Lratio = Lpost_new - Lpost;
    Ldraw = log(rand);
    fprintf('ep_new=%g\t ratio=%g\t draw=%g\n',ep_new,exp(Lratio),exp(Ldraw))
    if Ldraw<Lratio
        ephist(step) = ep_new;
        ep = ep_new;
        Lpost = Lpost_new;
        fprintf('\tnew step taken, step=%d\tep=%g\n',step,ep)
        num_accepted = num_accepted + 1;
    else
        ephist(step) = ep;
    end
end

num_accepted/Nsteps