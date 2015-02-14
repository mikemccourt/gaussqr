% FinanceMC.m
% This example experiments with Monte Carlo methods in finance
%
% Our motivating SPDE is 
%              dX = r*X*dt + sigma*X*dW
%   X is a random process
%     X_t is a random variable indexed at time t
%   r is the interest rate
%   sigma is the volatility
%   W is a Wiener process
%     dW is a Brownian motion
%
% We define a solution to this using Ito Calculus as
%         X_T = X_0*exp((r-sigma^2/2)*T + sigma*W(T))
%   T is the time when we can exercise our option
%   W(T) is a normal random variable with 0 mean and T variance
%
% If we choose to apply our call option to this security, our payoff is
%         C_e(X_T) = max(0,X-E)*exp(-r*T)
%   The exp(-r*T) term accounts for risk-free interest discounting
%   E is the strike price
%
% We define the fair price of this option to be
%         E[C(X)]
%   E[X] is the expected value of X
%   We compute this with an integral

% E is the exercise price
E = 1;

% T is the exercise time
T = 1;

% r is the risk-free interest rate
r = .05;

% s is the volatility of the asset
s = .3;

% mu is the drift of the asset, which we set to the interest rate
mu = r;

% payout(x) is the contract function which describes the payout
% This is also used to generate initial conditions
payout = @(x) exp(-r*T)*max(x-E,0);

% We require an initial condition
X_0 = 2;

% This is the computed expectation E[C(X)]
%   t is the time when the option can be exercised, which is T
d1_truesol = @(x,t) 1./(s*sqrt(t)).*(log(x/E)+(r+s^2/2)*t);
d2_truesol = @(x,t) d1_truesol(x,t) - s*sqrt(t);
C_truesol = @(x,t) normcdf(d1_truesol(x,t)).*x - E*normcdf(d2_truesol(x,t)).*exp(-r*t);

% We choose a certain number of random paths
Nvec = round(logspace(2,5,15));

MCvec = zeros(size(Nvec));
QMCvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Perform a standard Monte Carlo method for comparison
    W_T = sqrt(T)*randn(N,20);
    X_T = X_0*exp((r-s^2/2)*T + s*W_T);
    MCvec(k) = mean(abs(C_truesol(X_0,T) - mean(payout(X_T))));
    
    % Perform the quasi-Monte Carlo method
    u = pickpoints(0,1,N,'halt');
    Finvu = icdf('Lognormal',u,log(X_0)+(mu-s^2/2)*T,s*sqrt(T));
    QMCvec(k) = abs(C_truesol(X_0,T) - mean(payout(Finvu)));
    
    k = k + 1;
end

pMC = polyfit(log(Nvec),log(MCvec),1);
pQMC = polyfit(log(Nvec),log(QMCvec),1);

h = figure;
loglog(Nvec,MCvec,'--','linewidth',2)
hold on
% loglog(Nvec,1.5*exp(pMC(2))*Nvec.^pMC(1),':k')
loglog(Nvec,QMCvec,'linewidth',2)
% loglog(Nvec,.75*exp(pQMC(2))*Nvec.^pQMC(1),':k')
hold off
ylim([1e-4 1e-1])
xlabel('number of quadrature points')
ylabel('quadrature absolute error')
legend('standard Monte Carlo','quasi-Monte Carlo','location','southwest')