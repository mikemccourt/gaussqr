% FinanceMC.m
% This example experiments with Monte Carlo methods in finance
%
% Our motivating SPDE is 
%              dS = r*S*dt + sigma*S*dW
%   S is a random process
%     S(omega,t) is a random variable indexed at time t
%   r is the interest rate
%   sigma is the volatility
%   W is a Wiener process
%     dW is a Brownian motion
%
% We define a solution to this using Ito Calculus as
%         S(T) = S(0)*exp((r-sigma^2/2)*T + sigma*W(T))
%   T is the time until we can exercise our option
%   W(T) is a random variable 
%
% If we choose to apply our call option to this security, our payoff is
%         C(S) = max(0,S-K)*exp(-r*T)
%   The exp(-r*T) term accounts for risk-free interest discounting
%   K is the strike price
%
% We define the fair price of this option to be
%         E[C(S)]
%   E[X] is the expected value of X
%   We compute this with an integral

% K is the scaled exercise price, which we set as 1
K = 1;

% T is the exercise time, chosen to be 1
T = 1;

% d is the number of assets in the option
d = 1;

% sigma is the dxd volatility matrix (nonsingular)
% It has diagonal elements .3, off diagonal elements .05
sigma = .3*eye(d) + .05*(ones(d)-eye(d));

% r is the risk-free interest rate
r = .05;

% tol is our computational tolerance
tol = 1e-4;

% payout(x) is the contract function which describes the payout
% This is also used to generate initial conditions
payout = @(S) max(0,sum(S-K,2)/d)*exp(-r*T);

% We require an initial condition
S_0 = 2;

% This is the computed expectation E[C(S)]
%   t is the time when the option can be exercised, which is T
d1_truesol = @(x,t) 1./(sigma*sqrt(t)).*(log(x/K)+(r+sigma^2/2)*t);
d2_truesol = @(x,t) d1_truesol(x,t) - sigma*sqrt(t);
C_truesol = @(x,t) normcdf(d1_truesol(x,t)).*x - K*normcdf(d2_truesol(x,t)).*exp(-r*t);

% We choose a certain number of random paths
% N = 10000;

% We choose an integration technique
%   MC_scheme = 1 is for normal Monte Carlo
%             = 2 is for quasi Monte Carlo
MC_scheme = 2;

% Some options are only appropriate for Quasi Monte Carlo
%   steps is the number of time steps to take
%   N_pts_1D is the number of points per dimension for uniform points
steps = 1;
halton_setup = haltonset(steps);
dt = T/steps;
tvec = dt:dt:T;

% We evaluate our random walks
switch MC_scheme
    case 1
        rand_walks = randn(N,1);
    case 2
        % Pick some Halton points and make sure not to use 0
        X = net(halton_setup,N+1);
        Phi_X = norminv(X(2:end,:),0,1);
        % Covariance matrix must be formed and factored
        Sigma = min(repmat(tvec,steps,1),repmat(tvec',1,steps));
        A = chol(Sigma);
        % Compute our quasi-random walk
        rand_walks = Phi_X*A;
end

% We perform our summation
S_vals = S_0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*rand_walks);
MC_ans = mean(payout(S_vals));

% We evaluate the discrete mean of those outcomes
% This is then plugged into our payout function
% We display the error in distance to the true solution
MC_err = abs(MC_ans-C_truesol(S_0,T));
fprintf('N=%d\tMC=%g\terr=%g\n',N,MC_ans,MC_err)