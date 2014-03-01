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
d1_truesol = @(x,t) 1./(B*sqrt(t)).*(log(x/K)+(r+B^2/2)*t);
d2_truesol = @(x,t) d1_truesol(x,t) - B*sqrt(t);
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
steps = 3;
dt = T/steps;
tvec = dt:dt:T;
N_pts_1D = 20;
N_pts = N_pts_1D^steps;

% We evaluate our asset at each of these paths
switch MC_scheme
    case 1
        S_vals = S_0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*randn(N,1));
        MC_ans = mean(payout(S_vals));
    case 2
        % Choose low-discrepancy points in (steps) dimensions
        % Only in 3D right now
        x_1D = pickpoints(0.01,.99,N_pts_1D);
        X_1D = repmat(x_1D',N_pts_1D^2,1);
        Y_1D = repmat(x_1D',N_pts_1D,N_pts_1D);
        Z_1D = repmat(x_1D,N_pts_1D^2,1);
        % Apply the inverse CDF to draw instead from N(0,I)
        Phi_X = norminv([X_1D(:),Y_1D(:),Z_1D],0,1);
        % Covariance matrix must be formed and factored
        Sigma = min(repmat(tvec,steps,1),repmat(tvec',1,steps));
        A = chol(Sigma);
        integrand = (A'*Phi_X')';
        S_vals = integrand(:,steps);
        MC_ans = mean(payout(S_vals));
end

% We evaluate the discrete mean of those outcomes
% This is then plugged into our payout function
% We display the error in distance to the true solution
MC_err = abs(MC_ans-C_truesol(S_0,T));
fprintf('N=%d\tMC=%g\terr=%g\n',N,MC_ans,MC_err)