% FinanceEurCallPDE
% This example demonstrates the use of kernels to compute the solution to
% Black-Scholes type PDEs in only one dimension
%
% The problem is a parabolic problem,
%     du/dt = Lu
%     u(t,boundary) = bc(boundary)
%     u(0,x) = payout(x)
% on the domain 0<x<4K and 0<t<1
% K is the strike price, which implies payout max(x-K,0)
%
% Our differential operator for this problem has 3 pieces:
%      Lu = L1u + L2u + L3u
% L1u = r*x'*grad(u)
% L2u = (1/2)*x'*((S*S').*hess(u))*x
% L3u = -r*u
%
% The exact solution is the formula for the value of a European call option
%      C(x,t) = F_Z(d1)*u - F_Z(d2)*K*exp(-r*(T-t)) 
% where
%      F_Z(z) = P(Z<z) for Z~N(0,1)
%      d1(x) = 1/(S*sqrt(T-t))*(log(x/K)+(r+S^2/2)*(T-t))
%      d2(x) = d1(x) - S*sqrt(T-t)
% Our solution below uses t in place of T-t because we know the final value
% and not the initial value.  This means that when solving the problem we
% are really stepping backwards in time.

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define parameters of the PDE
K = 1;
T = 1;
S = .3;
r = .05;

% Define the payout, used as the initial condition
payout = @(x) max(x-K,0);

% Define the boundary conditions (works for x=0 and x=4)
bc = @(x,t) K*(4-exp(-r*t))*(x==4*K);

% Define the true solution
d1 = @(x,t) 1./(B*sqrt(t)).*(log(x/K)+(r+B^2/2)*t);
d2 = @(x,t) d1(x,t) - B*sqrt(t);
Ptrue = @(x,t) normcdf(d1(x,t)).*x - K*normcdf(d2(x,t)).*exp(-r*t);

% Define the possible kernels for this problem
rbfM2 = @(e,r) (1+e*r).*exp(-e*r);
rbfM2x = @(e,r,dx) -e^2*exp(-e*r).*dx;
rbfM2xx = @(e,r) e^2*exp(-e*r).*(e*r-1);
rbfM4 = @(e,r) (1+e*r+(e*r).^2/3).*exp(-e*r);
rbfM4x = @(e,r,dx) -e^2*exp(-e*r).*dx.*(1+e*r)/3;
rbfM4xx = @(e,r) e^2*exp(-e*r).*((e*r).^2-e*r-1)/3;
rbfM6 = @(e,r) (1+e*r+2/5*(e*r).^2+(e*r).^3/15).*exp(-e*r);
rbfM6x = @(e,r,dx) -e^2*exp(-e*r).*dx.*(3+3*e*r+(e*r).^2)/15;
rbfM6xx = @(e,r) e^2*exp(-e*r).*((e*r).^3-3*e*r-3)/15;
rbfG = @(e,r) exp(-(e*r).^2);
rbfGx = @(e,r,dx) -2*e^2*dx.*exp(-(e*r).^2);
rbfGxx = @(e,r) 2*e^2*(2*(e*r).^2-1).*exp(-(e*r).^2);

% Pick a specific kernel to solve with
% rbf = rbfM6;  rbfx = rbfM6x;  rbfxx = rbfM6xx;
rbf = rbfM4;  rbfx = rbfM4x;  rbfxx = rbfM4xx;
% rbf = rbfM2;  rbfx = rbfM2x;  rbfxx = rbfM2xx;
ep = 2;

% Choose a range of point values to consider for this problem
Nvec = 6:3:21;
pt_opt = 'cheb';

errvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Create points in both domains
    x1all = pickpoints(0,K,N,pt_opt);
    x1bc = x1all(x1all==0);  N1bc = length(x1bc);
    x1cc = x1all(x1all==K);  N1cc = length(x1cc);
    x1int = x1all(x1all~=0 & x1all~=K);  N1int = length(x1int);
    i1int = 1:N1int;
    i1bc = N1int+1:N1int+N1bc;
    i1cc = N1int+N1bc+1:N1int+N1bc+N1cc;
    x1 = [x1int;x1bc;x1cc];  N1 = length(x1);
    
    DM1 = DistanceMatrix(x1,x1);
    DM1int = DistanceMatrix(x1int,x1);
    DiffM1int = DifferenceMatrix(x1int,x1);
    DM1cc = DistanceMatrix(x1cc,x1);
    DiffM1cc = DifferenceMatrix(x1cc,x1);
    V1 = rbf(ep,DM1);
    Vx1int = rbfx(ep,DM1int,DiffM1int);
    Vxx1int = rbfxx(ep,DM1int);
    Vx1cc = rbfx(ep,DM1cc,DiffM1cc);
    Vxx1cc = rbfxx(ep,DM1cc);
    Vx1intV1inv = Vx1int/V1;
    Vxx1intV1inv = Vxx1int/V1;
    Vx1ccV1inv = Vx1cc/V1;
    
    x2all = pickpoints(K,4*K,3*N,pt_opt);
    x2bc = x2all(x2all==4*K);  N2bc = length(x2bc);
    x2cc = x2all(x2all==K);  N2cc = length(x2cc);
    x2int = x2all(x2all~=K & x2all~=4*K);  N2int = length(x2int);
    i2int = 1:N2int;
    i2bc = N2int+1:N2int+N2bc;
    i2cc = N2int+N2bc+1:N2int+N2bc+N2cc;
    x2 = [x2int;x2bc;x2cc];   N2 = length(x2);
    
    DM2 = DistanceMatrix(x2,x2);
    DM2int = DistanceMatrix(x2int,x2);
    DiffM2int = DifferenceMatrix(x2int,x2);
    DM2cc = DistanceMatrix(x2cc,x2);
    DiffM2cc = DifferenceMatrix(x2cc,x2);
    V2 = rbf(ep,DM2);
    Vx2int = rbfx(ep,DM2int,DiffM2int);
    Vxx2int = rbfxx(ep,DM2int);
    Vx2cc = rbfx(ep,DM2cc,DiffM2cc);
    Vxx2cc = rbfxx(ep,DM2cc);
    Vx2intV2inv = Vx2int/V2;
    Vxx2intV2inv = Vxx2int/V2;
    Vx2ccV2inv = Vx2cc/V2;
    
    % Form the full problem
    x = [x1;x2];
    i1 = 1:N1;
    i2 = N1+1:N1+N2;
    
    odeint1 = @(u) r*x1int.*(Vx1intV1inv*u(i1)) + ...
                   .5*S^2*x1int.^2.*(Vxx1intV1inv*u(i1)) - ...
                   r*u(i1int);
    odeint2 = @(u) r*x2int.*(Vx2intV2inv*u(i2)) + ...
                   .5*S^2*x2int.^2.*(Vxx2intV2inv*u(i2)) - ...
                   r*u(N1+i2int);
	odebc1  = @(t,u) u(i1bc)-bc(x1bc,t);
	odebc2  = @(t,u) u(N1+i2bc)-bc(x2bc,t);
    odecc1  = @(u) u(i1cc) - u(N1+i2cc);
    odecc2  = @(t,u) Vx1ccV1inv*u(i1) - Vx2ccV2inv*u(i2) + exp(-100000*t);
    odefun = @(t,u) [odeint1(u);odebc1(t,u);odecc1(u); ...
                     odeint2(u);odebc2(t,u);odecc2(t,u)];
    
    % Prepare the ODE solve
    Mass = sparse([[eye(N1int),zeros(N1int,N1bc+N1cc+N2);zeros(N1bc+N1cc,N1+N2)]; ...
                   zeros(N2,N1),[eye(N2int),zeros(N2int,N2bc+N2cc);zeros(N2bc+N2cc,N2)]]);
    odeopt = odeset('Mass',Mass,...
        'MStateDependence','none','MassSingular','yes');
    Pinit = payout(x);
    [tsol,Psol] = ode15s(odefun,[0,T],Pinit,odeopt);
    
    errvec(k) = errcompute(Psol(end,:)',Ptrue(x,T));
    k = k + 1;
end

h = figure;
loglog(4*Nvec,errvec)
xlabel('number of collocation points')
ylabel('absolute max norm error')

h_sol = figure;
iplot = [i1bc,i1int,i1cc,i2int+N1,i2bc+N1];
plot(x(iplot),Psol(1,iplot)','linewidth',2)
hold on
plot(x(iplot),Psol(end,iplot)','--','linewidth',2)
hold off
ylim([0,3.5])
legend('at expiry','at inception','location','southeast')
xlabel('spot price')
ylabel('expected payoff')

h_err = figure;
surf(tsol,x(iplot),abs(Psol(:,iplot)' - bsxfun(@(xe,te)Ptrue(xe,te),x(iplot),tsol')),'edgealpha',.2)
xlabel('time to expiry')
ylabel('spot price')
zlabel('option value error')
view([-70 6])
colormap(gray);C = colormap;colormap(flipud(C))

% For plotting multiple results on the same axes
% loglog(Nvec,errvecM2,'linewidth',2)
% hold on
% loglog(Nvec,errvecM4,'--','linewidth',2)
% loglog(Nvec,errvecM6,'-.','linewidth',2)
% hold off
% legend('C^2 Matern','C^4 Matern','C^6 Matern','location','southwest')
% xlabel('number of collocation points')
% ylabel('absolute max norm error')