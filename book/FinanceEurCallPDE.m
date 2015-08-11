% FinanceEurCallPDE
% This example demonstrates the use of kernels to compute the solution to
% Black-Scholes type PDEs in only one dimension
%
% The problem is a parabolic problem,
%     du/dt = Lu
%     u(t,boundary) = bc(boundary)
%     u(0,x) = payout(x)
% on the domain 0<x<4E and 0<t<1
% E is the strike price, which implies payout max(x-E,0)
%
% Our differential operator for this problem has 3 pieces:
%      Lu = L1u + L2u + L3u
% L1u = r*x'*grad(u)
% L2u = (1/2)*x'*((S*S').*hess(u))*x
% L3u = -r*u
%
% The exact solution is the formula for the value of a European call option
%      C(x,t) = F_Z(d1)*u - F_Z(d2)*E*exp(-r*(T-t)) 
% where
%      F_Z(z) = P(Z<z) for Z~N(0,1)
%      d1(x) = 1/(S*sqrt(T-t))*(log(x/E)+(r+S^2/2)*(T-t))
%      d2(x) = d1(x) - S*sqrt(T-t)
% Our solution below uses t in place of T-t because we know the final value
% and not the initial value.  This means that when solving the problem we
% are really stepping backwards in time.

global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% Define parameters of the PDE
E = 1;
T = 1;
B = .3;
r = .05;

% Define the payout, used as the initial condition
payout = @(x) max(x-E,0);

% Define the boundary conditions (works for x=0 and x=4)
bc = @(x,t) (x-E*exp(-r*t)).*(x==4*E);

% Define the true solution
d1 = @(x,t) 1./(B*sqrt(t)).*(log(x/E)+(r+B^2/2)*t);
d2 = @(x,t) d1(x,t) - B*sqrt(t);
Ptrue = @(x,t) normcdf(d1(x,t)).*x - E*normcdf(d2(x,t)).*exp(-r*t);

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
rbf = rbfM6;  rbfx = rbfM6x;  rbfxx = rbfM6xx;
% rbf = rbfM4;  rbfx = rbfM4x;  rbfxx = rbfM4xx;
% rbf = rbfM2;  rbfx = rbfM2x;  rbfxx = rbfM2xx;
ep = 2;

% Choose a range of point values to consider for this problem
Nvec = ceil(logspace(1,2,12));

errvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Create the points we want to work with
    xall = pickpoints(0,4*E,N);
    xbc = xall(xall==0 | xall==4*E); Nbc = length(xbc);
    xint = xall(xall~=0 & xall~=4*E); Nint = length(xint);
    iint = 1:Nint;
    ibc = Nint+1:Nint+Nbc;
    x = [xint;xbc];
    
    % Create the necessary differentiation matrices
    DM = DistanceMatrix(x,x);
    DMint = DistanceMatrix(xint,x);
    DiffMint = DifferenceMatrix(xint,x);
    V = rbf(ep,DM);
    Vxint = rbfx(ep,DMint,DiffMint);
    Vxxint = rbfxx(ep,DMint);
    VxintVinv = Vxint/V;
    VxxintVinv = Vxxint/V;
    
    % Form the functions for the ODE solver
    odeint = @(u) r*xint.*(VxintVinv*u) + ...
                  .5*B^2*xint.^2.*(VxxintVinv*u) - ...
                  r*u(iint);
	odebc  = @(t,u) u(ibc)-bc(xbc,t);
    odefun = @(t,u) [odeint(u);odebc(t,u)];
    
    % Prepare the ODE solve
    Mass = sparse([eye(Nint),zeros(Nint,Nbc);zeros(Nbc,Nint+Nbc)]);
    Jac = [.5*B^2*bsxfun(@times,xint.^2,VxxintVinv) + ...
                    r*bsxfun(@times,xint,VxintVinv) - ...
                       [r*eye(Nint),zeros(Nint,Nbc)]; ...
           [zeros(Nbc,Nint),eye(Nbc)]                      ];
    odeopt = odeset('Jacobian',Jac,'Mass',Mass,...
        'MStateDependence','none','MassSingular','yes');
    Pinit = payout(x);
    [tsol,Psol] = ode15s(odefun,[0,T],Pinit,odeopt);
    
    errvec(k) = errcompute(Psol(end,:)',Ptrue(x,T));
    k = k + 1;
end

h = figure;
loglog(Nvec,errvec)
xlabel('number of collocation points')
ylabel('absolute max norm error')

h_sol = figure;
iplot = [ibc(1),iint,ibc(end)];
plot(x(iplot),Psol(1,iplot)','linewidth',2)
hold on
plot(x(iplot),Psol(end,iplot)','--','linewidth',2)
hold off
ylim([0,3.5])
legend('at expiry','at inception','location','southeast')
xlabel('spot price')
ylabel('expected payoff')

h_err = figure;
surf(tsol,x(iplot),abs(Psol(:,iplot)' - bsxfun(@(xe,te)Ptrue(xe,te),x(iplot),tsol')))
xlabel('time to expiry')
ylabel('spot price')
zlabel('option value error')
view([-56 32])
colormap(gray);C = colormap;colormap(flipud(C)),colormap white
set(get(gca,'children'),'edgealpha',.8)
set(get(gca,'children'),'edgecolor',[1 1 1]*.5)
title(sprintf('error=%g',errvec(1)))

% For plotting multiple results on the same axes
% loglog(Nvec,errvecM2,'linewidth',2)
% hold on
% loglog(Nvec,errvecM4,'--','linewidth',2)
% loglog(Nvec,errvecM6,'-.','linewidth',2)
% hold off
% legend('C^2 Matern','C^4 Matern','C^6 Matern','location','southwest')
% xlabel('number of collocation points')
% ylabel('absolute max norm error')