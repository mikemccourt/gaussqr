% PDESpaceTimeTransport
% This script is going to consider the space-time kernel
% Now we want to work with a transport equation
%          u_t + move_const*u_x = 0
%          u(-1,t) = 0, u(x,0) = f(x)
% The true solution is u(x,t) = f(x-c*t)
%
% We approximate u(x,t) as
%       u(x,t) = sum_{i=1}^N c_i K([x,t],[x_i,t_i])
% where [x_i,t_i] are the points in the space-time domain
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

move_const = .25;
fic  = @(x) (512*(1-2*x(:,1)).^3.*(x(:,1)).^3).*(x(:,1)<=0.5);
true_sol  = @(x) (512*(1-2*(x(:,1)-move_const*x(:,2))).^3.* ...
                 (x(:,1)-move_const*x(:,2)).^3).*...
                 (abs(x(:,1)-move_const*x(:,2)-.25)<=.25);

% We call x both the time and space points
% Create points that are uniform in space and Chebyshev in time
% Make sure to dump the boundaries for the IBB kernel
Nx = 35;
Nt = 36;
tmax = 1;
x1Du = pickpoints(0,1,Nx+2); x1Du = x1Du(2:end-1);
x1Dc = pickpoints(0,tmax,Nt,'cheb');
[X,T] = meshgrid(x1Du,x1Dc);
xall = [X(:),T(:)];

% Come up with some evaluation points
Nevalx = 50;
Nevalt = tmax*Nevalx;
xeval = pick2Dpoints([0 0],[1 tmax],Nevalx);

% Choose to evaluate error at evaluation points or collocation points
%      1 - evaluation points
%      2 - collocation points
err_eval = 1;

% Choose a kernel for our problem
% We're making a crazy tensor product kernel up for this one
epvec = [10,.5];

% Analytic Chebyshev kernel
KCA = @(b,x,z) .5 + (1-b)* ...
         (b*(1-b^2) - 2*b*bsxfun(@plus,x.^2,z.^2') + (1+3*b^2)*x*z')./ ...
         ((1-b^2)^2 + 4*b*(b*bsxfun(@plus,x.^2,z.^2')-(1+b^2)*x*z'));
KCAdx = @(b,x,z) (1-b)* ...
        (((1-b^2)^2 + 4*b*(b*bsxfun(@plus,x.^2,z.^2')-(1+b^2)*x*z')).* ...
            bsxfun(@plus,        -4*b*x,        (1+3*b^2)*z') - ...
         (b*(1-b^2) - 2*b*bsxfun(@plus,x.^2,z.^2') + (1+3*b^2)*x*z').* ...
            (4*b).*bsxfun(@minus,  2*b*x,             (1+b^2)*z'))./ ...
         ((1-b^2)^2 + 4*b*(b*bsxfun(@plus,x.^2,z.^2')-(1+b^2)*x*z')).^2;
% Scale the Chebyshev kernels to the [-1,1] they love so much
% Don't forget the chain rule
KCAs = @(b,x,z) KCA(b,x*2/tmax-1,z*2/tmax-1);
KCAsdx = @(b,x,z) 2/tmax*KCAdx(b,x*2/tmax-1,z*2/tmax-1);
% C2 IBB Kernel
KI2 = @(e,x,z) ibb(x,z,e,2);
KI2dx = @(e,x,z) ibb(x,z,e,2,1);
% Form the tensor kernel
Kcell = {KI2,KCAs};
Kdxcell = {KI2dx,KCAs};
Kdtcell = {KI2,KCAsdx};
Kf = @(Kc,e,x,z) prod(cell2mat(reshape( ...
         cellfun(@(K,e,x,z) K(e,x,z), ...
         Kc,num2cell(epvec),num2cell(x,1), ...
         num2cell(z,1),'UniformOutput',0), ...
                            [1,1,length(epvec)])),3);
%%%%%%%%%%%%%%%%%%%%% Problem section
% There may be some issues with using below ... some sort of limit
% exists...
% This is the C2 Chebyshev kernel with a = 1
B4 = @(x) x.^4 - 2*x.^3 + x.^2 - 1/30;
B3 = @(x) x.^3 - 3/2*x.^2 + 1/2*x;
B2 = @(x) x.^2 - x + 1/6;
KC2 = @(x,z) -30*(B4(abs(bsxfun(@plus,acos(x),acos(z')))/(2*pi)) + ...
                    B4(abs(bsxfun(@minus,acos(x),acos(z')))/(2*pi)));
KC2dx = @(x,z) 60/pi*bsxfun(@ldivide,(x>-.995).*sqrt(1-x.^2)+eps, ...
    (sign(bsxfun(@plus,acos(x),acos(z'))).*B3(abs(bsxfun(@plus,acos(x),acos(z')))/(2*pi)) + ...
     sign(bsxfun(@minus,acos(x),acos(z'))).*B3(abs(bsxfun(@minus,acos(x),acos(z')))/(2*pi))) ...
                           ) + ...
                           ( ...
               40/pi*bsxfun(@ldivide,(B2(abs(bsxfun(@plus,acos(x),acos(z')))/(2*pi)) + ...
                      B2(abs(bsxfun(@minus,acos(x),acos(z')))/(2*pi)))) ...
                           );
%%%%%%%%%%%%%%%%%%%%%% End of problem section

% We organize the points so that:
%     [t=0]
%     [interior]
% Remember that the x boundary is removed because of the IBB kernel
xic  = xall(xall(:,2)==0,:);
xint = xall(xall(:,2)~=0,:);
x = [xic;xint];

% Evaluate the collocation matrix
Aic = Kf(Kcell,epvec,xic,x);
Aint = Kf(Kdtcell,epvec,xint,x) + move_const*Kf(Kdxcell,epvec,xint,x);
A = [Aic;Aint];

% Evaluate the RHS
rhsic  = fic(xic);
rhsint = zeros(size(xint,1),1);
rhs = [rhsic;rhsint];

% Solve the system to find the coefficients
coef = A\rhs;

% Evaluate at points in the domain
ueval = Kf(Kcell,epvec,xeval,x)*coef;

% Plot the results
h = figure;
subplot(1,2,1)
X = reshape(xeval(:,1),Nevalx,Nevalt);
T = reshape(xeval(:,2),Nevalx,Nevalt);
U = reshape(ueval,Nevalx,Nevalt);
surf(X,T,U,'edgecolor','none');
xlabel('x')
ylabel('t')
zlabel('u')

% Determine which points to check the error at
if err_eval==1
    xtest = xeval;
    Nxtest = Nevalx;
    Nttest = Nevalt;
else
    xtest = xall;
    Nxtest = Nx;
    Nttest = Nt;
end

% Compare this solution to the true solution at error test points
utest = Kf(Kcell,epvec,xtest,x)*coef;
utrue = true_sol(xtest);
title(sprintf('absolute max norm error %g',errcompute(utest,utrue)))

% Plot the difference between the collocation and
% Normalize it against the norm of the initial condition
subplot(1,2,2)
udiff = abs(utest-utrue)/norm(utrue(:,1));
X = reshape(xtest(:,1),Nxtest,Nttest);
T = reshape(xtest(:,2),Nxtest,Nttest);
Udiff = reshape(udiff,Nxtest,Nttest);
surf(X,T,Udiff,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('relative error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Separate interpolation test
% Compare to the best possible error that could occur, which would occur
% for the interpolation without any PDE stuff at all
uint = Kf(Kcell,epvec,xtest,x)*(Kf(Kcell,epvec,x,x)\true_sol(x));

h_int = figure;
subplot(1,2,1)
Uint = reshape(uint,Nxtest,Nttest);
surf(X,T,Uint,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('interpolation solution')
title(sprintf('absolute max norm error %g',errcompute(uint,utrue)))

subplot(1,2,2)
uintdiff = abs(uint-utrue)/norm(utrue(:,1));
Uintdiff = reshape(uintdiff,Nxtest,Nttest);
surf(X,T,Uintdiff,'edgecolor','none')
xlabel('x')
ylabel('t')
zlabel('interpolation error')