% disksourcetest_MFS.m
% This example studies the convergence rate for a 2D Laplace problem when
% the boundary conditions (or interface conditions) are tainted by noise.
% The problem we want to solve is an easy one:
%     Lap(u) = 0,     interior
%          u = v,     boundary

% The Dirichlet boundary condition function
%v = @(x) 1 + x(:,1) + x(:,2);
% v = @(x) x(:,1) + x(:,2) + 1./((x(:,1)-2).^2 + x(:,2).^2);
v = @(x) x(:,1) + x(:,2) + log((x(:,1)-2).^2 + x(:,2).^2);

% The fundamental solution
fs = @(x,z) log(DistanceMatrix(x,z));
fsd = @(x,z,dim) DifferenceMatrix(x(:,dim),z(:,dim))./(DistanceMatrix(x,z).^2);

% Error to introduce into the boundary condition (0 for no noise)
noise = 0.0001;

% Circle radii
% We're going to do a coupling between an inner circle and outer ring
% The inner circle will have radius ri and outer will have radius ro
% Fictitious boundary for the ring is Ro outside, Ri inside
% Fictitious boundary for the inner circle is the average of ro and ri
ri = .5;
ro = 1;
Ro = 1.5;
Ri = .3;
Rc = .5*(ri+ro);

% Number of collocation and source points and test points
NN = 100;
tt = pickpoints(0,2*pi,NN+1,'rand');tt = tt(1:NN);
xx = [cos(tt),sin(tt)];
yy = v(xx);

% Define the problem size vector
% This number of points will appear both on the outer boundary and on the
% interior boundary
Nvec = floor(logspace(1,3,20));

errvec = zeros(size(Nvec));
coefvec = zeros(size(Nvec));
condvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Define the number of points and the domains we are working on
    % xo is the ring collocation points, zo is the associated centers
    % xi is the inner circle collocation points, zi are the centers
    t = linspace(0,2*pi,N+1)';t = t(1:N);
    unit_circle = [cos(t),sin(t)];
    xo = ro*unit_circle;
    zo = [Ro*unit_circle(1:2:end,:);Ri*unit_circle(2:2:end,:)];
    xiI = ri*unit_circle(1:2:end,:);
    xiD = ri*unit_circle(2:2:end,:);
    zi = Rc*unit_circle;
%     zo = 3*unit_circle;
%     zi = 2.5*unit_circle;

    % Form the linear system, 6 blocks appear in this matrix
    % It takes the shape:
    %   [ K(xo,zo)        0    ] [ co ]   [ v ]
    %   [ K(xiI,zo)  -K(xiI,zi)] [    ] = [ 0 ]
    %   [GK(xiD,zo) -GK(xiK,zi)] [ ci ] = [ 0 ]
    % Also needed is the directional derivative in G
    H = [[fs(xo,zo),zeros(N)];
        [fs(xiI,zo),-fs(xiI,zi)];
        [diag(xiD(:,1)/ri)*fsd(xiD,zo,1)+diag(xiD(:,2)/ri)*fsd(xiD,zo,2),diag(xiD(:,1)/ri)*fsd(xiD,zi,1)+diag(xiD(:,2)/ri)*fsd(xiD,zo,2)]];

    K_eval = fs(xx,zo);
    rhs = [v(xo);zeros(length(xiI),1);zeros(length(xiD),1)];

    % Solve the system and compute the errors
    % Only using the coefficients for evaluating on the boundary 
    full_coef = H\rhs;
    outer_coef = full_coef(1:N);
    yp = K_eval*outer_coef;
    errvec(k) = errcompute(yp,yy);
    condvec(k) = cond(H);
    coefvec(k) = norm(outer_coef);
    k = k + 1;
end

h = figure;
loglog(Nvec,[errvec;condvec;coefvec],'linewidth',3)
xlabel('Boundary collocation points')
legend('Error','Cond','Coef Norm')
title('Coupling')