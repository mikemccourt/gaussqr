% disksourcetest_MFS.m
% This example studies the convergence rate for a 2D Laplace problem when
% the boundary conditions (or interface conditions) are tainted by noise.
% The problem we want to solve is an easy one:
%     Lap(u) = 0,     interior
%          u = v,     boundary

% The Dirichlet boundary condition function
% v = @(x) 1 + x(:,1) + x(:,2);
% v = @(x) x(:,1) + x(:,2) + 1./((x(:,1)-2).^2 + x(:,2).^2);
v = @(x) x(:,1) + x(:,2) + log((x(:,1)-2).^2 + x(:,2).^2);
vx = @(x) 1 + 2*(x(:,1)-2)./((x(:,1)-2).^2 + x(:,2).^2);
vy = @(x) 1 + 2*x(:,2)./((x(:,1)-2).^2 + x(:,2).^2);

% The fundamental solution
fs = @(x,z) log(DistanceMatrix(x,z));
fsd = @(x,z,dim) DifferenceMatrix(x(:,dim),z(:,dim))./(DistanceMatrix(x,z).^2);

% Choose to use the Neumann BC instead
neumannBC = 1;

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
oldtest = zeros(NN,1);
newtest = zeros(NN,1);
errdvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    % Define the number of points and the domains we are working on
    % xo is the ring collocation points, zo is the associated centers
    % xi is the inner circle collocation points, zi are the centers
    t = linspace(0,2*pi,N+1)';t = t(1:N);
    unit_circle = [cos(t),sin(t)];
    xo = ro*unit_circle;
    zo = [Ro*unit_circle(1:2:end,:);Ri*unit_circle(2:2:end,:)];
    xiD = ri*unit_circle(1:2:end,:);
    xiN = ri*unit_circle(2:2:end,:);
    zi = Rc*unit_circle;
%     zo = 1.7*unit_circle;
%     zi = 1.3*unit_circle;

    % Form the linear system, 6 blocks appear in this matrix
    % It takes the shape:
    %   [  K(xoD,zo)       0    ] [ co ]   [ v  ]
    %   [ GK(xoN,zo)       0    ] [    ]   [ Gv ]
    %   [  K(xiD,zo)  -K(xiD,zi)] [    ] = [ 0  ]
    %   [-GK(xiN,zo)  GK(xiN,zi)] [ ci ] = [ 0  ]
    % Also needed is the directional derivative in G'
    % Allow for Neumann BC if requested
    if neumannBC
        iN = 2:N;
    else
        iN = [];
    end
    iD = setdiff(1:N,iN);
    
    H = [[fs(xo(iD,:),zo),zeros(length(iD),N)];
        [diag(xo(iN,1)/ro)*fsd(xo(iN,:),zo,1)+diag(xo(iN,2)/ro)*fsd(xo(iN,:),zo,2),zeros(length(iN),N)];
        [fs(xiD,zo),-fs(xiD,zi)];
        [-diag(xiN(:,1)/ri)*fsd(xiN,zo,1)-diag(xiN(:,2)/ri)*fsd(xiN,zo,2),diag(xiN(:,1)/ri)*fsd(xiN,zi,1)+diag(xiN(:,2)/ri)*fsd(xiN,zi,2)]];
    rhs = [v(xo(iD,:));
           diag(xo(iN,1)/ro)*vx(xo(iN,:))+diag(xo(iN,2)/ro)*vy(xo(iN,:));
           zeros(length(xiD),1);
           zeros(length(xiN),1)];
    K_eval = fs(xx,zo);

    % Solve the system and compute the errors
    % Only using the coefficients for evaluating on the boundary 
    full_coef = H\rhs;
    outer_coef = full_coef(1:N);
    newtest = K_eval*outer_coef;
    errvec(k) = errcompute(newtest,yy);
    condvec(k) = cond(H);
    coefvec(k) = norm(outer_coef);
    if k>1
        errdvec(k) = errcompute(oldtest,newtest);
        if k==2
            errdvec(1) = errdvec(2);
        end
    end
    oldtest = newtest;
    k = k + 1;
end

h = figure;
loglog(Nvec,[errvec;condvec;coefvec;errdvec],'linewidth',3)
xlabel('Boundary collocation points')
legend('Error','Cond','Coef Norm','Test Diff')
title('Coupling')