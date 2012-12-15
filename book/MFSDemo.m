% MFSDemo.m
% MFS solver for a Helmholtz problem on a L-shaped domain
% Lap(u)+nu^2*u = 0                on interior, nu = 2*pi^2
% d/dn(u) = -pi*sin(pi*(x+y))    on y = 2
% u = cos(pi*(x+y))              on other boundary
% Calls on: DistanceMatrix, DifferenceMatrix
N = 64;nu = sqrt(2)*pi;
fundsol = @(r) 1i/4*besselh(0,2,nu*r);
fundsoldy = @(r,dy) -1i*nu/4*besselh(1,2,nu*r).*dy./r;
bc_DIR = @(x) cos(pi*(x(:,1)+x(:,2)));
bc_NEU = @(x) -pi*sin(pi*(x(:,1)+x(:,2)));
% Construct discretization of L-shaped domain
x1d = pickpoints(0,2,N/4+1);
x = [ [zeros(size(x1d)),x1d];
      [x1d,zeros(size(x1d))];
      [2*ones(size(x1d)),x1d];
      [x1d,2*ones(size(x1d))] ];
x = unique(1e-8*ceil(1e8*x),'rows'); % Dump duplicates
xLi = find(all(x<=1,2));
x(xLi,:) = ones(length(xLi),2) - x(xLi,:); % Build L-shape
xNEU = x(find(x(:,2)==2),:); % Neumann BC points
xDIR = setdiff(x,xNEU,'rows'); % Dirichlet BC points
% Choose Green's kernel source points
theta = pickpoints(0,2*pi,N);
z = ones(N,2) + 2*[cos(theta),sin(theta)];
% Set up the MFS system
DM_DIR = DistanceMatrix(xDIR,z); A_DIR = fundsol(DM_DIR);
DM_NEU = DistanceMatrix(xNEU,z);
DM_NEU_dy = DifferenceMatrix(xNEU(:,2),z(:,2));
A_NEU = fundsoldy(DM_NEU,DM_NEU_dy);
rhs_DIR = bc_DIR(xDIR); rhs_NEU = bc_NEU(xNEU);
A = [A_DIR;A_NEU]; rhs = [rhs_DIR;rhs_NEU];
% Solve the system, evaluate and plot the solution
coef = A\rhs;
xeval1 = pick2Dpoints([0,1],[1,2],10);
xe1x = unique(xeval1(:,1));xe1y = unique(xeval1(:,2));
xeval2 = pick2Dpoints([1,0],[2,1],10);
xe2x = unique(xeval2(:,1));xe2y = unique(xeval2(:,2));
xeval3 = pick2Dpoints([1,1],[2,2],10);
xe3x = unique(xeval3(:,1));xe3y = unique(xeval3(:,2));
A_eval = fundsol(DistanceMatrix([xeval1;xeval2;xeval3],z));
u = real(A_eval*coef); % Imaginary terms are noise
surf(xe1x,xe1y,reshape(u(1:100),10,10)), hold on
surf(xe2x,xe2y,reshape(u(101:200),10,10))
surf(xe3x,xe3y,reshape(u(201:300),10,10)), hold off
figure,plot(xDIR(:,1),xDIR(:,2),'or','linewidth',2), hold on
plot(xNEU(:,1),xNEU(:,2),'+m','linewidth',2)
plot(z(:,1),z(:,2),'xb','linewidth',2), hold off