% PDEEigenfunctionBasis3D
% This problem considers the solution to a 3D boundary value problem with
% the eigenfunction basis instead of the standard basis
% The BVP is
%       Lap(u) = 2                (x,y,z) in (-1,1)^3
%       u_y    = 2/3                 y = 1
%       u      = (x^2+y^2+z^2)/3    Boundary except y = 1

% Define the necessary functions
fint = @(x) 2*ones(size(x,1),1);
fneu = @(x) 2/3*ones(size(x,1),1);
fdir = @(x) (x(:,1).^2+x(:,2).^2+x(:,3).^2)/3;

% Define some points in the interior
% Scale them into [-1,1]^3
Nint = 300;
if exist('haltonset','file')
    point_generator = haltonset(3,'Skip',1);
    x3dhalton = net(point_generator,Nint);
else
    x3dhalton = haltonseq(Nint,3);
end
xint = 2*x3dhalton-1;

% We need boundary points, and some of the boundary points
% need to be associated with the Neumann condition
Nbc = 5;
x2d = pick2Dpoints(-1,1,Nbc);
xneu = [x2d(:,1),ones(Nbc^2,1),x2d(:,2)];
xdir = [x2d(:,1),-ones(Nbc^2,1),x2d(:,2);
        -ones(Nbc^2,1),x2d(:,1),x2d(:,2);ones(Nbc^2,1),x2d(:,1),x2d(:,2);
        x2d(:,1),x2d(:,2),-ones(Nbc^2,1);x2d(:,1),x2d(:,2),ones(Nbc^2,1)];
    
% Set up the collocation points
x = [xint;xneu;xdir];

% Set up the Gaussian eigenfunctions
GAUSSQR_PARAMETERS.DEFAULT_REGRESSION_FUNC = .3;
alpha = 1;
ep = .1;
GQR = gqr_solveprep(1,x,ep,alpha);

% Evaluate the components of the collocation problem
Phixx = gqr_phi(GQR,xint,[2 0 0]);
Phiyy = gqr_phi(GQR,xint,[0 2 0]);
Phizz = gqr_phi(GQR,xint,[0 0 2]);
Phiy  = gqr_phi(GQR,xneu,[0 1 0]);
Phi   = gqr_phi(GQR,xdir);

% Evaluate the right hand side function
rhsint = fint(xint);
rhsneu = fneu(xneu);
rhsdir = fdir(xdir);
rhs = [rhsint;rhsneu;rhsdir];

% Form the collocation matrix
A = [Phixx+Phiyy+Phizz;Phiy;Phi];

% Solve the system and store the coefficients
GQR.coef = A\rhs;

% Evaluate the error and plot the results
h = figure;
ucoll = gqr_eval(GQR,x);
utrue = fdir(x);
uplot = abs(ucoll - utrue);
scatter3(x(:,1),x(:,2),x(:,3),15,uplot,'filled')
xlabel('x'),ylabel('y'),zlabel('z')
h_color = colorbar;

% Do some more stupid plotting stuff
% For most of this stuff below, you should see the Matlab help file entry
% to understand what is happening - we basically copied out of that
h_contour = figure;
x1d = pickpoints(-1,1,15);
[X,Y,Z] = meshgrid(x1d,x1d,x1d);
x3v = [X(:),Y(:),Z(:)];
u3v = fdir(x3v);
u3coll = gqr_eval(GQR,x3v);
U = reshape(abs(u3v-u3coll),15,15,15);
contourslice(X,Y,Z,U,[],[-1,-.5,0,.5,1],[],[1e-4,1e-3,2e-3,3e-3,3.9e-3,5e-3])
xlabel('x'),ylabel('y'),zlabel('z')

% Weird slice plot stuff
h_slice = figure;
hslice = surf(linspace(-1,1,100),...
   linspace(-1,1,100),...
   zeros(100));

rotate(hslice,[-1,0,0],-45)
xd = get(hslice,'XData');
yd = get(hslice,'YData');
zd = get(hslice,'ZData');
delete(hslice)

h = slice(X,Y,Z,U,xd,yd,zd);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

hold on
hx = slice(X,Y,Z,U,1,[],[]);
hx.FaceColor = 'interp';
hx.EdgeColor = 'none';

hy = slice(X,Y,Z,U,[],1,[]);
hy.FaceColor = 'interp';
hy.EdgeColor = 'none';

hz = slice(X,Y,Z,U,[],[],-1);
hz.FaceColor = 'interp';
hz.EdgeColor = 'none';
hold off

h_color = colorbar;
xlabel('x'),ylabel('y'),zlabel('z')