function meshD = wamdisk(deg)

%produces a weakly admissible mesh in the unit disk

%INPUT:

%deg = interpolation degree

%OUTPUT:

%meshD = 2-columns array of mesh points coordinates


%FUNCTION BODY

%producing a weakly-admissible mesh
%discretization parameters
n1=deg;
if (mod(deg,2) == 0)
n2=deg+2;
else
n2=deg+1;
end;

%construction of the extraction symmetric polar grid 
j1=(0:1:n1); j2=(0:1:n2-1); [rho,theta]=meshgrid(cos(j1*pi/n1),j2*pi/n2);

meshD=[rho(:).*cos(theta(:)) rho(:).*sin(theta(:))];

% Dump all points at (0,0), which may be duplicates, and then append (0,0)
% back on to the set
meshD = meshD(sum(abs(meshD),2)>eps,:);
if mod(deg,2)==0
    meshD = [meshD;[0,0]];
end