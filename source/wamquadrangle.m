function [meshQ] = wamquadrangle(deg,Q);

%produces a (weakly) admissible mesh in a quadrangle

%INPUT:

%deg = interpolation degree

%Q = 2-columns array of the quadrangle vertices


%OUTPUT:

%meshQ = 2-columns array of mesh points coordinates


%FUNCTION BODY

%Chebyshev-Lobatto grid
j=(0:1:deg);
c=cos(j*pi/deg);

A=Q(1,:);
B=Q(2,:);
C=Q(3,:);
D=Q(4,:);

[u,v]=meshgrid(c,c);

meshQ=0.25*((1-u(:)).*(1-v(:))*A+(1+u(:)).*(1-v(:))*B);
meshQ=meshQ+0.25*((1+u(:)).*(1+v(:))*C+(1-u(:)).*(1+v(:))*D);

