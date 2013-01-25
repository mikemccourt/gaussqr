function h = SurfacePlot_dip( pnts, feval)
% function h = SurfacePlot_dip( pnts, feval)
%
% This function allows you to plot function values on the surface of a
% sphere.  I guess maybe it works for any manifold, but we only really use
% it for a sphere.
%
% 
% Input arguments
%   pnts - function locations
%   feval - function values
% Output arguments
%   h - plot handle

% Delaunay triangulation of the sphere
d_evalpnts = delaunayn(pnts);
tr = TriRep(d_evalpnts, pnts);
faces = freeBoundary(tr);

% figure('Color',[1 1 1]);
h = patch('Vertices',pnts,'Faces',faces,'FaceVertexCData',feval,...
      'FaceColor','interp','EdgeColor','none');
colorbar('FontSize',12);
axis('equal'); axis('square'); view(-30,30);
xlabel('x','FontWeight','bold','FontSize',14);
ylabel('y','FontWeight','bold','FontSize',14);
zlabel('z','FontWeight','bold','FontSize',14);
grid on;

if nargout==0 % Don't return anything if not requested
    clear h;
end

end

