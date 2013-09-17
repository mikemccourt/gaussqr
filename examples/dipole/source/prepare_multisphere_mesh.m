function [geometry, N_elements] = prepare_multisphere_mesh(N,R)
% Elements (triangles) are distributed on each surface by assuming the same
% ratio between number of elements and surface area and by considering that
% the number of elements on each surface is roughly twice the number of
% points on the same surface.
%
% Input:
% N          = Total number of elements;
% R          = array of spheres radii (in ascending order);
%
% Output:
% geometry   = structure containing points and faces defining the mesh;
% N_elements = actual number of elements.

% Compute the number of points to be distributed on the innermost sphere
N_pnt_innermost = N/(sum(R.^2,2)/R(1)^2);

N_layers = length(R);
N_elements = zeros(1,N_layers);
geometry = [];
for i=1:N_layers
    % Compute the number of points to be distributed on the i-th sphere
    N_pnt = ceil(0.5*N_pnt_innermost*(R(i)/R(1))^2);
    % Distribute points over the i-th sphere
    geometry.bnd(i).pnt = SphereSurfGoldPoints(N_pnt,R(i));
    % Delaunay triangulation of the i-th sphere
    d_pnts = delaunayn(geometry.bnd(i).pnt); 
    tr = TriRep(d_pnts, geometry.bnd(i).pnt);
    geometry.bnd(i).tri = freeBoundary(tr);
    % Flip normals direction for OpenMEEG convention (inwards oriented)
    geometry.bnd(i).tri = fliplr(geometry.bnd(i).tri);
    % Store the actual number of elements on the i-th sphere
    N_elements(i) = length(geometry.bnd(i).tri); % Store the actual number 
                                                 % of elements
end
% N_elements = sum(N_elements);
end

