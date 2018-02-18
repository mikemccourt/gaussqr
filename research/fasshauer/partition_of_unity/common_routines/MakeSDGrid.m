% gridpoints = MakeSDGrid(s,neval)
% Produces matrix of equally spaced points in s-dimensional unit cube
% (one point per row)
% Input
%   s:     space dimension
%   neval: number of points in each coordinate direction
% Output
%   gridpoints: neval^s-by-s matrix (one point per row,
%               d-th column contains d-th coordinate of point)
function gridpoints = MakeSDGrid(s,neval)
if (s==1)
    gridpoints = linspace(0,1,neval)';
    return;
end
% Mimic this statement for general s:
% [x1, x2] = ndgrid(linspace(0,1,neval));
outputarg = 'x1';
for d = 2:s
    outputarg = strcat(outputarg,',x',int2str(d));
end
makegrid = strcat('[',outputarg,'] = ndgrid(linspace(0,1,neval));');
eval(makegrid);
% Mimic this statement for general s:
% gridpoints = [x1(:) x2(:)];
gridpoints = zeros(neval^s,s);
for d = 1:s
    matrices = strcat('gridpoints(:,d) = x',int2str(d),'(:);');
    eval(matrices);
end
