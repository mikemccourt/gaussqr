function [ pnts, Npnts, NV, Dist ] = SphereRegGoldPoints3( R, Npnts_surf, solvertype )
% SphereRegGoldPoints creates a (almost) regular distribution of points 
% inside a multilayer sphere of given radii R starting receiving as input
% the desired number Npnts_surf of points equally distributed on the 
% outermost surface.
%
% Criterion:
% The distances between every interior point and the closest others are
% are related to the distance Dist between two neighboring points on the 
% surface by a factor b (see the first line of the function).
%
% Inputs:
% Npnts_surf    =    desired number of points on the surfaces given as 
%                    follows:
%                    [ Npnts_surf_innermost surface, ... , R_outermost Npnts_surf ]
%
% R             =    sphere radii row vector given as follows:
%                    [ R_innermost sphere, ... , R_outermost sphere ]
%
% Dist          =    distance between an internal point and the neighboring
%                    ones.
%
% solvertype    =    'kansa'  if data are needed to apply Kansa Method;
%                    'mfs'    if data are needed to apply Method of
%                             Fundamental Solutions
%                     	      If 'mfs' is selected, only boundary data will
%                             be provided
%
% Outputs:
% pnts         =    matrix of points' coordinates as the following 
%                   scheme:
%                   [ [innermost boundary points, first set];
%                     [innermost boundary points, second set];
%                     [innermost region points]; 
%                                ...
%                     [outermost region points]; 
%                     [outermost boundary region points] ]
%                   Note that the points of internal boundaries are divided
%                   in two sets, so that one set can be considered as
%                   belonging to the region inside the boundary and the
%                   other as belonging to the region outside the boundary.
%
% Npnts        =    vector of numbers of points as the previous scheme;
%
% NV           =    matrix of unit vectors (3 components) normal to 
%                   boundary surfaces as the following scheme:
%                   [ [innermost boundary normal unit vectors, first set];
%                     [innermost boundary normal unit vectors, second set];
%                                    ...
%                     [outermost boundary normal unit vectors] ].
%
% Dist         =    distance between points on the outermost surface 
%                   (approx);
%
% Calls on: 
%   SphereSurfGoldPoints.m
%   DistanceMatrix.m

% Arrays initialization
Npnts = [];
NV    = [];
pnts  = [];

R     = [0, R];
% Npnts_surf =  [0, Npnts_surf];
NSL = size(R,2);

switch lower(solvertype)
    % Geometry data to be employed by Kansa Method
    case 'kansa'
        for i = NSL:-1:2
            
            % Points on i-th boundary
            bdypnts = SphereSurfGoldPoints(Npnts_surf(i-1), R(i));
            
            DM = DistanceMatrix(bdypnts, bdypnts);
            dist = mean(min(DM+eye(size(DM))));
            Dist = dist;
            
            if i == NSL
                % The splitting of the set of point in two subsets is not
                % required for the outermost layer
                Npnts = [Npnts; size(bdypnts,1)];
                pnts =  [pnts; bdypnts];
                NV =    [NV; bdypnts/R(i)];
                clear bdypnts
            else
                A = bdypnts(1:2:end,:); % Pick points in even rows of bdypnts
                B = bdypnts(2:2:end,:); % Pick points in odd rows of bdypnts
                Npnts = [Npnts; size(A,1); size(B,1)];
                pnts =  [pnts; A; B];
                NV =    [NV; A/R(i); B/R(i)];
                clear bdypnts A B
            end
            
            intpnts = SphereLayerUnifPoints( Dist, R(i-1), R(i) );
            
            Npnts = [Npnts; size(intpnts,1)];
            pnts =  [pnts; intpnts];
            clear intpnts
                       
        end
        
    % Geometry data to be employed by Method of Fundamental Solutions
    case 'mfs'
         for i = NSL:-1:2
            
            % Points on i-th boundary
            bdypnts = SphereSurfGoldPoints(Npnts_surf(i), R(i));
                        
            if i == NSL
                % The splitting of the set of point in two subsets is not
                % required for the outermost layer
                Npnts = [Npnts; size(bdypnts,1)];
                pnts =  [pnts; bdypnts];
                NV =    [NV; bdypnts/R(i)];
                clear bdypnts
            else
                A = bdypnts(1:2:end,:); % Pick points in even rows of bdypnts
                B = bdypnts(2:2:end,:); % Pick points in odd rows of bdypnts
                Npnts = [Npnts; size(A,1); size(B,1)];
                pnts =  [pnts; A; B];
                NV =    [NV; A/R(i); B/R(i)];
                clear bdypnts A B
            end
           
        end
    otherwise
        error('Can only consider kansa or mfs: check input string')
end

% Rearrangements
Npnts = flipud(Npnts);
NV    = flipud(NV);
pnts  = flipud(pnts);

end

function [ pnts ] = SphereLayerUnifPoints( Dist, Rint, Rest )
%SphereLayerUnifPoints creates a uniform distribution of points inside a
% a sphere shell (layer) of given radii Rint and Rest starting from the 
% desired distance between points Dist.
%
% Inputs:
% Dist          =    desired distance between points;
% Rint          =    inner sphere radius.
% Rest          =    outer sphere radius.
%
% Outputs:
% pnts          =    matrix of points' coordinates.

[x, y, z] = meshgrid(-Rest:Dist:Rest, -Rest:Dist:Rest, -Rest:Dist:Rest);
ind = find((x(:).^2+y(:).^2+z(:).^2) < (Rest - Dist*0.5)^2 & (x(:).^2+y(:).^2+z(:).^2) > (Rint + Dist*0.5)^2);
n=length(ind);
pnts = zeros(n,3);
for i=1:n
  pnts(i,:) = [x(ind(i)), y(ind(i)), z(ind(i))];
end
end