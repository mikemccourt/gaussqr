function [ pnts, Npnts, NV, Dist ] = SphereRegGoldPoints2( Npnts_surf, R, solvertype )
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
% Npnts_surf    =    desired number of points on the outermost surface;
%
% R             =    sphere radii row vector given as follows:
%                    [ R_innermost sphere, ... , R_outermost sphere ]
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
%                   [ [innermost region points]; 
%                     [innermost interface points, first set];
%                     [innermost interface points, second set];
%                                ...
%                     [outermost region points]; 
%                     [outermost boundary points] ]
%                   Note that the points of internal interfaces are divided
%                   in two sets, so that one set can be considered as
%                   belonging to the region inside the boundary and the
%                   other as belonging to the region outside the interface.
%
% Npnts        =    vector of numbers of points as the previous scheme;
%
% NV           =    matrix of unit vectors (3 components) normal to 
%                   boundary surfaces as the following scheme:
%                   [ [innermost interface normal unit vectors, first set];
%                     [innermost interface normal unit vectors, second set];
%                                    ...
%                     [outermost boundary normal unit vectors] ].
%
% Dist         =    distance between points on the outermost boundary 
%                   (approx);
%
% Calls on: 
%   SphereSurfGoldPoints.m
%   DistanceMatrix.m

b = 1; % Distance factor (see help)

% Arrays initialization
Npnts = [];
NV    = [];
pnts  = [];
R     = [0, R];
% Number of layers
NSL = size(R,2);

switch lower(solvertype)
    % Geometry data to be employed by Kansa Method
    case 'kansa'
        for i = NSL:-1:2
            
            % Points on i-th boundary
            bdypnts = SphereSurfGoldPoints(Npnts_surf, R(i));
            
            % Evaluation of mean distance between neighboring points
            DM = DistanceMatrix(bdypnts, bdypnts);
            dist = mean(min(DM+eye(size(DM))));
            Dist = dist*b;
            clear DM
            
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
            
            NL = round( (R(i)-R(i-1))/Dist ); % Number of internal sublayers
                                              % for the i-th layer
            r = (R(i)-R(i-1))/NL;
            a = R(i)^2/Npnts_surf; % Number of points per unit of surface
                                   % over 4*pi
            
            Npnts_j = zeros(1,NL-1); R_j = zeros(1,NL-1);
            for j = 1:(NL-1)
                R_j(j) = R(i-1) + r * j; % radius of the i-th inner sphere
                Npnts_j(j) = round( R_j(j)^2/a );
            end
            
            % Interior points in the i-th layer
            intpnts = zeros(sum(Npnts_j),3);
            N = 0;
            for j = 1:(NL-1)
                intpnts(N+1:N+Npnts_j(j),:) = SphereSurfGoldPoints(Npnts_j(j), R_j(j));
                N = N + Npnts_j(j);
            end
            
            Npnts = [Npnts; size(intpnts,1)];
            pnts =  [pnts; intpnts];
            clear intpnts
            
            Npnts_surf = round( R(i-1)^2/a );            
        end
        
    % Geometry data to be employed by Method of Fundamental Solutions
    case 'mfs'
        for i = NSL:-1:2
            
            % Points on i-th boundary
            bdypnts = SphereSurfGoldPoints(Npnts_surf, R(i));
            
            % Evaluation of mean distance between neighboring points
            DM = DistanceMatrix(bdypnts, bdypnts);
            Dist = mean(min(DM+eye(size(DM))));
            clear DM
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
            
            a = R(i)^2/Npnts_surf; % Number of points per unit of surface
                                   % over 4*pi
            
            Npnts_surf = round( R(i-1)^2/a );
        end
    otherwise
        error('Can only consider kansa or mfs: check input string')
end

% Rearrangements
Npnts = flipud(Npnts);
NV    = flipud(NV);
pnts  = flipud(pnts);

end