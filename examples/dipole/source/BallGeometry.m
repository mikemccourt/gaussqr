function [ POINTS, NORMALS ] = BallGeometry( R, ...
                                              Npnts, ...
                                              solvertype, ...
                                              intptstype , ...
                                              bdyptstype)
% function [POINTS,NORMALS] = BallGeometry(R,Npnts,solvertype,ptstype)
% BALLGEOMETRY creates a distribution of points in a multilayer (N layers)
% ball.
% Input data are: spheres radii, desired number of interior points, desired
% number of boundary points (optional), meshfree solver type.
% Output data are points coordinates and unit vectors normal to boundary 
% spherical surfaces.
% Interior points are posed in a uniform grid. Boundary points are
% distributed by means of gold section spirals or Halton distributions.
%
% Inputs:
% R             =   sphere radii vector given as follows:
%                   [ R_innermost sphere, ... , R_outermost sphere ].
%
% Npnts         =   For 'kansa' solvertype, desired number of interior and 
%                   boundary points given as follows:
%                   [ N_interior points, N_points on the innermost sphere,
%                     ... , N_points on the outermost sphere ]
%                   If only the number of interior points is provided, the
%                   number of points on boundaries is automatically set.
%                   For 'mfs' solvertype, desired number of boundary points 
%                   given as follows:
%                   [ N_points on the innermost sphere, ...
%                     ... , N_points on the outermost sphere ].
%
% solvertype    =   'kansa'  if data are needed to apply Kansa Method;
%                   'mfs'    if data are needed to apply Method of
%                            Fundamental Solutions
%                   If 'mfs' is selected, only boundary data will be 
%                   provided as output.
%
% intptstype    =   Interior points distribution:
%                   'even'   Points uniformly distributed in the cube
%                   'halton' Points should have a Halton distribution
%                   'random' Points are uniform randomly chosen
%                   'cheb'   Points are Chebyshev spaced along radii
%                   ** the default choice is 'halton', which can be
%                   selected by passing []
%
% bdyptstype    =   Boundary points type
%                   'spiral' Points regularly distributed by means of a
%                            gold section spiral
%                   'halton' Points distributed by mapping Halton points
%                            from the unit square to the spherical surface
%                   ** the default choice is 'spiral', which can be
%                   selected by passing []
%
% avoidpts      =   Vector of points that there should be no points in the
%                   BallGeometry touching
%                   pass [] to include no such points
%
% avoidcushion  =   Size of the ball around the avoidpts which no points in
%                   the BallGeometry will be inside.
%                   This value can be as small as you would like, but the
%                   maximum allowable size is R, since that could negate
%                   the whole sphere.
%                   ** the default choice is R/5, which can be selected by
%                   passing []
%                   
%
% Outputs:
% POINTS        =   data structure containing points' coordinates as
%                   follows:
%                   POINTS.int1 = coordinates of innermost shell points;
%                   POINTS.bdy11 = coordinates of innermost boundary points, first set;
%                   POINTS.bdy12 = coordinates of innermost boundary points, second set;
%                   ...
%                   POINTS.intN = coordinates of outermost shell points;
%                   POINTS.bdyN = coordinates of outermost boundary points.
%
%                   Note that the points of internal interfaces are divided
%                   in two sets, so that one set can be considered as
%                   belonging to the region inside the boundary and the
%                   other as belonging to the region outside the interface.
%
% NORMALS      =    data structure containing interfaces' normal unit 
%                   vectors coordinates as follows:
%                   NORMALS.n11 = innermost interface normal unit vectors, first set;
%                   NORMALS.n12 = innermost interface normal unit vectors, second set;
%                   ...
%                   NORMALS.nN = outermost interface normal unit vectors.
%
% Required functions: 
%   SphereSurfGoldPoints.m
%   SphereSurfHaltonPoints.m

N = length(R); % Number of layers
lNpnts = length(Npnts);

% Check input data
if not(isvector(R)) || not(isvector(Npnts)) || not(ischar(solvertype))
    error('R and Npnts must be vectors, solvertype must be a string: check input data')
end
if N ~= 1 && not(issorted(R))
    error('R must be sorted in ascending order')
end
if any(Npnts<=0) || any(R<=0)
    error('All radii and point requests must be positive')
end
R = R(:)'; % Make sure it's a row vector

Rend = R(end);
intptstype_DEFAULT = 'halton';
bdyptstype_DEFAULT = 'spiral';
chebNumRegion = 12*lNpnts; % Maybe this should be passed or adaptive ...

% Run checks to make sure that the point distribution is okay
if nargin==3
    intptstype = intptstype_DEFAULT;
    bdyptstype = bdyptstype_DEFAULT;
else
    
    if isempty(intptstype)
        intptstype = intptstype_DEFAULT;
    else
        if ischar(intptstype)
            if ~any(strcmp(intptstype,{'even','halton','random','cheb'}))
                warning('intptstype=%s unacceptable, defaulting to %s',intptstype,intptstype_DEFAULT)
                intptstype = intptstype_DEFAULT;
            end
        else
            warning('intptstype=%g unacceptable, defaulting to %s',intptstype,intptstype_DEFAULT)
            intptstype = intptstype_DEFAULT;
        end
    end
    
    if isempty(bdyptstype)
        bdyptstype = bdyptstype_DEFAULT;
    else
        if ischar(bdyptstype)
            if ~any(strcmp(bdyptstype,{'spiral','halton'}))
                warning('bdyptstype=%s unacceptable, defaulting to %s',bdyptstype,bdyptstype_DEFAULT)
                bdyptstype = bdyptstype_DEFAULT;
            end
        else
            warning('bdyptstype=%g unacceptable, defaulting to %s',bdyptstype,bdyptstype_DEFAULT)
            bdyptstype = bdyptstype_DEFAULT;
        end
    end
    
end

% Actually determine the point distribution in and on the ball
switch lower(solvertype)
    case 'kansa'
        R1 = R;
        R = [0 R];
        
        if lNpnts == 1
            % Distribute the points on surfaces so that the surface point
            % density is equal in all surfaces and the number of interior 
            % points is twice the number of points on the outermost sphere 
            a = Npnts/(sum(R1.^2,2)+2*R1(end)^2);
            for i=1:N
                Npnts_bdy(i) = round(a*R1(i)^2);
            end
            Npnts_int = 2*Npnts_bdy(end);
        elseif lNpnts == N
            % Distribute as many points on each surface as the vector Npnts
            % indicates. The number of interior points will be twice the
            % number of the points on the outermost sphere
            Npnts_bdy = Npnts;
            Npnts_int = 2*Npnts_bdy(end);
        elseif lNpnts == N+1
            % Distribute as many points on the interior as the first 
            % elements in the vector Npnts indicates and as many points 
            % on each surface as the other elements in Npnts indicates.
            Npnts_bdy = Npnts(2:end);
            Npnts_int = Npnts(1);
        end
        
        if strcmp(intptstype,'even')
            d = 2*Rend / ( 6/pi * Npnts_int )^(1/3);
            x = -Rend:d:Rend;
            [x, y, z] = meshgrid(x);
            ptsvec = [x(:),y(:),z(:)];
        elseif strcmp(intptstype,'halton')
            ptsvec = Rend*(2*haltonseq(ceil(6/pi*Npnts_int),3)-1);
        elseif strcmp(intptstype,'random')
            ptsvec = Rend*(2*rand(ceil(6/pi*Npnts_int),3)-1);
        elseif strcmp(intptstype,'cheb')
            Nc = chebNumRegion;
            Nradii = Npnts_int/Nc;
            radii = SphereSurfGoldPoints(Nradii, Rend);
            temp = pickpoints(-R(2),R(2),Nc,'cheb');
            cfact = temp(ceil(Nc/2)+1:end-1);
            radMat = kron(radii,ones(length(cfact),1));
            cMat = repmat(cfact,length(radii),3);
            ptsvec = radMat.*cMat;
        else
            error('Unknown interior point distribution style')
        end
             
        ll = 0;
            
        int_DM = DistanceMatrix(ptsvec,[0,0,0]);
        for l = 1:N
            % Design necessary strings for final object
            lstring = num2str(l);
            lstring1 = num2str(l+1);
            name_int = strcat('int', lstring);
            name_bdy = strcat('bdy', lstring, lstring);
            name_bdy1 = strcat('bdy', lstring, lstring1);
            name_normals = strcat('n', lstring, lstring);
            name_normals1 = strcat('n', lstring, lstring1);
            
            % Interior points
            ind = logical( int_DM < R(l+1) & int_DM > R(l));
            POINTS.(name_int) = ptsvec(ind,:);
            
            % Boundary points
            if strcmp(bdyptstype,'spiral')
                POINTS.(name_bdy) = SphereSurfGoldPoints(Npnts_bdy(l), R(l+1));
                if l == N
                    % The splitting of the set of point in two subsets is not
                    % required for the outermost layer
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l+1);
                else
                    % Split bountary points into two sets
                    POINTS.(name_bdy1) = POINTS.(name_bdy)(1:2:end,:);
                    POINTS.(name_bdy) = POINTS.(name_bdy)(2:2:end,:);
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l+1);
                    NORMALS.(name_normals1) = POINTS.(name_bdy1) / R(l+1);
                end
            elseif strcmp(bdyptstype,'halton')
                POINTS.(name_bdy) = SphereSurfHaltonPoints(Npnts_bdy(l),R(l+1));
                if l == N
                    % The splitting of the set of point in two subsets is not
                    % required for the outermost layer
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l+1);
                else
                    % Split bountary points into two sets
                    POINTS.(name_bdy1) = POINTS.(name_bdy)(1:floor(Npnts_bdy(l)/2),:);
                    POINTS.(name_bdy) = POINTS.(name_bdy)(floor(Npnts_bdy(l)/2)+1:Npnts(l),:);
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l+1);
                    NORMALS.(name_normals1) = POINTS.(name_bdy1) / R(l+1);
                end
            else
                error('Unknown boundary point distribution style')
            end
        end
        
        
    case 'mfs'
        if lNpnts==1
            Npnts = Npnts*ones(N,1);
        end
%         if N ~= lNpnts
%             error('R and Npnts must have the same lenght for ''mfs'' solvertype')
%         end
        for l = 1:N
            lstring = num2str(l);
            lstring1 = num2str(l+1);
            name_int = strcat('int', lstring);
            name_bdy = strcat('bdy', lstring, lstring);
            name_bdy1 = strcat('bdy', lstring, lstring1);
            name_normals = strcat('n', lstring, lstring);
            name_normals1 = strcat('n', lstring, lstring1);
            % Interior points are not relevant
            POINTS.(name_int) = [];
            % Boundary points only
            if strcmp(bdyptstype,'spiral')
                POINTS.(name_bdy) = SphereSurfGoldPoints(Npnts(l), R(l));
                if l == N
                    % The splitting of the set of point in two subsets is not
                    % required for the outermost layer
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l);
                else
                    % Split bountary points into two sets
                    POINTS.(name_bdy1) = POINTS.(name_bdy)(1:2:end,:);
                    POINTS.(name_bdy) = POINTS.(name_bdy)(2:2:end,:);
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l);
                    NORMALS.(name_normals1) = POINTS.(name_bdy1) / R(l);
                end
            elseif strcmp(bdyptstype,'halton')
                POINTS.(name_bdy) = SphereSurfHaltonPoints(Npnts(l),R(l));
                if l == N
                    % The splitting of the set of point in two subsets is not
                    % required for the outermost layer
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l);
                else
                    % Split bountary points into two sets
                    POINTS.(name_bdy1) = POINTS.(name_bdy)(1:floor(Npnts(l)/2),:);
                    POINTS.(name_bdy) = POINTS.(name_bdy)(floor(Npnts(l)/2)+1:Npnts(l),:);
                    NORMALS.(name_normals) = POINTS.(name_bdy) / R(l);
                    NORMALS.(name_normals1) = POINTS.(name_bdy1) / R(l);
                end
            else
                error('Unknown boundary point distribution style')
            end
        end
    otherwise
        error('Method not recognized: check the input string')
end

end