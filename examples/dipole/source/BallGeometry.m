function [ POINTS, NORMALS ] = BallGeometry( R, ...
                                             Npnts, ...
                                             solvertype, ...
                                             ptstype , ...
                                             avoidpts, ...
                                             avoidcushion)
% function [POINTS,NORMALS] = BallGeometry(R,Npnts,solvertype,ptstype,avoidpts,avoidcushion)
% BALLGEOMETRY creates a distribution of points in a multilayer (N layers)
% ball.
% Input data are: spheres radii, desired number of interior points, desired
% number of boundary points (optional), meshfree solver type.
% Output data are points coordinates and unit vectors normal to boundary 
% spherical surfaces.
% Interior points are posed in a uniform grid. Boundary points are
% distributed in a regular manner by means of gold section spirals.
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
% ptstype       =   'even'   Points uniformly distributed in the cube
%                   'halton' Points should have a Halton distribution
%                   'random' Points are uniform randomly chosen
%                   ** the default choice is 'halton', which can be
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
%                   ** the default choice is R/30, which can be selected by
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
%
% To do list: Consider a Chebyshev points distribution

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
ptstype_DEFAULT = 'halton';
avoidcushion_DEFAULT = Rend/30;

% Run checks to make sure that the point distribution is okay
if nargin==3
    ptstype = ptstype_DEFAULT;
elseif isempty(ptstype)
    ptstype = ptstype_DEFAULT;
else
    if ischar(ptstype)
        if ~(strcmp(ptstype,'even') || strcmp(ptstype,'halton') || strcmp(ptstype,'random'))
            warning('ptstype=%s unacceptable, defaulting to %s',ptstype,ptstype_DEFAULT)
            ptstype = ptstype_DEFAULT;
        end
    else
        warning('ptstype=%g unacceptable, defaulting to %s',ptstype,ptstype_DEFAULT)
        ptstype = ptstype_DEFAULT;
    end
end

% Run checks to make sure the avoidance values are acceptable
if nargin<5
    avoidpts = [];
else
    if ~isempty(avoidpts)
        avoid_dim = size(avoidpts,2);
        if avoid_dim~=3
            error('Cushion points must be 3D, size(avoidpts,2)=%d',avoid_dim)
        end
        if sum(any(real(avoidpts)~=avoidpts))~=0
            error('Some cushion points are complex ... come on')
        elseif sum(any(sqrt(sum(avoidpts.^2,2))>Rend))~=0
            warning('Some of the cushion points are outside the ball')
        end
        avoid_on = 1;
        
        if not(exist('avoidcushion','var'))
            avoidcushion = avoidcusion_DEFAULT;
        else
            if abs(avoidcushion)~=avoidcushion || avoidcushion==0
                warning('Cushion = %g unacceptable, reverting to %g',avoidcushion,avoidcushion_DEFAULT)
                avoidcushion = avoidcusion_DEFAULT;
            elseif avoidcushion>Rend
                warning('Cushion = %g too large, reverting to %g',avoidcushion,avoidcushion_DEFAULT)
                avoidcushion = avoidcusion_DEFAULT;
            end
        end
    end
end

% Actually determine the point distribution in and on the ball
switch lower(solvertype)
    case 'kansa'
        d = 2*Rend / ( 6/pi * Npnts(1) )^(1/3);
        if strcmp(ptstype,'even')
            x = -Rend:d:Rend;
            [x, y, z] = meshgrid(x);
        elseif strcmp(ptstype,'halton')
            ptsvec = Rend*(2*haltonseq(ceil(6/pi*Npnts(1)),3)-1);
            x = ptsvec(:,1);
            y = ptsvec(:,2);
            z = ptsvec(:,3);
        elseif strcmp(ptstype,'random')
            ptsvec = Rend*(2*rand(ceil(6/pi*Npnts(1)),3)-1);
            x = ptsvec(:,1);
            y = ptsvec(:,2);
            z = ptsvec(:,3);
        else
            error('Unknown point distribution style')
        end
        R = [0 R];
        ll = 0;
        for l = 1:N
            lstring = num2str(l);
            lstring1 = num2str(l+1);
            name_int = strcat('int', lstring);
            name_bdy = strcat('bdy', lstring, lstring);
            name_bdy1 = strcat('bdy', lstring, lstring1);
            name_normals = strcat('n', lstring, lstring);
            name_normals1 = strcat('n', lstring, lstring1);
            % Interior points
            ind = find( (x(:).^2 + y(:).^2 + z(:).^2) < R(l+1)^2 & ...
                        (x(:).^2 + y(:).^2 + z(:).^2) > R(l)^2);
            n = length(ind);
            POINTS.(name_int) = zeros(n,3);
            for j = 1:n
                POINTS.(name_int)(j,:) = [x(ind(j)) y(ind(j)) z(ind(j))];
            end
            % Boundary points
            if lNpnts == 1
                ll = ll + length(POINTS.(name_int));
                Npnts_bdy = ll * .5; % TO BE IMPROVED ?!?!?!
            else
                Npnts_bdy = Npnts(l+1);
            end
            POINTS.(name_bdy) = SphereSurfGoldPoints(Npnts_bdy, R(l+1));
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
        end
    case 'mfs'
        if N ~= lNpnts
            error('R and Npnts must have the same lenght for ''mfs'' solvertype')
        end
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
        end
    otherwise
        error('Method not recognized: check the input string')
end

% Enforce the avoidance condition

end