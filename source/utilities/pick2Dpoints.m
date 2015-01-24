function [x,spacestr] = pick2Dpoints(a,b,N,spaceopt,ep)
% function [x,spacestr] = pick2Dpoints(a,b,N,spaceopt,ep)
% This function returns a matrix of 2D points, each row is a point
%
%   Inputs: a  - [dim1 lower bound,   dim2 lower bound]
%           b  - [dim1 upper bound,   dim2 upper bound]
%           N  - [num points in dim1, num points in dim2]
%                Note: N scalar is treated as [N,N]
%           spaceopt - (optional, default='even') choice of point distribution
%                      'even' - uniformly spaced points
%                      'cheb' - Chebyshev tensor grid
%                      'inner' - centrally clustered points
%                      'halton' - Halton quasi-random points
%                      'rand' - uniform random points
%                      'wam' - See wamdisk
%                      'wamq' - See wamquadrangle
%                               Only the first value in N is considered for
%                               the wam options
%           ep - (only for 'inner' points) shape parameter
%   Outputs: x - 2D points with spaceopt distribution
%            spacestr - string explaining your point choice
%
% N should be a 2-vector, but if N is scalar, we reset N=[N,N]
% Note that the Halton points only really need prod(N)
% Also, the Halton points do NOT include [0,0]
%
% If you pass scalar a & b, the same value is used in both dimensions:
%     passing a=1 & b=2 is equivalent to a=[1 1] & b=[2 2]
%
% PROGRAMMERS NOTE: More of this should be rewritten with bsxfun
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
stats_avail = GAUSSQR_PARAMETERS.STATISTICS_TOOLBOX_AVAILABLE;

if nargin<5
    ep=0;
    if nargin<4
        spaceopt='even';
    end
end

if length(N)==1
    N=[N,N];
end

% Orient the input boundaries the way we normally think of them
a = a(:)';
b = b(:)';
Na = length(a);
Nb = length(b);

% Allow the user to input scalars instead of arrays and use the same value
% for both points in the array
if Na==Nb
    if any(a>b)
        error('Lower bounds must be less than upper bounds')
    else
        if Na==1
            a = [a a];
            b = [b b];
        elseif Na~=2
            error('Bounds passed in more than 2 dimensions')
        end
    end
else
    error('Domain boundaries must be passed with same size')
end

switch lower(spaceopt)
    case {'even','unif'}
        x1 = repmat(linspace(a(1),b(1),N(1))',1,N(2));
        x2 = repmat(linspace(a(2),b(2),N(2)),N(1),1);
        x = [x1(:),x2(:)];
        spacestr=' evenly spaced';
    case 'inner'
        sigma = 1/(2*sqrt(ep^2+.5*sqrt(4*ep^4+1)));
        x1 = icdf('Normal',linspace(1/(2*N(1)-1),(2*N(1)-1)/(2*N(1)),N(1)),0,sigma)';
        xmax = max(x1);xmin = min(x1);xscale = (b(1)-a(1))/(xmax-xmin);xpt = a(1)-xscale*xmin;
        x1m = repmat(xpt + xscale*x1,1,N(2));
        x2 = icdf('Normal',linspace(1/(2*N(2)-1),(2*N(2)-1)/(2*N(2)),N(2)),0,sigma)';
        xmax = max(x2);xmin = min(x2);xscale = (b(2)-a(2))/(xmax-xmin);xpt = a(2)-xscale*xmin;
        x2m = repmat(xpt + xscale*x2',N(1),1);
        x = [x1m(:),x2m(:)];
        spacestr=' central cluster';
    case 'cheb'
        x1 = -.5*(b(1)-a(1))*cos(pi*(0:N(1)-1)/(N(1)-1))'+.5*(b(1)+a(1));
        x2 = -.5*(b(2)-a(2))*cos(pi*(0:N(2)-1)/(N(2)-1))'+.5*(b(2)+a(2));
        x1s = repmat(x1,1,N(2));
        x2s = repmat(x2',N(1),1);
        x = [x1s(:),x2s(:)];
        spacestr=' Chebyshev points';
    case {'halton','halt'}
        % Add the round in there to make sure an integer is used
        pN = round(prod(N));
        if stats_avail && exist('haltonset','file')
            point_generator = haltonset(2,'Skip',1);
            xh = net(point_generator,pN);
        else
            xh = haltonseq(pN,2);
        end
        x = bsxfun(@plus,bsxfun(@times,b-a,xh),a);
        spacestr=' Halton points';
    case {'random','rand'}
        pN = round(prod(N));
        x = bsxfun(@plus,bsxfun(@times,b-a,rand(pN,2)),a);
        spacestr=' uniform random points';
    case 'wamq'
        % Here, centered symmetric (actually Chebyshev points), but can be
        % used more flexibly, i.e., arbitrary quadrilateral
        x = wamquadrangle(N(1)-1,[a(1) a(2); b(1) a(2); b(1) b(2); a(1) b(2)]);
        spacestr = ' WAM quadrilateral';
    case 'wam'
        x = bsxfun(@plus,bsxfun(@times,(b-a)/2,wamdisk(N(1)-1)),a+(b-a)/2);
        spacestr = ' WAM disk';
    otherwise
        error('Unrecognized spaceopt=%s',spaceopt)
end