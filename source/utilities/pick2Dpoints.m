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
%                      'wam' - See wamquadrangle
%           ep - (only for 'inner' points) shape parameter
%   Outputs: x - 2D points with spaceopt distribution
%            spacestr - string explaining your point choice
%
% N should be a 2-vector, but if N is scalar, we reset N=[N,N]
% Note that the Halton points only really need prod(N)

if nargin<5
    ep=0;
    if nargin<4
        spaceopt='even';
    end
end

if length(N)==1
    N=[N,N];
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
        pN = prod(N);
        xh = haltonseq(pN,2);
        x = (repmat(b,pN,1)-repmat(a,pN,1)).*xh+repmat(a,pN,1);
        spacestr=' Halton points';
    case 'wam'
        % Here, centered symmetric (actually Chebyshev points), but can be
        % used more flexibly, i.e., arbitrary quadrilateral
        x = wamquadrangle(N(1)-1,[a(1) a(2); b(1) a(2); b(1) b(2); a(1) b(2)]);
        spacestr = ' WAM quadrilateral';
    otherwise
        error('Unrecognized spaceopt=%s',spaceopt)
end