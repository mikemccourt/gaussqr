% HSSVD_IBBSolve_Full
% This function allows the user to solve an interpolation problem with the
% Iterated Brownian Bridge kernel using the HSSVD basis

function yeval = HSSVD_IBBSolve_Full(ep,beta,x,y,xeval)
% Evaluate the HS-SVD basis IBB interpolant
%
% Requires a positive ep shape parameter and positive integer beta
% smoothness parameter
% Data locations x and data values y should be equal length column vectors
%
% function yeval = HSSVD_IBBSolve_Full(ep,beta,x,y,xeval)
%   Inputs:  ep    - shape parameter
%            beta  - smoothness parameter
%            x     - data locations
%            y     - data values
%            xeval - locations at which to evaluate the interpolant
%   Outputs: yeval - interpolant values at xeval
%
% Note that multiple sets of data values can be passed for simultaneous
% computation, e.g., y = [y1 y2 y3], seval = [seval1 seval2 seval3]
%
%
%
% This function can be called a second way, depending on if the user has
% already solved the interpolation problem once and just wants to evaluate
% at a different set of points
%
% function yeval = HSSVD_IBBSolve_Full(xeval)
%   Inputs:  xeval - locations at which to evaluate the interpolant
%   Outputs: yeval - interpolant values at xeval
%
% You can only call the function in this form if you have already performed
% a solve.  To solve the linear system but not perform evaluations:
%
% function HSSVD_IBBSolve_Full(ep,beta,x,y)

% Define an IBB data structure which will store the important solve content
% in it for use at a later time
% This is overwritten if the user passes new data to the problem
persistent IBB

% Perform some basic error testing to confirm things are acceptable
% Unpack the arguments passed to the function
% Solve for the IBB object if the appropriate data is passed
% Can also use varargin to manage input arguments if preferred
switch nargin
    case 1
        if isempty(IBB)
            error('Cannot evaluate because no coefficients have been stored')
        else
            % Rename the first info passed to the function
            xeval = ep;
        end
    case {4,5} 
        if nargin==4
            if nargout>0
                error('No evaluation points requested ... nothing can be returned')
            end
        end
        if numel(ep)>1 || numel(beta)>1
            error('Must have scalar ep and beta; check calling sequence')
        end
        if ep~=abs(ep)
            error('Shape parameter ep=%g unacceptable',ep)
        end
        if beta~=abs(floor(beta))
            error('Smoothness parameter beta=%g unacceptable',beta)
        end
        if size(x,1)~=size(y,1)
            error('Input data must be the same size, size(x,1)=%d, size(y,1)=%d',size(x,1),size(y,1))
        end
        if size(x,2)~=1
            error('The IBB kernel can only handle 1D problems')
        end
        
        % If we have reached this point, we can form the HSSVD basis
        % The try block and temporary storage prevents the IBB object from
        % being overwritten in the event of a faulty call
        try
            IBB_new = HSSVD_IBBForm(ep,beta,x,y);
        catch err
            rethrow(err);
        end
        IBB = IBB_new;
    otherwise
        error('Unacceptable inputs, nargin=%d',nargin)
end

% If there are points to be evaluated, do so
if exist('xeval','var')
    yeval = HSSVD_IBBEval(IBB,xeval);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These function can only be called within this Matlab file because they
% have been declared within a file of a different name

function IBB = HSSVD_IBBForm(ep,beta,x,y)
N = size(x,1);

% Define the eigenfunction and eigenvalue handles
phifunc = @(n,x) sqrt(2)*sin(pi*x*n);
lamfunc = @(b,e,n) ((pi*n).^2+e^2).^(-b);

% Determine how many eigenfunctions will be needed for accuracy
M = ceil(1/pi*sqrt(eps^(-1/beta)*(N^2*pi^2+ep^2)-ep^2));
narr = 1:M;

% Evaluate the Phi and lamvec
Phi1 = phifunc(narr(1:N),x);
Phi2 = phifunc(narr(N+1:end),x);
lamvec1 = lamfunc(beta,ep,narr(1:N));
lamvec2 = lamfunc(beta,ep,narr(N+1:end));

% Create the CbarT object
CbarT = bsxfun(@rdivide,lamvec2',lamvec1).*(Phi2'/Phi1');

% Form the Psi matrix and solve for the interpolation coefficients
Psi = Phi1 + Phi2*CbarT;
c = Psi\y;

% Pack the IBB object
IBB.CbarT = CbarT;
IBB.c = c;
IBB.narr = narr;
IBB.phifunc = phifunc;
end

function yeval = HSSVD_IBBEval(IBB,xeval)
% Unpack the IBB object
CbarT = IBB.CbarT;
c = IBB.c;
narr = IBB.narr;
phifunc = IBB.phifunc;
N = size(CbarT,2);

% Evaluate the Phieval components
Phieval1 = phifunc(narr(1:N),xeval);
Phieval2 = phifunc(narr(N+1:end),xeval);
yeval = Phieval1*c + Phieval2*(CbarT*c);
end