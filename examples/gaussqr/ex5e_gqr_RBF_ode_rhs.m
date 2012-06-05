function [fu,M] = ex5e_gqr_RBF_ode_rhs(x,u,t)
% This is the problem phrased as
%   Mass*dy/dt = ex5e_gqr_RBF_ode_rhs(x,u,t)
% The mass matrix is singular because of the BC
% Here I'm solving it with an RBF discretization using Gaussian
% eigenfunctions to avoid ill-conditioning
%
% You can call this in one of several ways:
% 
% [u,M] = ex5e_gqr_FD_ode_rhs(x)
% Returns the initial conditions of the system and the mass matrix
% Input  : x - grid discretization
%          M - mass matrix
% Output : u - Initial conditions [A;H]
%
% Jfu = ex5e_gqr_FD_ode_rhs(x,u)
% Inputs : x - grid discretization
%          u - state values for function evaluation
% Output : Jfu - Jacobian evaluation
% NOTE: Jacobian evaluation is not allowed right now
%
% fu = ex5e_gqr_FD_ode_rhs(x,u,t)
% Inputs : x - grid discretization
%          u - state values for function evaluation
%          t - time for function evaluation
% Output : du - derivative evaluation

persistent ix_dxBC ix_dyBC ix_INT A H AoH AH_dx AH_dy AH_L A_int H_int

% Size of the problem
N2 = size(x,1);
N = sqrt(N2);

% These are the simulation parameters
tau = .1;
ep = .04;
kappa = .0128;

% Setup the problem
if length(ix_dxBC) == 0
    % Determine where the BC are
    ix_dxBC = sort([find(x(:,1)==max(x(:,1)))',find(x(:,1)==min(x(:,1)))']);% Right or Left BC
    ix_dyBC = setdiff([find(x(:,2)==max(x(:,2)))',find(x(:,2)==min(x(:,2)))'],ix_dxBC);% Top or Bottom BC
    ix_INT = setdiff(1:N2,[ix_dxBC,ix_dyBC]); % Interior
    
    % Set up the RBF derivative function
    ApplyDerivative(x,ix_dxBC,ix_dyBC);
end

switch nargin
    case 1
        A = .5*(1+.001*sum(cos(pi/2*x(:,2)*(1:20)),2)).*sech(.5*sqrt(x(:,1).^2+x(:,2).^2)/ep).^2;
        H = cosh(1-sqrt(x(:,1).^2+x(:,2).^2))/(3*cosh(1));
        fu = [A;H];
        
        if nargout==2 % If the mass matrix is requested
            Mvec = zeros(2*N2,1);
            Mvec(ix_INT) = ones(length(ix_INT),1);
            Mvec(ix_INT+N2) = tau*ones(length(ix_INT),1);
            M = spdiags(Mvec,0,2*N2,2*N2);
        end
    case 2
        error('No Jacobian available yet')
    case 3
        A = u(1:N2);
        H = u(N2+1:end);
        A_int = A(ix_INT);
        H_int = H(ix_INT);
        
        % Take the derivatives of the data
        % This automatically satisfies the boundary conditions
        fu = ApplyDerivative(x,[A,H]);
        
        % We still must satisfy the interior conditions
        fu(ix_INT)    =  ep^2*fu(ix_INT)    - A_int + A_int.^2./H_int;
        fu(ix_INT+N2) = kappa*fu(ix_INT+N2) - H_int + 1/ep*A_int.^2;
    otherwise
        error('Illegal number of inputs: nargin=%d',nargin)
end

end % function ex5e_gqr_RBF_ode_rhs

function du = ApplyDerivative(x,ind_dxBC,ind_dyBC)
% This function serves as a differentiation matrix for the RBF system
% When evaluating the derivatives, pass values [A,H] but it returns values
% of the form [A;H]
%
% It evaluates the Laplacian on the interior points, and the appropriate
% derivatives on the boundaries
%
% It returns a vector such that
% du(ix_INT) - Laplacian of A
% du(ix_dxBC) - x derivative BC of A
% du(ix_dyBC) - y derivative BC of A
% du(ix_INT+N^2) - Laplacian of H
% du(ix_dxBC+N^2) - x derivative BC of H
% du(ix_dyBC+N^2) - y derivative BC of H
%
% There are two ways to call this function
% To initialize this function, you need to pass 5 inputs
% function ApplyDerivative(x,ind_dxBC,ind_dyBC)
% Inputs : x - computational grid
%          ind_dxBC - list of the dx boundary conditions
%          ind_dyBC - list of the dy boundary conditions
%
% function du = ApplyDerivative(x,u)
% Inputs : x - computational grid
%          u - function values in the form [A,H]
% Output : du - derivative of u as described above
    global GAUSSQR_PARAMETERS    
    ep = GAUSSQR_PARAMETERS.MY_EPSILON;
    alpha = GAUSSQR_PARAMETERS.MY_ALPHA;
    
    persistent phi_Q phi_R  phiDx  phiDy  phiL  ix_dxBC  ix_dyBC  ix_INT
    
    N2 = size(x,1);
    
    if nargin==3 & nargout==0
        % Set up Marr, and check if nonsensical ep or alpha passed
        [ep,alpha,Marr] = gqr_solveprep(1,x,ep,alpha);
        
        phi = gqr_phi(Marr,x,ep,alpha);
        [phi_Q,phi_R] = qr(phi,0);
        
        ix_dxBC = ind_dxBC;
        ix_dyBC = ind_dyBC;
        ix_INT = setdiff(1:N2,[ix_dxBC,ix_dyBC]);
        
        phiDx = gqr_phi(Marr,x(ix_dxBC,:),ep,alpha,[1 0]);
        phiDy = gqr_phi(Marr,x(ix_dyBC,:),ep,alpha,[0 1]);
        phiL = gqr_phi(Marr,x(ix_INT,:),ep,alpha,[2 0]) + gqr_phi(Marr,x(ix_INT,:),ep,alpha,[0 2]);
    elseif nargin==2 & nargout==1
        u = ind_dxBC; % Change the name of what is getting passed
        du = zeros(2*N2,1);
        
        opts.UT = true;
        beta = linsolve(phi_R,phi_Q'*u,opts);
        
        AH_dx = phiDx*beta;
        AH_dy = phiDy*beta;
        AH_L  =  phiL*beta;
        
        du(ix_INT) = AH_L(:,1);
        du(ix_dxBC) = AH_dx(:,1);
        du(ix_dyBC) = AH_dy(:,1);
        du(ix_INT+N2) = AH_L(:,2);
        du(ix_dxBC+N2) = AH_dx(:,2);
        du(ix_dyBC+N2) = AH_dy(:,2);
    else
        error('Unacceptable calling sequence, nargin=%d and nargout=%d',nargin,nargout)
    end
end