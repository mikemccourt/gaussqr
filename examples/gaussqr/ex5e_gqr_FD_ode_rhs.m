function [fu,M] = ex5e_gqr_FD_ode_rhs(x,u,t)
% This is the problem phrased as
%   Mass*dy/dt = ex5e_gqr_FD_ode_rhs(x,u,t)
% The mass matrix is singular because of the BC
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
%
% fu = ex5e_gqr_FD_ode_rhs(x,u,t)
% Inputs : x - grid discretization
%          u - state values for function evaluation
%          t - time for function evaluation
% Output : du - derivative evaluation

persistent ix_RBC ix_LBC ix_TBC ix_BBC ix_INT ix_BC D2 BC A H Aold Hold AoH I Jfu JA J1

% Size of the problem
N2 = size(x,1);
N = sqrt(N2);

% These are the solve parameters
tau = .1;
ep = .04;
kappa = .0128;

% Setup the problem
if length(D2) == 0
    % This is the 4th order finite difference stencil
    % There is a factor which accounts for the FD width
    stencil9 = (N-1)^2/4*[1/6 2/3 1/6  2/3 -10/3 2/3  1/6 2/3 1/6];
    stencil5 = (N-1)/2*[-25/12 4 -3 4/3 -1/4];

    % Determine where the BC are
    ix_RBC = find(x(:,1)==max(x(:,1)))'; % Right BC
    ix_LBC = find(x(:,1)==min(x(:,1)))'; % Left BC
    ix_TBC = setdiff(find(x(:,2)==max(x(:,2)))',[ix_RBC,ix_LBC]); % Top BC
    ix_BBC = setdiff(find(x(:,2)==min(x(:,2)))',[ix_RBC,ix_LBC]); % Bottom BC
    ix_BC = sort([ix_RBC,ix_LBC,ix_TBC,ix_BBC]);
    ix_INT = setdiff(1:N2,ix_BC); % Interior

    % Build the 4th order 2nd derivative matrix
    D2 = sparse(N2,N2);
    Dstenoff9 = [-N-1 -N -N+1  -1 0 1  N-1 N N+1];
    for k=ix_INT
        Dloc = Dstenoff9 + k;
        D2(k,Dloc) = stencil9;
    end

    % Build the boundary condition operator (4th order)
    % NOTE: size(BC) = [length(ix_BC) , N2] by the end of this
    BC = sparse(N2,N2);
    Bstenoff5 = [0 1 2 3 4];
    for k=ix_LBC % x=-1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = stencil5;
    end
    Bstenoff5 = [-4 -3 -2 -1 0];
    for k=ix_RBC % x=1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = -fliplr(stencil5);
    end
    Bstenoff5 = N*[0 1 2 3 4];
    for k=ix_BBC % y=-1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = stencil5;
    end
    Bstenoff5 = N*[-4 -3 -2 -1 0];
    for k=ix_TBC % y=1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = -fliplr(stencil5);
    end
    BC = BC(ix_BC,:); % Compress out the interior rows
    
    % Create an identity matrix, because I'll need it later
    % Also, define the Jacobian, and subJacobian for A
    I = speye(N2,N2);
    Jfu = sparse(2*N2,2*N2);
    JA = sparse(N2,N2);
    
    % Set the constant parts of the Jacobian here
    % First the interior
    Jfu(N2+1:end,N2+1:end) = kappa*D2 - I;
    % Now the boundary
    Jfu(ix_BC,1:N2) = BC;
    Jfu(N2+ix_BC,N2+1:end) = BC;
    % The part of JA that doesn't vary with u
    JA = ep^2*D2 - I;
    JA(ix_BC,:) = sparse(length(ix_BC),N2);
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
        % Make sure to zero out the BC points
        A([ix_BC,ix_INT]) = [zeros(length(ix_BC),1);u(ix_INT)];
        H = u(N2+1:end);
        AoH = A./H;
        
        % Handle the interior part of the Jacobian
        % Only the part that is non-constant is updated here
        % Note that the boundary is constant, and is not updated
        J1 = JA + spdiags(2*AoH,0,N2,N2);
        Jfu(ix_INT,1:N2) = J1(ix_INT,:);
        Jfu(1:N2,N2+1:end) = -spdiags(AoH.^2,0,N2,N2);
        % Only count the points on the interior
        Jfu(N2+1:end,1:N2) = 2/ep*spdiags(A,0,N2,N2);
        
        fu = Jfu; % Return the Jacobian
    case 3
        A = u(1:N2);
        H = u(N2+1:end);

        % Satisfy the interior condition
        fu = [ep^2*(D2*A) - A + A.^2./H;
              kappa*(D2*H) - H + 1/ep*A.^2];

        % Satisfy the boundary condition
        fu([ix_BC,N2+ix_BC]) = [BC*A;BC*H];
    otherwise
        error('Illegal number of inputs: nargin=%d',nargin)
end