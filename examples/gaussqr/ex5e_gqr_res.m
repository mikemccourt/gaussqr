function fu = ex5e_gqr_res(u,x,uold,dt,h,b)
% This evaluates the residual of the Geirer Meinhardt system
% The spatial discretization is finite differences
% The temporal discretization is Crank-Nicolson
%
% function fu = ex5e_gqr_res(x)
% Inputs : x - location of the FD discretization
% Outputs : fu - initial condition
%
% function fu = ex5e_gqr_res(u,x,uold,dt)
% Inputs : u - values of the activator and inhibitor
%          x - locations of the u values
%          uold - values at the previous time step
%          dt - time step
% Outputs : fu - residual
%
% function Jfu = ex5e_gqr_res(u,x,uold,dt,jac)
% Inputs : u - values of the activator and inhibitor
%          x - locations of the u values
%          uold - values at the previous time step
%          dt - time step
%          jac - pass 1 to return a preconditioning Jacobian
% Outputs : Jfu - the preconditioner (just a linearized operator)
%
% function Jfub = ex5e_gqr_res(u,x,uold,dt,h,b)
% Inputs : u - values of the activator and inhibitor
%          x - locations of the u values
%          uold - values at the previous time step
%          dt - time step
%          h - finite differencing parameter
%          b - vector for Jacobian-vector multiplication
% Outputs : Jfub - the approximate Jacobian-vector multiplication

persistent ix_RBC ix_LBC ix_TBC ix_BBC ix_INT ix_BC D2 BC A H Aold Hold AoH I

N2 = size(u,1)/2;
N = sqrt(N2);

% These are the solve parameters
tau = .1;
ep = .04;
kappa = .0128;

if length(D2) == 0 & nargin > 1
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
    BC = sparse(N2,N2);
    Bstenoff5 = [-4 -3 -2 -1 0];
    for k=ix_RBC % x=1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = -stencil5;
    end
    Bstenoff5 = [0 1 2 3 4];
    for k=ix_LBC % x=-1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = stencil5;
    end
    Bstenoff5 = N*[-4 -3 -2 -1 0];
    for k=ix_TBC % y=1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = stencil5;
    end
    Bstenoff5 = N*[0 1 2 3 4];
    for k=ix_BBC % y=-1 BC
        Bloc = Bstenoff5 + k;
        BC(k,Bloc) = stencil5;
    end
    
    % Create an identity matrix, because I'll need it later
    % Also, define the Jacobian if requested
    I = speye(N2,N2);
    Jfu = sparse(2*N2,2*N2);
end

switch nargin
    case 1
        A = .5*(1+.001*sum(cos(pi/2*u(:,2)*(1:20)),2)).*sech(.5*sqrt(u(:,1).^2+u(:,2).^2)/ep).^2;
        H = cosh(1-sqrt(u(:,1).^2+u(:,2).^2))/(3*cosh(1));
        fu = [A;H];
    case 4
        A = u(1:N2);
        H = u(N2+1:end);
        Aold = uold(1:N2);
        Hold = uold(N2+1:end);

        % Satisfy the interior condition
        fu = [1/dt*(A-Aold) + .5*(- ep^2*(D2*(A+Aold)) + (A+Aold) - (A.^2./H+Aold.^2./Hold));
            tau/dt*(H-Hold) + .5*(-kappa*(D2*(H+Hold)) + (H+Hold) - 1/ep*(A.^2+Aold.^2))];

        % Satisfy the boundary condition
        fu([ix_BC,N2+ix_BC]) = zeros(size([ix_BC,N2+ix_BC]));
        fu = fu + [BC*A;BC*H];
    case 5
        A = u(1:N2);
        H = u(N2+1:end);
        AoH = A./H;
        
        % Handle the interior part of the Jacobian
        Jfu(1:N2,1:N2) = (1/dt+.5)*I - spdiags(AoH,0,N2,N2) - ep^2/2*D2;
        Jfu(1:N2,N2+1:end) = .5*spdiags(AoH.^2,0,N2,N2);
        Jfu(N2+1:end,1:N2) = -1/ep*spdiags(A,0,N2,N2);
        Jfu(N2+1:end,N2+1:end) = (tau/dt+.5)*I - .5*kappa*D2;
        
        % Handle the boudary part of the Jacobian
        Jfu(ix_BC,:) = sparse(length(ix_BC),2*N2);
        Jfu(N2+ix_BC,:) = sparse(length(ix_BC),2*N2);
        Jfu(ix_BC,1:N2) = BC(ix_BC,:);
        Jfu(N2+ix_BC,N2+1:end) = BC(ix_BC,:);
        
        fu = Jfu; % Return the Jacobian
    case 6
        if size(b,2)==1
            fuh = ex5e_gqr_res(u,x,uold,dt);
            fuhb = ex5e_gqr_res(u+h*b,x,uold,dt);
            fu = 1/h*(fuhb - fuh);
        else
            fu = sparse(size(b,1),size(b,2));
            fuh = ex5e_gqr_res(u,x,uold,dt);

            k = 1;
            for bk=b
                fuhb = ex5e_gqr_res(u+h*bk,x,uold,dt);
                fu(:,k) = 1/h*(fuhb - fuh);
                k = k+1;
            end
        end
    otherwise
        error('Cannot call this with %d inputs',nargin)
end