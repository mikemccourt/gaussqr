function [V,D,A,H] = HSeigsolve(N,basis,rescale,n_eig_plot)
% This function approximates Hilbert-Schmidt eigenvalues
% of the Compact Matern kenrel with ep=0 and beta=1
%
% function PHI = HSeigsolve(N,basis,rescale,n_eig_plot)
% Inputs : N - number of points in the domain
%          basis - choice of approximating basis
%          rescale - <default=1> scale eigenfunctions to sqrt(2)
%          n_eig_plot - <optional> which eigenvalue(s) you want to plot
%                       Pass more than one as e.g., [1,4,6]
% Outputs : PHI - eigenfunction object (described below)
%
% function [V,D,A,H] = HSeigsolve(N,basis,rescale,n_eig_plot)
% Outputs : V - vectors of interpolating coefficients
%           D - Hilbert-Schmidt eigenvalues
%           A - integral collocation matrix
%           H - eigenfunction basis evaluation matrix
%
% The input basis can take the following values
%    1 - Standard polynomial basis, Chebyshev points
%    2 - CMatern basis (ep=0,beta=1), Uniform points
%
% In this function, we restrict L=1, although we could change that
%
% PHI Object details:
%
%    PHI.N         : Number of basis functions
%    PHI.basisName : Which basis is used for the approximation
%    PHI.centers   : What are the kernel centers for the basis
%    PHI.indices   : How are the kernels indexed
%    PHI.basisEval : Function for evaluating the basis
%                    calling sequence: p = PHI.basisEval(x,z,j)
%                    x - points at which to evaluate the basis
%                        x = [x1;x2;x3;...]
%                    z - centers for the basis
%                        z = [z1,z2,z3,...]
%                    j - which basis function to evaluate
%                        j = [j1,j2,j3,...]
%                    p - basis function evaluation
%                        p = [h1(x1),h2(x1),...;h1(x2),h2(x2),...]
%    PHI.eigvals   : The computed HS eigenvalues
%    PHI.reorient  : Function to make the eigenfunctions line up
%                    calling sequence: RV = PHI.basisEval(V)
%                    V - eigenfunction coefficients
%                        V = [v1,v2,...]
%                    RV - reoriented coefficients
%    PHI.coefs     : Coefficients for evaluating the eigenfunctions
%                    PHI.coefs(:,k) are for the kth eigenfunction
%   PHI.eigfuncEval: Function for evaluating the eigenfunctions
%                    calling sequence: f = PHI.eigfuncEval(x,k)
%                    x - points at which to evaluate the eigenfunction
%                        x = [x1;x2;x3;...]
%                    j - which eigenfunctions to evaluate
%                        k = [k1,k2,k3,...]
%                    f - eigenfunction evaluation
%                        f = [f1(x1),f2(x1),...;f1(x2),f2(x2),...]

L = 1;

% Account for inputs the user chooses
if nargin<3
    rescale = 1;
    if nargin<4
        n_eig_plot = 0;
    end
end
if length(rescale)==0
    rescale = 1;
end

% Determine if the user wants the PHI object back
returnPHI = 0;
if nargout==1
    returnPHI = 1;
end

% Create the PHI object, which we will use in this function
PHI.N = N;

% The Rescale function in here is used to make sure that the eigenfunctions
% are pointing in the proper setting.
switch basis
    case 1
        PHI.basisName = 'Standard Polynomial';
        
        Int_Kh = @(x,z,j) -1./(j.^2+j).*x.*(x.^j-1);
        H_mat = @(x,z,j) x.^(j-1);
        
        ptspace = 'cheb';
        x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
        X = repmat(x,1,N);
        z = [];
        Z = [];
        j = 1:N;
        J = repmat(j,N,1);
    case 2
        PHI.basisName = 'PP Spline Kernel';
        
        Int_Kh = @(x,z,j) (x<=z).*(1/3*(1-x).*(1-z).*x.^3 + 1/3*x.*z.*(1-z).^3 + x.*(1-z).*z.^2.*(1/2-z/3) - x.*(1-z).*x.^2.*(1/2-x/3))    +...
                          (x>z) .*(1/3*(1-x).*(1-z).*z.^3 + 1/3*x.*z.*(1-x).^3 + z.*(1-x).*x.^2.*(1/2-x/3) - z.*(1-x).*z.^2.*(1/2-z/3));
        % Below is the symmetric version of the A matrix evaluation
        Int_Kh = @(x,z,j) 1/3*(x.*z.*(1-max(x,z)).^3 + (1-x).*(1-z).*min(x,z).^3) +...
                          min(x,z).*(1-max(x,z)).*(max(x,z).^2.*(1/2-1/3*max(x,z)) - min(x,z).^2.*(1/2-1/3*min(x,z)));
        H_mat = @(x,z,j) min(x,z) - x.*z;
        
        ptspace = 'even';
        x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
        X = repmat(x,1,N);
        z = x';
        Z = repmat(z,N,1);
        j = [];
        J = [];
    otherwise
        error('Unacceptable basis=%e',basis)
end

% Store the data chosen by basis
PHI.centers   = z;
PHI.indices   = j;
PHI.basisEval = H_mat;

% Create the matrices needed for the eigenfunction problem
A = Int_Kh(X,Z,J);
H = H_mat(X,Z,J);

% Solve for the eigenfunction coefficients
[V,D] = eig(A,H);

% Sort the eigenvalues and reorder the matrices appropriately
[d,ix] = sort(diag(D),'descend');
V = V(:,ix);
D = D(ix,ix);

% Store the eigenvalues, now that they've been sorted
PHI.eigvals = diag(D);


% Now we need to do some stuff to assume the eigenfunctions can be compared
% to the analytically known eigenfunctions
% We need the eigenfunction to have positive derivative at 0
switch basis
    case 1
        PHI.reorient = @(V) V*diag(sign(V(2,:)));
    case 2
        PHI.reorient = @(V) V*diag(sign((1-PHI.centers)'*V));
end
    
% Reorient the eigenfunctions
V = PHI.reorient(V);

% Rescale the eigenfunctions to have max value sqrt(2)
if rescale~=0
    for k=1:N
        optimopts.Display = 'off';
        [xmax,fval,exitflag] = fminbnd(@(x)-H_mat(x*ones(1,N),z,j)*V(:,k),0,L,optimopts);
        if exitflag<1
            warning('Problems arose while rescaling eigenfunction %d',k)
        end
        V(:,k) = V(:,k)*sqrt(2)/(-fval);
    end
end

% Store the coefficients and a way of evaluating the eigenfunctions
PHI.coefs       = V;
PHI.eigfuncEval = @(x,k) PHI.basisEval(x*ones(1,PHI.N),repmat(PHI.centers,size(x,1),1),repmat(PHI.indices,size(x,1),1))*PHI.coefs(:,k);


% Now just some cleaning up of other stuff
% Plot any eigenfunctions the user wants
if n_eig_plot>0
    Neval = 300;
    xx = pickpoints(0,L,Neval);
    plot(xx,PHI.eigfuncEval(xx,n_eig_plot))
end

% Return the PHI object if the user wants it (which they should)
if returnPHI
    V = PHI;    
end