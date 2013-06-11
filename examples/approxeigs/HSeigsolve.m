function [V,D,A,H] = HSeigsolve(N,basis,n_eig_plot)
% function [V,D] = HSeigsolve(N,basis,n_eig_plot)
% This function approximates Hilbert-Schmidt eigenvalues
% of the Compact Matern kenrel with ep=0 and beta=1
%
% Inputs : N - number of points in the domain
%          basis - choice of approximating basis
%          n_eig_plot - <optional> which eigenvalue you want to plot
% Outputs : V - vectors of interpolating coefficients
%           D - Hilbert-Schmidt eigenvalues
%
% The input basis can take the following values
%    1 - Polynomial basis, Chebyshev points
%    2 - CMatern basis (ep=0,beta=1), Uniform points
%
% In this function, we restrict L=1, although we could change that

L = 1;

if nargin==2
    n_eig_plot = 0;
end

switch basis
    case 1
        Int_Kh = @(x,z,j) -1./(j.^2+j).*x.*(x.^j-1);
        H_mat = @(x,z,j) x.^(j-1);
        ptspace = 'cheb';

        x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
        X = repmat(x,1,N);
        Z = [];
        j = 1:N;
        J = repmat(j,N,1);
    case 2
        Int_Kh = @(x,z,j) (x<=z).*(1/3*(1-x).*(1-z).*x.^3 + 1/3*x.*z.*(1-z).^3 + x.*(1-z).*z.^2.*(1/2-z/3) - x.*(1-z).*x.^2.*(1/2-x/3))    +...
                          (x>z) .*(1/3*(1-x).*(1-z).*z.^3 + 1/3*x.*z.*(1-x).^3 + z.*(1-x).*x.^2.*(1/2-x/3) - z.*(1-x).*z.^2.*(1/2-z/3));
        % Below is the symmetric version of the A matrix evaluation
        Int_Kh = @(x,z,j) 1/3*(x.*z.*(1-max(x,z)).^3 + (1-x).*(1-z).*min(x,z).^3) +...
                          min(x,z).*(1-max(x,z)).*(max(x,z).^2.*(1/2-1/3*max(x,z)) - min(x,z).^2.*(1/2-1/3*min(x,z)));
        H_mat = @(x,z,j) min(x,z) - x.*z;
        ptspace = 'even';

        x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
        X = repmat(x,1,N);
        Z = repmat(x',N,1);
        j = [];
        J = [];
    otherwise
        error('Unacceptable basis=%e',basis)
end

A = Int_Kh(X,Z,J);
H = H_mat(X,Z,J);

[V,D] = eig(A,H);
[d,ix] = sort(diag(D),'descend');
V = V(:,ix);
D = D(ix,ix);

if n_eig_plot>0
    Neval = 300;
    c = V(:,n_eig_plot);
    
    xx = pickpoints(0,1,Neval);
    XX = repmat(xx,1,N);
    ZZ = repmat(x',Neval,1);
    JJ = repmat(j,Neval,1);
    HH = H_mat(XX,ZZ,JJ);
    P = HH*c;
    plot(xx,P)
end