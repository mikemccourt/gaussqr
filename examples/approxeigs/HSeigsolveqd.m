function PHI = HSeigsolveqd(N,kernel,B,M,epsilon,qdopts)
%This function approximates Hilbert-Schmidt eigenvalues. 
%The difference between this method and HSeigsolve is this function using
%quadrature to help us.
%
%function PHI = HSeigsolve(N,B,M,qdopts)
%
%Inputs :  N      - number of points in the domain
%          kernel - the kernel you want to use
%          B      - choice of approximating basis
%          M      - different quadrature method to compute the error
%          epsilon- value of epsilon
%          qdopts - quadrature related options
%          
%Outputs : PHI - eigenfunction object 
%         
% Should call qdopts = qdoptCHECK(M,qdopts)
%
%qdopts will have different structure for each choice of M
%            qdopts.npts        : The number of points to use for quadrature
%                                 methods(not for quadgk or quadl)
%            qdopts.ptspace     : The distribution of quadrature points
%            qdopts.RelTol      : The real tolerance for quadgk
%            qdopts.AbsTol      : the absolute tolerance for quadgk
%            qdopts.integralTOL : the tolerance for quadl if we don't have
%            quadgk
%
%PHI Object details:
%
%    PHI.N         : Number of basis functions
%    PHI.basisName : Which basis is used for the approximation
%    
%    PHI.eigvals   : The computed HSqd eigenvalues
%    PHI.quad      : The quadrature method that were used to compute
%                    eigenvalue
%    PHI.coefs     : Coefficients for evaluating the eigenfunctions
%                    PHI.coefs(:,k) are for the kth eigenfunction
%    

if nargin < 6
    qdopts.npts = N; % create the qdopts object if inputs don't have one.
    if nargin < 5
        epsilon = 1;
    end
end

quadgkEXISTS = 0;
if exist('quadgk')
    quadgkEXISTS = 1;
    RelTol = 1e-8;
    AbsTol = 1e-12;
    if isfield(qdopts,'RelTol')
        RelTol = qdopts.RelTol;
    end
    if isfield(qdopts,'AbsTol')
        AbsTol = qdopts.AbsTol;
    end
else
    integralTOL = 1e-4;
    if isfield(qdopts,'integralTol')
        integralTOL = qdopts.integralTOL;
    end
end


L=1;

PHI.N = N; %create object PHI
%pick kernel
switch kernel
    case 1 
         K_F= @(x,z,j) min(x,z)-x.*z;
    case 2
         K_F= @(x,z,j) sinh(epsilon.*min(x,z)).*sinh(epsilon.*(1-max(x,z)))./(epsilon.*sinh(epsilon));
end

%pick basis
switch B
    case 1
        PHI.basisName = 'Standard Polynomial';
        ptspace = 'cheb';

        H_mat = @(x,z,j) x.^(j-1);
        x = pickpoints(0,L,N+2,ptspace);
        x = x(2:end-1);
        j = 1:N;
        z =[];
        Z =[];
        X = repmat(x,1,N);
        J = repmat(j,N,1);

    case 2
        PHI.basisName = 'PP Spline Kernel';
        H_mat = @(x,z,j) min(x,z)-x.*z;
        ptspace = 'even';
        x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
        X = repmat(x,1,N);
        z = x';
        Z = repmat(z,N,1);
        j = [];
        J = [];
    case 3
        PHI.basisName = 'Chebyshev Polynomials';
        ptspace = 'cheb';

        H_mat = @(x,z,j) cos((j-1).*acos(2*x-1));
        x = pickpoints(0,L,N+2,ptspace);
        x = x(2:end-1);
        j = 1:N;
        z =[];
        Z =[];
        X = repmat(x,1,N);
        J = repmat(j,N,1);
    otherwise
        error('Unacceptable basis=%e',B)
end
%pick quadrature method
switch M
    case 1
        PHI.quad = 'left hand rule';
        Nqd = N;
        if(isfield(qdopts,'npts'))
            Nqd = qdopts.npts;
            if(isfield(qdopts,'ptspace'))
                qdspace = qdopts.ptspace;
            else
                qdspace = 'even';
            end
            v = pickpoints(0,L,Nqd+2,qdspace);
            v = v(2:end-1);
        else
            v = x;
        end
       
        X = repmat(x,1,Nqd);
        J = repmat(j,Nqd,1);
        V = repmat(v,1,N);
        VT = repmat(v',N,1);  
        v(2:Nqd) = v(2:Nqd)-v(1:Nqd-1);
        W = diag(v);
        K = K_F(X,VT,J);
        H = H_mat(V,Z,J);
        PHI.K = K;
        PHI.H = H;
        [eivec,eival] = eig(K*W);
       
    case 2
        PHI.quad = 'quadgk';
        for l = 1:N
            for i = 1:N
               if isempty(z) 
                    if quadgkEXISTS
                        K(l,i) = quadgk(@(p) H_mat(p,z,i).*K_F(x(l),p,i),0,1,'AbsTol',AbsTol,'RelTol',RelTol);
                    else
                        K(l,i) = quadl(@(p) H_mat(p,z,i).*K_F(x(l),p,i),0,1,integralTOL);
                    end
               else
                   if quadgkEXISTS
                        K(l,i) = quadgk(@(p) H_mat(p,z(i),i).*K_F(x(l),p,i),0,1,'AbsTol',AbsTol,'RelTol',RelTol);
                   else
                        K(l,i) = quadl(@(p) H_mat(p,z(i),i).*K_F(x(l),p,i),0,1,integralTOL);
                   end
               end
            end
        end
        H = H_mat(X,Z,J);
        [eivec,eival] = eig(K,H);
        PHI.K = K;
        PHI.H = H;
    otherwise
        error('Unacceptable quadrature method=%e',M);
end
 [esort,ix] = sort(diag(eival),'descend');
 eivec = eivec(:,ix);
 eival = eival(ix,ix);
 PHI.eigvals = diag(eival);
 PHI.coefs = eivec;
  











