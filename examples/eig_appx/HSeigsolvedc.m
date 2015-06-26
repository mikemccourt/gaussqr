function PHI = HSeigsolvedc(N,kernel,B,M,numqdp,epsilon)
%This function approximates Hilbert-Schmidt eigenvalues. 
%The difference between this method and HSeigsolve is this function using
%quadrature to help us.
%
%function PHI = HSeigsolve(N,B,M,qdopts)
%
%Inputs :  N      - number of points in the domain
%          kernel - the kernel you want to use
%          B      - choice of approximating basis
%          M      - different quadrature method to compute the eigenvalue
%          M:
%               1 - trapezoid rule, quadature point is same as collocation points
%             
%          epsilon- value of epsilon
%          qdopts - quadrature related options
%          
%Outputs : PHI - eigenfunction object 
%         
%
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
    if kernel ~= 1
        epsilon = 1;
    else
        epsilon = 0;
    end
    if nargin <5
        numqdp = 500; 
    end
end

L=1;
PHI.N = N; %create object PHI
PHI.epsilon = epsilon;
flag = 1;

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
                error('no polynomial basis this time');
           
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
                x = pickpoints(0,L,N+2,ptspace);x = x(2:end-1);
                X = repmat(x,1,N);
                j = 1:N;
                z =[];
                Z =[];
                J = repmat(j,N,1);
            otherwise
                error('Unacceptable basis=%e',B)
        end
        
       switch M
            case 1  
                PHI.quad = 'trapezoid rule with same collcation and quad points';
                if N <2
                    error('To use quadrature N need > 1');
                end
                XQ = repmat(x',N,1);
                w = (x(3:N) - x(1:N-2))/2;
                w = [(x(2)-x(1))/2 w' (x(N)-x(N-1))/2];
            case 2
                PHI.quad = 'simpson rule with same collcation and quad points';
               
                if (mod(N,2) == 1)
                    XQ = repmat(x',N,1);
                    w = zeros(1,N);
                    for i = 1 : 2 : N-2
                        [w1 w2 w3] = wsimpson(x(i),x(i+1),x(i+2));
                        w(i) = w(i) + w1;
                        w(i+1) = w(i+1) + w2;
                        w(i+2) = w(i+2) + w3;
                    end
                else
                    error('N must be odd');
                end
           case 3
                PHI.quad = 'clencurt with same collcation and quad points';
                [x,w] = clencurt(N+1);
                x = (x+1)/2;
                x = x(2:end-1);
                w = w(2:end-1);
                w = w/2;
                X = repmat(x,1,N);
                XQ = repmat(x',N,1);
                
           case 4
                PHI.quad = 'clencurt with different number of quadrature points';
                flag =0;
                [xq,w] = clencurt(numqdp+1);
                xq = (xq+1)/2;
                xq = xq(2:end-1);
                w = w(2:end-1);
                w = w/2;
                X = repmat(x,1,numqdp);
                XQ = repmat(xq',N,1);
                [x y]=size(XQ)
                if B==2
                    z = x';
                    Z = repmat(z,N,1);
                    j = [];
                    J = [];
                elseif B==3
                    j = 1:N;
                    J = repmat(j,numqdp,1);
                    [x y]= size(J)
                    z =[];
                    Z =[];
                    
                end
           otherwise
            error('Unacceptable quadrature method=%e',M);
       end
           
    if(flag)
         K = K_F(X,XQ);   
         H = H_mat(X,Z,J);
         W = diag(w);
         [eivec,eival] = eig(K*W);     
         [esort,ix] = sort(diag(eival),'descend');
         eivec = H\eivec;
         eivec = eivec(:,ix);
         eival = eival(ix,ix);
         PHI.eigvals = diag(eival);
         PHI.coefs = eivec;
    else
         K = K_F(X,XQ);
         H = H_mat(XQ',Z,J);
         W = diag(w); 
         [eivec,eival] = eig(H\K*W*H);     
         [esort,ix] = sort(diag(eival),'descend');
         eivec = eivec(:,ix);
         eival = eival(ix,ix);
         PHI.eigvals = diag(eival);
         PHI.coefs = eivec;
    end
  











