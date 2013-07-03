function PHI = HSeigsolveqd(N,B,M)
%This function approximates Hilbert-Schmidt eigenvalues. 
%The difference between this method and HSeigsolve is this function using
%quadrature to help us.
%
%function PHI = HSeigsolve(N,B,M,qdopts)
%
%Inputs :  N - number of points in the domain
%          B - choice of approximating basis
%          M - different quadrature method to compute the error
%          qdopts - quadrature related options
%Outputs : PHI - eigenfunction object 
%
%
%
% Should call qdopts = qdoptCHECK(M,qdopts)
%   qdopts will have different structure for each choice of M
%      M=1 - qdopts.npts is the number of points to use
%            qdopts.ptspace is the distribution of quadrature points

L=1;

PHI.N = N; %create object PHI
%pick basis
switch B
    case 1
        PHI.basisName = 'Standard Polynomial';
        ptspace = 'cheb';
        K_F= @(x,z,j) min(x,z)-x.*z;
                     
        H_mat = @(x,z,j) x.^(j-1);
        x = pickpoints(0,L,N+2,ptspace);
        x = x(2:end-1);
        j = 1:N;
        z =[];
        Z =[];
        %x =0:1/(N+1):1;
        %x = x(2:end-1)';
       
       
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
        [eival,eivec] = eig(K*W);
       
    case 2
        PHI.quad = 'quadgk';
        for l = 1:N
            for i = 1:N
                K(l,i) = quadgk(@(p) H_mat(p,z,i).*K_F(x(l),p,i),0,1);
            end
        end
        [eival,eivec] = eig(K);
    otherwise
        error('Unacceptable quadrature method=%e',M);
end
 [~,ix] = sort(diag(eival),'descend');
 eivec = eivec(:,ix);
 eival = eival(ix,ix);
 eivec = H\eivec;
 PHI.eigvals = diag(eival);
 PHI.coefs = eivec;
% SiK = size(K);
% SiW = size(W);
% if SiK(1) == SiW(2)
    
% else
%     L = chol(H,'lower');
%     [eival,eivec] = eig(L\K*W*L);
%     [~,ix] = sort(diag(eival),'descend');
%     eivec = eivec(:,ix);
%     eival = eival(ix,ix);
%     %eivec = L'\eivec;
%     PHI.eigvals = diag(eival);
%     PHI.coefs = eivec;
% end





