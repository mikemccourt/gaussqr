function [V,D]= EFunction(N,nplot)
    if nargin == 1
        nplot = 0;
    end
    Int_Kh = @(x,j) -1./(j.^2+j).*x.*(x.^j-1);
    x = cos((0:(N+1))/(N+1)*pi)/2+.5;
    j = 1:N;
    J = repmat(j,N,1);
    x = x(2:end-1)';
    X = repmat(x,1,N);
    A = Int_Kh(X,J);
    H_mat = @(z,j) z.^(j-1);
    H = H_mat(X,J);
    HinvA = H\A;
    [V,d] = eig(HinvA);
    [D,index] = sort(diag(d),'descend');
    if nplot ~= 0 
        xx = linspace(0,1,300)';
        XX = repmat(xx,1,10);
        JJ = repmat(j,300,1);
        HH = H_mat(XX,JJ);
        c = V(:,index(nplot));
        Z = HH*c;
        plot(xx,Z);
    end
    
    

