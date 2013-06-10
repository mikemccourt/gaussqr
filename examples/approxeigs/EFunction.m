function [V,D]= EFunction(N)
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
    D = sort(diag(d),'descend');
    

