% This was the history from Tuesday (6/3) work session
% Turn this into a real program
% The program should take in N, the number of points in the domain
% It will return the eigenvalues and eigenfunction coefficients
%   associated with the collocation approximation
% Many of the things I ahve done may not make sense at first blush
% Ask me if you have questions

x = linspace(0,1,10)';
DM = abs(repmat(x,1,10)-repmat(x',10,1))
ep = 1;
K = exp(-ep^2*DM.^2)
eig(K)
plot(sort(eig(K)))
semilogy(sort(eig(K),'descend'))
ep = .1;
K = exp(-ep^2*DM.^2)
semilogy(sort(eig(K),'descend'))
x
X = repmat(x,1,10)
j = 1:10
J = repmat(j,10,1)
Int_Kh = @(x,j) -1./(j.^2+j).*x.*(x.^j-1)
A = Int_Kh(X,J)
H_mat = @(z,j) z.^(j-1);
H = H_mat(X,J)
HinvA = H\A;
[V,D] = eig(HinvA)
A
HinvA
format short
HinvA
HinvA.*(abs(HinvA)>1e-10)
x = cos((0:11)/12*pi)
x = cos((0:11)/11*pi)
x = cos((0:11)/11*pi)/2+.5
x = cos((0:11)/11*pi)/2+.5; x = x(2:end-1)'
X = repmat(x,1,10)
Int_Kh = @(x,j) -1./(j.^2+j).*x.*(x.^j-1)
A = Int_Kh(X,J)
H_mat = @(z,j) z.^(j-1);
H = H_mat(X,J)
HinvA = H\A;
HinvA
eig(H)
eig(A)
log(abs(eig(A)))
log10(abs(eig(A)))
log10(abs(eig(HinvA)))
eig(HinvA)
[V,D] = eig(A)
[V,D] = eig(HinvA)
H_mat
c = V(:,1);
xx = linspace(0,1,300)';
XX = repmat(xx,1,10);
JJ = repmat(j,300,1);
HH = H_mat(XX,JJ);
Z = HH*c;
plot(xx,Z)
c = V(:,2);
c = V(:,2);Z = HH*c;plot(xx,Z)
c = V(:,3);Z = HH*c;plot(xx,Z)
[ix,d] = sort(diag(D))
[ix,d] = sort(diag(D),'ascend')
[ix,d] = sort(diag(D),'descend')
n_eig = 1;c = V(:,d(n_eig));Z = HH*c;plot(xx,Z)
n_eig = 2;c = V(:,d(n_eig));Z = HH*c;plot(xx,Z)
n_eig = 3;c = V(:,d(n_eig));Z = HH*c;plot(xx,Z)
n_eig = 4;c = V(:,d(n_eig));Z = HH*c;plot(xx,Z)
n_eig = 5;c = V(:,d(n_eig));Z = HH*c;plot(xx,Z)
ZZ = HH*V(:,d);
plot(xx,ZZ(:,1:3))
