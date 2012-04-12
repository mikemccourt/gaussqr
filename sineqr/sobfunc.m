% Sobolev function
function s = sobfunc(x,z,L,sigma)

sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L)

beta = 2;
lamfunc = @(n,L,sigma,beta) ((pi*n/L).^2+sigma^2).^(-beta)

% ci = zeros(size(x));
 t = 6; % 20 does not work; take t <= 11 
 M = 10;
 
 MM = 1/pi*sqrt(10^(t/beta)*((pi*M)^2+(sigma*L)^2)-(sigma*L)^2)
 %ceil(MM), floor(MM)
 
 Xmat = sinfunc(1:MM,L,x);
 Zmat = sinfunc(1:MM,L,z);
 Lmat = diag(lamfunc(1:MM,L,sigma,beta)); size(Lmat) % take t <= 6
 Smat = (Xmat.*Zmat)*Lmat;
 
 s = sum(Smat,2);

