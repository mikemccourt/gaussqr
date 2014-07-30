% GSSEx2.m
f = @(x) sin(2*pi*x) + x;
Nvec = floor(logspace(1,3,25));
NN = 150;
xx = pickpoints(0,1,NN);
yy = f(xx);

MM = @(x,z) {repmat(x,1,size(z,1)),repmat(z',size(x,1),1)};
K = @(x,z) min(x,z) + 1;

errvec = zeros(size(Nvec));
k = 1;
for N=Nvec
    x = pickpoints(0,1,N);
    y = f(x);
    tmp = MM(x,x);
    X = tmp{1};
    Z = tmp{2};
    tmp = MM(xx,x);
    XX = tmp{1};
    ZZ = tmp{2};
    
    K_int = K(X,Z);
    K_eval = K(XX,ZZ);
    yp = K_eval*(K_int\y);
    errvec(k) = errcompute(yp,yy);
    k = k + 1;
end

h = figure;
loglog(Nvec,errvec,'linewidth',3)
xlabel('\epsilon')
ylabel('error')