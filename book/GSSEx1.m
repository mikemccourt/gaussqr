f = @(x) sin(2*pi*x) + x;
N = 10;
x = pickpoints(0,1,N);
x = pickpoints(0,1,N+2);x = x(2:end-1);
y = f(x);
NN = 150;
xx = pickpoints(0,1,NN);
yy = f(xx);

MM = @(x,z) {repmat(x,1,size(z,1)),repmat(z',size(x,1),1)};
tmp = MM(x,x);
X = tmp{1};
Z = tmp{2};
tmp = MM(xx,x);
XX = tmp{1};
ZZ = tmp{2};

KK{1} = @(x,z) min(x,z) - .129*x.*z   - .3261*(x + z)  + .2656;
KK{2} = @(x,z) min(x,z) - 0.3598*x.*z - 0.2055*(x + z) + 0.1376;
KK{3} = @(x,z) min(x,z) + 3.3171*x.*z - 2.1692*(x + z) + 1.1006;
KK{4} = @(x,z) min(x,z) + .998*x.*z   - .9988*(x + z)  + .9999;

h = figure;
for k=1:4
    Kf = KK{k};
    K_int = Kf(X,Z);
    K_eval = Kf(XX,ZZ);
    yp = K_eval*(K_int\y);
    
    subplot(2,2,k)
    plot(xx,yp,'linewidth',3)
    title(sprintf('Error = %g',errcompute(yp,yy)))
end