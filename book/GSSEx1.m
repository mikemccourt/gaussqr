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

pause
%%%%%% Something below regarding reproducing property - more to do
f = @(x) x.^2;
fp = @(x) 2*x;
Ap = @(p) [0,p(2)-1,p(2)+p(1);0,p(3)+1,p(2)+p(3);p(2)-1,p(2)+p(1),0;p(3)+1,p(2)+p(3),0];
rf = @(p) Ap(p)\[1;0;-p(2);-p(3)];
p = [7,1,2];r=rf(p),norm(Ap(p)*r-[1;0;-p(2);-p(3)])
Kr = @(x,z,r) min(x,z) + r(1)*x.*z + r(2)*(x+z) + r(3);
x = pickpoints(0,1,500);p = [7,1,2];r=rf(p);figure,plot(x,[Kr(x,0,r),Kr(x,.2,r),Kr(x,.4,r),Kr(x,.6,r),Kr(x,.8,r),Kr(x,1,r)],'linewidth',2)
Krp = @(x,z,r) (x<z) + r(1)*z + r(2);
ipP = @(fp,gp) quadl(@(x)fp(x).*gp(x),0,1);
ipP_known = @(f,z,r) f(z) - f(0) + (f(1)-f(0))*(r(1)*z+r(2));
ipB = @(f,g,p) p(1)*f(0)*g(0) + p(2)*(f(0)*g(1) + f(1)*g(0)) + p(3)*f(1)*g(1);
ipB_known = @(f,z,r,p) p(1)*f(0)*(r(2)*z+r(3)) + p(2)*(f(0)*(z*(1+r(1)+r(2))+r(2)+r(3)) + f(1)*(r(2)*z+r(3))) + p(3)*f(1)*(z*(1+r(1)+r(2))+r(2)+r(3))
z = .5;
p = [7,1,2];
r = rf(p);
ipP(fp,@(x)Krp(x,z,r))
ipP_known(f,z,r)
ipB(f,@(x)Kr(x,z,r),p)
ipB_known(f,z,r,p)
zz = pickpoints(0,1,200)';
fz = [];for z=zz fz = [fz ipP_known(f,z,r) + ipB_known(f,z,rf(p),p)];end
figure,plot(zz,fz)

fr = @(x) [-1+x(3)+x(4)+2*x(1)*x(2)*(x(4)-x(3));x(1)*x(2)*(x(3)-x(4))-x(3)+x(1)^2*(x(3)-x(4));x(3)+x(1)^2*(x(4)-x(3));x(1)^2+x(2)^2];
v = fsolve(@(x)fr(x)-[r;1],[1;1;1;1])
Kr_test = @(x,z,v) min(x,z) + x.*z*(-1+v(3)+v(4)+2*v(1)*v(2)*(v(4)-v(3))) + (x+z)*(v(1)*v(2)*(v(3)-v(4))-v(3)+v(1)^2*(v(3)-v(4))) + v(3)+v(1)^2*(v(4)-v(3));
x = pickpoints(0,1,500);figure,plot(x,[Kr_test(x,0,v),Kr_test(x,.2,v),Kr_test(x,.4,v),Kr_test(x,.6,v),Kr_test(x,.8,v),Kr_test(x,1,v)],'linewidth',2)
fp = @(x) [(1-x(1)^2)/x(3)+x(1)^2/x(4)-1;x(1)*x(2)*(1/x(3)-1/x(4))+1;x(1)^2/x(3)+(1-x(1)^2)/x(4)-1;x(1)^2+x(2)^2]
v = fsolve(@(x)fp(x)-[p;1],[1;1;1;1])