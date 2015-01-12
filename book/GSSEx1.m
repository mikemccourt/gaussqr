% A demonstration of the computation of GSS kernels for a given inner
% product.  This requires either the solution of a nonlinear or linear
% system and there may be multiple solutions that produce the same kernel.
%
% There may also be some discussion of the reproducing property but at the
% moment, that portion of the code below has been commented out.

% The values in this vector allow are the coefficients of the inner
% product which are chosen to define a specific kernel
% I make this a column vector for later computation
p = [101;100;102];

% This test confirms that the values in this inner product actually produce
% a valid semi-norm
if not((p(1)+p(3)>=0) && (p(2)^2-p(1)*p(3)<=0))
    error('The proposed inner product is invalid:\np(1)+p(3)=%g, p(2)^2-p(1)*p(3)=%g',p(1)+p(3),p(2)^2-p(1)*p(3))
end

% Form the linear system which defines the kernel coefficients
% Note that this can be done linearly only because this particular problem
% happens to be simple, and in general will require a nonlinear solve
Ap = @(p) [0,p(2)-1,p(2)+p(1);...
           0,p(3)+1,p(2)+p(3);...
           p(2)-1,p(2)+p(1),0;...
           p(3)+1,p(2)+p(3),0];

% Using the linear system, solve for the kernel coefficients
% This is a function that can be used to solve for any given p
rf = @(p) Ap(p)\[1;0;-p(2);-p(3)];

% Evaluate the kernel using known coefficients
Kr = @(x,z,r) min(x,z) + r(1)*x.*z + r(2)*(x+z) + r(3);

% Plot the kernel with a variety of centers
x = pickpoints(0,1,500);
r = rf(p);
plot(x,[Kr(x,0,r),Kr(x,.2,r),Kr(x,.4,r),Kr(x,.6,r),Kr(x,.8,r),Kr(x,1,r)],'linewidth',3)
fprintf('Inner product coefficients p=[%g %g %g]\n',p)
fprintf('Kernel coefficients r=[%g %g %g]\n',r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below, we test if the same solution can be achieved with the nonlinear
% solution strategy as with the linear solution strategy above.
% Define the nonlinear components of the problem, noting that there are 4
% here: ca = [c_1 c_2 a_1 a_2]
ca_comp = @(ca) [(1-ca(1)^2)/ca(3) + ca(1)^2/ca(4) - 1;...
                 ca(1)*ca(2)*(1/ca(3)-1/ca(4)) + 1;...
                 ca(1)^2/ca(3) + (1-ca(1)^2)/ca(4) - 1;...
                 ca(1)^2+ca(2)^2];

% Define the full nonlinear function, including the nonhomogeneous RHS
% This is where we use the fact that p is a column vector
ca_solve = @(ca,p) ca_comp(ca) - [p;1];

% Solve the nonlinear system, using a dumb initial guess
init_guess = [2;2;2;2];
fsolve_opt = optimoptions('fsolve','Display','off','MaxFunEvals',10000);
ca = fsolve(@(ca)ca_solve(ca,p),init_guess,fsolve_opt);
fprintf('Kernel c & a (nonlinear) coefficients c=[%g %g] a=[%g %g]\n',ca)

% Define the kernel in terms of the ca values instead of the p values
Kca = @(x,z,ca) min(x,z) + ...
                x.*z*(-1+ca(3)+ca(4)+2*ca(1)*ca(2)*(ca(4)-ca(3))) + ...
                (x+z)*(ca(1)*ca(2)*(ca(3)-ca(4))-ca(3)+ca(1)^2*(ca(3)-ca(4))) + ...
                ca(3)+ca(1)^2*(ca(4)-ca(3));

% Test if the kernel evaluated with this strategy and with the linear
% strategy are the same kernels
norm_diff = errcompute(Kca(x,.23456,ca),Kr(x,.23456,r));
fprintf('Norm of difference between linear and nonlinear solves: %g\n',norm_diff)
            
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Something below regarding reproducing property - NOT WORKING
Ap = @(p) [0,p(2)-1,p(2)+p(1);0,p(3)+1,p(2)+p(3);p(2)-1,p(2)+p(1),0;p(3)+1,p(2)+p(3),0];
rf = @(p) Ap(p)\[1;0;-p(2);-p(3)];
p = [7,1,2];r=rf(p),norm(Ap(p)*r-[1;0;-p(2);-p(3)])
f = @(x) x.^2;
fp = @(x) 2*x;
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

pause

fr = @(x) [-1+x(3)+x(4)+2*x(1)*x(2)*(x(4)-x(3));x(1)*x(2)*(x(3)-x(4))-x(3)+x(1)^2*(x(3)-x(4));x(3)+x(1)^2*(x(4)-x(3));x(1)^2+x(2)^2];
v = fsolve(@(x)fr(x)-[r;1],[1;1;1;1])
Kr_test = @(x,z,v) min(x,z) + x.*z*(-1+v(3)+v(4)+2*v(1)*v(2)*(v(4)-v(3))) + (x+z)*(v(1)*v(2)*(v(3)-v(4))-v(3)+v(1)^2*(v(3)-v(4))) + v(3)+v(1)^2*(v(4)-v(3));
x = pickpoints(0,1,500);figure,plot(x,[Kr_test(x,0,v),Kr_test(x,.2,v),Kr_test(x,.4,v),Kr_test(x,.6,v),Kr_test(x,.8,v),Kr_test(x,1,v)],'linewidth',2)
fp = @(x) [(1-x(1)^2)/x(3)+x(1)^2/x(4)-1;x(1)*x(2)*(1/x(3)-1/x(4))+1;x(1)^2/x(3)+(1-x(1)^2)/x(4)-1;x(1)^2+x(2)^2]
v = fsolve(@(x)fp(x)-[p;1],[1;1;1;1])