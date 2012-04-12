N = 10;
L = 1;
sigma = 1;

[x,spacestr] = pickpoints(0,L,N,'cheb',sigma);
z = x;

s = sobfunc(x,z,L,sigma)
