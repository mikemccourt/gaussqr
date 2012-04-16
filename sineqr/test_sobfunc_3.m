% Script for: test_sobolev_3.m

N = 10; % number of points
L = 1;

sigma = 1.0; % [1e-2,1e-1,1,10,100]
beta = 3;  % [1,2,3]

t = 5;

[x,spacestr] = pickpoints(0,L,N,'cheb');
z = x;

[MM,s] = sobfunc_3(x,z,L,sigma,beta,t)
