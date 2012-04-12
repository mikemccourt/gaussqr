% Script for: test_sobolev_2.m

N = 10; % number of points
L = 5;

sigma = 10; % [1e-2,1e-1,1,10,100]
beta = 3;  % [1,2,3]

M = 1;
t = 16;

[x,spacestr] = pickpoints(0,L,N,'cheb');
z = x;

s = sobfunc_2(x,z,L,sigma,beta,M,t)
