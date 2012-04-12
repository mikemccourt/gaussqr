clear all;

t = 16; % 10^(-t)
tvec = linspace(1,t,100);

sigma = 10.;
beta = 3;
L = 1;

M = 1; N = 10;
x = pickpoints(0,L,N);
z = x;

 for k = 2:length(tvec) 
     t1 = tvec(k); t0=tvec(k-1); 
     s1 = sobfunc_2(x,z,L,sigma,beta,M,t1);
     s0 = sobfunc_2(x,z,L,sigma,beta,M,t0); 
     fprintf('%e\t%e\t%e\n', t1, 10^(-t1), norm(s1-s0)/norm(s1)); 
 end
 
%plot(x,sobfunc_2(x,z,L,sigma,beta,M,t1))
