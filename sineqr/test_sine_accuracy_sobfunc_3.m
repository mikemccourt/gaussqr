clear all;

t = 20;
tvec = linspace(1,t,20);

sigma = .1;
beta = 3;
L = 5;

N = 10;
x = pickpoints(0,L,N);
z = x;

fprintf('\t\tM0\t\t    M1\t\t\t   t0\t\t\t   t1\t     norm(s0)\t     norm(s1)   norm(s1-s0)/norm(s0)\n');
 for k = 2:length(tvec) 
     t1 = tvec(k); t0=tvec(k-1); 
     [MM1,s1] = sobfunc_3(x,z,L,sigma,beta,t1);
     [MM0,s0] = sobfunc_3(x,z,L,sigma,beta,t0); 
     fprintf('%10d\t%10d\t%e\t%e\t%e\t%e\t%e\n', MM0, MM1, t0, t1, norm(s0), norm(s1), norm(s1-s0)/norm(s0)); 
 end
 
%plot(x,sobfunc_2(x,z,L,sigma,beta,M,t1))
