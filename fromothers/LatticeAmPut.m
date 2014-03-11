function [price, lattice] = LatticeAmPut(S0,X,r,T,sigma,N)
% Code from Fig. 1.24 in Brandimarte (2003) to price American put options
% via a binomial lattice model.
% As an example, use S0=50, X=50, r=0.1, T=5/12, sigma=0.4, N=5, i.e.,
% LatticeAmPut(50,50,0.1,5/12,0.4,5)
deltaT = T/N;
u = exp(sigma*sqrt(deltaT));
d = 1/u;
p = (exp(r*deltaT) - d)/(u-d);
lattice = zeros(N+1,N+1);
for j=0:N
    lattice(N+1,j+1) = max(0, X-S0*(u^j)*(d^(N-j)));
end
for i=N-1:-1:0
    for j=0:i
        lattice(i+1,j+1) = max(X-S0*u^j*d^(i-j), ...
            exp(-r*deltaT)*(p*lattice(i+2,j+2) + (1-p)*lattice(i+2,j+1)));
    end
end
price = lattice(1,1);