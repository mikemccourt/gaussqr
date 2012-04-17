%ex3_SQR
% This examples considers the required length of the series for the Sobolev
% function to be accurate for use in RBF-QR style error analysis.
% We can do this in 1D because we have a closed form for the function, and
% then use the known relationship for the change in beta.
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 4; % Relative RMS

% This is the closed form solution in 1D for beta = 1
sinhfunc = @(x,z,sigma,L) sinh(sigma*min(x,z)).*sinh(sigma*(L-max(x,z)))./(sigma*sinh(L*sigma));

L = 1;
sigma = 1;
beta = 1;

N = 100;
x = pickpoints(0,L,N);
z = .4*L*ones(size(x));

truesf = sinhfunc(x,z,sigma,L);

Mvec = 1000*2.^(0:7);
errvec = zeros(size(Mvec));
k = 1;
fprintf('Length\tRelative Err\tExec Time\tEval Ratio\n')
for M=Mvec
    tic
    compsf = sobfunc(x,z,L,sigma,beta,M);
    tcomp = toc;
    errvec(k) = errcompute(compsf,truesf);
    evalrat = (pi^2+sigma^2*L^2)/(pi^2*M^2+sigma^2*L^2);
    fprintf('%d\t%g\t%8g\t%g\n',M,errvec(k),tcomp,evalrat^beta)
    k = k + 1;
end

b = [ones(size(Mvec))',log(Mvec)']\log(errvec)';
g = ceil(exp(1/b(2)*(log(eps)-b(1))));
fprintf('For beta=%d, L=%g, sigma=%g\n',beta,L,sigma)
fprintf('According to these results, M=%d is required for machine precision\n',g)

beta = 2;

% Create a "true" answer by using a very large M.  I expect this series to
% converge twice as fast as the one above.  That means when the error above
% is 1e-5, this error should be 1e-10
tic
truesf = sobfunc(x,z,sigma,L,beta,100000);
tcomp = toc;
fprintf('\n%g seconds spent computing "true" function for beta = %d\n',tcomp,beta)

Mvec = 100*2.^(0:7);
errvec = zeros(size(Mvec));
k = 1;
fprintf('Length\tRelative Err\tExec Time\tEval Ratio\n')
for M=Mvec
    tic
    compsf = sobfunc(x,z,L,sigma,beta,M);
    tcomp = toc;
    errvec(k) = errcompute(compsf,truesf);
    evalrat = (pi^2+sigma^2*L^2)/(pi^2*M^2+sigma^2*L^2);
    fprintf('%d\t%g\t%8g\t%g\n',M,errvec(k),tcomp,evalrat^beta)
    k = k + 1;
end

b = [ones(size(Mvec))',log(Mvec)']\log(errvec)';
g = ceil(exp(1/b(2)*(log(eps)-b(1))));
fprintf('For beta=%d, L=%g, sigma=%g\n',beta,L,sigma)
fprintf('According to these results, M=%d is required for machine precision\n',g)

beta = 3;

% Create a "true" answer by using a very large M.  I expect this series to
% converge twice as fast as the one above.  That means when the error above
% is 1e-5, this error should be 1e-10
tic
truesf = sobfunc(x,z,sigma,L,beta,10000);
tcomp = toc;
fprintf('\n%g seconds spent computing "true" function for beta = %d\n',tcomp,beta)

Mvec = 10*2.^(0:7);
errvec = zeros(size(Mvec));
k = 1;
fprintf('Length\tRelative Err\tExec Time\tEval Ratio\n')
for M=Mvec
    tic
    compsf = sobfunc(x,z,L,sigma,beta,M);
    tcomp = toc;
    errvec(k) = errcompute(compsf,truesf);
    evalrat = (pi^2+sigma^2*L^2)/(pi^2*M^2+sigma^2*L^2);
    fprintf('%d\t%g\t%8g\t%g\n',M,errvec(k),tcomp,evalrat^beta)
    k = k + 1;
end

b = [ones(size(Mvec))',log(Mvec)']\log(errvec)';
g = ceil(exp(1/b(2)*(log(eps)-b(1))));
fprintf('For beta=%d, L=%g, sigma=%g\n',beta,L,sigma)
fprintf('According to these results, M=%d is required for machine precision\n',g)


N = 100;
L = 1;
x = pickpoints(0,L,N);
sigma = 10;
beta = 3;
Vk = zeros(N,9);
for k=1:9
    Vk(:,k) = sobfunc(x,L*k/10,L,sigma,beta);
end
plot(x,Vk)