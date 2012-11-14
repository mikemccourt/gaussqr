% Test studying the leave half out cross-validation
% Uses compact Matern kernels
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

f = @(x) sin(2*pi*x);

L = 1;

N = 24;
x = pickpoints(0,L,N+2);x = x(2:end-1);
y = f(x);

NN = 100;
xx = pickpoints(0,L,NN);
yy = f(xx);

beta = 5;

h1 = 1:2:N;
h2 = setdiff(1:N,h1);

t1 = 1:3:N;
t2 = 2:3:N;
t3 = setdiff(1:N,[t1,t2]);

DM = DistanceMatrix(x,x);

epvec = [logspace(-2,-1,5),logspace(-1,.3,33),logspace(.3,1,6)];
halfvec = [];
thirdvec = [];
loocvvec = [];
errvec = [];
k = 1;
for ep=epvec
    %%%%%%%%%%%%%%%%%%%%%%%%
    % First the leave half out
    x_train = x(h1);
    y_train = f(x_train);
    x_valid = x(h2);
    y_valid = f(x_valid);
    MQR = mqr_solve(x_train,y_train,L,ep,beta);
    yp = mqr_eval(MQR,x_valid);
    halfvec(k) = errcompute(yp,y_valid);
    
    x_train = x(h2);
    y_train = f(x_train);
    x_valid = x(h1);
    y_valid = f(x_valid);
    
    MQR = mqr_solve(x_train,y_train,L,ep,beta);
    yp = mqr_eval(MQR,x_valid);
    halfvec(k) = halfvec(k) + errcompute(yp,y_valid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Then the leave 1/3 out
    x_train = x([t1,t2]);
    y_train = f(x_train);
    x_valid = x(t3);
    y_valid = f(x_valid);
    
    MQR = mqr_solve(x_train,y_train,L,ep,beta);
    yp = mqr_eval(MQR,x_valid);
    thirdvec(k) = errcompute(yp,y_valid);
    
    x_train = x([t1,t3]);
    y_train = f(x_train);
    x_valid = x(t2);
    y_valid = f(x_valid);
    
    MQR = mqr_solve(x_train,y_train,L,ep,beta);
    yp = mqr_eval(MQR,x_valid);
    thirdvec(k) = thirdvec(k) + errcompute(yp,y_valid);
    
    x_train = x([t2,t3]);
    y_train = f(x_train);
    x_valid = x(t1);
    y_valid = f(x_valid);
    
    MQR = mqr_solve(x_train,y_train,L,ep,beta);
    yp = mqr_eval(MQR,x_valid);
    thirdvec(k) = thirdvec(k) + errcompute(yp,y_valid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now the actual solution
    MQR = mqr_solve(x,y,L,ep,beta);
    yp = mqr_eval(MQR,xx);
    errvec(k) = errcompute(yp,yy);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now an attempt at LOOCV
    [L,ep,beta,M] = mqr_solveprep(x,L,ep,beta);
    Phi = mqr_phi(1:M,x,L);
    Phi1 = Phi(:,1:N);
    Psi = Phi*[eye(N);MQR.Rbar];
    invPsi = pinv(Psi);
    invPhi1 = pinv(Phi1');
    Lambda1 = diag(((pi*(1:N)/L).^2+ep^2).^(-beta));
    invLambda1 = pinv(Lambda1);
    invA = invPhi1*invLambda1*invPsi;
    EF = (invA*y)./diag(invA);
    loocvvec(k) = norm(EF,1);
    
    
    fprintf('k=%d \t ep=%g \n',k,ep)
    k = k + 1;
end

loglog(epvec,loocvvec,'--k','linewidth',3), hold on
loglog(epvec,thirdvec,'r','linewidth',3)
loglog(epvec,halfvec,'g','linewidth',3)
loglog(epvec,errvec,'linewidth',3), hold off
xlabel('\epsilon')
ylabel('error')
legend('Leave one out','Leave 1/2 out','Leave 1/3 out','True solution')
