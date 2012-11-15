% Test studying the leave half out cross-validation
% Uses compact Matern kernels
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

close all

%f = @(x) 30*x.^2.*(L-x).^2.*sin(2*pi*x/L).^4;
%f = @(x) 30*x.^2.*(L-x).^2.*sin(10*pi*x/L).^4;
%f = @(x) 1./(1+(x/L).^2)-(1-.5*(x/L));
%f = @(x) sinh(3/L*x)./(1+cosh(3/L*x));    % Interesting 1/3 out predicts clear minimum which doesn't seem to exist
%f = @(x) cos(x)+exp(-(x-1).^2)-exp(-(x+1).^2);
%f = @(x) (x.^2.*(1-x).^2).*(cos(10*x)+exp(-(x-1).^2)-exp(-(x+1).^2));
%testfuncN = 14; f = @(x) 10^(testfuncN+1).*(max(0,x-(1/4))).^testfuncN.*(max(0,(3/4)-x)).^testfuncN;
%f = @(x) exp(x) - (1-x) - x*exp(1);
f = @(x) exp(x) - 1 + (1/6).*(x.^3 - 3.*x.^2 + 8.*x - exp(1).*(x.^3+5.*x));
%f = @(x) min(abs(x-pi/6)-.25,0);
%f = @(x) 1./(1+(2*x-1).^2);
              
L = 1;

N = 50;
x = pickpoints(0,L,N+2);x = x(2:end-1);
%x = pickpoints(0,L,N+2,'cheb');x = x(2:end-1);
y = f(x);

NN = 200;
xx = pickpoints(0,L,NN);
yy = f(xx);

h1 = 1:2:N;
h2 = setdiff(1:N,h1);

t1 = 1:3:N;
t2 = 2:3:N;
t3 = setdiff(1:N,[t1,t2]);

DM = DistanceMatrix(x,x);

epvec = unique([logspace(-1,0,5),logspace(0,1,10),logspace(1,3,30)]);
betavec = [1:8];
halfvec = [];
thirdvec = [];
loocvvec = [];
errvec = [];
l = 1;
for beta=betavec
    k = 1;
    fprintf('l=%d \t beta=%g \n',l,beta)
    for ep=epvec
        %%%%%%%%%%%%%%%%%%%%%%%%
        % First the leave half out
        x_train = x(h1);
        y_train = f(x_train);
        x_valid = x(h2);
        y_valid = f(x_valid);
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        halfvec(k,l) = errcompute(yp,y_valid);
        
        x_train = x(h2);
        y_train = f(x_train);
        x_valid = x(h1);
        y_valid = f(x_valid);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        halfvec(k,l) = halfvec(k,l) + errcompute(yp,y_valid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Then the leave 1/3 out
        x_train = x([t1,t2]);
        y_train = f(x_train);
        x_valid = x(t3);
        y_valid = f(x_valid);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        thirdvec(k,l) = errcompute(yp,y_valid);
        
        x_train = x([t1,t3]);
        y_train = f(x_train);
        x_valid = x(t2);
        y_valid = f(x_valid);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        thirdvec(k,l) = thirdvec(k,l) + errcompute(yp,y_valid);
        
        x_train = x([t2,t3]);
        y_train = f(x_train);
        x_valid = x(t1);
        y_valid = f(x_valid);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        thirdvec(k,l) = thirdvec(k,l) + errcompute(yp,y_valid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now the actual solution
        MQR = mqr_solve(x,y,L,ep,beta);
        yp = mqr_eval(MQR,xx);
        errvec(k,l) = errcompute(yp,yy);
        
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
        loocvvec(k,l) = norm(EF,1);
        
        fprintf('k=%d \t ep=%g \n',k,ep)
        k = k + 1;
    end
    l = l + 1;
end

%%%%%%%%%%%%%%%%%
% Finally, plot all epsilon-beta error surfaces
[X,Y] = meshgrid(epvec,betavec);
surf(X,Y,errvec'), title('True solution')
set(gca,'XScale','log')
set(gca,'ZScale','log')
xlabel('\epsilon')
ylabel('\beta')
zlabel('error')
[i,j]=find(errvec==min(min(errvec)));
fprintf('True solution: optimal epsilon=%f, optimal beta=%d\n',epvec(i),betavec(j))
figure, surf(X,Y,halfvec'), title('Leave 1/2 out')
set(gca,'XScale','log')
set(gca,'ZScale','log')
xlabel('\epsilon')
ylabel('\beta')
zlabel('error')
figure, surf(X,Y,thirdvec'), title('Leave 1/3 out')
set(gca,'XScale','log')
set(gca,'ZScale','log')
xlabel('\epsilon')
ylabel('\beta')
zlabel('error')
[i,j]=find(thirdvec==min(min(thirdvec)));
fprintf('Third out CV: optimal epsilon=%f, optimal beta=%d\n',epvec(i),betavec(j))
figure, surf(X,Y,loocvvec'), title('LOOCV')
set(gca,'XScale','log')
set(gca,'ZScale','log')
xlabel('\epsilon')
ylabel('\beta')
zlabel('error')
