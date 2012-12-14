% Test studying the leave half out cross-validation
% Uses compact Matern kernels and titanium test data from Curve Fitting
% Toolbox (see "doc titanium" for original reference)
rbfsetup
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;

close all

[xx,yy] = titanium;
NN = length(xx);
xx = xx'; yy = yy';
xx = (xx-min(xx))/(max(xx)-min(xx));
yy = yy - min(yy);
pick = [2 5 11 21 27 29 31 33 35 40 45 48];
N = length(pick);
x = xx(pick);
y = yy(pick);
L = 1;

h1 = 1:2:N;
h2 = setdiff(1:N,h1);

t1 = 1:3:N;
t2 = 2:3:N;
t3 = setdiff(1:N,[t1,t2]);

DM = DistanceMatrix(x,x);

epvec = unique([logspace(-1,0,5),logspace(0,1,5),logspace(1,2,40)]);
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
        y_train = y(h1);
        x_valid = x(h2);
        y_valid = y(h2);
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        halfvec(k,l) = errcompute(yp,y_valid);
        
        x_train = x(h2);
        y_train = y(h2);
        x_valid = x(h1);
        y_valid = y(h1);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        halfvec(k,l) = halfvec(k,l) + errcompute(yp,y_valid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Then the leave 1/3 out
        x_train = x([t1,t2]);
        y_train = y([t1,t2]);
        x_valid = x(t3);
        y_valid = y(t3);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        thirdvec(k,l) = errcompute(yp,y_valid);
        
        x_train = x([t1,t3]);
        y_train = y([t1,t3]);
        x_valid = x(t2);
        y_valid = y(t2);
        
        MQR = mqr_solve(x_train,y_train,L,ep,beta);
        yp = mqr_eval(MQR,x_valid);
        thirdvec(k,l) = thirdvec(k,l) + errcompute(yp,y_valid);
        
        x_train = x([t2,t3]);
        y_train = y([t2,t3]);
        x_valid = x(t1);
        y_valid = y(t1);
        
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
[i,j]=find(halfvec==min(min(halfvec)));
fprintf('Half out CV: optimal epsilon=%f, optimal beta=%d\n',epvec(i),betavec(j))
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
[i,j]=find(loocvvec==min(min(loocvvec)));
fprintf('LOOCV: optimal epsilon=%f, optimal beta=%d\n',epvec(i),betavec(j))
