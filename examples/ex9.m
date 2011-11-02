% This example will help determine what the cost of the RBF-QR computation
% is from a time perspective
%
% In order to get useful times, this example may take a long time to run

clear all
rbfsetup

Nvec = [20 40 80 160 320 640 1280 2560 5120 10240];
epvec = [.005 .02 .08 .32 1.28];
times1d = zeros(length(Nvec),length(epvec));
teval1d = zeros(length(Nvec),length(epvec));
Msize1d = zeros(length(Nvec),length(epvec));

yf = @(x) 1+tanh(x); % Irrelevant
aa = -4;bb = 4;

i = 1;
for N=Nvec
    x = pickpoints(aa,bb,N);
    y = yf(x);
    NN = 10*N;
    xx = pickpoints(aa,bb,NN);
    yy = yf(xx);
    j = 1;
    for ep=epvec
        clear rbfqrOBJ yp
        rbfsetup
        alpha = rbfalphasearch(ep,aa,bb);
        tic
        rbfqrOBJ = rbfqr_solve(x,y,ep,alpha);
        times1d(i,j) = toc;
        tic
        yp = rbfqr_eval(rbfqrOBJ,xx);
        teval1d(i,j) = toc;
        Msize1d(i,j) = size(rbfqrOBJ.Marr,2);
        j = j+1;
        fprintf(' %d ',j)
    end
    i = i+1;
    fprintf(' %d \n',i)
end

Nvec = [3,5,7,9,11,13;3,5,7,9,11,13];
epvec = [.005 .01 .02 .04 .08 .16 .32];
times2d = zeros(size(Nvec,2),length(epvec));
teval2d = zeros(size(Nvec,2),length(epvec));
Msize2d = zeros(size(Nvec,2),length(epvec));

yf = @(x) 1+tanh(x(:,1)+x(:,2));
aa = [-4,-4];bb = [4,4];

i = 1;
for N=Nvec
    x = pick2Dpoints(aa,bb,N);
    y = yf(x);
    NN = 3*N;
    xx = pick2Dpoints(aa,bb,NN);
    yy = yf(xx);
    j = 1;
    for ep=epvec
        clear rbfqrOBJ yp
        rbfsetup
        alpha = rbfalphasearch(ep,aa,bb);
        tic
        rbfqrOBJ = rbfqr_solve(x,y,ep,alpha);
        times2d(i,j) = toc;
        tic
        yp = rbfqr_eval(rbfqrOBJ,xx);
        teval2d(i,j) = toc;
        Msize2d(i,j) = size(rbfqrOBJ.Marr,2);
        j = j+1;
        fprintf(' %d ',j)
    end
    i = i+1;
    fprintf(' %d \n',i)
end

save timeruns.mat times1d Msize1d teval1d times2d teval2d Msize2d
