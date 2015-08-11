% ex1
% Tests the HS-SVD to show that it can avoid the ill-conditioning
% associated with a small shape parameter.  The polynomial interpolant on
% the Chebyshev points is also computed as a demonstration of the ep->0
% limit correspondence.  Note that this matching limit is best demonstrated
% on the Chebyshev points because of better stability.
global GAUSSQR_PARAMETERS

% Define the closed Gaussian form of the kernel
rbf = @(e,r) exp(-(e*r).^2);

% Different base test options, although others are possible
%     1 - sin function on [-4,4]
%     2 - sinc function on [-1,1]
test_case = 1;

switch test_case
    case 1
        epvecd = logspace(-2,1,40);
        epvec = logspace(-2,0.3,100);
        Nvec = [10,20,30];
        NN = 100;
        spaceopt = 'cheb';
        fopt = 'sin';
        aa = -4;bb = 4;
        gqr_alpha = .65;
        GAUSSQR_PARAMETERS.NORM_TYPE = 2;
    case 2
        epvecd = logspace(-2,1.5,41);
        epvec = logspace(-2,1.5,81);
        Nvec = [3,5,9,17];
        NN = 80;
        spaceopt = 'even';
        fopt = 'sinc';
        aa = -1;bb = 1;
        gqr_alpha = 1;
        GAUSSQR_PARAMETERS.NORM_TYPE = inf;
    otherwise
        error('No such test exists')
end

[yf,fstr] = pickfunc(fopt,1);
xx = pickpoints(aa,bb,NN);
yy = yf(xx);
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvecd));
errpoly = zeros(length(Nvec),1);

% Check the error of the various interpolants
n = 1;
for N=Nvec
    [x,spacestr] = pickpoints(aa,bb,N,spaceopt);
    y = yf(x);
    
        % The HS-SVD Gaussian interpolant
    k = 1;
    for ep=epvec
        GQR = gqr_solve(x,y,ep,gqr_alpha);
        yp = gqr_eval(GQR,xx);
        errvec(n,k) = errcompute(yp,yy);
        k = k + 1;
    end
    
    % The direct Gaussian interpolant
    k = 1;
    for ep=epvecd
        K = rbf(ep,DistanceMatrix(x,x));
        warning('off','MATLAB:nearlySingularMatrix')
        beta = K\y;
        warning('on','MATLAB:nearlySingularMatrix')
        yp = rbf(ep,DistanceMatrix(xx,x))*beta;
        errvecd(n,k) = errcompute(yp,yy);
        k = k + 1;
    end
    
    % Create the polynomial interpolants
    warning('off','MATLAB:polyfit:RepeatedPoints')
    [p,S,mu] = polyfit(x,y,N-1);
    warning('on','MATLAB:polyfit:RepeatedPoints')
    yp = polyval(p,xx,S,mu);
    errpoly(n) = errcompute(yp,yy);
    
    n = n+1;
end

switch test_case
    case 1
        loglog(epvecd,errvecd(1,:),'-bx')
        hold on
        loglog(epvecd,errvecd(2,:),'-g+')
        loglog(epvecd,errvecd(3,:),'-r^')
        loglog(epvec,errvec(1,:),'b','LineWidth',3)
        loglog(epvec,errvec(2,:),'g','LineWidth',3)
        loglog(epvec,errvec(3,:),'r','LineWidth',3)
        loglog(epvecd,errpoly(1)*ones(size(epvecd)),'--b')
        loglog(epvecd,errpoly(2)*ones(size(epvecd)),'--g')
        loglog(epvecd,errpoly(3)*ones(size(epvecd)),'--r')
        hold off
        xlabel('\epsilon')
        ylabel('average error')
        ylim([10^-15 10])
        ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
        title(strcat(fstr,ptsstr,spacestr))
        legend('N=10 (Direct)','N=20 (Direct)','N=30 (Direct)','N=10 (QR)','N=20 (QR)','N=30 (QR)','Location','SouthEast')
    case 2
        loglog(epvec,errvec(1,:),'r','LineWidth',3)
        hold on
        loglog(epvec,errvec(2,:),'g','LineWidth',3)
        loglog(epvec,errvec(3,:),'b','LineWidth',3)
        loglog(epvec,errvec(4,:),'c','LineWidth',3)
        loglog(epvecd,errvecd(1,:),'-rx')
        loglog(epvecd,errvecd(2,:),'-g+')
        loglog(epvecd,errvecd(3,:),'-b^')
        loglog(epvecd,errvecd(4,:),'-c^')
        loglog(epvecd,errpoly(1)*ones(size(epvecd)),'--r')
        loglog(epvecd,errpoly(2)*ones(size(epvecd)),'--g')
        loglog(epvecd,errpoly(3)*ones(size(epvecd)),'--b')
        loglog(epvecd,errpoly(4)*ones(size(epvecd)),'--c')
        hold off
        xlabel('\epsilon')
        ylabel('Error')
        ylim([10^-15 10])
        xlim([min([epvec,epvecd]),max([epvec,epvecd])])
        ptsstr=strcat(', x\in[',num2str(aa),',',num2str(bb),'],');
        title(strcat(fstr,ptsstr,spacestr))
        legend('N=3','N=5','N=9','N=17','Location','SouthEast')
    otherwise
        error('You already know you messed up if you are here')
end