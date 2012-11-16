% ex6a_mqr
% Tests RBF-QR to show that it can avoid the ill-conditioning
% associated with a small shape parameter
rbfsetup
global GAUSSQR_PARAMETERS
useSplines = or(GAUSSQR_PARAMETERS.SPLINE_TOOLBOX_AVAILABLE,GAUSSQR_PARAMETERS.CURVEFIT_TOOLBOX_AVAILABLE);
GAUSSQR_PARAMETERS.ERROR_STYLE = 2; % Absolute error
GAUSSQR_PARAMETERS.NORM_TYPE = inf; % Sup norm
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = 500; % Need more eigenfunctions for N=40, 
                                          % but then the direct method for beta=3 is VERY slow

% The range of N values to consider
Nvec = [12,24,48];
% The spacing choice for the points
spaceopt = 'even';
% The order (smoothness) of the kernel
beta = 2;
% The  range of kernel shape parameters to consider
epvec =  logspace(0,2,31);
epvec = [logspace(0,1,10),logspace(1,2,35)]; % Focus on the right region
% The length of the domain
L = 1;
% The embedding width for nonhomogeneous functions (must be <.5)
embed_cushion = .1;
% The number of evenly spaced points at which to sample error
NN = 397;

% This is the function we are interested in considering
% Depending on which function consider, it will choose embedding
fopt = 9;
fpar = [3,.0567];
switch fopt
    case 1
        yf = @(x) sin(2*pi*x/L) + 1;
        fstr = 'u(x) = sin(2\pi{x}/L)+1';
        embed = embed_cushion;
    case 2
        yf = @(x) sin(2*pi*x/L);
        fstr = 'u(x) = sin(2\pi{x}/L)';
        embed = 0;
    case 3
        yf = @(x) 30*x.^2.*(L-x).^2.*sin(2*pi*x/L).^4;
        fstr = 'u(x) = 30x^2(L-x)^2sin(2\pi{x}/L)^4';
        embed = 0;
    case 4
        yf = @(x) 1./(1+(x/L).^2);
        yf = @(x) 1./(1+(x/L).^2) - ((L-x)/L + 1/2*x/L);
        fstr = 'u(x) = 1/(1+(x/L)^2)';
        embed = embed_cushion;
        embed = 0;
    case 5
        yf = @(x) 1./(1+(x/L).^2)-(1-.5*(x/L));
        fstr = 'u(x) = 1/(1+(x/L)^2)+.5(x/L)-1';
        embed = 0;
    case 6
        yf = @(x) sinh(3/L*x)./(1+cosh(3/L*x)) - sinh(3)/(1+cosh(3))*x/L;
        fstr = 'u(x) = sinh(3x/L)./(1+cosh(3x/L)) - (sinh(3)x)/((1+cosh(3))L)';
        embed = embed_cushion;
        embed = 0;
    case 7
        fstr = 'y(x)=cos(x)+e^{-(x-1)^2}-e^{-(x+1)^2}';
        yf = @(x) cos(x)+exp(-(x-1).^2)-exp(-(x+1).^2) - ((L-x)/L + (cos(L)+exp(-(L-1).^2)-exp(-(L+1).^2))*x/L);
        embed = embed_cushion;
        embed = 0;
    case 8 % this family is from the Hubbert-Muller paper on thin plate spline interpolation
        HMN = fpar(1);
        gamval = fpar(2);
        CC = 1/(.5-gamval)^2;
        yf = @(x) CC^(HMN).*(max(0,x-gamval)).^HMN.*(max(0,(1-gamval)-x)).^HMN;
        fstr = ['y(x) = 10^{',num2str(HMN+1),'} (max(0,x-(1/4)))^{',num2str(HMN),'}(max(0,(3/4)-x)^{',num2str(HMN),'}'];
        embed = 0;
    case 9 % this function is good for beta=2, with fpar = [3,.0567]
        HMN = fpar(1);
        gamval = fpar(2);
        CC = 1/(.5-gamval)^2;
        %yf = @(x) CC^(HMN).*(max(0,x-gamval)).^HMN.*(max(0,(1-gamval)-x)).^HMN.*(2*(.5-x)-3*exp(-30^2*(x-.4).^2));
        yf = @(x) CC^(HMN).*(max(0,x-gamval)).^HMN.*(max(0,(1-gamval)-x)).^HMN.*exp(-25^2*(x-.4).^2);
        fstr = ['y(x) = 10^{',num2str(HMN+1),'} (max(0,x-(1/4)))^{',num2str(HMN),'}(max(0,(3/4)-x)^{',num2str(HMN),'}'];
        embed = 0;
    otherwise
        error('This function does not exist')
end

% Selecting some points in the domain, not including the boundary
aa = embed*L; bb = (1-embed)*L;
xx = pickpoints(aa,bb,NN);
yy = yf(xx);

% Set up the error vectors to store the results
errvec = zeros(length(Nvec),length(epvec));
errvecd = zeros(length(Nvec),length(epvec));
errvecp = zeros(length(Nvec),1);
errvecs = zeros(length(Nvec),1);


i = 1;
for N=Nvec
    K_solve = zeros(N);
    K_eval = zeros(NN,N);
    
    if embed==0
    % Don't consider the homogeneous end points because they are fixed
        [x,spacestr] = pickpoints(aa,bb,N+2,spaceopt);
        x = x(2:end-1);
    else
        [x,spacestr] = pickpoints(aa,bb,N,spaceopt);
    end
    y = yf(x);
    
    fprintf('N=%d\t',N)
    k = 1;
    for ep=epvec
        MQR = mqr_solve(x,y,L,ep,beta);
        yp = mqr_eval(MQR,xx);
        errvec(i,k) = errcompute(yp,yy);
        
        for j=1:N
            K_solve(:,j) = cmatern(x,x(j),L,ep,beta);
            K_eval(:,j) = cmatern(xx,x(j),L,ep,beta);
        end
        
        % Ill-conditioning will show up in the graph
        warning off
        [b,r] = linsolve(K_solve,y);
        warning on
        
        yp = K_eval*b;
        errvecd(i,k) = errcompute(yp,yy);
        fprintf('%d ',k)
        k = k+1;
    end
    fprintf('\n')
    
    for j=1:N
        K_solve(:,j) = ppsplinekernel(x,x(j),L,beta);
        K_eval(:,j) = ppsplinekernel(xx,x(j),L,beta);
    end
    [b,r] = linsolve(K_solve,y);
    r
    yp = K_eval*b;
    errvecp(i) = errcompute(yp,yy);
    
    if embed==0
        x = [0; x; L];
        y = yf(x);
    end
    switch beta
        case 1
            yp = interp1(x,y,xx);
            errvecs(i) = errcompute(yp,yy);
        case 2
            % Custom routine modeled after Moler's code from NCM
            % yp = splinetx_natural(x,y,xx);
            % Routine from curve fitting toolbox (also only exists for beta=2)
            pp = csape(x,y,'variational');
            yp = ppval(pp,xx);
            errvecs(i) = errcompute(yp,yy);
        otherwise
            fprintf('Cannot do exact test with beta>2\n')
    end
    
    i = i+1;
end

clf reset
% Direct method plots:
% hdLowN = loglog(epvec,errvecd(1,:),'-bx');
% hdMedN = loglog(epvec,errvecd(2,:),'-g+');
% hdHighN = loglog(epvec,errvecd(3,:),'-r^');
% MaternQR plots:
hqLowN = loglog(epvec,errvec(1,:),'b','LineWidth',3);
hold on
hqMedN = loglog(epvec,errvec(2,:),'g','LineWidth',3);
hqHighN = loglog(epvec,errvec(3,:),'r','LineWidth',3);

% PP Spline plots:
% hpLowN = loglog(epvec,errvecp(1)*ones(size(epvec)),'-.b','LineWidth',2);
% hpMedN = loglog(epvec,errvecp(2)*ones(size(epvec)),'-.g','LineWidth',2);
% hpHighN = loglog(epvec,errvecp(3)*ones(size(epvec)),'-.r','LineWidth',2);

% Natural Spline plots:
if beta==1 || beta==2
    hsLowN = loglog(epvec,errvecs(1)*ones(size(epvec)),'--b','LineWidth',2);
    hsMedN = loglog(epvec,errvecs(2)*ones(size(epvec)),'--g','LineWidth',2);
    hsHighN = loglog(epvec,errvecs(3)*ones(size(epvec)),'--r','LineWidth',2);
end
hold off

xlabel('\epsilon')
ylabel('absolute error')

% Plot labels:
hdHighNLabel = ['Direct (N = ',num2str(Nvec(3)),')'];
hqHighNLabel = ['MaternQR (N = ',num2str(Nvec(3)),')'];
hpHighNLabel = ['PP Spline (N = ',num2str(Nvec(3)),')'];
hsHighNLabel = ['Natural Spline (N = ',num2str(Nvec(3)),')'];
hdMedNLabel = ['Direct (N = ',num2str(Nvec(2)),')'];
hqMedNLabel = ['MaternQR (N = ',num2str(Nvec(2)),')'];
hpMedNLabel = ['PP Spline (N = ',num2str(Nvec(2)),')'];
hsMedNLabel = ['Natural Spline (N = ',num2str(Nvec(2)),')'];
hdLowNLabel = ['Direct (N = ',num2str(Nvec(1)),')'];
hqLowNLabel = ['MaternQR (N = ',num2str(Nvec(1)),')'];
hpLowNLabel = ['PP Spline (N = ',num2str(Nvec(1)),')'];
hsLowNLabel = ['Natural Spline (N = ',num2str(Nvec(1)),')'];

ah1 = gca;
ah2=axes('position',get(gca,'position'), 'visible','off');
legend(ah1,[hsLowN,hsMedN,hsHighN],'N=12, Spline','N=24, Spline','N=48, Spline',1)
legend(ah2,[hqLowN,hqMedN,hqHighN],'N=12, CMatern','N=24, CMatern','N=48, CMatern',2)

% if beta==1 || beta==2
%     legend([ hqLowN hdLowN hpLowN hsLowN ], hqLowNLabel, hdLowNLabel, hpLowNLabel, hsLowNLabel, 'location','Best')
% else
%     legend([ hqLowN hdLowN hpLowN ], hqLowNLabel, hdLowNLabel, hpLowNLabel, 'location','Best')
% end



