% ex6_mqr
% This puts together convergence order pictures for MaternQR
% We consider two choices:
%     General functions with arbitrary "BC" embedded in [0,L]
%     Functions which are homogeneous on [0,L]
% These tests will fix L, ep and beta and consider an increase in N
% Points will be evenly spaced at first, althought we could consider other
%   choices of point distributions
% For the embedded setting, the embedding cushion will be fixed.  The
%   effect of the embedding cushion still needs to be considered.
% I guess we will need to consider both the kernel version and the RBF-QR
%   version of the solve.  Need to think about that because one may be
%   better/faster than the other in different conditions
rbfsetup
global GAUSSQR_PARAMETERS
symavail = GAUSSQR_PARAMETERS.SYMBOLIC_TOOLBOX_AVAILABLE;
GAUSSQR_PARAMETERS.ERROR_STYLE = 2;
GAUSSQR_PARAMETERS.NORM_TYPE = inf;

% The range of N values to consider
Nvec = [8,16,32,64,128,256];
% The orders (smoothness) of the kernel to consider
betavec = 1:4;
% The kernel shape parameter
ep = 1;
% The length of the domain
L = 1;
% The embedding width for nonhomogeneous functions
embed_cushion = .1;
% The number of evenly spaced points at which to sample error
NN = 397;
% Choice of interpolation function and associated parameters
fopt = 8; fopt = 15;
fpar = [8,.0567]; fpar = 0;

% This determines how many extra basis functions should be added to the
% RBF-QR evaluation to get the necessary accuracy: M = Mfactor*N
% Actually picking a good value for this may be difficult
% I guess the minimum should be something like 1.1
GAUSSQR_PARAMETERS.MAX_EXTRA_EFUNC = -700;

% This is the function we are interested in considering
% Depending on which function consider, it will choose embedding
switch fopt
    case 1
        yf = @(x) sin(2*pi*x/L) + 1;
        fstr = 'u(x) = sin(2\pi{x}/L)+1';
        embed = embed_cushion;
    case 2
        yf = @(x) sin(2*pi*x/L);
        fstr = 'u(x) = sin(2\pi{x}/L)';
        embed = 0;
    case 3     %optimal beta = 4?
        yf = @(x) 30*x.^2.*(L-x).^2.*sin(2*pi*x/L).^4;
        fstr = 'u(x) = 30x^2(L-x)^2sin(2\pi{x}/L)^4';
        embed = 0;
    case 4
        yf = @(x) 1./(1+(x/L).^2);
        fstr = 'u(x) = 1/(1+(x/L)^2)';
        embed = embed_cushion;
    case 5 %optimal beta = 3?
        yf = @(x) 1./(1+(x/L).^2)-(1-.5*(x/L));
        fstr = 'u(x) = 1/(1+(x/L)^2)+.5(x/L)-1';
        embed = 0;
    case 6
        yf = @(x) sinh(3/L*x)./(1+cosh(3/L*x));
        fstr = 'u(x) = sinh(3x/L)./(1+cosh(3x/L))';
        embed = embed_cushion;
    case 7
        fstr = 'y(x)=cos(x)+e^{-(x-1)^2}-e^{-(x+1)^2}';
        yf = @(x) cos(x)+exp(-(x-1).^2)-exp(-(x+1).^2);
        embed = embed_cushion;
 %--------------------------------------------------------------
 % Casey and Will's test functions
    case 8 % this family is from the Hubbert-Muller paper on thin plate spline interpolation
        % meant to be used on the unit interval
        %which beta is optimal?
        HMN = fpar(1);
%        yf = @(x) 10^(HMN+1).*(max(0,x-(1/4))).^HMN.*(max(0,(3/4)-x)).^HMN;
        gamval = fpar(2);
        CC = 1/(.5-gamval)^2;
        yf = @(x) CC^(HMN).*(max(0,x-gamval)).^HMN.*(max(0,(1-gamval)-x)).^HMN;
        fstr = ['y(x) = 10^{',num2str(HMN+1),'} (max(0,x-(1/4)))^{',num2str(HMN),'}(max(0,(3/4)-x)^{',num2str(HMN),'}'];
        embed = 0;
    case 9
        fstr = 'y(x) = 1-x';
        yf = @(x) ones(size(x));
        embed = 0;
    case 10
        fstr = 'y(x) = -(x-.5)^2+.25';
        yf = @(x) -(x-.5).^2+.25;
        embed = 0;
    case 11
        fstr = 'y(x) = x-2x^3+x^4';
        yf = @(x) x - 2*x.^3 + x.^4;
        embed = 0;
    case 12
        fstr = 'example with first deriv'
        aa = (1/9-8/105)*L^(5/2);
        yf = @(x) 8/105*x.^(7/2) - 1/9*sqrt(L)*x.^3 + aa*x;
        embed = 0;
    case 13
        yf = @(x) x.^5.*(x-L).^5;
        embed = 0;
    case 14
        yf = @(x) -sin(3*pi*x);
        embed = 0;
    case 15 % Choosing C=384/77 makes this a Type IV L-spline
        C = fpar(1);
        fstr = sprintf('D5 broken, C=%g',C);
        yf = @(x) 1/27720*((1155*C-5888)*x + (14080-2310*C)*x.^3 + 1155*C*x.^4-8192*x.^(15/4));
        embed = 0;
    case 16
        fstr = 'Does not satisfy beta=2 conditions';
        yf = @(x) 1/21*(x-1/2).*(-64*abs(x-1/2).^(3/4)+2^(1/4)*(28*(x-1/2).^2+25));
        embed = 0;
    otherwise
        error('This function does not exist')
end

% At some point I'll move the solves into a separate function and consider
% warnings there.  For now I'm tired of seeing the warnings pop up.
warning off

errvec = zeros(length(betavec),length(Nvec));
k = 1;
j = 1;

for N=Nvec
    if embed==0
    % Don't consider the homogeneous end points because they are fixed
        [x,spacestr] = pickpoints(0,L,N+2);
        x = x(2:end-1);
%         if fopt==8 % Include the discontinuities
%             [h,i] = min(abs(x-.5));
%             x(i) = .5;
%             [h,i] = min(abs(x-log(1.25)));
%             x(i) = log(1.25);
%             [h,i] = min(abs(x-(1-log(1.25))));
%             x(i) = 1-log(1.25);
%         end
    else
        [x,spacestr] = pickpoints(embed*L,(1-embed)*L,N);
    end
    y = yf(x);
    xx = pickpoints(embed*L,(1-embed)*L,NN);
    yy = yf(xx);
    
    j = 1;
    for beta=betavec
        if beta==1% || beta==2 % Work with the kernel form
            K_solve = zeros(N);
            K_eval = zeros(NN,N);
            for n=1:N
                K_solve(:,n) = cmatern(x,x(n),L,ep,beta);
                K_eval(:,n) = cmatern(xx,x(n),L,ep,beta);
            end
            b = K_solve\y;
            yp = K_eval*b;
        else % Work with the series form
            MQR = mqr_solve(x,y,L,ep,beta);
            yp = mqr_eval(MQR,xx);
        end

        errvec(j,k) = errcompute(yp,yy);
        j = j + 1;
    end
    
    k = k + 1;
end
warning on

% Plot RMS error as beta and N vary:
% close all
% errorplot = figure('NumberTitle','off','Name',num2str(fopt));
loglog(Nvec,errvec,'linewidth',2);
% semilogy(Nvec,errvec,'linewidth',2)
xlabel('input points N');
ylabel('absolute error');
title(fstr)
%title(strcat(sprintf('x \\in [0,%g] ',L),sprintf('     \\epsilon = %g',ep)));

% Build and create the legend
legendvals = [];
for beta=betavec
    legendvals = [legendvals;sprintf('beta=%d',beta)];
end
legend(cellstr(legendvals),'location','southwest')
%legend('\beta = 1','\beta = 2','\beta = 3','\beta = 4','\beta = 5','location','best');