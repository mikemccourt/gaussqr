% ex5_SQR
% This puts together convergence order pictures for SineQR
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

% The range of N values to consider
Nvec = 10:5:100;
% The orders (smoothness) of the kernel to consider
betavec = 1:5;
% The kernel shape parameter
ep = .1;
% The length of the domain
L = 1;
% The embedding width for nonhomogeneous functions
embed_cushion = .1;
% The number of evenly spaced points at which to sample error
NN = 397;

% This determines how many extra basis functions should be added to the
% RBF-QR evaluation to get the necessary accuracy: M = Mfactor*N
% Actually picking a good value for this may be difficult
% I guess the minimum should be something like 1.1
Mfactor = 4.5;

% Define the eigenfunctions and eigenvalues
sinfunc = @(n,L,x) sqrt(2/L)*sin(pi*x*n/L);
lamfunc = @(n,L,ep,beta) ((pi*n/L).^2+ep^2).^(-beta);

% This is the function we are interested in considering
% Depending on which function consider, it will choose embedding
fopt = 3;
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
        fstr = 'u(x) = 1/(1+(x/L)^2)';
        embed = embed_cushion;
    case 5
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
        testfuncN = 10;
        yf = @(x) 10^(testfuncN+1).*(max(0,x-(1/4))).^testfuncN.*(max(0,(3/4)-x)).^testfuncN;
        fstr = ['y(x) = 10^{',num2str(testfuncN+1),'} (max(0,x-(1/4)))^{',num2str(testfuncN),'}(max(0,(3/4)-x)^{',num2str(testfuncN),'}'];
        embed = 0;
    case 9 % ADDITIVE boundary condition forcing
        if symavail % requires Symbolic Math Toolbox :(
            bSatTestFunc = sym(franke(sym('x'),0.5)); % this must be a symbolic expression
            
            % the test function will satisfy boundary conditions for all *even*
            % derivatives up to the (2*bSatDegree)th derivative:
            bSatDegree = 5; % a value of -1 here doesn't force any conditions to be satisfied
            % specify desired boundary values here:
            lbval = 0;
            ubval = 0;
            
            symyf = forceBCsatADD(bSatTestFunc,bSatDegree,0,L,lbval,ubval);
            symf = symf / 10^(2*bSatDegree); % we're adding polynomials on the order of 2*bSatDegree, so we should normalize by that much
            yf = matlabFunction(symyf);
            fstr = char(simplify(symyf));
            embed = 0;
        else
            error('Cannot call this function without the symbolic toolkit available')
        end
    case 10 % MULTIPLICATIVE boundary condition forcing
        if symavail % requires Symbolic Math Toolbox :(
            bSatTestFunc = sym(franke(sym('x'),0.5)); % this must be a symbolic expression
            
            % the test function will satisfy boundary conditions for *all*
            % derivatives up to the (2*bSatDegree)th derivative:
            bSatDegree = 5; % a value of -1 here doesn't force any conditions to be satisfied
            
            symyf = forceBCsatMULT(bSatTestFunc,bSatDegree,0,L);
            yf = matlabFunction(symyf);
            fstr = char(simplify(symyf));
            embed = 0;
        else
            error('Cannot call this function without the symbolic toolkit available')
        end
 %--------------------------------------------------------------
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
    else
        [x,spacestr] = pickpoints(embed*L,(1-embed)*L,N);
    end
    y = yf(x);
    xx = pickpoints(embed*L,(1-embed)*L,NN);
    yy = yf(xx);
    
    j = 1;
    for beta=betavec
        if beta==1 % Work with the kernel form
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

semilogy(Nvec,errvec,'linewidth',2)
xlabel('input points N')
ylabel('RMS relative error')
title(strcat(fstr,sprintf('x\\in[%g,%g], ',embed*L,(1-embed)*L),sprintf('\\epsilon=%g',ep)))
legend('\beta=1','\beta=2','\beta=3','\beta=4','\beta=5','location','northeast')