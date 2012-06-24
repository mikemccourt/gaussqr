function alpha = gqr_alphasearch(ep,a,b)
% function alpha = gqr_alphasearch(ep,a,b)
% ep - shape parameter for the Gaussians
% a - minimum bound(s) of domain
% b - maximum bound(s) of domain
%   In multiple dimensions these should be a=[low_1 ... low_d]
%                                          b=[high_1 ... high_d]
%   These vectors must be the same size, obviously
%
% This function tries to return the smallest alpha for which orthonormality
% is maintained throughout the domain until a set eigenfunction.
%
% You can set the top eigenfunction for which orthonornmality should be
% maintained in rbfsetup with the value
%      GAUSSQR_PARAMETERS.DEFAULT_ORTH_REQUESTED
%
% If this function cannot find an alpha up to what was requested, it will
% try to get as high as possible.  The index up to which orthogonality was
% maintained is returned as k.
%
% Should change the way this is called, to a symmetric interval

% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
orthmax = GAUSSQR_PARAMETERS.ORTH_INDEX_REQUESTED;
alphamin = GAUSSQR_PARAMETERS.ORTH_MINIMUM_ALPHA;
alphamax = GAUSSQR_PARAMETERS.ORTH_MAXIMUM_ALPHA;
options.TolX = GAUSSQR_PARAMETERS.ORTH_SEARCH_ACCURACY;
inttol = GAUSSQR_PARAMETERS.ORTH_INTEGRATION_TOL;

k = orthmax;

a = a(:);
b = b(:);

if length(a)~=length(b)
    error('Integration bounds have mismatched dimension: len(a)=%d, len(b)=%d',length(a),length(b))
end

% This is a little dummy iteration to walk through the region which is
% specified by the alphamin and alphamax values.  The first value it finds
% which produces good orthonormality is used as an initial guess for the
% optimization program.  If it gets through the walk without finding a spot
% it will take smaller steps or eventually give up.
goodAreaFound = 0;
intfactor = 5;
while intfactor>1.1 && goodAreaFound==0
    intval = gqr_orthintegral_PRIVATE(k,ep,alphamin,a,b);
    alphaguess = alphamin;
    while intval==0 && alphaguess<alphamax
        alphaguess = alphaguess*intfactor;
        intval = gqr_orthintegral_PRIVATE(k,ep,alphaguess,a,b,inttol);
        if intval<0 % Integral didn't converge because alpha too high
            alphaguess = alphamax;
        end
    end
    if alphaguess>=alphamax
        intfactor = intfactor/2;
    else
        goodAreaFound = 1;
    end
end

% Now we solve the actual minimization problem after getting a good bound
% from earlier.
if goodAreaFound == 1
    [alpha,fminval,exitflag] = fminbnd(@(alpha)-gqr_orthintegral_PRIVATE(k,ep,alpha,a,b),alphaguess/intfactor,alphaguess,options);
    if exitflag<1
        warning('Solution not found, consider reducing orthonormality tolerance \n or increasing integration tolerance')
    end
else
    alpha = sqrt(alphamin*alphamax);
    warning('Failed to find an acceptable alpha in [%g,%g]',alphamin,alphamax);
end

end


% This function is only an auxilliary to be used by rbfalphasearch
% It computes the integral of phi_k^2 rho which should equal 1 if
% orthonormality is preserved.
% a and b are the bounds of the integral
% inttol is the accuracy we demand of the quadrature.  If this value is not
% passed, it is set to 1e-6, the default from Matlab.
% orthtol is a value which determines to what point we accept
% computational orthonormality.  If this value is not passed, it is set to
% the default value from rbfsetup.
%
% If we accept the computational orthogonality of the function then we can
% return 1/alpha, otherwise we return 0.
% These values are needed to do the optimization in rbfalphasearch.
%
% Note that for the quadrature, I need to be able to pass row vectors, for
% some reason (can't imagine why they couldn't be column vectors).  That's
% why you see the weird transpose in/transpose out thing.
function integral = gqr_orthintegral_PRIVATE(k,ep,alpha,a,b,inttol,orthtol)
% Import global parameters from rbfsetup
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
if not(exist('inttol'))
    inttol = 1e-6;
elseif length(inttol)==0 % So the user can pass only inttol
    inttol = 1e-6;
elseif inttol<0 || inttol>1
    warning('inttol=%g is irrational, reset to %g',inttol,1e-6)
    inttol = 1e-6;
end
if not(exist('orthtol'))
    orthtol = GAUSSQR_PARAMETERS.DEFAULT_ORTH_TOLERANCE;
elseif orthtol<0 || orthtol>1
    warning('orthtol=%g is irrational, reset to %g',orthtol,GAUSSQR_PARAMETERS.DEFAULT_ORTH_TOLERANCE)
    orthtol = GAUSSQR_PARAMETERS.DEFAULT_ORTH_TOLERANCE;
end
weight = @(a,x) a/sqrt(pi)*exp(-(a*x).^2);

if alpha<=0 % The optimization could pass whatever
    integral = 0;
else
    intappx = 1;
    lastwarn('','')
    
    % The multidimensional integral is actually just a bunch of 1D
    % integrals because of the tensor product structure of the
    % eigenfunctions.
    warning off MATLAB:quadl:ImproperFcnValue
    if exist('quadgk')
        for j=1:length(a)
            intappx = intappx*quadgk(@(x)gqr_phi(k,x',ep,alpha)'.^2.*weight(alpha,x),a(j),b(j),'RelTol',inttol);
        end
    else
        for j=1:length(a)
            intappx = intappx*quadl(@(x)gqr_phi(k,x',ep,alpha)'.^2.*weight(alpha,x),a(j),b(j),inttol);
        end
    end
    warning on MATLAB:quadl:ImproperFcnValue
    
    % Check in case the integral returned a crap value
    % If it did, return -1 to indicate a problem
    % Otherwise return the optimization value
    [matlabwarnMSG,matlabwarnID] = lastwarn;
    if strcmp(matlabwarnID,'MATLAB:quadl:ImproperFcnValue')
        integral = -1;
    else
        integral = (abs(1-intappx)<orthtol)*(1/alpha);
    end
end

end

% Developers note: I need to figure out a way to easily test lower function
% indexes if the highest one the user asks for fails
%
% Eventually, I'd like to switch all these functions to using options
% structs instead of passing individual parameters