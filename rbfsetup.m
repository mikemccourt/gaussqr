% rbfsetup.m
% This file puts the appropriate directories in your path
% This is called a function because I don't want the user to
%   see these internal variables after this is called
function rbfsetup
P = path;
thisDir = pwd;

if(length(strfind(thisDir,'\'))>0) % We are in Windows
    sourceDir = strcat(thisDir,'\source');
    exampleDir = strcat(thisDir,'\examples');
elseif(length(strfind(thisDir,'/'))>0) % We are in Unix
    sourceDir = strcat(thisDir,'/source');
    exampleDir = strcat(thisDir,'/examples');
end
path(P,sourceDir);
P = path;
path(P,exampleDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup global constants and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GAUSSQR_PARAMETERS

% At what point should the asymptotic approximation to Hermite be used
% Anything between 35-60 you shouldn't be able to tell the difference
% Beyond 70 it's pretty necessary
GAUSSQR_PARAMETERS.HERMITE_ASYMPTOTIC_MIN_INDEX = 40;

% Sets up the polynomial coefficients for later use
GAUSSQR_PARAMETERS.HERMITE_COEFFICIENTS = cell(GAUSSQR_PARAMETERS.HERMITE_ASYMPTOTIC_MIN_INDEX,1);
for k=1:GAUSSQR_PARAMETERS.HERMITE_ASYMPTOTIC_MIN_INDEX
    GAUSSQR_PARAMETERS.HERMITE_COEFFICIENTS{k} = HermitePoly(k-1);
end

% Use logarithms when computing rbfphi
% Shouldn't be an issue unless you have M>200 or x>20
% In general, you should use logs
% Set to 0 to turn off, 1 to turn on
GAUSSQR_PARAMETERS.RBFPHI_WITH_LOGS = 1;

% Tolerance for using asymptotic approximation of exponential within rbfphi
% The term sqrt(1+(2ep/a)^2)-1 pops up, and you can get cancelation
% This will cause the switch to sqrt(1+(2ep/a)^2)-1 = 2(ep/a)^2-2(ep/a)^4
% The switch occurs when (1+(2ep/a)^2)^(1/4)-1<tol
GAUSSQR_PARAMETERS.RBFPHI_EXP_TOL = 1e-4;

end
