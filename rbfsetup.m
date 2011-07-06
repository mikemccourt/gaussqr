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
