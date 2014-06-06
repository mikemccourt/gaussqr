function savefig(h,fname,style,newBaseDir)
% function savefig(h,fname,style,newBaseDir)
%
% This function saves the plot handle h with the standard book format
% Input: h     - plot handle (use gcf for current figure)
%        fname - string for name of the file
%        style - which plotting style to use (see bottom) <default=0>
%                0 - No Input (use whatever is currently in the figure)
%                1 - Black and White
%                2 - Color
%                3 - Grayscale
%         newBaseDir - <optional> directory for saving figures
%
% The idea of newBaseDir is to define the directory for your pitures once
% and then just pass in different names for the figures you want to save.
% You can set the default directory in rbfsetup.m with the variable
%      GAUSSQR_PARAMETERS.FIGURE_DIRECTORY
% If you don't want to set the directory in rbfestup.m, you can
% pass the directory in newBaseDir.  That variable is set as persistent, so
% if you pass it once at the beginning of a session you won't have to pass
% it again.  If you reset your Matlab session, you will have to start over.
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
defaultDir = GAUSSQR_PARAMETERS.FIGURE_DIRECTORY;

% So that you can pass newBaseDir once at the start of a session
persistent baseDir

% This switch statement assigns the appropriate value to baseDir
% Either the persistent value, the default value or the new value
switch nargin
    case {2,3} % Handle default baseDir choice
        if nargin==2 % Default style value
            style = 0;
        end

        % Checks to make sure baseDir exists
        if exist(defaultDir,'dir')
            if isempty(baseDir) || not(exist(baseDir,'dir'))
                baseDir = defaultDir;
                if not(exist(baseDir,'dir'))
                    warning('Persistent directory baseDir=%s does not exist; deleted since first call?\n Resetting to %s',baseDir,defaultDir)
                end
            end
        else
            if exist(baseDir,'dir')
                warning('Default figure directory %s does not exist; deleted after calling rbfsetup?\n Using persistent directory %s',defaultDir,baseDir)
            else
                error('Default figure directory %s does not exist; deleted after calling rbfsetup?\n Pass a newBaseDir to this function (pwd for current directory)',defaultDir)
            end
        end
    case 4
        if exist(newBaseDir,'dir')
            baseDir = newBaseDir;
        else
            error('Requested directory %s does not exist',newBaseDir)
        end
    otherwise
        error('Incorrent calling sequence')
end

%
% Black & White
%   Rendering
%      Colorspace black and white
%      Custom color w
%      Resolution 600
%      Keep axis limits
%      Show uicontrols
%   Fonts
%      Custom - fixed font size 14
%      Custom name - Helvetica ?
%   Lines
%      Fixed width 2 points
% Color
%   Rendering
%      Colorspace CMYK
%      Custom color w
%      Resolution 300
%      Keep axis limits
%      Show uicontrols
%   Fonts
%      Custom - fixed font size 14
%      Custom name - Helvetica ?
%   Lines
%      Fixed width 2 points
% Grayscale
%   Rendering
%      Colorspace grayscale
%      Custom color w
%      Resolution 300
%      Keep axis limits
%      Show uicontrols
%   Fonts
%      Custom - fixed font size 14
%      Custom name - Helvetica ?
%   Lines
%      Fixed width 2 points