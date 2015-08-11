function gqr_savefig(h,fname,style,newBaseDir)
% function gqr_savefig(h,fname,style,newBaseDir)
%
% This function saves the plot handle h with the standard book format
% Input: h - figure number (handle) that you want plotted
%        fname - string for name of the file (no .eps or .fig)
%        style - which plotting style to use (see bottom) <default=0>
%                0 - No Input (use whatever is currently in the figure)
%                1 - Black and White
%                2 - Color
%                3 - Grayscale
%         newBaseDir - <optional> directory for saving figures
%
% You must pass a figure handle, not an axes handle, for this to work
% Example: h = figure;plot(x,y) --- not --- h = plot(x,y)
% I honestly can't figure out why Matlab can't make it work with just the
% axis handle, but it can't and I can't figure out how to extract the
% necessary information.
% To pass the current figure, use h = gcf
%
% Example: gqr_savefig(gcf,'book\happy')
%      Saves in C:\Users\ironmike\Documents\fasshauer\gaussqr\book\happy.eps
%        because the default base directory is the gaussqr main directory
% Example: gqr_savefig(gcf,'happy',2,'/home/mccomic/Documents')
%      Saves a color eps (and figure) to /home/mccomic/Documents/happy.eps
% Example: gqr_savefig(gcf,'Documents/happy',2,'/home/mccomic')
%      Same output as above
% Example: gqr_savefig(h,'happy',2,'/home/mccomic/Documents')
%          gqr_savefig(h,'sad',3)
%      Saves color /home/mccomic/Documents/happy.eps and
%      grayscale /home/mccomic/Documents/sad.eps
%
% The idea of newBaseDir is to define the directory for your pitures once
% and then just pass in different names for the figures you want to save.
% You can set the default directory in rbfsetup.m with the variable
%      GAUSSQR_PARAMETERS.FIGURE_DIRECTORY
% If you don't want to set the directory in rbfestup.m, you can
% pass the directory in newBaseDir.  That variable is set as persistent, so
% if you pass it once at the beginning of a session you won't have to pass
% it again.  If you reset your Matlab session, you will have to start over.
%
% DEVELOPER'S NOTE: Better work needs to be done to handle surface plots
% with colorbars and colors on the surface: example
% a = get(gcf,'Children')
% a(3)
% set(a(3),'fontsize',14)

global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
defaultDir = GAUSSQR_PARAMETERS.FIGURE_DIRECTORY;
dirSlash = GAUSSQR_PARAMETERS.DIRECTORY_SLASH;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
defaultFS = GAUSSQR_PARAMETERS.FIGURE_FONTSIZE;

% So that you can pass newBaseDir once at the start of a session
persistent baseDir

% This switch statement assigns the appropriate value to baseDir
% Either the persistent value, the default value or the new value
switch nargin
    case {2,3} % Handle default baseDir choice
        fs = defaultFS;
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
            if alertuser && exist(baseDir,'dir')
                warning('Default figure directory %s does not exist; deleted after calling rbfsetup?\n Using persistent directory %s',defaultDir,baseDir)
            else
                error('Default figure directory %s does not exist; deleted after calling rbfsetup?\n Pass a newBaseDir to this function (pwd for current directory)',defaultDir)
            end
        end
    case 4
        if ~ischar(newBaseDir)
            fs = newBaseDir;
            baseDir = defaultDir;
        else
            if exist(newBaseDir,'dir')
                baseDir = newBaseDir;
            else
                error('Requested directory %s does not exist',newBaseDir)
            end
            fs = defaultFS;
        end
    otherwise
        error('Incorrent calling sequence')
end

% Form the directory and file name to save to
if baseDir(end)~=dirSlash
    baseDir = strcat(baseDir,dirSlash);
end
saveDirBase = strcat(baseDir,fname);
epsDir = strcat(saveDirBase,'.eps');

% Check to make sure that you can write to the directory you want to
[dirstat,dirinfo] = fileattrib(baseDir);
if dirstat~=1 || ~dirinfo.UserWrite
    error('You do not appear to have permission to write to %s\nYou may not have permission, or that may not be a valid directory',baseDir)
end

% Use the appropriate color plotting pattern as requested by the user
switch style
    case 0 % Just save the fig/eps of the handle
        saveas(h,saveDirBase,'fig')
        print(h,'-depsc2',epsDir)
        epsDir = strcat(saveDirBase,'.png');
        print(h,'-dpng',epsDir)
    case {1,2,3}
        % Font Size 14        
        % This command should change the axis labels, title
        set(findall(findobj(h),'Type','text'),'FontSize',fs)
        % This command should change the xticklabels and legend
        c = get(h,'Children');
%         set(c(isgraphics(c,'axes')),'FontSize',fs)
        set(c(isprop(c,'fontsize')),'FontSize',fs)

        % Keep axis limits
        % I think that this is equivalent to "keep axis limits" in the sense that
        % the limits will not change unless the user changes them
        % This I am not sure about though
        set(get(h,'CurrentAxes'),'XLimMode','manual')
        set(get(h,'CurrentAxes'),'YLimMode','manual')
        set(get(h,'CurrentAxes'),'ZLimMode','manual')

        % Custom background color 'w'
        % This should set the background color to white, which is the default
        % There may be some difference between this and transparent that I don't
        % yet fully appreciate
        set(h,'Color','w')

        % Show uicontrols
        % I'm not entirely sure what this means, but I don't think it has anything
        % to do with us because we don't have any GUIs.  Therefore I am not
        % implementing this right now

        % Lines width 2 points
        % Changes all line plots to have width 2
        % The way this is written should omit any non line plots, but we might want
        % something comparable for surfact or bar plots
        set(findobj(get(get(h,'CurrentAxes'),'Children'),'type','line'),'linewidth',2)
        
        % Save the figure
        saveas(h,saveDirBase,'fig')
        % Save the eps in the desired format
        switch style
            case 1 % Black and White
                % Must change image colors to black lines
                % Not sure how to handle this for surface plots
                % Note this won't change the figure, in case you want the
                % figure to still be in color
                set(findobj(h,'type', 'line'),'color','k');
                set(findobj(findobj(h,'type', 'line'),'-not','MarkerFaceColor','none'),'MarkerFaceColor','k')
                print(h,'-deps2','-r600',epsDir)
            case 2 % Color
                print(h,'-depsc2','-r300','-cmyk',epsDir)
            case 3 % Grayscale
                print(h,'-deps2','-r300',epsDir)
            otherwise
                error('How on Earth did you reach this point?')
        end
    otherwise
        error('Unacceptable plotting pattern style=%g',style)
end

% Standards for all plots:
%   Custom background color 'w'
%   Keep axis limits
%   Show uicontrols
%   Font size 14
%   Lines width 2 points
%
% Black & White
%   Rendering
%      Colorspace black and white
%      Resolution 600
% Color
%   Rendering
%      Colorspace CMYK
%      Resolution 300
% Grayscale
%   Rendering
%      Colorspace grayscale
%      Resolution 300