function gaussqr_list = gqr_downloaddata(filestr,webbase,saveloc)
% This function tries to download data from the gaussqr website, or from
% some other location if one is desired
% It will store the data in the gaussqr data directory as defined in
% rbfsetup
% NOTE: This will not overwrite existing files in your directory, so if you
% want to download an updated version of a file, you must manually delete
% the existing file
%
%    gqr_downloaddata(filestr)
% Downloads filestr from the GaussQR website into the data directory
% defined in rbfsetup
%    Input: filestr - string of file to download
%                     Ex: 'sphereMDpts.mat'
%
%    gqr_downloaddata(filestr,webbase,saveloc)
%    Inputs: webbase - (optional) website from which to download
%                      <default=http://math.iit.edu/~mccomic/gaussqr>
%            saveloc - (optional) directory to save to
%                      <default=GAUSSQR_PARAMETERS.DATA_DIRECTORY>
%
%    gaussqr_list = gqr_downloaddata()
%    Output: gaussqr_list - returns a list of all the data files
%                           currently recognized by GaussQR
%
% When downloading, this function does not have any return values, it
% simply throws an error if the desired operation cannot be completed.
% I will reconsider this at some point, if needed, and perhaps return
% different values to gracefully exit the program
%
% Developer Note: Allow for cell array filestr so that multiple files can
% be downloaded from the same place with one function call
%
% Developer Note: I'm pretty sure that newer Matlab wants use to use
% websave instead of urlwrite.  I'll have to look into that at some point
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
dataDir = GAUSSQR_PARAMETERS.DATA_DIRECTORY;
dirslash = GAUSSQR_PARAMETERS.DIRECTORY_SLASH;
gqrwebloc = GAUSSQR_PARAMETERS.WEB_DIRECTORY;

% Below is the list of files that are currently recognized as GaussQR
% supported.  This could become out of date, so please alert me if you
% notice inconsistencies
gaussqr_list = {'sphereMDpts_data.mat'};

if nargout==0
    % If the user did not pass webbase, the file to be downloaded must be
    % in gaussqr_list or this will throw an error
    if not(exist('webbase','var'))
        if not(any(cellfun(@(b)strcmp(b,filestr),gaussqr_list)))
            error('You are downloading a file which is not registered with GaussQR, %s\nCheck the spelling, but this may be our error; please alert us if needed.',filestr)
        else
            webbase = strcat(gqrwebloc,'/data');
        end
    end
    % Create the full web address
    if not(strcmp(webbase(end),'/'))
        webbase = strcat(webbase,'/');
    end
    webaddress = strcat(webbase,filestr);
    
    % Determine if the user passed a viable save location
    if not(exist('saveloc','var'))
        saveloc = dataDir;
    end
    [dirstat,dirinfo] = fileattrib(saveloc);
    if dirstat~=1 || ~dirinfo.UserWrite
        error('You do not appear to have permission to write to %s\nYou may not have permission, or that may not be a valid directory',saveloc)
    end
    % Append a slash as needed
    if not(strcmp(saveloc(end),dirslash))
        saveloc = strcat(saveloc,dirslash);
    end
    % Form the full file name to check if such a file already exists
    fullfilename = strcat(saveloc,filestr);
    
    % Check if the file already exists in that directory
    % If it does, do nothing
    if not(exist(fullfilename,'file'))
        [filestr,downloadSuccessful] = urlwrite(webaddress,fullfilename,'Timeout',20);
        if ~downloadSuccessful
            error('Download of data failed or could not be written to %s',filestr)
        end
    end
    
    % To prevent the code from returning a value when one is not desired
    clear gaussqr_list
end