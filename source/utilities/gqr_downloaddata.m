function retval = gqr_downloaddata(filestr,webbase,saveloc)
% This function tries to download data from the gaussqr website, or from
% some other location if one is desired
% It will store the data in the gaussqr data directory as defined in
% rbfsetup
% NOTE: This will not overwrite existing files in your directory, so if you
% want to download an updated version of a file, you must manually delete
% the existing file
%
%    download_occurred = gqr_downloaddata(filestr)
% Downloads filestr from the GaussQR website into the data directory
% defined in rbfsetup
%    Input: filestr - string of file to download
%                     Ex: 'sphereMDpts.mat'
%    Output: download_occurred - 1 if the file was downloaded, 0 else
%
%    download_occurred = gqr_downloaddata(filestr,webbase,saveloc)
%    Inputs: webbase - (optional) website from which to download
%                      <default=http://math.iit.edu/~mccomic/gaussqr>
%            saveloc - (optional) directory to save to
%                      <default=GAUSSQR_PARAMETERS.DATA_DIRECTORY>
%
%    gaussqr_list = gqr_downloaddata()
%    Output: gaussqr_list - returns a list of all the data files
%                           currently recognized by GaussQR
% This may now be obsolete since this list is available in
% GAUSSQR_PARAMETERS, but I'm leaving this here because it doesn't hurt
% anything.
%
% When downloading, this function does not have any return values, it
% simply throws an error if the desired operation cannot be completed.
% I will reconsider this at some point, if needed, and perhaps return
% different values to gracefully exit the program
%
% Developer Note: Allow for cell array filestr so that multiple files can
% be downloaded from the same place with one function call
%
% Developer Note: websave has been implemented in favor of urlwrite, but it
% should probably be implemented in a try/catch block to better notice user
% errors and return helpful error messages.
global GAUSSQR_PARAMETERS
if ~isstruct(GAUSSQR_PARAMETERS)
    error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
end
dataDir = GAUSSQR_PARAMETERS.DATA_DIRECTORY;
dirslash = GAUSSQR_PARAMETERS.DIRECTORY_SLASH;
gqrwebloc = GAUSSQR_PARAMETERS.WEB_DIRECTORY;
alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;
gaussqr_list = GAUSSQR_PARAMETERS.AVAILABLE_DATA;

if nargin>0
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
    % If it does, do nothing; otherwise, attempt to download it
    % There is a rehash here to make sure Matlab knows this new file is in
    % the path
    download_occurred = 0;
    if not(exist(fullfilename,'file'))
        websave(fullfilename,webaddress,'Timeout',20);
        if alertuser
            fprintf('File %s was downloaded to %s\n',filestr,saveloc)
        end
        rehash
    end
    
    % Return a value in case the user wants to see but wants alerts off
    if nargout==1
        retval = download_occurred;
    end
else
    retval = gaussqr_list;
end