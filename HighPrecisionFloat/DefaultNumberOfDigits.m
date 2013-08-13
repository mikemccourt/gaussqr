function NDig = DefaultNumberOfDigits(NDig,varargin)
% Allows the user to set the default value for the number of digits for a new HPF object.
% Usage: DefaultNumberOfDigits(Ndig);  % Resets the default to a new value
% Usage: Ndig = DefaultNumberOfDigits; % Returns the current default
% Usage: DefaultNumberOfDigits reset   % Resets the default to 64 2, the original default
%
% DefaultNumberOfDigits does NOT change the number of digits
% stored in any existing hpf numbers, only in those numbers
% created in the future. Also the user can always override the
% default for any specific number created, by specifying a
% specific second argument in the initial call to hpf. hpf
% numbers can also have this property modified by later calls
% to the augmentdigits method.
%
% Arguments: (input)
%  Ndig - scalar positive integer, or a vector of length 2.
%        If a vector of length 2, then the first element defines
%        the number of decimal digits that will be displayed.
%        (Trailing zeros will generally be dropped.) The second
%        element will define the number of spare or shadow digits
%        to be carried. Since the last digit or so of any floating
%        point number should never be fully trusted, this allows
%        those last few digits to be hidden from view.
%
%        If NDig is a scalar value, then the number of shadow
%        digits will remain unchanged.
%
%        If DefaultNumberOfDigits is called as a command at
%        the MATLAB command line then the argument will be
%        converted from character to numeric form by str2num.
%
%        If Ndig is the character string 'reset', then the default
%        is reset to the original default of [64 2].
%
%        Default: [64 4]
%
%        The number of digits reported must always be at least 1,
%        while the number of shadow digits may be as low as zero.
%
% permanence - (OPTIONAL) character string - defines whether the
%        default you supply here will be remembered the next time
%        you start up matlab, or after a "clear functions" command.
%
%        permanence == 'permanent' --> Always remember this as
%             the number of digits to be carried, even after
%             matlab is restarted, or your computer is rebooted.
%             (Only deleting/overwriting the associated preference
%             will affect this property.)
%
%        permanence == 'session' --> the default number of
%             digits will apply until matlab is restarted,
%             or until a "clear functions" command is executed.
%             At that point, the permanent DefaultNumberOfDigits
%             will then apply.
%
%        Default value: 'permanent'
%
%        permanence may be shortened down as far as single
%        characters. Thus even 's' or 'p' may be supplied.
%        Capitalization is ignored.
% 
% Arguments: (output)
%  NDig - The current number of digits stored by default in a
%        newly constructed HPF number, along with the number of
%        shadow digits retained but not displayed.
%
% Examples:
% % Set the default number of digits to be used in an HPF number
% % to 25, with the default value of 2 shadow digits carried.
% DefaultNumberOfDigits(25)
% F = hpf('pi')
% F =
%     3.141592653589793238462643
%
% % In command form, set the default number of digits to be used,
% % to 50, with a default value of 3 shadow digits carried.
% DefaultNumberOfDigits 50 3
%
% % Interrogate the current number of default number of digits
% % to be used in all new HPF numbers.
% DefaultNumberOfDigits
% ans =
%     50 3
%
% % Set the number of digits to be used in an HPF number,
% % disregarding the default value for that property as previously
% % specified. This will not change the default number of digits
% % for future hpf numbers.
% F = hpf('pi',137)
% F =
%     3.141592653589793238462643383279502884197169399375105820974944592307
% 8164062862089986280348253421170679821480865132823066470938446095505822  
%
% % Set the number of digits for future hpf numbers, but only during
% % this current MATLAB session. The next time MATLAB is restarted,
% % the permanently stored value will be used.
% DefaultNumberOfDigits([64 3],'session')
%
%
% See also: hpf
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com

% Is this the first time executed after a restart or a clear functions?
% by setting the value as a persistent variable, this is much faster
% than repeated calls to getpref, which must do a disk read.
persistent NumberOfDigits
if isempty(NumberOfDigits)
  % We can get the Default from a preference, IF that pref has been set.
  if ispref('hpf','DefaultNumberOfDigits')
    % yes, it has been set already
    NumberOfDigits = getpref('hpf','DefaultNumberOfDigits');
  else
    % there is no preference that has been set so far, so choose
    % the root default value of [64 4].
    NumberOfDigits = [64 4];
    % Now we need to add this as a pref
    addpref('hpf','DefaultNumberOfDigits',NumberOfDigits)
  end
end

% The default permanence is 'permanent';
permanence = 'permanent';

% in what form was the call in?
if nargin < 1
  % If no arguments, then it must be a request for the current
  % number of digits to be used.
  NDig = NumberOfDigits;
  
  % nothing was set, so there is no need to change a preference
  permanence = 'session';
  
elseif nargin == 1
  % it MUST be the number of digits has been supplied, or the word 'reset'
  
  % since no permanence was specified, then the change will
  % also go into the prefs
  permanence = 'permanent';
  
  % Was this in character form?
  if ~isnumeric(NDig)
    % a character form was supplied
    
    % was it 'reset'?
    if strcmpi('reset',NDig)
      NDig = [64 4];
    else
      NDig = str2num(NDig); %#ok
      if ~ismember(numel(NDig),[1 2])
        error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
          'Number of digits (and shadow digits) must be positive integers')
      end
    end
  end
  
  % a numeric value was supplied, or character, and was
  % then converted to numeric
  
  % was it a scalar?
  if numel(NDig) == 1
    if (NDig == round(NDig)) && (NDig > 0)
      % re-use the previous number of shadow digits
      NDig = [NDig , NumberOfDigits(2)];
    else
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of digits must be positive integer')
    end
  elseif numel(NDig) == 2
    if all(NDig == round(NDig)) && all(NDig >= 0) && (NDig(1) > 0)
      % NDig is ok as it is
    else
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of digits (and shadow digits) must be positive integers')
    end
  else
    error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
      'Number of digits (and shadow digits) must be positive integers')
  end
elseif nargin > 3
  % exclude this as an error
  error('DEFAULTNUMBEROFDIGITS:invaliddigitspec','No more than 3 arguments are allowed')
elseif nargin == 3
  SDig = varargin{1};
  permanence = varargin{2};
  % if 3 arguments, then the last argument will be 'permanent' or 'session'
  % or a simple shortening thereof.
  valid = {'permanent' 'session'};
  k = find(strncmpi(permanence,valid,length(permanence)));
  if isempty(k)
    error('DEFAULTNUMBEROFDIGITS:invalidpermanence', ...
      'Permanence must be either ''permanent'' or ''session''')
  else
    permanence = valid{k};
  end
  
  % Was NDig in character form?
  if ~isnumeric(NDig)
    % a character form was supplied
    NDig = str2num(NDig); %#ok
    if (numel(NDig) ~= 1) || (NDig <= 0) || (NDig ~= round(NDig))
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of digits must be a positive integer')
    end
  end
  
  % SDig must be the number of shadow digits carried.
  % Was SDig in character form?
  if ~isnumeric(SDig)
    % a character form was supplied
    SDig = str2num(SDig); %#ok
    if (numel(SDig) ~= 1) || (SDig < 0) || (SDig ~= round(SDig))
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of shadow digits must be a non-negative integer')
    end
  end
  
  % combine the number of and shadow digits
  NDig = [NDig,SDig];
  
elseif nargin == 2
  % two arguments supplied. Was the second argument permanence, or
  % was it the number of shadow digits?
  arg2 = varargin{1};
  
  % the default for shadow digits is 2
  SDig = 2;
  if ischar(arg2)
    % Was it 'permanent' or 'session' or a simple shortening thereof?
    valid = {'permanent' 'session'};
    k = find(strncmpi(arg2,valid,length(permanence)));
    if ~isempty(k)
      permanence = valid{k};
      SDig = [];
    else
      % was it a scalar numeric form, but as characters because of
      % execution as a command?
      SDig = str2num(arg2); %#ok
      if (numel(SDig) ~= 1) || (SDig < 0) || (SDig ~= round(SDig))
        error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
          'Number of shadow digits must be a non-negative integer')
      end
    end
  else
    % the second argument must have been numeric, so the
    % number of shadow digits
    SDig = arg2;
  end
  
  % Was NDig in character form?
  if ~isnumeric(NDig)
    % a character form was supplied
    NDig = str2num(NDig); %#ok
    if (numel(NDig) == 2) && ~isempty(SDig)
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Only a number of reported and shadow digits may be provided')
    end
    
    if (NDig(1) <= 0) || any(NDig ~= round(NDig))
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of digits must be positive integer')
    end
  end
  
  % SDig must be the number of shadow digits carried.
  % Was SDig in character form?
  if ~isnumeric(SDig)
    % a character form was supplied
    SDig = str2num(SDig); %#ok
    if (numel(SDig) ~= 1) || (SDig < 0) || (SDig ~= round(SDig))
      error('DEFAULTNUMBEROFDIGITS:invaliddigitspec', ...
        'Number of shadow digits must be a non-negative integer')
    end
  end
  
  % combine the number of and shadow digits
  NDig = [NDig,SDig];
  if numel(NDig) == 1
    NDig = [NDig,2];
  end
  
end

% was this a permanent change, to be set into the prefs?
NumberOfDigits = NDig;
if strncmpi(permanence,'permanent',length(permanence))
  % and stuff the new value into the prefs for hpf
  setpref('hpf','DefaultNumberOfDigits',NumberOfDigits)
end

% do we need to return anything?
if (nargout < 1) && (nargin > 0)
  clear NDig
end


