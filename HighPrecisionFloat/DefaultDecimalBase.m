function DBase = DefaultDecimalBase(DBaseNew)
% Allows the user to set the default value for the decimal base for a new HPF object.
% Usage: DefaultDecimalBase(DBaseNew); % Resets the default to a new value
% Usage: DBase = DefaultDecimalBase; % Returns the current default
%
% DefaultDecimalBase can be used in either command or functional mode.
%
% DefaultDecimalBase does NOT change the base of the mantissa
% (migits) stored in any existing hpf numbers, only in those numbers
% created in the future. This change will be permanent, remembered
% even after matlab is restarted, or your computer is rebooted.
% (Only deleting/overwriting the associated preference will affect
% this property.)
%
% Arguments: (input)
%  DBaseNew - scalar positive integer, from the set {1,2,3,4,5,6}
%        or scalar character from the set {'1','2','3','4','5','6'}.
%        The word 'reset' is also accepted, in which case the
%        installation default is to use 4-migits.
%
%        This controls the number of decimal digits carried in
%        each migit of the mantissa.
%
%             DBaseNew = 1 --> Base 10 digits
%             DBaseNew = 2 --> Base 100 migits
%             DBaseNew = 3 --> Base 1000 migits
%             DBaseNew = 4 --> Base 10000 migits
%             DBaseNew = 5 --> Base 100000 migits
%             DBaseNew = 6 --> Base 1000000 migits
%
%        Different values of DBase will affect various properties
%        of HPF numbers. Thus...
%
%          Memory requirements are inversely proportional to DBase.
%
%          Addition (or subtraction) speed is proportional to DBase
%
%          Multiplication (or division) speed is proportional to DBase^2
%
%          The maximum theoretical number of decimal digits allowed
%          is also related to DBase. Thus, as long as you know you
%          will never wish to do arithmetic on numbers with more than
%          36000 decimal digits, then HPF will be most efficient
%          with a DefaultDecimalBase of 6.
%
%             DBase = 1 --> 3.6e14 decimal digits
%             DBase = 2 --> 3.6e12 decimal digits
%             DBase = 3 --> 3.6e10 decimal digits
%             DBase = 4 --> 3.6e8 decimal digits
%             DBase = 5 --> 3.6e6 decimal digits
%             DBase = 6 --> 36000 decimal digits
%    
%        You should see that unless you have a huge amount of RAM
%        as well as an incredibly fast processor, even DBase of 4
%        will allow arithmetic on numbers with 360,000,000 decimal
%        digits. While I've allowed the user to change that value
%        to smaller values of DBase in case someone actually wanted
%        to do arithmetic with billions of decimal digits, I expect
%        that computers will need great increases in power before
%        anyone ever seriously wants this capability.
%
%        If you are looking for yet more speed, then moving to a
%        higher decimal base will give that. Thus, for large numbers,
%        a decimal base of 6 will roughly double the speed of a
%        multiplication.
%
%        Default: 4
%
% Arguments: (output)
%  DBase - The current default value for the decimal base to be used
%        whenever new HPF numbers are created.
%
%
% Examples:
% % Set the Default decimal base for future HPF numbers as a command
% DefaultDecimalBase 5
%
% % In function form, set the default decimal base
% DefaultNumberOfDigits(4)
%
% % Extract the current default decimal base that is used for HPF numbers
% DBase = DefaultDecimalBase
% DBase =
%    4
%
% See also: hpf, DefaultNumberOfDigits
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com

% Is this the first time executed after a restart or a clear functions?
% by setting the value as a persistent variable, this is much faster
% than repeated calls to getpref, which must do a disk read.
persistent CurrentDecimalBase
if isempty(CurrentDecimalBase)
  % We can get the default from a preference, IF that pref has been set.
  if ispref('hpf','DefaultDecimalBase')
    % yes, it has been set already
    CurrentDecimalBase = getpref('hpf','DefaultDecimalBase');
  else
    % there is no preference that has been set so far, so choose
    % the root default value of 4.
    CurrentDecimalBase = 4;
    % Now we need to add this as a pref
    addpref('hpf','DefaultDecimalBase',CurrentDecimalBase)
  end
end

% in what form was the call made?

% If no arguments, then it must be a request for the current
% number of digits to be used. DBase already exists as the
% current value we have stored, so we are done in that case.
if nargin == 1
  % it MUST be a new value for the decimal base
  
  % Was this in character form or a number?
  if isnumeric(DBaseNew)
    % a number was supplied.
    if ismember(DBaseNew,[1 2 3 4 5 6])
      CurrentDecimalBase = DBaseNew;
    else
      error('DEFAULTDECIMALBASE:invaliddecimalbase', ...
        'Decimal base must be one of [1 2 3 4 5 6]')
    end
  elseif ischar(DBaseNew)
    % a character form was supplied. A switch will do it easily enough,
    % doing the necessary validity checks for me.
    switch lower(DBaseNew)
      case '1'
        CurrentDecimalBase = 1;
      case '2'
        CurrentDecimalBase = 2;
      case '3'
        CurrentDecimalBase = 3;
      case '4'
        CurrentDecimalBase = 4;
      case '5'
        CurrentDecimalBase = 5;
      case '6'
        CurrentDecimalBase = 6;
      case {'reset' 'rese' 'res' 're' 'r'}
        % resets to the initial default of 4
        CurrentDecimalBase = 4;
      otherwise
        error('DEFAULTDECIMALBASE:invaliddigitspec', ...
          'Decimal Base must be an integer in the set [1,2,3,4,5,6], or the word ''reset''')
    end
  end
  
  % and stuff the new value into the prefs for hpf
  setpref('hpf','DefaultDecimalBase',CurrentDecimalBase)
elseif nargin > 1
  error('DEFAULTDECIMALBASE:arguments', ...
    'Only one or zero arguments are allowed')
end

% do we need to return anything?
if (nargout > 0) || (nargin == 0)
  DBase = CurrentDecimalBase;
end


