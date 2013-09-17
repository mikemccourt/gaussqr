function [Mantissa,Exponent] = powerssum2dec(pow2list,NDig)
% Convert a list of binary exponents into a decimal Mantissa + Exponent 
% usage: [Mantissa,Exponent] = powerssum2dec(pow2list,NDig)
% 
% Computes the number sum(2.^pow2list)
% 
% Arguments: (input)
%  pow2list - (nonempty) vector of integers between
%       -1023 and 1024
%
%  NDig - Scalar integer - defines the number of digits
%       in the final Mantissa to be returned
%
% Arguments: (output)
%  Mantissa - vector of double integers (0 - 9) defining
%       the decimal mantissa of the sum
%
%  Exponents - (Double) Integer that defines the power
%       of 10 of the sum.

% we need to bring in this array only once
persistent powersof2  %#ok
if isempty(powersof2)
  load _powersof2_
end

% set up Mantissa to accumulate into
Mantissa = [];

% get the set of non-negative powers from the list
pos = sort(pow2list(pow2list >= 0),'descend');

if ~isempty(pos)
  Mantissa = double(powersof2.digits{1024 + pos(1)});
  nmant = numel(Mantissa);
  for i = 2:numel(pos)
    nextpow = double(powersof2.digits{1024 + pos(i)});
    nextpow = padz(nextpow,nmant,'left');
    Mantissa = Mantissa + nextpow;
  end
  
  % do we need to perform any carries?
  if numel(pos) > 1
    % we may need to do so
    carryindex = find(Mantissa > 9);
    while ~isempty(carryindex)
      % if we must do a carry past the highest
      % order digit in the loop, then add a zero
      % digit on the left for that carry to flow
      % into
      if carryindex(1) == 1
        carryindex = carryindex + 1;
        Mantissa = [0,Mantissa]; %#ok
      end
      remainder = rem(Mantissa(carryindex),10);
      carry = (Mantissa(carryindex) - remainder)/10;
      Mantissa(carryindex) = remainder;
      carryindex = carryindex - 1;
      % do the carry itself
      Mantissa(carryindex) = Mantissa(carryindex) + carry;
      
      % Of the elements carried into, which ones are
      % still a problem?
      carryindex(Mantissa(carryindex) < 10) = [];
    end
    
  end
  
  % The final exponent on the result is given by the
  % number of decimal digits in Mantissa (less 1)
  Exponent = numel(Mantissa) - 1;
end

% is there a fractional part here?
fraction = [];

% get the set of negative powers from the list
neg = sort(pow2list(pow2list < 0),'ascend');
if ~isempty(neg)
  % there is at least ONE negative power of 2 in that list
  % start with the most negative one
  fdigits = double(powersof2.digits{1024 + neg(1)});
  fraction = padz(fdigits,numel(fdigits) + powersof2.leadingzeros(1024 + neg(1)),'left');
  
  nfract = numel(fraction);
  for i = 2:numel(neg)
    fdigits = double(powersof2.digits{1024 + neg(i)});
    fdigits = padz(fdigits,numel(fdigits) + powersof2.leadingzeros(1024 + neg(i)),'left');
    fdigits = padz(fdigits,nfract,'right');
    fraction = fraction + fdigits;
  end
  
  % do we need to perform any carries on the fractional part?
  % note that we KNOW this sum can never exceed 1, so we will
  % never need to add a new digit on the left end
  if numel(neg) > 1
    % we possibly need to do a carry
    carryindex = find(fraction > 9);
    while ~isempty(carryindex)
      remainder = rem(fraction(carryindex),10);
      carry = (fraction(carryindex) - remainder)/10;
      fraction(carryindex) = remainder;
      carryindex = carryindex - 1;
      % do the carry itself
      fraction(carryindex) = fraction(carryindex) + carry;
      
      % Of the elements carried into, which ones are
      % still a problem?
      carryindex(fraction(carryindex) < 10) = [];
    end
    
  end
  
end

% we need to combine the integer and fractional parts together
% at least one of these must be non-empty if pow2list was non-empty
%
% if isempty(fraction)
%   all the powers would have been positive, so both Mantissa and Exponent
%   are correct, so this would be a no-op case.
if isempty(Mantissa)
  % all powers were negative, so we must take Exponent from the
  % number of leading zeros in the result.
  k = find(fraction ~= 0,1,'first');
  Exponent = -k;
  Mantissa = fraction(k:end);
elseif ~isempty(Mantissa) && ~isempty(fraction)
  % mixed integer and fractional part
  % Exponent has already been set correctly.
  Mantissa = [Mantissa,fraction];
end

% finally, do we need to truncate the digits, or must we
% pad zeros on the end for the target number of digits?
if numel(Mantissa) ~= NDig
  Mantissa = padz(Mantissa,NDig,'right');
end






