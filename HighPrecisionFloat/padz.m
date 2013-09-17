function vec = padz(vec,finallength,paddirection)
% pads a vector to have a desired length, by appending zeros at either end
% usage: vec = padz(vec,finallength,paddirection)
%
% If the initial vector is longer than the goal, then digits from the
% right end will be truncated. Of course, if the two lengths are the
% same, this is a no-op
%
% arguments: (input)
%  vec - initial vector
%
%  finallength - integer scalar, denoting the desired length of the vector
%
%  paddirection - character flag
%        'right' --> pad on the right end
%        'left'  --> pad on the left end
%
%        Note that if finallength is less than numel(vec), 
%        then the truncation will always be taken from the
%        right end of the digit string, irregardless of the
%        choice of paddirection.

% is this a truncation?
nvec = numel(vec);
if nvec > finallength
  % a truncation operation
  vec = vec(1:finallength);
elseif nvec < finallength
  if strcmpi(paddirection,'left')
    vec = [zeros(1,finallength - nvec),vec];
  else
    vec(finallength) = 0;
  end
end


