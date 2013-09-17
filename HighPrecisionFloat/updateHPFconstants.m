function updateHPFconstants(maxdigits,constantname)
% updateHPFconstants: generates new values for the HPF stored constants: e, pi, and the natural log of 10
% usage: updateHPFconstants(maxdigits,constantname)
%
% When I wrote the HPF tools, I generated those important constants to
% what seemed like a sufficient number of digits for any reasonable
% computation. However, Moore's law being what it is, perhaps 100,000
% digits or so will be too small next year. This tool allows matlab to
% update those constants programmatically, then saving them out for future
% use.
%
% Recognize that updateHPFconstants will be a cpu hog, taking a serious
% amount of time for large numbers of digits. This is NOT a tool that you
% will use on any daily basis, since those constants are already stored
% internally to plenty of digits for any sane, rational person. In fact,
% expect the running time for this code to be in the hours, or even days
% to update these constants beyond a million digits or so.
%
% arguments: (input)
%  maxdigits - scalar integer - defines a new maximum number of stored
%      digits for the three constants usd internally, e, pi, and the
%      natural log of 10.
%
%      if maxdigits is less than the currently stored number of digits
%      for those constants, then nothing is done.
%
%  constantname - string containing the name of the constant to be updated
%      in value. If the new value has more digits in it for this constant
%      than are currently stored by HPF, then an update is done to the
%      file on disk for that constant.
%
%      constant name may be a string containing one constant name, or a
%      cell array of strings.
%
%      Currently, the list of constants in HPF is {'pi', 'e', 'ln10'}.

% get the path to _special_numbers_.mat
snpath = which('_special_numbers_.mat');

% bring in the current values in _special_numbers_.mat
load(snpath)

% 50 spare digits should be entirely adequate
sparedigits = 50;
NDig = [maxdigits,sparedigits];

if (nargin < 2) || isempty(constantname)
  constantname = {'pi' 'e', 'ln10'};
else
  constantname = cellstr(constantname);
  if any(~ismember(constantname,{'pi', 'e', 'ln10'}))
    error('Legal constants that may be updated are: ''pi'', ''e'', ''ln10''')
  end
end

% do we need to update e?
if ismember('e',constantname) && (maxdigits > numel(edigits))
  % we need to update e. do so by computing exp(2^-20), then raise that to
  % the power 1048576. I'll keep plenty of spare digits on the side, so there
  % is no fear of a loss of precision here.
  z = hpf('0.00000095367431640625',NDig); % 1/1048576 = 2^-20
  
  % How many terms in the Taylor series do we need to compute?
  % Since the exponential series has as its last term
  % z^m/factorial(m), we can quit when the natural log of that
  % expression is less than the log of our exponential goal.
  % Since log(exp(z)) is never less than -0.5, we have a simple
  % cutoff, based on the number of significant digits in our
  % hpf number.
  %
  %  -0.5 + log(10^-NDig) = log(z^m/m!) = m log(z) - log(m!)
  %
  % the logs here are natural logs of course. Using Stirling's
  % approximation for m!, we get
  %
  %  -0.5 - NDig*log(10) = m*log(z) - [1/2*log(2*pi) + 1/2*log(m) + m*log(m) - m]
  %
  % Combining terms, this reduces to
  %
  %  1/2*(log(2*pi)-1) - NDig*log(10) = m*(log(z) + 1) - log(m)*(m + 1/2)
  %
  % Bump up NDig by 1, just to ensure convergence for the full digits.
  % fun is a decreasing function. The zero crossing defines the number
  % of the last term we need to compute.
  fun = @(m) m.*(log(1/1048576) + 1) - log(m).*(m + 1/2) + ...
    (sum(NDig)+1)*log(10) - 1/2*(log(2*pi) - 1);
  
  mterms = termsbisector(fun);
  
  % run the Taylor series in reverse now. This is why I wanted to know how
  % many terms would be necessary. By running the loop backwards, I
  % avoid divisions, and a divide is more expensive than a multiply for
  % an hpf number.
  E1 = hpf('1',NDig);
  Fact = E1;
  H = waitbar(0,['Computing exp(1) to ',num2str(maxdigits),' digits']);
  for m = mterms:-1:1
    % update a waitbar
    if mod(m,100) == 0
      waitbar((mterms - m)/mterms,H)
    end
    
    Fact = Fact*m;
    E1 = z*E1 + Fact;
  end
  delete(H)
  
  % because we ran the loop backwards, in the end we need to
  % divide by factorial(mterms)
  E1 = E1./Fact;
  
  % raise E1 to the 1048576'th power, to get exp(1). Squaring E1 20 times
  % will suffice. Along the way, compute a good approximation to
  % exp(log(10)) in binary fractional powers of e
  E1 = E1.*E1;  % 2^-19'th power
  E1 = E1.*E1;  approx10 = E1; % 18
  E1 = E1.*E1;  % 17
  E1 = E1.*E1;  % 16
  E1 = E1.*E1;  approx10 = approx10.*E1; % 15
  E1 = E1.*E1;  approx10 = approx10.*E1; % 14
  E1 = E1.*E1;  % 13
  E1 = E1.*E1;  approx10 = approx10.*E1; % 12
  E1 = E1.*E1;  approx10 = approx10.*E1; % 11
  E1 = E1.*E1;  approx10 = approx10.*E1; % 10
  E1 = E1.*E1;  % 9
  E1 = E1.*E1;  approx10 = approx10.*E1; % 8
  E1 = E1.*E1;  % 7
  E1 = E1.*E1;  approx10 = approx10.*E1; % 6
  E1 = E1.*E1;  approx10 = approx10.*E1; % 5
  E1 = E1.*E1;  % 4
  E1 = E1.*E1;  % 3
  E1 = E1.*E1;  approx10 = approx10.*E1; % 2
  E1 = E1.*E1;  % 1
  E1 = E1.*E1;  approx10 = approx10.*E1.*E1;% 0
  
  % extract the digits as a uint8 vector to be saved back into
  % the special numbers mat file.
  E1dig = mantissa(E1);
  edigits = uint8(E1dig(1:maxdigits));
  
end

% update pi
if ismember('pi',constantname) && (maxdigits > numel(pidigits))
  NDig = [maxdigits,50];
  H = waitbar(0,['Computing pi to ',num2str(maxdigits),' digits (quadratic convergence)']);
  
  a0 = hpf('1',NDig);
  b0 = sqrt(hpf('0.5',NDig));
  t0 = hpf('0.25',NDig);
  p0 = hpf('1',NDig);
  
  notdone = true;
  while notdone
    a1 = (a0 + b0)/2;
    b1 = sqrt(a0.*b0);
    t1 = t0 - p0.*(a0 - a1).^2;
    p1 = 2 .*p0;
    
    % has our estimate converged?
    k = find(t0.Migits ~= t1.Migits,1,'first');
    if isempty(k)
      % t has converged to all of its digits, so we can generate
      % the final estimate for pi, and terminate the loop
      piest = (a1 + b1).^2./(4 .*t1);
      
      notdone = false;
    else
      % update the waitbar...
      waitbar(k/maxdigits,H)
      
      a0 = a1;
      b0 = b1;
      t0 = t1;
      p0 = p1;
    end
  end
  % zap the waitbar
  delete(H)
  piestdig = mantissa(piest);
  pidigits = uint8(piestdig(1:maxdigits));
  
end

% update log(10)
if ismember('ln10',constantname) && (maxdigits > numel(ln10digits))
  % use the current value for log(10) as an offset to make the
  % natural log series converge like a bat out of a very warm place.
  ln10 = hpf('ln10',[numel(ln10digits),0]);
  NDig = [maxdigits,50];
  ln10 = augmentdigits(ln10,NDig);
  
  % first compute exp(ln10) for the old approximation, to the
  % desired number of digits, so it will not be exactly 10.
  % I can't use exp to do so directly, since it uses the old
  % value of ln10 inside and exp uses ln10.
  z = ln10 - 2 - 1/4 - 1/32 - 1/64 - 1/256 - 1/1024 - 1/2048 - 1/4096 - 1/16384 - 1/32768 - 1/262144;
  
  fun = @(m) m.*(log(double(z)) + 1) - log(m).*(m + 1/2) + ...
    (sum(NDig)+1)*log(10) - 1/2*(log(2*pi) - 1);
  
  mterms = termsbisector(fun);
  
  % again, run the exponential Taylor series in reverse.
  Ez = hpf('1',NDig);
  Fact = Ez;
  H = waitbar(0,['Computing log(10) to ',num2str(maxdigits),' digits']);
  for m = mterms:-1:1
    % update a waitbar
    if mod(m,100) == 0
      waitbar((mterms - m)/mterms,H)
    end
    
    Fact = Fact*m;
    Ez = z*Ez + Fact;
  end
  delete(H)
  
  % because we ran the loop backwards, in the end we need to
  % divide by factorial(mterms)
  Ez = Ez./Fact;
  
  % multiply by our approximation to exp(log(10)) from before
  Ez = Ez.*approx10;
  
  % Ez is now exp(ln10). We can use this to bootstrap our computation
  % log(10) to the desired number of digits. Thus divide 10 by Ez. That
  % will be a number fairly close to 1, but slightly above or below 1.
  z = reciprocal(Ez).*10;
  
  % transform z. this transformation creates y VERY near zero, since z
  % was incredibly close to 1.
  y = (z - 1)./(z + 1);
  
  % we can find the log(z) with a series in y, that will take vanishingly
  % few terms, since we already had a very good approximation for log(x).
  logz = 2.*y;
  ysq = y.*y;
  yterm = logz;
  notdone = true;
  n = 1;
  while notdone
    n = n + 2;
    yterm = yterm.*ysq;
    logz = logz + yterm./n;
    if yterm.Exponent < -sum(NDig)
      notdone = false;
    end
  end
  
  % we have found an approximation for ln10 that is accurate to the
  % requested number of digits
  ln10 = ln10 + logz;
  
  % extract those digits as uint8
  ln10dig = mantissa(ln10);
  ln10digits = uint8(ln10dig(1:maxdigits));
  
end

% save out the new file to be used for future computations
save(snpath,'edigits','pidigits','ln10digits')

% we need to cause the persistent variables in hpf to be reinitialized
clear hpf

% =============================================================
%      end mainline hpf, begin subfunctions
% =============================================================

function mterms = termsbisector(fun)
% mterms is the last term in the Taylor series that we will
% need to get a good approximation from that series. Always take
% at least 5 terms in the series.

if fun(4) <= 0
  % always take at least 5 terms in the series
  mterms = 5;
else
  % bracket the root
  ab = [4 8];
  fab = fun(ab);
  while fab(2) > 0
    ab = [ab(2),2*ab(2)];
    fab = [fab(2),fun(ab(2))];
  end
  % simple bisection step
  while diff(ab) > 1
    c = mean(ab);
    fc = fun(c);
    if fc <= 0
      ab = [ab(1),c];
      fab = [fab(1),fc];
    else
      ab = [c,ab(2)];
      fab = [fc,fab(2)];
    end
  end
  
  % pick the upper end point.
  mterms = ab(2);
end % function termsbisector

