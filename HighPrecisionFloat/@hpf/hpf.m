classdef  (InferiorClasses = {?vpi}) hpf
  % HPF - High Precision Floating point format
  %
  % HPF is a full long precision form, allowing the user to specify any
  % fixed number of significant (decimal) digits,
  %
  % HPF is not a variable precision format, in the sense that the number
  % of digits carried always remains fixed at the level set by the user.
  % Nor is HPF a fixed point format, where the number of digits to the
  % right of the decimal point stay fixed. HPF is a FLOATING point
  % arithmetic where the number of digits is arbitrary, and can be
  % specified by the user.
  %
  % HPF stores the digits of your number in a base 10^k format, where
  % each element encodes k decimal digits of the number. Thus the total
  % number of true decimal digits stored must always be a multiple of k.
  % Since the user can control the value of k as a default parameter,
  % you can control the speed of computation. Larger values of k here
  % will scale the speed of multiplies and divides by k^2, while
  % additions and subtractions will be scaled in speed by a factor of k.
  % Furthermore, the memory requirements for an HPF number will scale
  % in a negative linear relationship with the value of k,
  %
  % The nominal default value of k is chosen to be 4, so all internal
  % arithmetic is effectively carried out in base 10000. Note that the
  % base also effectively implies a limit on the total number of digits
  % that can be stored, in the sense that overflows in multiplication
  % (using convolution) can become a limiting factor. More information
  % on the impact of the exponent k can be found in the HPF.pdf file,
  % as well as in the help for DefaultDecimalBase.m and in the demo
  % files.
  %
  % For example, the number 123.4567 will be stored (given a decimal
  % base of 4) with migits of [1234 5670], a sign of +1, and an
  % exponent of 3. You can think of it in scientific notation as
  % 0.1234567 x 10^3, or as 0.1234567e3.
  %
  % In general, you should never completely trust the least significant
  % digits in a floating point number, as there will always be round-off
  % errors. So use more digits than you absolutely need. In fact, this
  % capability is built in to HPF in the use of shadow digits, which
  % are provided to absorb the precision loss in computations. To be
  % conservative, one should always allow some shadow digits to avoid
  % round-off issues. Some computations may require a large number of
  % shadow digits. A good understanding of floating point arithmetic will
  % always help you in your work.
  %
  %    Total number of digits = displayed digits + shadow digits
  %
  % Finally, HPF is a tool written for fun by me. You will always be
  % better off using good practices of numerical analysis to improve
  % your results rather than resorting to the lazy solution of increased
  % precision. I do foresee a few interesting uses for HPF:
  %
  %  - You may enjoy seeing how I chose to tease dozens or even many
  %    thousands of digits of precision for any given computation.
  %
  %  - There are a few innovative ideas in HPF to apply to your
  %    own tools, for example, the DefaultNumberOfDigits tool.
  %
  %  - An example of use for classdef.
  %
  %  - There will be a few (rare) people who will find a way to use
  %    this tool for some useful purpose.
  %
  %  - For fun, perhaps solving the (classic) railroad rail problem
  %    using brute force, or computing a few million digits of pi.
  %
  %  - A teaching tool to learn about floating point arithemetic
  %    in general.
  %
  %  - A floating point simulation tool to work in a non-standard
  %    number of digits of precision.
  %
  % For more information about the HPF format, read the HPF.pdf file,
  % which provides many details about the code. You should also read the
  % published hpf_demo.html file for many examples of use.
  %
  % Author: John D'Errico
  % e-mail: woodchips@rochester.rr.com
  
  properties (SetAccess = private)
    % NumberOfDigits sets the number of digits carried within the
    % hpf representation for this number. This parameter has a default
    % value set by the function DefaultNumberOfDigits. It contains a
    % vector of length 2 - the first element will be the number of decimal
    % digits shown by disp. The second number is the number of shadow
    % (or guard) digits, i.e., additional digits that are used and
    % carried in the computation, yet are not displayed in the result.
    NumberOfDigits
    
    % The number of decimal digits carried in each migit of the Mantissa.
    % This is controlled by the function DefaultDecimalBase, and should
    % generally not be varied within a session.
    DecimalBase
    
    % Base = 10.^DecimalBase. Stored here for efficiency.
    Base
    
    % Even with a huge number of digits available to you, infs and NaNs
    % are still necessary values. Therefore I need to include a flag to
    % indicate these special numbers. Numeric may be any of: 0, inf,
    % inf, NaN. All normal numbers are identified as such by a 0 in
    % this field.
    Numeric
    
    % Sign is +1 for positive numbers, -1 for negative, 0 for 0
    Sign
    
    % The power of 10 to be applied, in a scientific notation. Exponents
    % that exceed +/- 2^53-1 will overflow the exponent, as the Exponent
    % is a double precision number. Those are really vast numbers, but
    % beware that even HPF has limits.
    Exponent
    
    % The actual decimal migits of the number, to be interpreted as
    % having a decimal point immediately preceding the first migit.
    % Additional powers of 10 as indicated by the Exponent field will
    % shift the decimal point as necessary to the right or to the left.
    % Note I have called the elements of the mantissa "migits", as
    % they will generally be 4 digit integers, thereby encoding
    % Multiple decimal digits into each migit. The base of these
    % migits is controlled by the DecimalBase property.
    Migits
    
  end
  
  methods (Static)
    % Static methods for some basic array generation tools
    mat0 = zeros(varargin)
    mat1 = ones(varargin)
    mati = eye(varargin)
    mat10 = ten(varargin)
  end
  
  methods
    function F = abs(F)
      % absolute value for a hpf number
      %
      % abs(F) returns F when F >= 0, -F when F < 0
      %
      %  See also: uminus
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % This works for any case, even NaNs, infs and empty
      for i = 1:numel(F)
        % leave the sign alone for NaNs
        if ~isnan(F(i).Numeric)
          F(i).Sign = abs(F(i).Sign);
        end
      end
    end % F = abs(F)
    
    function F = acos(X)
      % evaluates acos(X) (the inverse cosine function) for a hpf number X
      
      % Use a simple scheme - just a transformation from asin
      F = X;
      for i = 1:numel(X)
        F(i) = hpf('pi',X(i).NumberOfDigits)./2 - asin(X(i));
      end
      
    end % function F = acos(X)
    
    function F = acosd(X)
      % evaluates acos(X) for a hpf number X (returning a result in degrees)
      
      % simple conversion of degrees to radians
      F = X;
      for i = 1:numel(X)
        F(i) = acos(X(i)).*(180./hpf('pi',X(i).NumberOfDigits));
      end
      
    end % function F = acosd(X)
    
    function F = acosh(X)
      % evaluates acosh(X) (the inverse hyperbolic cosine) for an hpf number X
      
      F = X;
      for i = 1:numel(X)
        if isnan(X(i).Numeric) || (X(i).Sign <= 0) || (X(i) < 1)
          F(i) = hpf('NaN',X(i).NumberOfDigits);
        elseif isinf(X(i).Numeric)
          F(i) = hpf('inf',X(i).NumberOfDigits);
        elseif X(i) == 1
          F(i) = hpf('0',X(i).NumberOfDigits);
        else
          % simple inverse
          F(i) = log(X(i) + sqrt(X(i).*X(i) - 1));
        end
      end
      
    end % function F = acosh(X)
    
    function F = acot(X)
      % evaluates acot(X) (the inverse cotangent function) for a hpf number X
      
      % Use a simple scheme - just a transformation from atan
      F = X;
      for i = 1:numel(X)
        if isinf(X(i).Numeric)
          % +/-inf both map to zero
          F(i) = hpf('0',X(i).NumberOfDigits);
        elseif isnan(X(i).Numeric)
          F(i) = hpf('NaN',X(i).NumberOfDigits);
        elseif X(i).Sign == 0
          % X(i) was a zero, so no reason to do any work here
          % actually, acot(0) should be undefined, probably NaN, since
          % acot(-eps) would be -pi/2. But MATLAB has acot(0) as pi/2,
          % so I've chosen to be consistent.
          F(i) = hpf('pi',X(i).NumberOfDigits)./2;
        else
          % any other x
          F(i) = atan(reciprocal(X(i)));
        end
      end
    end % function F = acot(X)
    
    function F = acotd(X)
      % evaluates acotd(X) (the inverse cotangent function, in degrees) for a hpf number X
      
      % Use a simple scheme - just a transformation from atand
      F = X;
      for i = 1:numel(X)
        F(i) = acot(X(i)) .*180 ./hpf('pi',X(i).NumberOfDigits);
      end
      
    end % function F = acotd(X)
    
    function F = acsc(X)
      % evaluates acsc(X) (the inverse cosecant function) for a hpf number X
      
      % Use a simple scheme - just a transformation from asin
      F = X;
      for i = 1:numel(X)
        if isinf(X(i))
          F(i) = hpf('0',X(i).NumberOfDigits);
        else
          F(i) = asin(reciprocal(X(i)));
        end
      end
    end % function F = acsc(X)
    
    function F = acscd(X)
      % evaluates acscd(X) (the inverse cosecant function, in degrees) for a hpf number X
      
      % Use a simple scheme - just a transformation from asind
      F = X;
      for i = 1:numel(X)
        if isinf(X(i))
          F(i) = hpf('0',X.NumberOfDigits);
        else
          F(i) = asind(reciprocal(X(i)));
        end
      end
      
    end % function F = acscd(X)
    
    function F = adjustdecimalbase(F,N)
      % allows the user to modify the decimal base of an HPF number
      % usage: F = augmentdigits(F,N)
      %
      % N must be a scalar, positive integer from the set {1,2,3,4,5,6}
      %
      % this tool would rarely be used, only in the rare case between
      % a dyadic operation with two HPF numbers having different decimal
      % bases.
      %
      % Note that when the decimal base changes, there will often be
      % the need to append shadow digits, because the total number of
      % digits must always be a multiple of the base. Existing significant
      % digits will not be truncated in this adjustment.
      if (nargin ~= 2)
        error('HPF:adjustdecimalbase','Must provide two arguments')
      elseif (numel(N) ~= 1) || ~ismember(N,[1 2 3 4 5 6]) 
        error('HPF:adjustdecimalbase', ...
          'N must be scalar integer from the set {1,2,3,4,5,6}')
      end
      
      % IS F an array?
      if numel(F) > 1
        for i = 1:numel(F)
          F(i) = adjustdecimalbase(F(i),N);
        end
        return
      end
      
      % its a no-op if the two bases are the same
      if N ~= F.DecimalBase
        D = m2d(F.Migits,F.DecimalBase);
        
        NDigLimit = [3.6e14 3.6e12 3.6e10 3.6e8 3.6e6 36000];
        if NDigLimit(N) <= sum(F.NumberOfDigits)
          switch N
            case 1
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 3.6e14 total digits for a DecimalBase of 1')
            case 2
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 3.6e12 total digits for a DecimalBase of 2')
            case 3
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 3.6e10 total digits for a DecimalBase of 3')
            case 4
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 3.6e8 total digits for a DecimalBase of 4')
            case 5
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 3.6e6 total digits for a DecimalBase of 5')
            case 6
              error('HPF:adjustdecimalbase', ...
                'Cannot exceed 36000 total digits for a DecimalBase of 6')
          end
        end
        
        % do we need to add a few shadow digits?
        r = mod(sum(F.NumberOfDigits),N);
        if r ~= 0
          F.NumberOfDigits(2) = F.NumberOfDigits(2) + (N - r);
        end
        
        % expand D, right paddding with zeros
        if numel(D) ~= sum(F.NumberOfDigits)
          D = padz(D,sum(F.NumberOfDigits),'right');
        end
        
        % fix F up with the new base
        F.DecimalBase = N;
        F.Base = 10.^N;
        F.Migits = d2m(D,F.DecimalBase);
        
      end
      
    end % function adjustdecimalbase
      
    function F = asec(X)
      % evaluates asec(X) (the inverse secant function) for a hpf number X
      
      % Use a simple scheme - just a transformation from acos
      F = X;
      for i = 1:numel(X)
        if isinf(X(i))
          F(i) = hpf('pi',X(i).NumberOfDigits)./2;
        else
          F(i) = acos(reciprocal(X(i)));
        end
      end
      
    end % function F = asec(X)
    
    function F = asecd(X)
      % evaluates asecd(X) (the inverse secant function, in degrees) for a hpf number X
      
      % Use a simple scheme - just a transformation from acosd
      F = X;
      for i = 1:numel(X)
        if isinf(X(i))
          F(i) = hpf('90',X.NumberOfDigits);
        else
          F(i) = acos(reciprocal(X(i)));
        end
      end
      
    end % function F = asecd(X)
    
    function F = asin(X)
      % evaluates asin(X) (the inverse sine function) for a hpf number X
      %
      % The method used is the standard trigonometric series for asin,
      % but with a twist. If we have a good estimate of asin(X), then
      % can transform the problem using a simple trig identity.
      %
      % sin(u + v) = sin(u)*cos(v) + cos(u)*sin(v)
      %
      % So if v is our approximate estimate for asin(X), then we wish to
      % compute the value of u, such that sin(u+v) = sin(X). Along the way,
      % we will need to do a little algebraic manipulation on that identity
      % above. Starting from the point
      %
      %   X = sin(u)*cos(v) + cos(u)*sin(v)
      %
      % we transform that into
      %
      %   sin(u) = X*cos(v) - sin(v)*sqrt(1-X.*X)
      %
      % v is known, and is typically within 1 part in 10^16 of the sum (u+v)
      % where it is u+v that we wish to compute.
      %
      % Once we have done the above transformation, we use a series for
      % asin that will yield roughly an additional 32 digits of accuracy
      % per term, just the standard asin Taylor series, but for very asmall
      % numbers, it will be efficient, and go like blazes.
      
      % is X a scalar?
      if isempty(X)
        F = X;
        return
      elseif numel(X) > 1
        % vector or array
        F = X;
        for i = 1:numel(X)
          F(i) = asin(X(i));
        end
        return
      end
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if ~isfinite(X.Numeric)
        % asin(NaN) = asin(inf) = asin(-inf) = NaN
        F = hpf('NaN',NDig);
        return
      elseif X.Sign == 0
        % asin(0) == 0
        F = hpf('0',NDig);
        return
      else
        % compare abs(X) to 1. If it is greater than 1, then we
        % return NaN. If they are equal, then we return +/- pi/2
        Xminus1 = abs(X) - 1;
        if Xminus1 > 0
          F = hpf('NaN',NDig);
          return
        elseif Xminus1 == 0
          % asin(1) = pi/2, Asin(-1) = -pi/2
          F = hpf('pi',NDig)./2;
          if X.Sign < 0
            F.Sign = -1;
          end
          return
        end
      end
      % X must now lie in the open interval (-1,1).
      
      % get a double precision version of X
      dx = double(X);
      
      % compute asin(dx). This gives us a pretty good estimate of the
      % arcsine, accurate to roughly 16 digits or so. Convert that
      % number back into an hpf form, with the desired number of digits.
      v = hpf(asin(dx),NDig);
      
      % we will need to compute sin(v) and cos(v)
      sinv = sin(v);
      one = hpf('1',NDig);
      %      cosv = sqrt(one - sinv.*sinv); % this introduces too much error
      cosv = cos(v);
      
      % now, do the transformation to get sin(u)
      z = X.*cosv - sinv.*sqrt(one-X.*X);
      
      % sin(u) is a very small number. Typically on the order of 1e-16,
      % but we are working with hpf numbers, so that is not a dauntingly
      % small number.
      u = z;
      term = z;
      z2 = z.*z; % this is a series that goes up with z^2
      % we need to know when to stop the series.
      logz = log10(double(z));
      logterm = logz;
      n = 1;
      while -logterm <= sum(NDig)
        term = term.*z2.*(n.*n)./((n+1).*(n+2));
        logterm = logterm + 2*logz + log10((n.*n)./((n+1).*(n+2)));
        u = u + term;
        n = n + 2;
      end
      
      % combine the two parts of our solution
      F = u + v;
      
    end % function F = asin(X)
    
    function F = asind(X)
      % evaluates asin(X) for a hpf number X (returning a result in degrees)
      
      % simple conversion of degrees to radians
      F = X;
      for i = 1:numel(X)
        F(i) = asin(X(i)).*180./hpf('pi',X(i).NumberOfDigits);
      end
      
    end % function F = asind(X)
    
    function F = asinh(X)
      % evaluates asinh(X) (the inverse hyperbolic sine) for an hpf number X
      
      F = X;
      for i = 1:numel(X)
        if isnan(X(i).Numeric)
          F(i) = hpf('NaN',X(i).NumberOfDigits);
        elseif isinf(X(i).Numeric)
          F(i) = hpf('inf',X(i).NumberOfDigits);
          F(i).Sign = X(i).Sign;
        elseif X(i).Sign == 0
          F(i) = hpf('0',X(i).NumberOfDigits);
        else
          % simple inverse for normal numbers
          F(i) = log(X(i) + sqrt(X(i).*X(i) + 1));
        end
      end
      
    end % function F = asinh(X)
    
    function F = atan(X)
      % evaluates atan(X) (the inverse tangent function with one argument) for a hpf number X
      %
      % The method used is the standard trigonometric series for atan,
      % but with a twist. If we have a good estimate of atan(X), then
      % can transform the problem using a simple trig identity.
      %
      %   tan(u + v) = (tan(u) + tan(v))/(1 - tan(u)*tan(v))
      %
      % So if v is our approximate estimate for atan(X), then we wish to
      % compute the value of u, such that sin(u+v) = sin(X). Along the way,
      % we will need to do a little algebraic manipulation on that identity
      % above. Starting from the point
      %
      %   X = (tan(u) + tan(v))/(1 - tan(u)*tan(v))
      %
      % we transform that into
      %
      %   tan(u) = (tan(v) - X)/(1 - X*tan(v))
      %
      % v is chosen to make u small, and is typically within 1 part in
      % 10^16 of the sum (u+v) where it is u+v that we wish to compute.
      % This is done by converting X into a double precision form, then
      % using the double version of atan to give us a leg up.
      %
      % Once we have done the above transformation, we use a series for
      % atan that will yield roughly an additional 32 digits of accuracy
      % per term, It is just the standard atan Taylor series, but for
      % very small numbers, it will be tremendously efficient.
      % Unfortunately, that series still requires one divide for each
      % term.
      
      % is X a scalar?
      if isempty(X)
        F = X;
        return
      elseif numel(X) > 1
        % vector or array
        F = X;
        for i = 1:numel(X)
          F(i) = atan(X(i));
        end
        return
      end
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if isnan(X.Numeric)
        % asin(NaN) = NaN
        F = hpf('NaN',NDig);
        return
      elseif isinf(X.Numeric)
        % inf or -inf mapps to +/- pi/2
        F = hpf('pi',NDig)/2;
        if X.Sign < 0
          F.Sign = -1;
        end
        return
      elseif X.Sign == 0
        % asin(0) == 0, so atan(0) == 0 too
        F = hpf('0',NDig);
        return
      end
      
      % get a double precision version of X
      dx = double(X);
      
      % compute atan(dx). This gives us a pretty good estimate of the
      % arctan, accurate to roughly 16 digits or so. Convert that
      % number back into an hpf form, with the desired number of digits.
      v = hpf(atan(dx),NDig);
      
      % we will need to compute tan(v)
      tanv = tan(v);
      
      % now, do the transformation to get tan(u)
      z = (X - tanv)./(1 + X.*tanv);
      
      if z ~= 0 %#ok
        % tan(u) is a very small number. Typically on the order of 1e-16,
        % but we are working with hpf numbers, so that is not a dauntingly
        % small number.
        u = z;
        term = z;
        % logterm is an approximation to the negative base 10 log of each
        % term since it won't include the 1/n in it, this is a conservative
        % approximation.
        logterm = log10(abs(double(z)));
        logz2 = 2*logterm;
        z2 = z.*z; % this is a series that goes up with z^2
        n = 1;
        while logterm >= -sum(NDig)
          term = uminus(term.*z2);
          logterm = logterm + logz2;
          n = n + 2;
          u = u + term./n;
        end
        % combine the two parts of our solution
        F = u + v;
      else
        % we got lucky and nailed v exactly with the double computation
        F = v;
      end
      
    end % F = atan(X)
    
    function F = atand(X)
      % evaluates atan(X) for a hpf number X (returning a result in degrees)
      
      % simple conversion of degrees to radians
      F = X;
      for i = 1:numel(X)
        F(i) = atan(X(i)).*180./hpf('pi',X(i).NumberOfDigits);
      end
      
    end % function F = atand(X)
    
    function F = atanh(X)
      % evaluates atanh(X) (the inverse hyperbolic tangent) for an hpf number X
      F = X;
      for i = 1:numel(X)
        if isnan(X(i).Numeric)
          F(i) = hpf('NaN',X(i).NumberOfDigits);
        elseif X(i).Sign == 0
          F(i) = hpf('0',X(i).NumberOfDigits);
        else
          F(i) = log((1+X(i))./(1-X(i)))./2;
        end
      end
    end % function F = atanh(X)
    
    function F = augmentdigits(F,N)
      % allows the user to increase or decrease the number of significant digits carried
      % usage: F = augmentdigits(F,N)
      %
      % N must be a scalar, positive integer or integer vector
      %    of length 2. If N is provided as a scalar, then the
      %    number of hidden (shadow) digits will be unchanged.
      if (nargin ~= 2)
        error('HPF:augmentdigits','Must provide two arguments')
      elseif ~ismember(numel(N),[1 2]) || any(N < 0) || (N(1) == 0) || any(N ~= round(N))
        error('HPF:augmentdigits','N must be integer, N(1) > 0, N(2) >= 0')
      end
      
      % IS F an array?
      if numel(F) > 1
        for i = 1:numel(F)
          F(i) = augmentdigits(F(i),N);
        end
        return
      end
      
      % was only the first element supplied?
      if numel(N) == 1
        N(2) = F.NumberOfDigits(2);
      end
      
      % remember though, that the total digits must be divisible by
      % the DecimalBase, so we may need to increase the number of shadow
      % digits.
      r = mod(sum(N),F.DecimalBase);
      if r > 0
        N(2) = N(2) + (F.DecimalBase - r);
      end
      
      NDigLimit = [3.6e14 3.6e12 3.6e10 3.6e8 3.6e6 36000];
      if NDigLimit(F.DecimalBase) <= sum(N)
        switch F.DecimalBase
          case 1
            error('HPF:augmentdigits', ...
              'Cannot exceed 3.6e14 total digits for a DecimalBase of 1')
          case 2
            error('HPF:augmentdigits', ...
              'Use a smaller base. Cannot exceed 3.6e12 total digits for a DecimalBase of 2')
          case 3
            error('HPF:augmentdigits', ...
              'Use a smaller base. Cannot exceed 3.6e10 total digits for a DecimalBase of 3')
          case 4
            error('HPF:augmentdigits', ...
              'Use a smaller base. Cannot exceed 3.6e8 total digits for a DecimalBase of 4')
          case 5
            error('HPF:augmentdigits', ...
              'Use a smaller base. Cannot exceed 3.6e6 total digits for a DecimalBase of 5')
          case 6
            error('HPF:augmentdigits', ...
              'Use a smaller base. Cannot exceed 36000 total digits for a DecimalBase of 6')
        end
      end
      
      % do we increase or decrease the number, or just leave it alone?
      if sum(N) > sum(F.NumberOfDigits)
        % it is an increase
        % padding as necessary will be done with zero migits on the right
        F.Migits = padz(F.Migits,sum(N)/F.DecimalBase,'right');
        
      elseif sum(N) < sum(F.NumberOfDigits)
        % we have been asked to decrease the number of digits
        % since we know that both sum(N) and sum(F.NumberOfDigits)
        % are divisible by the DecimalBase, then we are actually
        % dropping off one or more migits here, so we need to do
        % rounding using that truncated migit. The change will always
        % be at a migit boundary, because the total number of digits
        % in all cases is maintained as a multiple of the DecimalBase.
        k = sum(N)/F.DecimalBase; % an integer
        
        truncatedmigit = F.Migits(k+1);
        if truncatedmigit < (F.Base/2)
          % just a chop
          F.Migits = F.Migits(1:k);
        else
          % a round up
          F.Migits = F.Migits(1:k);
          F.Migits(end) = F.Migits(end) + 1;
          
          % did this cause a carry?
          if F.Migits(end) >= F.Base
            [mant,exponentshift] = carryop(F.Migits,F.DecimalBase,F.Base);
            F.Migits = mant;
            F.Exponent = F.Exponent + exponentshift;
          end
        end
        
      else
        % the total number of digits are unchanged, but
        % this may be a change in the number of shadow digits.
        % so just drop down.
      end
      F.NumberOfDigits = N;
      
    end % function F = augmentdigits(F,N)
    
    function F = ceil(F)
      % A round towards plus inf
      %
      % Equivalent to fix(F) for negative F
      for i = 1:numel(F)
        if isfinite(F(i).Numeric)
          if F(i).Sign < 0
            F(i) = fix(F(i));
          elseif F(i).Sign > 0
            % positive numbers round up, so there may be a carry
            % to be chased down.
            db = F(i).DecimalBase;
            
            % Where does the decimal place fall?
            decind = F(i).Exponent;
            
            % if decind is larger than the total number of digits in the
            % number, then the ceil is a no-op.
            if decind <= 0
              % ceil rounds up for positive F
              F(i) = hpf('1',F(i).NumberOfDigits);
            elseif decind < sum(F(i).NumberOfDigits)
              % decind is in the middle somewhere, so the ceil operation
              % must be done, but on exactly what migit will that happen,
              % or does it fall between a pair of migits?
              mind = ceil(decind./db);
              mindi = mod(decind,db);
              if mindi == 0
                if any(F(i).Migits(mind:end) > 0)
                  F(i).Migits((mind+1):end) = 0;
                  F(i).Migits(mind) = F(i).Migits(mind) + 1;
                  if F(i).Migits(mind) >= F(i).Base
                    [F(i).Migits,exponentshift] = carryop(F(i).Migits,F(i).DecimalBase,F(i).Base);
                    F(i).Exponent = F(i).Exponent + exponentshift;
                  end
                end
              else
                % splitting a migit
                d = m2d(F(i).Migits(mind),db);
                if any(d((mindi+1):end) > 0) || any(F(i).Migits(mind:end) > 0)
                  F(i).Migits((mind+1):end) = 0;
                  d(mindi) = d(mindi) + 1;
                  d((mindi+1):end) = 0;
                  m = d*10.^((db-1):-1:0)';
                  F(i).Migits(mind) = m;
                  if F(i).Migits(mind) >= F(i).Base
                    [F(i).Migits,exponentshift] = carryop(F(i).Migits,F(i).DecimalBase,F(i).Base);
                    F(i).Exponent = F(i).Exponent + exponentshift;
                  end
                end
              end
            end
          end
        end
      end

    end % F = ceil(F)
    
    function F = cos(X)
      % evaluates cos(X) for a hpf number X (in radians)
      
      % is X a scalar?
      if isempty(X)
        F = X;
        return
      elseif numel(X) > 1
        % vector or array
        F = X;
        for i = 1:numel(X)
          F(i) = cos(X(i));
        end
        return
      end
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if ~isfinite(X.Numeric)
        % cos(NaN) = cos(inf) = cos(-inf) = NaN
        F = hpf('NaN',NDig);
        return
      elseif X.Sign == 0
        % trap exact zeros.
        F = hpf('1',NDig);
        return
      end
      
      % get pi for that number of digits because cos(X) is a periodic
      % function. This allows us to do range reduction.
      pie = hpf('pi',NDig);
      twopie = 2.*pie;
      
      % We need to reduce X into the interval [-pi,pi].
      k = X./twopie;
      kfrac = fractionalpart(k);
      X = X - floor(k)*twopie;
      
      % convert kfrac to a double, to most easily determine the subinterval
      kfrac= double(kfrac);
      if kfrac < 0
        kfrac = 1 + kfrac;
      end
      
      % what sub-interval did k live in?
      % the break points for our bins are {1/8, 3/8, 5/8, 7/8}
      % so one trick here might have been to multiply by 8, then just take
      % the first digit pf the result. But it will probably be just as fast
      % to convert to double.
      bins = [1 3 5 7]/8;
      [~,kindex] = histc(kfrac,bins);
      switch kindex
        case 1
          % X was in the interval 2*pi*(n + [1/8 3/8]). This will be most
          % efficiently approximated using a shift, then a call to coscore
          F = uminus(sincore(X-pie./2));
        case 2
          % X was in the interval 2*pi*(n + [3/8 5/8]). This will be most
          % efficiently approximated using a shift, then a call to sincore
          F = uminus(coscore(X-pie));
        case 3
          % X was in the interval 2*pi*(n + [1/8 3/8]). This will be most
          % efficiently approximated using a shift, then a call to coscore
          % Note that 3/4 is exactly converted to a hpf number here.
          F = sincore(X-pie.*(3/2));
        otherwise
          % X was in the interval 2*pi*(n + [-1/8 1/8]). This will be
          % simply done as a call to sincore.
          F = coscore(X);
      end
      
    end % function F = cos(X)
    
    function F = cosd(X)
      % evaluates cos(X) for a hpf number X (in degrees)
      
      % simple conversion of degrees to radians
      F = X;
      for i = 1:numel(X)
        F(i) = cos(X(i).*hpf('pi',X(i).NumberOfDigits)./180);
      end
      
    end % function F = cosd(X)
    
    function F = cosh(X)
      % evaluates cosh(X) (the hyperbolic cosine) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology.
      F = (exp(X) + exp(-X))./2;
      
    end % function F = sinh(X)
    
    function F = coth(X)
      % evaluates coth(X) (the hyperbolic cotangent) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology.
      % A loop is necessary however, to catch any +infs, which would
      % otherwise end up as inf/inf, therefore NaNs. To be consistent with
      % MATLAB, coth(inf) = 1,
      F = X;
      for i = 1:numel(X)
        if isinf(X(i).Numeric) && (X(i).Sign > 0)
          F(i) = hpf('1',X(i).NumberOfDigits);
        else
          Fi = exp(2*X(i));
          F(i) = (Fi+1)./(Fi-1);
        end
      end
      
    end % function F = coth(X)
    
    function F = csc(X)
      % evaluates csc(X) (csc(X) = 1/sin(X)) for an hpf number X in radians
      
      % a simple scheme suffices for now, costing only a spare divide
      F = reciprocal(sin(X));
      
    end % function F = csc(X)
    
    function F = cscd(X)
      % cscd(X) (cosecant function, where X is in degrees) for an hpf number X
      
      % a simple scheme suffices for now, costing only a spare divide
      F = reciprocal(sind(X));
      
    end % function F = cscd(X)
    
    function F = csch(X)
      % evaluates csch(X) (the hyperbolic cosecant) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology.
      F = exp(X);
      F = (2.*F)./(F.*F - 1);
      
    end % function F = csch(X)
    
    
    function xroot = cubrt(x)
      % Computes the cube root of x, where x is an hpf number.
      % usage: xroot = cubrt(x)
      %
      % The method employed uses only multiplications to compute the cube
      % root, in a quadratically convergent algorithm, i.e., Newton's
      % method, applied to the inverse cube root.
      %
      % A derivation for this scheme can be found online, at:
      %  http://www.azillionmonkeys.com/qed/sqroot.html
      
      % is x a scalar or something else?
      xroot = x;
      if isempty(x)
        return
      elseif numel(x) > 1
        % vector or array
        for i = 1:numel(x)
          xroot(i) = cubrt(x(i));
        end
        return
      elseif ~isfinite(x)
        % catch the nans or infs.
        return
      end
      
      % check for zero or negative
      xSign = x.Sign;
      if xSign == 0
        % a no-op
        return
      end
      
      % how many digits do we need?
      NDig = x.NumberOfDigits;
      
      % reduce the problem to bring x into the interval [0.1,10)
      xroot = hpf('1',NDig);
      p = floor(x.Exponent/3);
      xroot.Exponent = p + 1;
      x.Exponent = x.Exponent - 3*p;
      % we have now reduced x to a reasonable number in case it would
      % have otherwise exceeded the dynamic range of a double. x must
      % now lie in the half open interval [.1,1).
      
      % get a (very) good estimate for the inverse cube root of x
      % using the double precision cube root.
      xest = hpf(1./nthroot(double(x),3),NDig);
      
      % Newton on the inverse cube root, a divide-free method. This
      % will double the number of correct digits for each iteration
      % of the process, and it is a divide-free iteration so very fast.
      % note that the division by 3 is done as a multiply.
      threeinv = hpf('1',NDig)./3;
      four = hpf('4',NDig);
      flag = true;
      cycle = 0;
      while flag
        xprior = xest.Migits;
        xest = threeinv.*xest.*(four - x.*xest.*xest.*xest);
        
        k = find(xest.Migits ~= xprior,1,'first');
        if isempty(k) || (cycle > 2)
          flag = false;
        else
          if k == numel(xprior)
            cycle = cycle + 1;
          end
        end
      end
      % merge the pieces together.
      % so this recovers the desired root with nary a hpf divide required.
      xroot = xroot.*x.*xest.*xest;
      
    end % function xroot = cubrt(x)
    
    function Asum = cumsum(A,dim)
      % cumsum: cumulative sum of an HPF array
      % usage: Asum = cumsum(A,dim)
      %
      % Arguments:
      %  A - an HPF object array
      %
      %  dim - (optional) dimension of A to sum over
      %      DEFAULT: dim is the first non-singleton
      %      dimension of A.
      %
      % Arguments: (output)
      %  Asum - the cumulative sum HPF object array
      %
      % Example:
      %  A = HPF(-3:2:25);
      %  double(cumsum(A))
      %
      %  ans =
      %    -3  -4  -3  0  5  12  21  32  45  60  77  96 117 140 165
      %
      %  See also: prod, sum, cumprod
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com

      if (nargin<1) || (nargin>2)
        error('cumsum takes one or two arguments')
      end
      
      if numel(A) == 1
        % the sum of a scalar is a no-op
        Asum = A;
      else
        % a vector or array
        
        % default for dim?
        if (nargin==1) || isempty(dim)
          dim = find(size(A)>1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        
        % sum over the dimension dim. Do so by
        % a permutation, then a sum and an
        % inverse permute.
        P = 1:length(size(A));
        P([1,dim]) = [dim,1];
        
        A = permute(A,P);
        NAsum = size(A);
        N1 = NAsum(1);
        Asum = reshape(A,N1,[]);
        for j = 1:size(Asum,2)
          for i = 2:N1
            Asum(i,j) = Asum(i-1,j) + Asum(i,j);
          end
        end
        
        % do the inverse permutation
        Asum = ipermute(reshape(Asum,NAsum),P);
        
      end
      
    end % function Asum = cumsum(A,dim)
    
    function Aprod = cumprod(A,dim)
      % cumprod: cumulative product of an HPF array
      % usage: Aprod = cumprod(A,dim);
      %
      % Arguments:
      %  A - an HPF object array
      %
      %  dim - (optional) dimension of A to prod over
      %      DEFAULT: dim is the first non-singleton
      %      dimension of A.
      %
      % Arguments: (output)
      %  Aprod - the product HPF object array
      %
      % Example:
      %  A = hpf(-3:2:12);
      %  double(cumprod(A))
      %
      % ans =
      % HPF element: (1,1)
      % -3
      % HPF element: (1,2)
      % 3
      % HPF element: (1,3)
      % 3
      % HPF element: (1,4)
      % 9
      % HPF element: (1,5)
      % 45
      % HPF element: (1,6)
      % 315
      % HPF element: (1,7)
      % 2835
      % HPF element: (1,8)
      % 31185
      %
      %  See also: cumsum, prod, sum
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      
      if (nargin<1) || (nargin>2)
        error('prod takes one or two arguments')
      end
      
      if numel(A) == 1
        % the product of a scalar is a no-op
        Aprod = A;
      else
        % a vector or array
        
        % default for dim?
        if (nargin==1) || isempty(dim)
          dim = find(size(A)>1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        
        % product over the dimension dim. Do so by
        % a permutation, then a prod and an
        % inverse permute.
        P = 1:length(size(A));
        P([1,dim]) = [dim,1];
        
        A = permute(A,P);
        NAprod = size(A);
        N1 = NAprod(1);
        Aprod = reshape(A,N1,[]);
        for j = 1:size(Aprod,2)
          for i = 2:N1
            Aprod(i,j) = Aprod(i-1,j).*Aprod(i,j);
          end
        end
        
        % do the inverse permutation
        Aprod = ipermute(reshape(Aprod,NAprod),P);
        
      end

    end %     function Aprod = cumprod(A,dim)
    
    function Fchar = disp(F)
      % hpf/disp: Forms a character representation of a hpf object
      %
      %  See also: hpf/display, disp
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % Is F empty, or an array or vector?
      if isempty(F)
        Fchar = '  []';
        if nargout == 0
          disp(Fchar)
          clear Fchar
        end
        return
      elseif numel(F) > 1
        % an array or vector
        
        Nf = numel(F);
        Fsize = size(F);
        
        % Must be an array or vector with 2 or more elements.
        % An array will be dumped out in order, one element at a time.
        sub = cell(1,numel(Fsize));
        for i = 1:Nf
          % make a tag that identifies each element
          tag = 'HPF element: (';
          [sub{:}] = ind2sub(Fsize,i);
          for j = 1:length(Fsize)
            if j == 1
              tag = [tag,num2str(sub{j})]; %#ok
            else
              tag = [tag,',',num2str(sub{j})]; %#ok
            end
          end
          tag = [tag,')']; %#ok
          disp(tag)
          disp(F(i))
        end
        return
      else
        % disp for a scalar hpf number
        
        % is it a finite number?
        if isnan(F.Numeric)
          % its a NaN
          Fchar = '  NaN';
        elseif isinf(F.Numeric)
          % its an inf
          if F.Sign > 0
            Fchar = '  Inf';
          else
            Fchar = ' -Inf';
          end
        else
          % a simple number, but we need to dump it out in character form
          
          % do we need a minus sign in front?
          if F.Sign == 0
            % exact zero
            Fchar = '    0';
            if nargout == 0
              disp(Fchar)
              clear Fchar
            end
            
            % we are all done
            return
          elseif F.Sign>0
            % positive
            Fchar = '';
          else
            % negative, so stick in the -
            Fchar = '-';
          end
        
          % What is the base for each migit element?
          DBase = F.DecimalBase;
          
          % and the number of digits to be displayed
          NDig = F.NumberOfDigits;
          
          if NDig(2) > 0
            % rounding any shadow digits away is necessary
            % to display the correctly rounded value.
            F = roundn(F,F.Exponent - F.NumberOfDigits(1));
          end
          NDig = NDig(1);
          
          % extract the digits themselves, as a chararacter string
          dig = m2d(F.Migits,DBase);
          cdig = char(dig(1:NDig) + '0');
          
          % decide if scientific notation is merited for large or small
          % exponents, or if not, then where does the decimal point go
          if F.Exponent == 0
            nnzlast = find(cdig ~= '0',1,'last');
            if isempty(nnzlast)
              nnzlast = 1;
            elseif nnzlast >= (NDig*0.75)
              nnzlast = NDig;
            end
            cdig = ['0.',cdig(1:nnzlast)];
            
          elseif (F.Exponent < 0) && ((-F.Exponent*4) <= NDig) && (F.Exponent > -50)
            % small negative exponents compared to the size of the number,
            % best to just add a few zeros in front
            nnzlast = find(cdig ~= '0',1,'last');
            if isempty(nnzlast)
              nnzlast = 1;
            elseif nnzlast >= (NDig*0.75)
              nnzlast = NDig;
            end
            cdig = ['0.',repmat('0',1,abs(F.Exponent)),cdig(1:nnzlast)];
            
          elseif (F.Exponent > 0) && ((F.Exponent+1) < numel(cdig))
            % we can just insert a decimal point here
            cdig = [cdig(1:(F.Exponent)),'.',cdig((F.Exponent+1):end)];
            nnzlast = find(cdig ~= '0',1,'last');
            if cdig(nnzlast) == '.'
              nnzlast = nnzlast - 1;
            elseif nnzlast >= (NDig*0.75)
              nnzlast = NDig + 1 + (F.Sign < 0);
            end
            cdig = cdig(1:min(nnzlast,length(cdig)));
          else
            % use scientific notation
            nnzlast = find(cdig ~= '0',1,'last');
            if isempty(nnzlast)
              nnzlast = 1;
            elseif nnzlast >= (NDig*0.75)
              nnzlast = NDig;
            end
            cdig = [cdig(1),'.',cdig(2:nnzlast)];
            cdig = [cdig,'e',num2str(F.Exponent - 1)];
          end
          
          Fchar = [Fchar,cdig];
        end
        
        if nargout == 0
          disp(Fchar)
          clear Fchar
        end
      end
      
    end % Fchar = disp(F)
    
    function display(F)
      % hpf/display: Display a hpf object, calls disp
      %
      %  See also: hpf/disp, display
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      name = inputname(1);
      if ~isempty(name)
        disp([name,' ='])
      else
        disp('ans =')
      end
      disp(F)
      
    end % display(F)
    
    function F = double(X)
      % Convert an hpf number into a double precision approximation.
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      if isempty(X)
        % empty begets empty
        F = [];
      elseif numel(X) > 1
        % F is an array
        F = zeros(size(X));
        for i = 1:numel(X)
          F(i) = double(X(i));
        end
      else
        % X is scalar. any special cases?
        if ~isfinite(X.Numeric)
          F = X.Numeric.*X.Sign;
        else
          % This is a scalar real number, although we may get an overflow.
          % It should also map a flint into a flint, as long as it
          % fits into the 2^53-1 limits imposed by a double.
          
          % Extract the mantissa, but only take the first 20 or so digits
          % at a max. We have no need for any lower migits than this.
          Mmax = min(ceil(20./X.DecimalBase),numel(X.Migits));
          highmigits = X.Migits(1:Mmax);
          
          % convert the high order digits to simple digits
          highdigits = m2d(highmigits,X.DecimalBase);
          
          F = 0;
          E = 10.^(X.Exponent - 1);
          for ind = 1:numel(highdigits)
            F = F + highdigits(ind).*E;
            E = E ./ 10;
          end
          if X.Sign < 0
            F = -F;
          end
        end
      end
    end % F = double(X)
    
    function D = eps(X)
      % D = EPS(X), is the positive distance from ABS(X) to the next larger
      % number in magnitude of the floating point HPF number.
      %
      % Example:
      %  X = hpf(1,[20 0])
      % % X =
      % %   1.
      %
      %  X + eps(X)
      % % ans =
      % %     1.0000000000000000001
      %
      %  X + eps(X)*0.9999999
      % % ans =
      % %     1.
      %
      %  See also: eps
      %
      %  Author: John D'Erricoj
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      if isempty(X)
        % empty begets empty
        D = X;
      elseif numel(X) > 1
        % X is an array
        D = X;
        for i = 1:numel(X)
          D(I) = eps(X(i));
        end
      else
        % X is scalar
        
        % eps (inf) = eps(nan) = nan
        if ~isfinite(X.Numeric)
          D = hpf('NaN',X.NumberOfDigits);
        else
          D = hpf('1',X.NumberOfDigits);
          D.Exponent = X.Exponent + 1 - sum(X.NumberOfDigits);
        end
      end
    end % D = eps(X)
    
    function result = eq(x,y)
      % Test for numerical equality between numbers, returning a boolean result
      % usage: result = eq(x,y);
      % usage: result = x == y;
      
      % is either x or y an array? If so, we need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x) || isnan(y)
          result = false;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        %
        % Note that comparing an HPF number to a real number might be
        % a problem if the HPF number has less than 55 digits or so,
        % AND if the real number if not representable exactly.
        if ~isa(x,'hpf')
          if (x == 0)
            result = (y.Sign == 0);
            return
          end
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          if (y == 0)
            result = (x.Sign == 0);
            return
          end
          y = hpf(y,x.NumberOfDigits);
        end
        
        % was one or both numbers an inf?
        if isinf(x.Numeric)
          if isinf(y.Numeric)
            result = (x.Sign == y.Sign);
          else
            result = false;
          end
        elseif isinf(y.Numeric)
          % since x was not an inf
          result = false;
        else
          
          % We need to ensure that both are hpf numbers to compare them.
          % At least ONE of the numbers must have been an hpf though.
          if ~isa(x,'hpf')
            x = hpf(x,y.NumberOfDigits);
          end
          if ~isa(y,'hpf')
            y = hpf(y,x.NumberOfDigits);
          end
          
          % do they have the same decimal base? choose the smaller one
          % if they differ.
          if x.DecimalBase ~= y.DecimalBase
            if x.DecimalBase < y.DecimalBase
              y = adjustdecimalbase(y,x.DecimalBase);
            else
              x = adjustdecimalbase(x,y.DecimalBase);
            end
          end
          
          % augment a shorter number so that both are the same lengths
          % spare zero digits carried are not relevant to equality.
          if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
            x = augmentdigits(x,y.NumberOfDigits);
          elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
            y = augmentdigits(y,x.NumberOfDigits);
          end
          
          % now that we have made x and y comparable (and excluded the
          % special numbers, so the comparison itself is easy
          result = isequal(x,y);
          
        end
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = false(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = eq(x,y(i));
        end
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = false(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = eq(x(i),y);
        end
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:EQ:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = eq(x(i),y(i));
        end
      end
      
    end % result = eq(x,y)
    
    function F = erf(X)
      % erf(X) for an hpf number X
      
      % special cases
      if isempty(X)
        % empty begets empty
        F = X;
      elseif numel(X) > 1
        % just do scalar computation
        F = X;
        for i = 1:numel(X)
          F(i) = erf(X(i));
        end
      elseif isnan(X)
        % nan begets nan
        F = X;
      elseif isinf(X)
        %
        if X.Sign > 0
          F = hpf('1',X.NumberOfDigits);
        else
          F = hpf('-1',X.NumberOfDigits);
        end
      elseif X.Sign == 0
        % zero goes to zero
        F = hpf('0',X.NumberOfDigits);
      else
        % get the sign of X. We know it is not zero.
        Xsign = sign(X);
        
        % just make it positive
        X = abs(X);
        
        % how many digits do we need here?
        NDig = X.NumberOfDigits;
        
        % convert X into a double to make the decision about the
        % method to be employed
        dx = double(X);
        
        % How many terms will it take to get the required number of digits
        % of precision from the simple erf series? I've used the simple
        % stirling's approximation for n! here, which tends to underpredict
        % n! by a small amount. Consequently, that will be slightly
        % conservative in the number of terms necessary.
        T = @(N) (2*N*log(dx) - log(2*N+1) - N.*log(N)- log(2*pi*N)/2 + ...
          N + log(2/sqrt(pi)))/log(10) + sum(NDig);
        
        if T(4) < 0
          mterms = 4;
        else
          mterms = termsbisector(T);
        end
        
        % how many digits can we get from the asymptotic approximation
        % for erfc for this value of X? The terms stop getting smaller when
        % 2*x^2 = 2*n-1, therefore at n = (2*x^2+1)/2
        m = floor((2*dx^2+1)/2);
        
        % what is the magnitude of the asymptotic approximation term
        % at term m? We can do no better than that.
        amag = (-dx.^2 - (2*m+1)*log(dx) - log(sqrt(pi)) - ...
          m*log(2) + gammaln(2*m + 1) - gammaln(m + 1))/log(10);
        
        if (amag < -sum(NDig)) && m < mterms
          % use the asymptotic expansion for erfc(X)
          X2 = X.*X;
          term = exp(-X2)./X./sqrt(hpf('pi',NDig));
          F = term;
          n = 1;
          coef = hpf('0.5',NDig).*reciprocal(X2);
          while (term.Exponent > -sum(NDig)) && (n <= m)
            term = -term.*coef;
            F = F + term;
          end
          
          % when all done, erf(X) = 1 - erfc(X)
          F = 1 - F;
          
        else
          % use the direct Taylor series expansion
          
          % If X <= 1, then the terms in the Taylor series expansion have
          % their maximum at the very first term. For X greater than 1
          % however, they may grow significantly larger than the first term.
          % All of this is important, because terms that grow in size (and
          % alternate in sign) then cause subtractive cancellation
          % problems. In fact, the index of the term of maximum value
          % occurs at roughly the solution of the equation
          %  2*log(X) - log(N) - 2/(2*N+1) = 0
          % A quick and dirty approximate solution to that equation is
          %  N = X^2.
          % At that value of N, what is the magnitude of that term?
          % Get the common log of that term
          nmax = max(1,floor(dx^2) + [-10 10]);
          nmax = nmax(1):nmax(2);
          t = ((2*nmax + 1).*log(dx) - gammaln(nmax + 1) - log(2*nmax + 1) + ...
            log(2/sqrt(pi)))/log(10);
          tmax = max(t);
          
          % the number of spare digits we will need to carry is given by
          % ceil(tmax).
          erfspares = max(0,1 + ceil(tmax));
          X = augmentdigits(X,NDig + [erfspares 0]);
          
          term = X.*2./sqrt(hpf('pi',NDig + [erfspares 0]));
          F = term;
          X2 = X.*X;
          for n = 1:mterms
            term = -term.*X2./n;
            F = F + term./(2*n+1);
          end
          
          % undo those spare digits
          F = augmentdigits(F,NDig);
        end
        
        % restore the sign
        F.Sign = Xsign;
        
      end
    end % function F = erf(X)
    
    function F = exp(X)
      % exp(X) for a hpf number X
      %
      % Running the Taylor series backwards does two things. First, it is
      % more accurate, since we add the smallest terms into the sum first.
      % Secondly, it allows us to avoid almost all of the necessary divides,
      % replacing them by multiplies, which are far faster for hpf mumbers.
      
      % Is this an array?
      if isempty(X)
        % empty begets empty
        F = X;
        return
      elseif numel(X) > 1
        F = X;
        for i = 1:numel(X)
          F(i) = exp(X(i));
        end
        return
      end
      % if we drop through, then X was scalar
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if isnan(X)
        % exp(NaN) = NaN
        F = hpf('NaN',NDig);
        return
      elseif X.Sign == 0
        % exp(0) == 1
        F = hpf('1',NDig);
        return
      elseif (X.Numeric == inf) && (X.Sign < 0)
        % exp(-inf) == 0
        F = hpf('0',NDig);
        return
      elseif X.Numeric == inf
        % exp(inf) == inf
        F = hpf('inf',NDig);
        return
      end
      
      % find the nearest integer to X
      Xi = round(X);
      
      % z is the fractional part that we have rounded off
      z = X - Xi;
      
      % exp(0) == 1, so initialize it that way
      F = hpf('1',NDig);
      
      % we might be able to gain a bump in convergence if we can
      % choose an appropriate multiple of log(10) (natural log) to
      % shift X. The goal here is to choose integer K such that
      % (X + K*log(10)) has a minimal difference from an integer,
      % AND where the integer part of the sum is a positive integer.
      %   dx = double(X);
      dz = double(z);
      ztol = 0.01;
      Kbest = 0;
      if abs(dz) > ztol
        dx = double(X);
        
        zbest = dz;
        if dx <= 0
          Kdir = 1;
        else
          Kdir = -1;
        end
        K = 0;
        Kmax = 100;
        Kflag = true;
        ln10 = log(10);
        % a simple loop on K will suffice
        while Kflag
          K = K + Kdir;
          znew = dz + K*ln10;
          znew = znew - round(znew);
          xnew = dx + K*ln10;
          if abs(znew) < abs(zbest)
            zbest = znew;
            Kbest = K;
            if abs(zbest) < ztol
              % this is good enough. so we can stop here.
              Kflag = false;
            end
          end
          if (abs(K) >= Kmax) || (abs(xnew) > abs(dx))
            Kflag = false;
          end
        end
        % when we drop out, we have either found a value of K such that
        % the series will converge as quickly as possible, or there is no
        % such value that does not cause a loss of precision in the result.
        
        if Kbest ~= 0
          ln10 = hpf('ln10',NDig);
          X = X + Kbest*ln10;
          
          % recompute z as an hpf, and Xi.
          Xi = round(X);
          % z is the fractional part that we have rounded off
          z = X - Xi;
        end
      end
      
      % If X overflows uint64, then we have problems.
      Xi64 = uint64(abs(Xi));
      
      % first compute the exponential of the integer part. get the
      % binary representation of Xi. Recall that Xi is an hpf integer.
      S = sign(Xi);
      if S ~= 0
        % we will need the value of E = exp(1) for this part
        E = hpf('e',NDig);
        
        % X has an integer part, so we start by computing exp(Xi64)
        % We will convert Xi to binary, one bit at a time.
        % Find the least significant bit of Xi, one bit at a time.
        % There MUST be at least one non-zero digit, since Xi was
        % non-zero.
        while Xi64 ~= 0
          if mod(Xi64,2) == 1
            F = F.*E;
            
            % drop off that last bit of Xi
            Xi64 = Xi64 - 1;
          end
          Xi64 = Xi64/2;
          
          if Xi64 ~= 0
            E = E.*E;
          end
        end
      end
      % Xi will now have been bit-reduced to zero, while at the
      % same time, we have built up F as exp(abs(Xi)). If the sign
      % of Xi is less than zero, then we need to invert F now.
      if S < 0
        F = reciprocal(F);
      end
      
      % shift the exponent of F to reflect the number of
      % powers of 10 in Kbest
      if Kbest ~= 0
        F.Exponent = F.Exponent - Kbest;
      end
      
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
      fun = @(m) m.*(log(abs(double(z))) + 1) - log(m).*(m + 1/2) + ...
        (sum(NDig)+1)*log(10) - 1/2*(log(2*pi) - 1);
      
      % we use many series approximations in hpf, so a bisection
      % scheme will be used in sin, cos, and exp, (to name a few functions)
      % to determine the number of terms required for the desired accuracy
      % in this series. That step is needed because we can't simply test
      % the individual terms to know when convergence has occurred. Why
      % not? Because we run the sum for the series in reverse, starting
      % with the very last term! This avoids almost all possible divides,
      % expensive ops in hpf.
      mterms = termsbisector(fun);
      
      % run that loop in reverse now. This is why I wanted to know how
      % many terms would be necessary. By running the loop backwards, I
      % avoid divisions, and a divide is more expensive than a multiply for
      % an hpf number.
      Fz = hpf('1',NDig);
      Fact = Fz;
      for m = mterms:-1:1
        Fact = Fact*m;
        Fz = z*Fz + Fact;
      end
      % because we ran the loop backwards, in the end we need to
      % divide by factorial(mterms)
      Fz = Fz./Fact;
      
      % and multiply Fz into F
      F = F.*Fz;
      
    end % function F = exp(X)
    
    function F = factorial(X)
      % evaluates factorial(X) for non-negative integer HPF arguments,
      % in the interval [0,1e6].
      
      % is X a scalar or empty?
      if isempty(X)
        % empty begets empty
        F = X;
      else
        % copy X into F to make it the proper size
        F = X;
        NDig = X(1).NumberOfDigits;
        
        % there is no need to worry about X being larger than 2^53 anyway,
        % since factorial would be impossible to compute there, at least
        % before the computer melts into a molten mass of silicon.
        X = double(X(:));
        % a realistic maximum is 1e6.
        Xmax = 1e6;
        
        if any(X ~= round(X)) || any(X < 0)
          error('HPF:FACTORIAL:invalid','Positive integer arguments required')
        end
        
        % first, sort the elements of X in ascending order.
        if numel(X) > 1
          [X,tags] = sort(X,'ascend');
        else
          tags = 1;
        end
        
        X0 = 1;
        Fact = hpf('1',NDig);
        for i = 1:numel(X)
          X1 = X(i);
          if X1 == X0
            % we already did this one
            F(tags(i)) = Fact;
          elseif X1 <= Xmax
            % we only need to compute a few products
            for n = (X0+1):X1
              Fact = Fact*n;
            end
            F(tags(i)) = Fact;
          else
            % just too big.
            F(tags(i)) = hpf('inf',NDig);
          end
          
          X0 = X1;
        end
      end
    end % function F = factorial(N)
    
    function F = fix(F)
      % Round towards zero - Truncates that part of A below the decimal point
      
      for i = 1:numel(F)
        if isfinite(F(i).Numeric)
          % if F was not finite, then fix is a no-op
          db = F(i).DecimalBase;
          
          % Where does the decimal place fall?
          decind = F(i).Exponent;
          
          % if decind is larger than the total number of digits in the
          % number, then the fix is a no-op.
          if decind <= 0
            % the result is zero for either sign, because fix always
            % rounds towards zero.
            F(i) = hpf('0',F(i).NumberOfDigits);
          elseif decind < sum(F(i).NumberOfDigits)
            % decind is in the middle somewhere, so the fix operation
            % must be done, but on exactly what migit will that happen,
            % or does it fall between a pair of migits?
            mind = ceil(decind./db);
            mindi = mod(decind,db);
            if mindi > 0
              % splitting a migit
              d = m2d(F(i).Migits(mind),db);
              d((mindi+1):end) = 0;
              F(i).Migits(mind) = d2m(d,db);
            end
            % anything beyond that point goes to zero
            F(i).Migits((mind + 1):end) = 0;
          end
        end
      end
      
    end % F = fix(F)
    
    function F = floor(F)
      % Round towards minus inf
      %
      % Equivalent to fix(F) for positive F, -ceil(-F) for negative F.
      for i = 1:numel(F)
        if isfinite(F(i).Numeric)
          if F(i).Sign > 0
            F(i) = fix(F(i));
          elseif F(i).Sign < 0
            F(i) = ceil(-F(i));
            F(i).Sign = -1;
          end
        end
      end
    end % F = floor(F)
    
    function F = fractionalpart(F)
      % Extracts the fractional part of a number F, thus F - fix(F)
      % Note that this will be negative for negative F.
      F = F - fix(F);
    end % F = fix(F)
    
    function result = ge(x,y)
      % Test for numerical inequality (greater than or equal to) between numbers, returning a boolean result
      % usage: result = ge(x,y);
      % usage: result = x >= y;
      
      % is either x or y an array? If so, we need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans, nor can they be compared to other numbers.
        if isnan(x) || isnan(y)
          result = false;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        if ~isa(x,'hpf')
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          y = hpf(y,x.NumberOfDigits);
        end
        
        % do they have the same decimal base? choose the smaller one
        % if they differ.
        if x.DecimalBase ~= y.DecimalBase
          if x.DecimalBase < y.DecimalBase
            y = adjustdecimalbase(y,x.DecimalBase);
          else
              x = adjustdecimalbase(x,y.DecimalBase);
          end
        end
        
        % augment a shorter number so that both are the same lengths
        % spare zero digits carried are not relevant to equality.
        if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
          x = augmentdigits(x,y.NumberOfDigits);
        elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
          y = augmentdigits(y,x.NumberOfDigits);
        end
        
        % was one or both numbers an inf?
        if isinf(x.Numeric) || isinf(y.Numeric)
          % this works in any case
          result = (x.Numeric*x.Sign >= y.Numeric*y.Sign);
        elseif x.Sign < y.Sign
          % if there is a sign differential, then we are done
          result = false;
        elseif x.Sign > y.Sign
          % if there is a sign differential, then we are done
          result = true;
        elseif x.Sign == 0
          % if we skip past the previous tests, then the signs must be
          % equal. If they are both zero, then we are done.
          result = true;
        elseif x.Sign > 0
          % positive x and y
          if x.Exponent > y.Exponent
            result = true;
          elseif x.Exponent < y.Exponent
            result = false;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = true;
            elseif x.Migits(k) > y.Migits(k)
              result = true;
            else
              result = false;
            end
          end
        else
          % negative x and y
          if x.Exponent > y.Exponent
            result = false;
          elseif x.Exponent < y.Exponent
            result = true;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = true;
            elseif x.Migits(k) > y.Migits(k)
              result = false;
            else
              result = true;
            end
          end
        end
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = false(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = ge(x,y(i));
        end
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = false(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = ge(x(i),y);
        end
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:GE:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = ge(x(i),y(i));
        end
      end
      
    end % result = ge(x,y)
    
    function result = gt(x,y)
      % Test for strict numerical inequality (greater than) between numbers, returning a boolean result
      % usage: result = gt(x,y);
      % usage: result = x > y;
      
      % is either x or y an array? If so, we need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans, nor can they be compared to other numbers.
        if isnan(x) || isnan(y)
          result = false;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        if ~isa(x,'hpf')
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          y = hpf(y,x.NumberOfDigits);
        end
        
        % do they have the same decimal base? choose the smaller one
        % if they differ.
        if x.DecimalBase ~= y.DecimalBase
          if x.DecimalBase < y.DecimalBase
            y = adjustdecimalbase(y,x.DecimalBase);
          else
            x = adjustdecimalbase(x,y.DecimalBase);
          end
        end
        
        % augment a shorter number so that both are the same lengths
        % spare zero digits carried are not relevant to equality.
        if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
          x = augmentdigits(x,y.NumberOfDigits);
        elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
          y = augmentdigits(y,x.NumberOfDigits);
        end
        
        % was one or both numbers an inf?
        if isinf(x.Numeric) || isinf(y.Numeric)
          % this works in any case
          result = (x.Numeric*x.Sign > y.Numeric*y.Sign);
        elseif x.Sign > y.Sign
          % if there is a sign differential, then we are done
          result = true;
        elseif x.Sign < y.Sign
          % if there is a sign differential, then we are done
          result = false;
        elseif x.Sign == 0
          % if we skip past the previous tests, then the signs must be
          % equal. If they are both zero, then we are done.
          result = false;
        elseif x.Sign > 0
          % positive x and y
          if x.Exponent > y.Exponent
            result = true;
          elseif x.Exponent < y.Exponent
            result = false;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = false;
            elseif x.Migits(k) > y.Migits(k)
              result = true;
            else
              result = false;
            end
          end
        else
          % negative x and y
          if x.Exponent > y.Exponent
            result = false;
          elseif x.Exponent < y.Exponent
            result = true;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = false;
            elseif x.Migits(k) > y.Migits(k)
              result = false;
            else
              result = true;
            end
          end
        end
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = false(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = gt(x,y(i));
        end
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = false(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = gt(x(i),y);
        end
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:LT:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = gt(x(i),y(i));
        end
      end
      
    end % result = gt(x,y)
    
    function F = hpf(N,NDig)
      % constructor for a hpf number from the scalar value N
      % usage: F = hpf(N,NDig)
      %
      % this construction will, given sufficient digits, be an
      % exact decimal form.
      %
      % arguments: (input)
      %  N - any scalar number in matlab of class: double, single,
      %      logical, uint8, int8, uint16, int16, uint32, int32, etc.
      %
      %      complex numbers are not yet supported in HPF form.
      %
      %      If N is a numeric array, then it will be converted into
      %      an array of HPF numbers.
      %
      %      If N is a character vector, then it will be presumed to
      %      indicate a scalar HPF number. 'pi', 'e' and 'ln10'
      %      are special numbers that HPF can supply with a high
      %      number of digits (up to 500,000 digits are stored for
      %      these special numbers.) Other numbers in character form,
      %      such as '1.2' or '6.0221417930e23' will be parsed into
      %      their exact decimal representation, subject to the number
      %      of significant digits specified, either by default or
      %      actual specification by the user.
      %
      %  NDig - (OPTIONAL) the number of digits to be stored for the
      %      resulting hpf number. The first element is the number
      %      of digits reported to you, whereas the second number
      %      is the number of additional, hidden (shadowed) digits
      %      carried to avoid issues of round-off.
      %
      %      DEFAULT: 64 4
      %
      %      This default can be reset if you prefer. For example,
      %      to change that default in your MATLAB session to 100
      %      decimal digits, you would execute the following command:
      %
      % % A permanent change to 100 digits of precision reported, without
      % % changing the current number of digits shadowed
      %       DefaultNumberOfDigits 100
      %
      % % A permanent change to 50 digits of precision reported, with
      % % 3 more digits shadowed. DefaultNumberOfDigits can be executed
      % % as a function or a command.
      %       DefaultNumberOfDigits(50,3)
      %
      % % A temporary change to 2000 digits of precision reported, with
      % % 5 more digits shadowed. Note that a "clear functions" call
      % % would override a change made only for the current matlab session.
      %       DefaultNumberOfDigits 2000 5 session
      
      % There are some useful numbers we might need to load in. make
      % them persistent in case we need to use them more than once.
      persistent specialnumbers
      
      % was hpf called with no args? If so, then it is the same as a
      % call to hpf with N == 0.
      if nargin == 0
        N = 0;
      end
      
      % hpf is a no-op on an existing hpf, unless there was
      % a change to the number of digits.
      if isa(N,'hpf')
        % was there a number of digits specified?
        if nargin > 1
          % if so, then this is effectively just a change in the number
          % of digits, which augmentdigits is designed to do.
          F = augmentdigits(N,NDig);
        else
          F = N;
        end
        return
      end
      
      % The decimal base used will come always from the default value.
      % It cannot be overridden in the call to hpf, but only changed by
      % adjustdecimalbase.
      F.DecimalBase = DefaultDecimalBase;
      F.Base = 10.^F.DecimalBase;
      
      % The default number of decimal digits to be used is given by
      % the function DefaultNumberOfDigits. The default value there
      % will be 64 digits.
      if nargin < 2
        NDig = DefaultNumberOfDigits;
        F.NumberOfDigits = NDig;
      elseif (numel(NDig) == 1)
        % only one number was supplied, so it must be the number of
        % reported digits.
        if (NDig < 1) || (NDig ~= round(NDig))
          error('HPF:HPF:improperdigitrequest', ...
            'The number of digits must be integer, with at least 1 digit required')
        end
        
        % get the number of shadow digits from the current stored
        % default
        defaultdigits = DefaultNumberOfDigits;
        NDig = [NDig,defaultdigits(2)];
        F.NumberOfDigits = NDig;
        
      elseif (numel(NDig) == 2)
        if (NDig(1) < 1) || (NDig(2) < 0) || any(NDig ~= round(NDig))
          error('HPF:HPF:improperdigitrequest', ...
            'The number of digits must be positive integer')
        end
        
        F.NumberOfDigits = NDig;
      elseif (numel(NDig) > 2)
        error('HPF:HPF:improperdigitrequest', ...
          'The number of digits specified must be a scalar or a vector of length 2')
      end
      
      NDigLimit = [3.6e14 3.6e12 3.6e10 3.6e8 3.6e6 36000];
      if NDigLimit(F.DecimalBase) <= sum(F.NumberOfDigits)
        switch F.DecimalBase
          case 1
            error('HPF:hpf', ...
              'Cannot exceed 3.6e14 total digits for a DecimalBase of 1')
          case 2
            error('HPF:hpf', ...
              'Cannot exceed 3.6e12 total digits for a DecimalBase of 2')
          case 3
            error('HPF:hpf', ...
              'Cannot exceed 3.6e10 total digits for a DecimalBase of 3')
          case 4
            error('HPF:hpf', ...
              'Cannot exceed 3.6e8 total digits for a DecimalBase of 4')
          case 5
            error('HPF:hpf', ...
              'Cannot exceed 3.6e6 total digits for a DecimalBase of 5')
          case 6
            error('HPF:hpf', ...
              'Cannot exceed 36000 total digits for a DecimalBase of 6')
        end
      end
      
      res = mod(sum(F.NumberOfDigits),F.DecimalBase);
      if res > 0
        % we need to augment the number of shadow digits
        F.NumberOfDigits(2) = F.NumberOfDigits(2) + (F.DecimalBase - res);
        NDig = F.NumberOfDigits;
      end
      
      % Initialize the remaining properties
      F.Numeric = 0;
      F.Sign = 0;
      F.Exponent = 0;
      F.Migits = zeros(1,ceil(sum(NDig)./F.DecimalBase));
      
      if isempty(N)
        % return an empty HPF object.
        F(1) = [];
        F = reshape(F,size(N));
        return
      elseif isnumeric(N) && (numel(N) > 1)
        % F will be an HPF array of the same shape as N
        F = repmat(F,size(N));
        for i = 1:numel(N)
          F(i) = hpf(N(i),NDig);
        end
        return
      end
      
      % we can trip out now if the number was identically zero
      % or some special value
      if isequal(N,0) || (nargin == 0)
        return
      elseif isnan(N)
        % a Nan has a flag of its own
        F.Numeric = NaN;
        F.Sign = 1;
        return
      elseif isinf(N)
        % infs also go here
        F.Numeric = abs(N);
        F.Sign = sign(N);
        return
      end
      
      % HPF does not admit complex numbers at this time
      if ~isreal(N) && ~isa(N,'vpi')
        error('HPF:complex','HPF does not currently allow complex numbers')
      end
      
      % a zero hpf to start with
      Nclass = class(N);
      switch Nclass
        case {'single', 'double'}
          % extract the hex digits for N from the ieee representation
          
          % nans and infs are already gone
          
          % in the case of a single, just convert it to a double.
          % that will be exact, with no loss of precision, and
          % cause no significant additional work later on. I
          % MAY fix this later on to work directly on singles
          % from the ieee form, but the conversion to double and
          % then the rest of the processing is quite fast anyway.
          if isa(N,'single')
            N = double(N);
          end
          
          % convert the double precision number to an exact
          % hex form. remember that the special case of 0
          % has already been removed.
          hex = sprintf('%bx',N);
          
          % convert the hex digits into decimal form
          dec = hex2dec(hex');
          
          % now convert to binary form, 4 bits for each hex digit
          bin = dec2bin(dec,4);
          
          % extract the sign. we have already ruled out 0
          F.Sign = 1 - 2*(bin(1) == '1');
          
          % and convert to a single binary string
          bitstr = reshape(bin',[1 64]);
          
          % the exponent here is the power of 2 to multiply by...
          Expon2 = bin2dec(bitstr(2:12)) - 1023;
          
          % the non-zero bits of the mantissa are...
          Mantissabits = find(bitstr(64:-1:13) == '1');
          
          % we can write N in the form
          %
          %   N = Sign*(1+Mantissa*2^(-52))*(2^Expon2)
          % or,
          %   N = Sign*(2^52+Mantissa)*(2^(Expon2 - 52))
          %
          % This means we can write the number as the sum of
          % this list of powers of 2 (ignoring the sign)
          pow2list = [Mantissabits-53 , 0] + Expon2;
          
          % pow2list will be non-empty, and sorted in increasing order.
          [D,exponent] = powerssum2dec(pow2list,sum(F.NumberOfDigits));
          F.Exponent = exponent + 1;
          
          % convert to migits
          F.Migits = d2m(D,F.DecimalBase);
          
        case {'uint8', 'uint16', 'uint32'}
          % convert uint8 numbers to hpf form
          % Sign = 1 by default here, since 0 has already been removed
          F.Sign = 1;
          D = dec2base(N,10) - '0';
          F.Exponent = numel(D);
          D = [D,wrepvec(0,sum(F.NumberOfDigits)-numel(D))];
          F.Migits = d2m(D,F.DecimalBase);
          
        case {'int8', 'int16', 'int32'}
          % its integer, but watch out for negatives
          F.Sign = double(sign(N));
          N = abs(N);
          
          % dec2base does all the work here
          D = dec2base(N,10) - '0';
          F.Exponent = numel(D);
          D = [D,wrepvec(0,sum(F.NumberOfDigits)-numel(D))];
          F.Migits = d2m(D,F.DecimalBase);
          
        case 'logical'
          % convert logical numbers to hpf form
          % In this case, the Sign bit is the same as the number
          % itself.
          if N
            F.Migits(1) = F.Base/10;
            F.Sign = 1;
            F.Exponent = 1;
          else
            F.Migits(1) = 0;
            F.Sign = 0;
            F.Exponent = 0;
          end
          
        case {'int64' 'uint64'}
          % Convert int64 and uint64 numbers into hpf format.
          
          % We need to do this the hard way, taking the number
          % apart one digit at a time. first extract the sign.
          F.Sign = double(sign(N));
          N = abs(N);
          
          % we already know the number is not 0
          D = [];
          F.Exponent = 0;
          while N ~= 0
            lastdigit = mod(N,10);
            N = (N - lastdigit)/10;
            % this will not grow for more than 16 digits or so,
            % so not a preallocation problem
            D = [lastdigit,D]; %#ok
            F.Exponent = F.Exponent + 1;
          end
          D = [D,wrepvec(0,sum(F.NumberOfDigits)-numel(D))];
          F.Migits = d2m(D,F.DecimalBase);
          
        case 'vpi'
          % converting a vpi to a hpf number is easy
          D = digits(N);
          
          if (numel(D) == 1) && (D == 0)
            % this was a vpi zero
            return
          end
          
          % get the sign from vpi.sign
          if sign(N) < 0
            F.Sign = -1;
          else
            F.Sign = 1;
          end
          
          % The exponent is easy, since this is an integer
          F.Exponent = numel(D);
                    
          % We may end up truncating the digits or appending
          % zeros to the right of the decimal point. All depends
          % on the number of digits the HPF number carries.
          D = padz(D,sum(F.NumberOfDigits),'right');
          
          % the mantissa is easy too
          F.Migits = d2m(D,F.DecimalBase);
          
        case 'char'
          % the number was represented as a character. There are a few
          % special cases, 'pi' being the first.
          
          % just in case, strip out all leading and trailing blanks
          N = strtrim(N);
          
          switch lower(N)
            case 'pi'
              % we need to bring in the digits of pi. this is not the
              % same as if pi itself had been passed in, as then MATLAB
              % has already converted it to a double.
              if isempty(specialnumbers)
                specialnumbers = load('_special_numbers_.mat');
              end
              
              if sum(F.NumberOfDigits) > numel(specialnumbers.pidigits)
                warning('HPF:HPF:TooManyDigitsRequested', ...
                  ['A maximum of ',num2str(numel(specialnumbers.pidigits)),' digits are available for pi'])
                F.NumberOfDigits(1) = numel(specialnumbers.pidigits) - F.NumberOfDigits(2);
              end
              
              % pi is positive
              F.Sign = 1;
              
              % and since it will be stored as .314159...
              F.Exponent = 1;
              
              % convert to migits
              F.Migits = d2m(specialnumbers.pidigits(1:sum(F.NumberOfDigits)),F.DecimalBase);
              
            case 'e'
              % the digits of e - another number worth storing
              if isempty(specialnumbers)
                specialnumbers = load('_special_numbers_');
              end
              if sum(F.NumberOfDigits) > numel(specialnumbers.edigits)
                warning('HPF:HPF:TooManyDigitsRequested', ...
                  ['A maximum of ',num2str(numel(specialnumbers.edigits)),' digits are available for e'])
                F.NumberOfDigits(1) = numel(specialnumbers.edigits) - F.NumberOfDigits(2);
              end

              % e is positive
              F.Sign = 1;
              
              % and since it will be stored as .2718...
              F.Exponent = 1;
              
              % convert to migits
              F.Migits = d2m(specialnumbers.edigits(1:sum(F.NumberOfDigits)),F.DecimalBase);
              
            case 'ln10'
              % the digits of log(10), the natural log of 10.
              if isempty(specialnumbers)
                specialnumbers = load('_special_numbers_');
              end
              
              if sum(F.NumberOfDigits) > numel(specialnumbers.ln10digits)
                warning('HPF:HPF:TooManyDigitsRequested', ...
                  ['A maximum of ',num2str(numel(specialnumbers.ln10digits)),' digits are available for log(10)'])
                F.NumberOfDigits(1) = numel(specialnumbers.ln10digits) - F.NumberOfDigits(2);
              end
              
              % log(10) is positive
              F.Sign = 1;
              
              % and since it will be stored as ..2308...
              F.Exponent = 1;
              
              % convert to migits
              F.Migits = d2m(specialnumbers.ln10digits(1:sum(F.NumberOfDigits)),F.DecimalBase);
              
            case {'1' 'one'}
              % one is the loneliest number
              F.Sign = 1;
              F.Migits(1) = F.Base/10;
              F.Exponent = 1;
              
            case {'2' 'two'}
              % two can be as bad as one
              F.Sign = 1;
              F.Migits(1) = 2*F.Base/10;
              F.Exponent = 1;
              
            case {'3' 'three'}
              % 's company
              F.Sign = 1;
              F.Migits(1) = 3*F.Base/10;
              F.Exponent = 1;
              
            case {'0' 'zero'}
              % zero (is already covered, so just drop out here)
              
            case 'nan'
              % a NaN
              F.Numeric = NaN;
              % set the sign flag to 1, since we test for zero against
              % the sign flag.
              F.Sign = 1;
              
            case '-inf'
              % an inf
              F.Numeric = inf;
              F.Sign = -1;
              
            case 'inf'
              % an inf
              F.Numeric = inf;
              F.Sign = 1;
              
            otherwise
              % the number should be represented as a character string
              
              % make it all lower case.
              N = lower(N);
              
              % Determine if the number is a valid scalar form
              % both sscanf and str2num must succeed here and result
              % in a scalar if the result is deemed valid.
              doubleres1 = sscanf(N,'%g');
              doubleres2 = str2num(N); %#ok
              if (numel(doubleres1) ~= 1) || (numel(doubleres2) ~= 1) || ...
                  ~all(ismember(N,'0123456789+-ed.')) || (sum(ismember(N,'ed')) > 1)
                error('HPF:HPF:characterinput', ...
                  'N must be a valid scalar numerical form if supplied as character')
              end
              
              % strip out blanks
              N(N == ' ') = '';
              
              % replace d characters with e
              N(N == 'd') = 'e';
              
              % is the first character a '-' ?
              F.Sign = 1;
              if N(1) == '-'
                F.Sign = -1;
                N(1) = '';
              elseif N(1) == '+'
                N(1) = '';
              end
              
              % we now know the sign, if any was given.
              % filter out the mantissa next
              num1 = regexp(N,'\d+[.]\d+','match');
              num2 = regexp(N,'[.]\d+[e]','match');
              num3 = regexp(N,'\d+[.e]','match');
              num4 = regexp(N,'[.]\d+','match');
              num5 = regexp(N,'\d+','match');
              if ~isempty(num1)
                % we have a mantissa of the form nnn.nnn
                % where is the decimal point?
                num1 = num1{1};
                k = find(num1 == '.');
                F.Exponent = k - 2;
                num1(k) = '';
                
                % how many leading zeros were there?
                [D,nlz] = parsenumericdigits(num1,0);
                F.Exponent = F.Exponent - nlz + 1;
                
                % do we need to pad zeros onto the end?
                % I should probably round the trailing digits if too many
                % were supplied. padz will simply do a truncate.
                D = padz(D,sum(F.NumberOfDigits),'right');
                
                % store as migits
                F.Migits = d2m(D,F.DecimalBase);
                if all(F.Migits == 0)
                  F.Sign = 0;
                  F.Exponent = 0;
                end
                
              elseif ~isempty(num2)
                % of the form .nnne
                % strip off the trailing '.' and 'e'
                num2 = regexp(num2{1},'\d+','match');
                num2 = num2{1};
                [D,nlz] = parsenumericdigits(num2,0);
                
                % We do care about leading zeros in this case,
                % Strip them off, adjusting the exponent.
                % Trailing zeros are unimportant.
                F.Exponent = F.Exponent - nlz;
                k = find(D ~= 0,'1','first');
                D = D(k:end);
                
                D = padz(D,sum(F.NumberOfDigits),'right');
                F.Migits = d2m(D,F.DecimalBase);
                if all(D == 0)
                  F.Sign = 0;
                  F.Exponent = 0;
                end
                
              elseif ~isempty(num3)
                % of the form nnn. or nnne or nnn.e
                % strip off the trailing '.' or 'e'
                num3 = regexp(num3{1},'\d+','match');
                num3 = num3{1};
                D = parsenumericdigits(num3,0);
                
                % We don't care about leading zeros in this case,
                % Strip them off, but don't adjust the exponent.
                % Trailing zeros are important, as they add to the
                % exponent.
                k = find(D ~= 0,'1','last');
                if k < numel(D)
                  F.Exponent = F.Exponent + (numel(D) - k) + 1;
                else
                  F.Exponent = 1;
                end
                
                D = padz(D,sum(F.NumberOfDigits),'right');
                F.Migits = d2m(D,F.DecimalBase);
                if all(D == 0)
                  F.Sign = 0;
                  F.Exponent = 0;
                end
              elseif ~isempty(num4)
                % of the form .nnn
                % strip off the leading '.'
                num4 = num4{1};
                num4(1) = '';
                F.Exponent = -1;
                % how many leading zeros were there?
                [D,nlz] = parsenumericdigits(num4,0);
                D = padz(D,sum(F.NumberOfDigits),'right');
                F.Migits = d2m(D,F.DecimalBase);
                F.Exponent = F.Exponent - nlz + 1;
                if all(D == 0)
                  F.Sign = 0;
                  F.Exponent = 0;
                end
              elseif ~isempty(num5)
                % Of the form nnnnnn , so an integer part, with no decimal point
                % strip off the leading zeros, but we can ignore them
                D = parsenumericdigits(num5{1},0);
                F.Exponent = numel(D);
                D = padz(D,sum(F.NumberOfDigits),'right');
                F.Migits = d2m(D,F.DecimalBase);
                if all(D == 0)
                  F.Sign = 0;
                  F.Exponent = 0;
                end
              end
              
              % next, look for a power of 10 on there. It will
              % be of one of the general forms e-nnn or e+nnn or ennn.
              % this regular expression will catch an e folowed by any
              % set of characters.
              num = regexp(N,'[e].+','match');
              if ~isempty(num)
                % we have an exponent of the desired form
                num = num{1};
                
                % strip off the e
                num(1) = '';
                
                exponentsign = 1;
                % is the next character a + or a -?
                if num(1) == '-'
                  exponentsign = -1;
                  num(1) = '';
                elseif num(1) == '+'
                  num(1) = '';
                end
                
                % The remaining portion of num must be pure numeric,
                % as an integer. Just convert it to a double
                F.Exponent = F.Exponent + exponentsign*str2double(num);
                if (round(F.Exponent) ~= F.Exponent)
                  error('HPF:HPF:characterinput', ...
                    'N must be a valid scalar numerical form if supplied as character')
                elseif (F.Exponent >= 2e53)
                  error('HPF:HPF:characterinput', ...
                    'Sorry, but this exponent is itself too large to be represented as a double.')
                end
              end
              
              % make sure that for a zero, we did not get a non-zero
              % exponent
              if F.Sign == 0
                F.Exponent = 0;
              end
              
          end % switch lower(N)
      end % switch Nclass
      
    end % F = hpf(...)
    
    function F = int8(X)
      % Convert an hpf number into its int8 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = int8(int64(X));
      
    end % F = int8(X)
    
    function F = int16(X)
      % Convert an hpf number into its int16 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = int16(int64(X));
      
    end % F = int32(X)
    
    function F = int32(X)
      % Convert an hpf number into its int32 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = int32(int64(X));
      
    end % F = int32(X)
    
    function F = int64(X)
      % Convert an hpf number into its int64 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      if isempty(X)
        % empty begets empty
        F = [];
      elseif numel(X) > 1
        % F is an array
        F = zeros(size(X));
        for i = 1:numel(X)
          F(i) = int64(X(i));
        end
      else
        % X is scalar. any special cases?
        if ~isfinite(X.Numeric)
          % inf and -inf will map to +/- 2^63. NaN will map to 0,
          % which is consistent with MATLAB, as int64(NaN) = 0.
          F = int64(X.Numeric);
        elseif X.Sign == 0
          F = int64(0);
        elseif X.Exponent > 19
          % we can be sure this is an overflow.
          F = int64(inf)*X.Sign;
        else
          % something in the middle, although it may still
          % be an overflow for int64
          X = floor(X);
          
          db = X.DecimalBase;
          E = X.Exponent;
          base = int64(X.Base);
          migits = int64(X.Migits);
          mlast = find(migits ~= 0,1,'last');
          F = int64(0);
          for i = 1:(mlast-1)
            F = F.*base + migits(i);
            E = E - db;
          end
          
          % extract the high non-zero digits from the last non-zero migit
          Dlast = m2d(double(migits(mlast)),db);
          for i = 1:db
            F = F.*int64(10) + Dlast(i);
            E = E - 1;
            if E <= 0
              break
            end
          end
          
          if E > 0
            % there were still a few powers of 10 in there
            F = F.*int64(10).^E;
          end
          
          if X.Sign < 0
            F = -F;
          end
        end
      end
    end % F = int64(X)
    
    function Boolean = isinf(F)
      % test for an inf hpf number
      % Boolean = isinf(F)
      %
      %  See also: isnan, isfinite
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      Boolean = false(size(F));
      for i = 1:numel(F)
        Boolean(i) = isinf(F(i).Numeric);
      end
    end % F = isinf(F)
    
    function Boolean = isfinite(F)
      % test for a finite hpf number
      % Boolean = isfinite(F)
      %
      %  See also: isinf, isnan
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      Boolean = true(size(F));
      for i = 1:numel(F)
        Boolean(i) = isfinite(F(i).Numeric);
      end
    end % F = isfinite(F)
    
    function Boolean = isnan(F)
      % test for a NaN hpf number
      % Boolean = isnan(F)
      %
      %  See also: isinf, isfinite
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      Boolean = false(size(F));
      for i = 1:numel(F)
        Boolean(i) = isnan(F(i).Numeric);
      end
    end % F = isnan(F)
    
    function Boolean = isnumeric(F)
      % test that returns true for an HPF number
      % Boolean = isnumeric(F)
      %
      %  See also: isinf, isfinite
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      
      % trivially true if we got into this routine
      Boolean = true;
      
    end % Boolean = isnumeric(F)
    
    function Boolean = iszero(F)
      % quick test for a zero hpf number
      % Boolean = iszero(F)
      %
      %  See also: isnan, isfinite, isinf
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      Boolean = false(size(F));
      for i = 1:numel(F)
        Boolean(i) = F(i).Sign == 0;
      end
    end % F = iszero(F)
    
    function result = le(x,y)
      % Test for numerical inequality (less than or equal to) between numbers, returning a boolean result
      % usage: result = le(x,y);
      % usage: result = x <= y;
      
      % is either x or y an array? If so, we need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans, nor can they be compared to other numbers.
        if isnan(x) || isnan(y)
          result = false;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        if ~isa(x,'hpf')
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          y = hpf(y,x.NumberOfDigits);
        end
        
        % do they have the same decimal base? choose the smaller one
        % if they differ.
        if x.DecimalBase ~= y.DecimalBase
          if x.DecimalBase < y.DecimalBase
            y = adjustdecimalbase(y,x.DecimalBase);
          else
            x = adjustdecimalbase(x,y.DecimalBase);
          end
        end
          
        % augment a shorter number so that both are the same lengths
        % spare zero digits carried are not relevant to equality.
        if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
          x = augmentdigits(x,y.NumberOfDigits);
          elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
            y = augmentdigits(y,x.NumberOfDigits);
        end
        
        % was one or both numbers an inf?
        if isinf(x.Numeric) || isinf(y.Numeric)
          % this works in any case
          result = (x.Numeric*x.Sign <= y.Numeric*y.Sign);
        elseif x.Sign < y.Sign
          % if there is a sign differential, then we are done
          result = true;
        elseif x.Sign > y.Sign
          % if there is a sign differential, then we are done
          result = false;
        elseif x.Sign == 0
          % if we skip past the previous tests, then the signs must be
          % equal. If they are both zero, then we are done.
          result = true;
        elseif x.Sign > 0
          % positive x and y
          if x.Exponent > y.Exponent
            result = false;
          elseif x.Exponent < y.Exponent
            result = true;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = true;
            elseif x.Migits(k) > y.Migits(k)
              result = false;
            else
              result = true;
            end
          end
        else
          % negative x and y
          if x.Exponent > y.Exponent
            result = true;
          elseif x.Exponent < y.Exponent
            result = false;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = true;
            elseif x.Migits(k) > y.Migits(k)
              result = true;
            else
              result = false;
            end
          end
        end
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = false(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = le(x,y(i));
        end
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = false(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = le(x(i),y);
        end
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:LE:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = le(x(i),y(i));
        end
      end
      
    end % result = le(x,y)
    
    function y = linspace(d1, d2, n)
      % LINSPACE Linearly spaced vector.
      %   LINSPACE(X1, X2) generates a row vector of 100 linearly
      %   equally spaced points between X1 and X2.
      %
      %   LINSPACE(X1, X2, N) generates N points between X1 and X2.
      %   For N = 1, LINSPACE returns X2.
      %
      %   See also LOGSPACE, COLON.
      
      if nargin == 2
        n = 100;
      end
      n = double(n);
      
      if ~isa(d1,'hpf')
        d1 = hpf(d1,d2.NumberOfDigits);
      elseif ~isa(d2,'hpf')
        d2 = hpf(d2,d1.NumberOfDigits);
      end
      
      if n < 2
        y = zeros(1,floor(n)) + d2;
      else
        % at least two end points
        n1 = floor(n)-1;
        c = (d2 - d1).*(n1-1); % opposite signs may cause overflow
        if isinf(c)
          y = d1 + (d2./n1).*(0:n1) - (d1./n1).*(0:n1);
        else
          y = d1 + (0:n1).*((d2 - d1)./n1);
        end
        y(1) = d1;
        y(end) = d2;
      end
      
    end % y = linspace(d1, d2, n)
    
    function F = log(X)
      % evaluates log(X) (the natural log of x) for a hpf number X
      
      % Is this an array?
      if isempty(X)
        % empty begets empty
        F = X;
        return
      elseif numel(X) > 1
        F = X;
        for i = 1:numel(X)
          F(i) = log(X(i));
        end
        return
      end
      % if we drop through, then X was scalar
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases. was X zero, negative or inf? If so, we can trip out
      if X.Sign == 0
        F = hpf('-inf',NDig);
        return
      elseif (X.Sign < 0) || isnan(X.Numeric)
        F = hpf('NaN',NDig);
        return
      elseif isinf(X.Numeric)
        F = hpf('inf',NDig);
        return
      end
      
      % Bring in log(10) for the necessary number of digits
      ln10 = hpf('ln10',NDig);
      
      % we can get the first part of our result from the exponent
      % on 10 in our hpf number X
      F = X.Exponent*ln10;
      
      % The remainder comes from the Mantissa. X will now live in the
      % bounded interval [1, 10)
      X.Exponent = 0;
      
      % Now we use a simple trick. If we knew a number VERY close to
      % log(X), then we could shift the problem so the log series I will
      % use will become very efficient. So convert to this remainder in
      % X to a double precision number.
      dx = double(X);
      
      % If X is identically 1, then there will be a subtle problem,
      % causing the loop at the end to be infinite. catch those cases
      if dx == 1
        dx = dx + eps;
      end
      
      % compute the log as a double, but then revert to hpf, with the
      % requisite number of digits.
      lndx = hpf(log(dx),NDig);
      
      % and exponentiate lndx. negate lndx to avoid a later divide.
      dxinv = exp(uminus(lndx));
      
      % shift X by the exact log we just computed. X will be VERY close to
      % 1 at this point, in fact, within about 16 digits of accuracy of
      % being exactly 1.
      X = X.*dxinv;
      
      % as well as our current approximation for F ~ log(X)
      F = F + lndx;
      
      % transform X, so that I can use a more rapidly convergent series
      % than the standard series for log(X).
      %
      %   Y = (X-1)/(X+1)
      %
      % which is equivalent to
      %   X = (1+Y)/1-Y)
      %
      % This transformation provides roughly 32 additional digits of
      % precision per term of the series that follows. We could in fact
      % use a second iteration of this trick to gain hundreds or even
      % thousands of correct digits. This is actually how I computed
      % the 500,000 digits of log(10) used in several places in hpf.
      %
      % http://en.wikipedia.org/wiki/Natural_logarithm
      one = hpf('1',NDig);
      Y = (X-one)/(X+one);
      
      % is Y identically zero? If so, then we are done.
      if Y.Sign == 0
        return
      end
      
      Ysq = Y.*Y;
      Fy = one;
      term = one;
      n = 1;
      while (-term.Exponent) < sum(NDig)
        term = term.*Ysq;
        n = n + 2;
        % there is a divide I cannot avoid inside this loop
        Fy = Fy + term./n;
      end
      F = F + 2.*Y.*Fy;
      
    end % function F = log(X)
    
    function F = log10(X)
      % evaluates log10(X) (the common log of x) for a hpf number X
      
      % Is this an array?
      if isempty(X)
        % empty begets empty
        F = X;
        return
      elseif numel(X) > 1
        F = X;
        for i = 1:numel(X)
          F(i) = log10(X(i));
        end
        return
      end
      % if we drop through, then X was scalar
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % strip out the powers of 10
      F = hpf(X.Exponent,NDig);
      X.Exponent = 0;
      % then compute the log10 of the mantissa.
      F = F + log(X)./hpf('ln10',NDig);
      
    end % function F = log10(X)
    
    function F = log2(X)
      % evaluates log2(X) (log to the base 2 of x) for a hpf number X
      
      % Is this an array?
      if isempty(X)
        % empty begets empty
        F = X;
        return
      elseif numel(X) > 1
        F = X;
        for i = 1:numel(X)
          F(i) = log2(X(i));
        end
        return
      end
      % if we drop through, then X was scalar
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      if isnan(X.Numeric)
        % NaN begets NaN
        F = hpf('NaN',NDig);
      elseif isinf(X.Numeric) && (X.Sign < 0)
        % -inf --> 0
        F = hpf('0',NDig);
      elseif isinf(X.Numeric) && (X.Sign > 0)
        % inf --> inf
        F = hpf('inf',NDig);
      elseif (X.Sign == 0)
        % 0 --> 1
        F = hpf('1',NDig);
      else
        % simply a re-scaling of the natural log
        F = log(X)./log(hpf(uint8(2),NDig));
      end
      
    end % function F = log2(X)
    
    function result = lt(x,y)
      % Test for strict numerical inequality (less than) between numbers, returning a boolean result
      % usage: result = lt(x,y);
      % usage: result = x < y;
      
      % is either x or y an array? If so, we need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans, nor can they be compared to other numbers.
        if isnan(x) || isnan(y)
          result = false;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        if ~isa(x,'hpf')
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          y = hpf(y,x.NumberOfDigits);
        end
        
        % do they have the same decimal base? choose the smaller one
        % if they differ.
        if x.DecimalBase ~= y.DecimalBase
          if x.DecimalBase < y.DecimalBase
            y = adjustdecimalbase(y,x.DecimalBase);
          else
            x = adjustdecimalbase(x,y.DecimalBase);
          end
        end
        
        % augment a shorter number so that both are the same lengths
        % spare zero digits carried are not relevant to equality.
        if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
          x = augmentdigits(x,y.NumberOfDigits);
        elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
          y = augmentdigits(y,x.NumberOfDigits);
        end
        
        % was one or both numbers an inf?
        if isinf(x.Numeric) || isinf(y.Numeric)
          % this works in any case
          result = (x.Numeric*x.Sign < y.Numeric*y.Sign);
        elseif x.Sign < y.Sign
          % if there is a sign differential, then we are done
          result = true;
        elseif x.Sign > y.Sign
          % if there is a sign differential, then we are done
          result = false;
        elseif x.Sign == 0
          % if we skip past the previous tests, then the signs must be
          % equal. If they are both zero, then we are done.
          result = false;
        elseif x.Sign > 0
          % positive x and y
          if x.Exponent > y.Exponent
            result = false;
          elseif x.Exponent < y.Exponent
            result = true;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = false;
            elseif x.Migits(k) > y.Migits(k)
              result = false;
            else
              result = true;
            end
          end
        else
          % negative x and y
          if x.Exponent > y.Exponent
            result = true;
          elseif x.Exponent < y.Exponent
            result = false;
          else
            % equal exponents, so now we need to test the migits
            k = find(x.Migits ~= y.Migits,1,'first');
            if isempty(k)
              % all migits were identical
              result = false;
            elseif x.Migits(k) > y.Migits(k)
              result = true;
            else
              result = false;
            end
          end
        end
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = false(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = lt(x,y(i));
        end
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = false(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = lt(x(i),y);
        end
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:LT:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = lt(x(i),y(i));
        end
      end
      
    end % result = lt(x,y)
    
    function digits = mantissa(X)
      % Extract the mantissa itself as a string of decimal digits
      if X.DecimalBase == 1
        % need do nothing special here, as migits and digits are
        % synonymous for this base.
        digits = X.Migits;
      else
        % don't forget that we will need any leading zeros
        digits = dec2base(X.Migits,10,X.DecimalBase) - '0';
        
        % return a row vector of pure digits
        digits = reshape(digits.',1,[]);
      end
      
    end % digits = mantissa(X)
    
    function [res,ind] = max(X,Y,dim)
      % Maximum of an array or between a pair of arrays (for HPF numbers)
      %
      % For vectors, MAX(X) is the largest element in X. For matrices,
      %    MAX(X) is a row vector containing the maximum element from each
      %    column. For N-D arrays, max(X) operates along the first
      %    non-singleton dimension.
      %
      %    [Y,I] = MAX(X) returns the indices of the maximum values in vector I.
      %    If the values along the first non-singleton dimension contain more
      %    than one maximal element, the index of the first one is returned.
      %
      %    MAX(X,Y) returns an array the same size as X and Y with the
      %    largest elements taken from X or Y. Either one can be a scalar.
      %
      %    [Y,I] = MAX(X,[],dim) operates along the dimension dim.
      
      % which mode are we in?
      if nargin == 0
        error('HPF:MAX:improperarguments','Insufficient input arguments')
      elseif nargin > 3
        error('HPF:MAX:improperarguments','No more than 3 input arguments are allowed')
      elseif nargin == 2
        if numel(X) == 1
          % scalar expansion for X
          res = hpf(Y);
          res(X > Y) = X;
        elseif numel(Y) == 1
          % scalar expansion for Y
          res = hpf(X);
          res(Y > X) = Y;
        elseif ~isequal(size(X),size(Y))
          % sizes are inconsistent
          error('HPF:MAX:unmatchedsizes','Matrix dimensions must agree')
        else
          % two compatibly sized arrays
          res = hpf(X);
          ind = Y > X;
          res(ind) = Y(ind);
        end
      else
        % single argument, or 3 arguments
        
        if nargin ==1
          % determine the dimension to work on, as the first non-singleton
          % dimension
          dim = find(size(X) > 1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        % dim is now defined
        
        % check for problems with dim
        Sx = size(X);
        if (dim < 0) || (dim ~= round(dim))
          error('HPF:MAX:improperdim','dim must be a positive integer')
        elseif (dim > numel(Sx)) || (Sx(dim) == 1)
          % this is a no-op
          res = X;
          ind = ones(Sx);
          return
        else
          % compute the max along dimension dim
          
          % create result and ind as the correct sizes
          Sr = Sx;
          Sr(dim) = 1;
          res = repmat(X(1),Sr);
          ind = ones(Sr);
          
          % permute and reshape X into a 2-d array
          Ix = 1:numel(Sx);
          Ix(dim) = [];
          Ix = [dim,Ix];
          X = permute(X,Ix);
          X = reshape(X,Sx(dim),[]);
          
          % now loop over the columns of X
          ndim = Sx(dim);
          for cind = 1:size(X,2)
            res(cind) = X(1,cind);
            for rind = 2:ndim
              if res(cind) < X(rind,cind)
                ind(cind) = rind;
                res(cind) = X(rind,cind);
              end
            end
          end
        end
      end
    end
    
    function [res,ind] = min(X,Y,dim)
      % Maximum of an array or between a pair of arrays (for HPF numbers)
      %
      % For vectors, MIN(X) is the largest element in X. For matrices,
      %    MIN(X) is a row vector containing the minimum element from each
      %    column. For N-D arrays, MIN(X) operates along the first
      %    non-singleton dimension.
      %
      %    [Y,I] = MIN(X) returns the indices of the minimum values in vector I.
      %    If the values along the first non-singleton dimension contain more
      %    than one minimal element, the index of the first one is returned.
      %
      %    MIN(X,Y) returns an array the same size as X and Y with the
      %    largest elements taken from X or Y. Either one can be a scalar.
      %
      %    [Y,I] = MIN(X,[],dim) operates along the dimension dim.
      
      % which mode are we in?
      if nargin == 0
        error('HPF:MIN:improperarguments','Insufficient input arguments')
      elseif nargin > 3
        error('HPF:MIN:improperarguments','No more than 3 input arguments are allowed')
      elseif nargin == 2
        if numel(X) == 1
          % scalar expansion for X
          res = Y;
          res(X < Y) = X;
        elseif numel(Y) == 1
          % scalar expansion for Y
          res = X;
          res(Y < X) = Y;
        elseif ~isequal(size(X),size(Y))
          % sizes are inconsistent
          error('HPF:MIN:unmatchedsizes','Matrix dimensions must agree')
        else
          % two compatibly sized arrays
          res = X;
          i = Y < X;
          res(i) = Y(i);
        end
      else
        % single argument, or 3 arguments
        
        if nargin ==1
          % determine the dimension to work on, as the first non-singleton
          % dimension
          dim = find(size(X) > 1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        % dim is now defined
        
        % check for problems with dim
        Sx = size(X);
        if (dim < 0) || (dim ~= round(dim))
          error('HPF:MIN:improperdim','dim must be a positive integer')
        elseif (dim > numel(Sx)) || (Sx(dim) == 1)
          % this is a no-op
          res = X;
          ind = ones(Sx);
          return
        else
          % compute the min along dimension dim
          
          % create result and ind as the correct sizes
          Sr = Sx;
          Sr(dim) = 1;
          res = repmat(X(1),Sr);
          ind = ones(Sr);
          
          % permute and reshape X into a 2-d array
          Ix = 1:numel(Sx);
          Ix(dim) = [];
          Ix = [dim,Ix];
          X = permute(X,Ix);
          X = reshape(X,Sx(dim),[]);
          
          % now loop over the columns of X
          ndim = Sx(dim);
          for i = 1:size(X,2)
            res(i) = X(1,i);
            for j = 2:ndim
              if res(i) > X(j,i)
                ind(i) = j;
                res(i) = X(j,i);
              end
            end
          end
        end
      end
    end
    
    function result = minus(A,B)
      % Subtracts a hpf number from another, or a hpf and another numeric class
      % Simplest is just to change the sign on B, then do an add.
      %
      % This lets plus do all of the work, and avoids the need for
      % duplicated code. Note that if B is an integer form, like uint8,
      % this may yield strange results, but then why would you do that
      % operation anyway?
      result = A + uminus(B);
    end % result = minus(A,B)
    
    function result = mod(A,B)
      % Modulus (signed remainder after division)
      % usage: result = mod(A,B)
      Q = floor(A./B);
      result = A - Q.*B;
      
    end % result = mod(A,B)
    
    function result = mpower(x,y)
      % x^y in MATLAB when either of the operands is a hpf number.
      %
      % If x is a square array, then y must be a non-negative scalar
      % integer, no larger than 2^53 - 1 as a double, or 2^64-1 as a
      % uint64 or a larger integer if hpf or vpi.
      
      
      % at least one of these numbers must be a hpf, so if not num, then den
      % how many digits will we have in the result?
      if isa(x,'hpf')
        if isa(y,'hpf')
          % both are hpf
          NDig = combineNDig(x.NumberOfDigits,y.NumberOfDigits);
        else
          % only x is hpf
          NDig = x.NumberOfDigits;
        end
      else
        % only y is hpf
        NDig = y.NumberOfDigits;
        x = hpf(x,NDig);
      end
      
      % check for the special cases
      if (numel(x) == 1) && (numel(y) == 1)
        result = x.^y;
      elseif isempty(x) || isempty(y)
        % empty begets empty
        result = hpf([],NDig);
      elseif numel(y) > 1
        error('HPF:mpower:notimplemented', ...
          'Sorry, but mpower is not implemented where the exponent is an array')
      elseif (size(x,1) ~=size(x,2)) || (numel(size(x)) > 2)
        error('HPF:mpower:notimplemented', ...
          'Sorry, but mpower is not implemented where the base is not a 2-d square matrix')
      elseif (y ~= round(y)) || (y < 0)
        error('HPF:mpower:notimplemented', ...
          'Sorry, but mpower is only implemented where the exponent is a non-negative integer')
      elseif y == 0
        % an array to the zero power will be the identity matrix
        result = hpf(eye(size(x)),NDig);
      elseif y == 1
        % unit power, an identity op
        result = hpf(f,NDig);
      else
        % x must be a square array, at least 2x2, and y a scalar integer.
        
        % how large is y? can it be converted into a uint64 integer?
        if 0
          
          
          
          
          
          
          
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        end
        
      end
      
    end % function mpower
    
    function result = mrdivide(num,den)
      % Arithmetic operation num/den, where one or both are hpf
      % usage: result = mrdivide(num,den)
      % usage: result = num/den
      %
      % Where the two numbers are represented with a different number of
      % digits, use the smaller of the two lengths.
      
      % are either num or den vectors/arrays or empty?
      if isempty(num)
        if (numel(den) > 1) || (isempty(den) && ~isequal(size(num),size(den)))
          error('HPF:MRDIVIDE:impropersize','Matrix dimensions must agree.')
        else
          result = num;
        end
      elseif isempty(den)
        if (numel(num) > 1) || (isempty(num) && ~isequal(size(num),size(den)))
          error('HPF:MRDIVIDE:impropersize','Matrix dimensions must agree.')
        else
          result = den;
        end
      elseif (numel(num) >= 1) && (numel(den) == 1)
        % scalar expansion of den. ./ does it.
        result = num./den;
      elseif (numel(num) >= 1) && (numel(den) > 1)
        % both num and den are arrays, or only den
        error('HPF:MRDIVIDE:undefined','Sorry, but mrdivide is not yet implemented for array inputs.')
      end
      
    end % function mrdivide
    
    function result = mtimes(A,B)
      % Matrix multiplication of two hpf numbers or arrays of numbers,
      % or an hpf and another numeric class
      
      if (numel(A) == 1) || (numel(B) == 1)
        % propagate scalar multiplication
        result = A.*B;
      else
        % both A and B must have more then one element.
        % do they conform in size?
        sizea = size(A);
        sizeb = size(B);
        if (numel(sizea) > 2) || (numel(sizeb) > 2)
          error('HPF:MTIMES:noncomformingarrays', ...
            'Higher dimensional array multiplication undefined in MATLAB')
        elseif sizea(2) ~= sizeb(1)
          error('HPF:MTIMES:noncomformingarrays', ...
            'Matrix multiplication incompatibility - inner dimensions do not conform')
        end
        
        % they do conform
        % preallocate a hpf result. Don't worry about the
        % number of digits, times will take care of that.
        result = repmat(hpf('0'),[sizea(1),sizeb(2)]);
        for i = 1:sizea(1)
          for j = 1:sizeb(2)
            % plus will turn this into an hpf number,
            % also worrying about the digits.
            accum = 0;
            for k = 1:sizea(2)
              accum = accum + A(i,k).*B(k,j);
            end
            result(i,j) = accum;
          end
        end
      end
    end % result = mtimes(A,B)
    
    function result = ne(x,y)
      % Test for numerical inequality between numbers, returning a boolean result
      % usage: result = ne(x,y);
      % usage: result = x ~= y;
      
      % is either x or y an array? If so, we may need to do a scalar expansion
      % for the comparison.
      if isempty(x) || isempty(y)
        % empty begets empty
        result = [];
        
      elseif (numel(x) == 1) && (numel(y) == 1)
        % both are scalars
        
        % if either one is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x) || isnan(y)
          result = true;
          return
        end
        
        % We need to ensure that both are hpf numbers to compare them.
        % At least ONE of the numbers must have been an hpf though.
        if ~isa(x,'hpf')
          x = hpf(x,y.NumberOfDigits);
        end
        if ~isa(y,'hpf')
          y = hpf(y,x.NumberOfDigits);
        end
        
        % do they have the same decimal base? choose the smaller one
        % if they differ.
        if x.DecimalBase ~= y.DecimalBase
          if x.DecimalBase < y.DecimalBase
            y = adjustdecimalbase(y,x.DecimalBase);
          else
            x = adjustdecimalbase(x,y.DecimalBase);
          end
        end
        
        % augment a shorter number so that both are the same lengths
        % spare zero digits carried are not relevant to equality.
        if sum(x.NumberOfDigits) < sum(y.NumberOfDigits)
          x = augmentdigits(x,y.NumberOfDigits);
        elseif sum(x.NumberOfDigits) > sum(y.NumberOfDigits)
          y = augmentdigits(y,x.NumberOfDigits);
        end
        
        % now that we have made x and y comparable (and excluded the
        % NaNs, so the comparison itself is easy
        result = ~isequal(x,y);
        
      elseif (numel(x) == 1) && (numel(y) > 1)
        % scalar expansion for x
        
        % if x is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(x)
          result = true(size(y));
          return
        end
        
        result = true(size(y));
        for i = 1:numel(y)
          result(i) = ne(x,y(i));
        end
        
      elseif (numel(y) == 1) && (numel(x) > 1)
        % scalar expansion for y
        
        % if y is a nan, then we are done. nans equal nothing,
        % including other nans.
        if isnan(y)
          result = true(size(x));
          return
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = ne(x(i),y);
        end
        
      elseif (numel(y) > 1) && (numel(x) > 1)
        % both are arrays/vectors
        if ~isequal(size(x),size(y))
          error('HPF:NE:unmatchedarrays','x and y did not match in size/shape')
        end
        
        result = true(size(x));
        for i = 1:numel(x)
          result(i) = ne(x(i),y(i));
        end
      end
      
    end % result = ne(x,y)
    
    function result = nthroot(x,n)
      % Forms the nth root of x, for x an HPF number, and n a positive integer
      %
      % usage: result = nthroot(x,n)
      
      % x or n must be scalar, or both the same size
      if nargin < 2
        error('HPF:NTHROOT:incorrectarguments','exactly 2 arguments must be provided')
      end
      
      % n must be integer, and not too large, so make sure it is just a
      % double
      if isa(n,'hpf')
        NDig = n(1).NumberOfDigits;
        n = double(n);
      end
      
      % one or the other of x and n must have been an HPF number
      if ~isa(x,'hpf')
        x = hpf(x,NDig);
      end
      
      % if n was scalar, and x was not, make n the same size as x for
      % simplicity later
      if isscalar(n) && ~isscalar(x)
        n = repmat(n,size(x));
      end
      
      % how many digits?
      if isa(x,'hpf')
        NDig = x(1).NumberOfDigits;
      else
        % only n was an hpf
        NDig = n.NumberOfDigits;
        x = hpf(x,NDig);
      end
      
      % preallocate result to be the size of x (or n, if n was an array)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    end % function nthroot
    
    function result = plus(A,B)
      % Adds two HPF numbers, or a hpf and another numeric class
      % Where two HPF numbers are represented with a different number of
      % digits, use the smaller of the two lengths, as no significant
      % digits beyond that point will generally be meaningful.
      %
      % Adding a double or an integer of any class to an HPF number
      % converts that number to HPF, then does the addition. Where
      % the number in question is stored as a double, the conversion
      % will use the IEEE representation of the double to create the HPF
      % number.
      
      if ~isa(A,'hpf')
        % B must be a hpf
        A = hpf(A,B(1).NumberOfDigits);
      end
      if ~isa(B,'hpf')
        % A must be a hpf
        B = hpf(B,A(1).NumberOfDigits);
      end
      % both numbers must now be hpf numbers
      
      % are either A or B arrays of any size, or empty?
      if isempty(A)
        if (numel(B) > 1) || (isempty(B) && ~isequal(size(A),size(B)))
          error('HPF:PLUS:impropersize','Matrix dimensions must agree.')
        else
          result = A;
        end
        return
      elseif isempty(B)
        if (numel(A) > 1) || (isempty(A) && ~isequal(size(A),size(B)))
          error('HPF:PLUS:impropersize','Matrix dimensions must agree.')
        else
          result = B;
        end
        return
      end
      
      % final number of digits
      NDig = combineNDig(A(1).NumberOfDigits,B(1).NumberOfDigits);
        
      if (numel(A) > 1) && (numel(B) == 1)
        % scalar expansion of B
        result = A;
        for i = 1:numel(A)
          result(i) = result(i) + B;
        end
        return
      elseif (numel(B) > 1) && (numel(A) == 1)
        % scalar expansion of A
        result = B;
        for i = 1:numel(B)
          result(i) = result(i) + A;
        end
        return
      elseif (numel(A) > 1) && (numel(B) > 1)
        % both A and B are arrays
        if ~isequal(size(A),size(B))
          error('HPF:PLUS:impropersize','Matrix dimensions must agree.')
        end
        
        result = A;
        for i = 1:numel(B)
          result(i) = A(i) + B(i);
        end
        return
      end
      % if we drop down, then both A and B are scalars
      
      % do the numbers have the same decimal base? choose the smaller one
      % if they differ.
      if A.DecimalBase ~= B.DecimalBase
        if A.DecimalBase < B.DecimalBase
          B = adjustdecimalbase(B,A.DecimalBase);
        else
          A = adjustdecimalbase(A,B.DecimalBase);
        end
      end
      
      % check for nan or infs or a zero
      if isnan(A.Numeric) || isnan(B.Numeric)
        % the result will be a NaN, as NaN + anything is NaN
        result = hpf('nan',NDig);
        return
      elseif isinf(A.Numeric)
        % inf + anything is inf, unless B was also an inf, of opposite sign
        if isinf(B.Numeric)
          % both are infs. do they have the same signs?
          if (A.Sign*B.Sign) < 0
            % opposite signs, so a NaN results
            result = hpf('NaN',NDig);
            return
          else
            % same signs, so an inf results
            result = augmentdigits(A,NDig);
            return
          end
        else
          % inf + anything else is inf
          result = augmentdigits(A,NDig);
          return
        end
      elseif isinf(B.Numeric)
        % inf + anything else is inf
        result = augmentdigits(B,NDig);
        return
      elseif A.Sign == 0
        % 0 + B = B
        result = B;
        return
      elseif B.Sign == 0
        % A + 0 = A
        result = A;
        return
      end
      % if we drop down here, they are both finite and non-zero
      
      % First, just copy A into result to set up result as a hpf
      result = A;
      
      % The end result will have no more significant digits than
      % the smaller of the two terms
      result.NumberOfDigits = NDig;
      
      % make sure that A and B conform in the number of digits
      if sum(NDig) ~= sum(A.NumberOfDigits)
        A = augmentdigits(A,NDig);
      end
      if sum(NDig) ~= sum(B.NumberOfDigits)
        B = augmentdigits(B,NDig);
      end
      
      % zero has already been excluded
      if (A.Sign*B.Sign) > 0
        % both had the same signs, so the result will have the
        % same sign as either input.
        result.Sign = A.Sign;
        A.Sign = 1;
        B.Sign = 1;
      end
      
      % we can determine the sign easily enough in most cases
      % if there is a difference in the exponents
      ediff = A.Exponent - B.Exponent;
      
      % do we need to shift the digits of A or B?
      if ediff == 0
        % no shifts required. extract the migits for the add
        Am = A.Migits;
        Bm = B.Migits;
      elseif ediff < 0
        % A is the smaller one
        Am = A.Migits;
        Bm = B.Migits;
        
        % make sure we take the exponent from B
        result.Exponent = B.Exponent;
        
        Ashift = mod(-ediff,A.DecimalBase);
        if Ashift > 0
          % we only get in here if the decimal base is greater
          % than 1, AND if the two numbers had an offset in the
          % exponents. We need to shift A by a factor of Ashift
          % powers of 10, to the right.
          A = shiftmig(A,Ashift);
          Am = A.Migits;
          ediff = ediff + Ashift;
        end
        
        % how many migits from A do we need, and how many zero
        % migits will there be?
        nA0 = -ediff/result.DecimalBase; % this will be an integer
        nA = max(0,numel(A.Migits) - nA0);
        if nA == 0
          % A was so small compared to B, that A has no contribution
          result.Migits = B.Migits;
          return
        end
        Am = [zeros(1,nA0),Am(1:nA)];
        
        % do we need to round the last migit up by 1?
        if (nA0 > 0) && (A.Migits(nA+1) > (.5*A.Base))
          Am(end) = Am(end) + 1;
        end
      elseif ediff > 0
        % B is the smaller one
        Am = A.Migits;
        Bm = B.Migits;
        
        % make sure we take the exponent from A
        result.Exponent = A.Exponent;
        
        Bshift = mod(ediff,B.DecimalBase);
        if Bshift > 0
          % we only get in here if the decimal base is greater
          % than 1, AND if the two numbers had an offset in the
          % exponents. We need to shift B by a factor of Bshift
          % powers of 10, to the right.
          B = shiftmig(B,Bshift);
          Bm = B.Migits;
          ediff = ediff - Bshift;
        end
        
        % how many migits from B do we need, and how many zero
        % migits will there be?
        nB0 = ediff/result.DecimalBase; % this will be an integer
        nB = max(0,numel(B.Migits) - nB0);
        if nB == 0
          % B was so small compared to A, that B has no contribution
          result.Migits = A.Migits;
          return
        end
        Bm = [zeros(1,nB0),Bm(1:nB)];
        
        % do we need to round the last migit up by 1?
        if (nB0 > 0) && (B.Migits(nB+1) > (.5*B.Base))
          Bm(end) = Bm(end) + 1;
        end
      end
      
      % add (or subtract) the migits
      if A.Sign < 0
        result.Migits = Bm - Am;
        % am unsure about the sign bit at this point, so set it to 1
        % the test below will correct this if I am wrong.
        result.Sign = 1;
      elseif B.Sign < 0
        result.Migits = Am - Bm;
        % am unsure about the sign bit at this point, so set it to 1
        % the test below will correct this if I am wrong.
        result.Sign = 1;
      else
        % both were positive or negative
        result.Migits = Am + Bm;
        % the sign bit was set correctly before, so leave it alone
      end
      
      % did we get the wrong sign before? If so, then the
      % highest order migit will be negative
      if result.Migits(1) < 0
        result.Migits = -result.Migits;
        result.Sign = -1;
      end
      
      % resolve any carries
      [result.Migits,exponentshift,signflag] = carryop(result.Migits,result.DecimalBase,result.Base);
      result.Exponent = result.Exponent + exponentshift;
      result.Sign = result.Sign*signflag;
      % just in case we left any exponent on a zero result, catch that now.
      if result.Sign == 0
        result.Exponent = 0;
      end
      
      % =====================================================
      % nested function, only used by plus.m
      function F = shiftmig(F,shift)
        % shifts the digints in the migits of F to the right by shift places.
        %
        % Positive (rightward) shifts are used to align the exponents of
        % the two terms in an addition operation. Reverse carries are
        % employed for this operation.
        F.Exponent = F.Exponent + shift;
        
        m = 10.^shift;
        
        % do rightward reverse carries
        carry = mod(F.Migits,m);
        
        % this will result in an exact positive integer,
        % so floor will suffice
        F.Migits = floor(F.Migits/m);
        
        nM = numel(F.Migits);
        if nM > 1
          k = 2:nM;
          L = F.Base./m;
          F.Migits(k) = F.Migits(k) + L*carry(1:(end-1));
        end
        
        % there will be one or more digits lost at the bottom end,
        % so do rounding at that point. In some circumstances,
        % this may result in the last migit element rounding up
        % to an "overflow". This does not matter, since the
        % number in F will be used for an addition, and then
        % a final carry operation will be done, so any potential
        % migit overflow will be resolved in that next step.
        if carry(end) >= m/2
          F.Migits(end) = F.Migits(end) + 1;
        end
      end % function shiftmig
      % =====================================================
        
    end % result = plus(A,B)
    
    function result = power(x,y)
      % Raises a HPF number to a real power, or a hpf and another numeric class.
      % The functional equivalent of x.^y in MATLAB when either of the
      % operands is a HPF number.
      %
      % usage: result = power(x,y)
      % usage: result = x.^y
      %
      % When the two numbers are represented with a different number of
      % digits, use the smaller of the two lengths for result.
      
      % are x or y arrays?
      
      
      
      
      
      
      
      
      
      
      % how many digits?
      if isa(x,'hpf')
        if isa(y,'hpf')
          NDig = combineNDig(x.NumberOfDigits,y.NumberOfDigits);
        else
          % only x was an hpf
          NDig = x.NumberOfDigits;
        end
      else
        % only y was an hpf
        NDig = y.NumberOfDigits;
        x = hpf(x,NDig);
      end
      
      % some special cases first
      if y == 0
        if isnan(x)
          result = hpf('NaN',NDig);
        else
          % anything else to the 0 power is 1, including inf
          % even 0^0 is 1 in matlab
          result = hpf('1',NDig);
        end
        return
      elseif y == 1
        % anything to the first power is itself
        result = hpf(x,NDig);
        return
      elseif isnan(y)
        result = hpf('NaN',NDig);
        return
      elseif isinf(y)
        % any number to an infinite power will be 0, inf, or nan
        if abs(x) == 1
          % an indeterminate case is at 1.
          result = hpf('NaN',NDig);
          return
        elseif x < 0
          % -N^
          result = hpf('0',NDig);
          return
        end
        
        % is y +inf or -inf?
        if y > 0
          % +inf
          if x > 1
            result = hpf('inf',NDig);
          else
            % x must be positive, < 1
            result = hpf('0',NDig);
          end
        else
          % -inf
          if abs(x) > 1
            result = hpf('0',NDig);
          else
            % x must be positive, < 1
            result = hpf('inf',NDig);
          end
        end
        return
      end
      
      % simple test for zero on a resultant hpf number
      iszero = @(F) F.Sign == 0;
      
      if isa(y,'hpf') && iszero(fractionalpart(y))
        if y.Exponent <= 52
          % y can be safely converted to a flint
          y = double(y);
        end
      end
      
      if isnumeric(y) && (y == round(y))
        % y is a numeric integer, but not zero or 1 or an inf or nan
        
        % convert y to binary, forming the powers of x by repeated squaring
        ybin = dec2bin(y) == '1';
        
        if ybin(end)
          result = hpf(x,NDig);
        else
          result = hpf('1',NDig);
        end
        xsq = x.*x;
        for i = (numel(ybin)-1):-1:1
          if ybin(i)
            result = result.*xsq;
          end
          if i > 1
            xsq = xsq.*xsq;
          end
        end
        
        return
      end
      
      % having dropped down to here, we know that if y is numeric, it is
      % either not a flint, or it is too large to fit into a double flint.
      % so make sure that both x and y are hpf numbers, and do the power
      % computation the hard way. There are some special cases to check for
      % first based on the value of x. We have already excluded the
      % special cases in y.
      if x == 0
        % we have already excluded y == 0
        if y > 0
          % zero to any positive power is zero.
          result = hpf('0',NDig);
        else
          % zero to any negative power is inf.
          result = hpf('inf',NDig);
        end
        return
      elseif x == 1 % 1 to any power is 1
        result = hpf('1',NDig);
        return
      elseif isnan(x) % nan to any power is nan
        result = hpf('NaN',NDig);
        return
      end
      
      % nothing left but this
      result = exp(y.*log(x));
      
    end % function power
    
    function Aprod = prod(A,dim)
      % Aprod: product of a HPF array
      % usage: Aprod = prod(A,dim);
      %
      % Arguments:
      %  A - an HPF object array
      %
      %  dim - (optional) dimension of A to product over
      %      DEFAULT: dim is the first non-singleton
      %      dimension of A.
      %
      % Arguments: (output)
      %  Aprod - the product HPF object array
      %
      % Example:
      %  A = hpf(-3:2:25);
      %  prod(A)
      %
      %  ans =
      %      165
      %
      %  See also: sum, cumsum, cumprod
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com

      if (nargin<1) || (nargin>2)
        error('prod takes one or two arguments')
      end
      
      if numel(A) == 1
        % the product of a scalar is a no-op
        Aprod = A;
      else
        % a vector or array
        
        % default for dim?
        if (nargin==1) || isempty(dim)
          dim = find(size(A)>1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        
        % product over the dimension dim. Do so by
        % a permutation, then a prod and an
        % inverse permute.
        P = 1:length(size(A));
        P([1,dim]) = [dim,1];
        
        A = permute(A,P);
        NAprod = size(A);
        N1 = NAprod(1);
        NAprod(1) = 1;
        Aprod = repmat(hpf(1,A(1).NumberOfDigits),NAprod);
        for j = 1:prod(NAprod(2:end))
          for i = 1:N1
            k = i + (j-1)*N1;
            Aprod(j) = Aprod(j).*A(k);
          end
        end
        
        % do the inverse permutation
        Aprod = ipermute(Aprod,P);
        
      end

    end % function Aprod = prod(A,dim)
    
    function result = rdivide(num,den)
      % Arithmetic operation num./den, where one or both are hpf
      % usage: result = rdivide(num,den)
      % usage: result = num./den
      %
      % Where the two numbers are represented with a different number of
      % digits, use the smaller of the two lengths.
      
      % are either num or den vectors/arrays or empty?
      if isempty(num)
        if (numel(den) > 1) || (isempty(den) && ~isequal(size(num),size(den)))
          error('HPF:RDIVIDE:impropersize','Matrix dimensions must agree.')
        else
          result = num;
        end
        return
      elseif isempty(den)
        if (numel(num) > 1) || (isempty(num) && ~isequal(size(num),size(den)))
          error('HPF:RDIVIDE:impropersize','Matrix dimensions must agree.')
        else
          result = den;
        end
        return
      elseif (numel(num) > 1) && (numel(den) == 1)
        % scalar expansion of den
        if isa(num,'hpf')
          result = num;
        else
          result = hpf(num,den.NumberOfDigits);
        end
        recip = reciprocal(hpf(den,num(1).NumberOfDigits));
        for i = 1:numel(num)
          result(i) = num(i).*recip;
        end
        return
      elseif (numel(den) > 1) && (numel(num) == 1)
        % scalar expansion of num
        if isa(den,'hpf')
          result = den;
        else
          result = hpf(den,num.NumberOfDigits);
        end
        for i = 1:numel(den)
          result(i) = num./den(i);
        end
        return
      elseif (numel(num) > 1) && (numel(den) > 1)
        % both num and den are arrays
        if ~isequal(size(num),size(den))
          error('HPF:RDIVIDE:impropersize','Matrix dimensions must agree.')
        end
        
        if isa(num,'hpf')
          result = num;
        else
          result = hpf(num,den.NumberOfDigits);
        end
        for i = 1:numel(den)
          result(i) = num(i)./den(i);
        end
        return
      end
      % if we drop down here, then both num and den are scalars
      
      % at least one of these numbers must be a hpf, so if not num, then den
      % how many digits will we have in the result?
      if isa(num,'hpf')
        if isa(den,'hpf')
          % both are hpf
          NDig = combineNDig(num.NumberOfDigits,den.NumberOfDigits);
        else
          % only num is hpf
          NDig = num.NumberOfDigits;
        end
      else
        % only den is hpf
        NDig = den.NumberOfDigits;
      end
      
      % test for special cases in either num or den
      if isnan(num) || isnan(den)
        % nan begets nan
        result = hpf('nan',NDig);
        return
      elseif num == 0
        % 0/den is always zero, unless den == 0
        if den == 0
          % 0/0 is undefined, so nan
          result = hpf('nan',NDig);
        else
          result = hpf('0',NDig);
        end
        return
      elseif den == 0
        % num/0 is always +/- inf, unless num was zero,
        % or a nan is in there, but we have already trapped
        % out the case 0/0, and nan/0.
        result = hpf('inf',NDig);
        result.Sign = sign(num);
        return
      elseif isinf(den)
        % N/inf is always zero, and signs are not relevant.
        % unless of course, N is also inf. inf/inf = nan.
        if isinf(num)
          result = hpf('nan',NDig);
        else
          result = hpf('0',NDig);
        end
        return
      elseif isinf(num)
        % inf/den is always inf. We have excluded any other
        % special cases for den.
        result = hpf('inf',NDig);
        result.Sign = sign(num).*sign(den);
        return
      end
      % If we drop through here, then num and den are both scalar,
      % and neither one is any of {0, +/- inf, nan}.
      
      % if den is not an hpf number, then num MUST be one to have
      % gotten into here at all. We may be able to do something efficient.
      if isa(den,'hpf')
        % den is an hpf already, but not zero or nan or inf. I won't
        % check for things like a small integer as an hpf number, as
        % that will take more effort than it is perhaps worth.
        result = num.*reciprocal(den);
      else
        % den is not hpf, so num must be hpf.
        
        % Just copy num into result, to set up result as a hpf
        result = num;
        
        % set the known sign of the result
        result.Sign = num.Sign * ((den>=0)*2 - 1);
        % we have taken care of the sign of den, so make it positive
        den = abs(den);
        
        % the number of digits in the result will be the same as for num
        
        % check for the special cases of small integers. That will make
        % other operations more efficient, such as Taylor series expansions
        switch den
          case 1
            % a no-op, with result already correct;
            
          case 10
            % just shift the exponent
            result.Exponent = result.Exponent - 1;
            
          case 100
            % just shift the exponent
            result.Exponent = result.Exponent - 2;
            
          case 1000
            % just shift the exponent
            result.Exponent = result.Exponent - 3;
            
          case 10000
            % just shift the exponent
            result.Exponent = result.Exponent - 4;
            
          case 100000
            % just shift the exponent
            result.Exponent = result.Exponent - 5;
            
          case 1000000
            % just shift the exponent
            result.Exponent = result.Exponent - 6;
            
          case 2
            % divide by 2 is best done as a multiply by 5 and a divide by 10
            % the sign is already set. Note, while I recognize that I could
            % accomplish the divide by 2 fairly efficiently otherwise,
            % the multiply by 5 and exponent shift is still faster.
            result.Exponent = result.Exponent - 1;
            result = result.*hpf(uint8(5),result.NumberOfDigits);
            
          case 4
            % divide by 4 is best done as a multiply by 25 and a divide by 100
            % the sign is already set
            result.Exponent = result.Exponent - 2;
            result = result.*hpf(uint8(25),result.NumberOfDigits);
            
          case 8
            % divide by 8 is best done as a multiply by 125 and a divide by 1000
            % the sign is already set
            result.Exponent = result.Exponent - 3;
            result = result.*hpf(uint8(125),result.NumberOfDigits);
            
          case 16
            % divide by 16 is best done as a multiply by 0.0625
            % the sign is already set
            result = result.*hpf('0.0625',result.NumberOfDigits);
            
          case 32
            % divide by 32 is best done as a multiply by 0.03125
            % the sign is already set
            result = result.*hpf('0.03125',result.NumberOfDigits);
            
          case 64
            % divide by 64 is best done as a multiply by 0.015625
            % the sign is already set
            result = result.*hpf('0.015625',result.NumberOfDigits);
            
          case 5
            % divide by 5 is most easily done as a multiply by 2 and a divide by 10
            % the sign is already set
            result.Exponent = result.Exponent - 1;
            result = result.*hpf('2',result.NumberOfDigits);
            
          case 25
            % divide by 25 is most easily done as a multiply by 4 and a divide by 100
            % the sign is already set
            result.Exponent = result.Exponent - 2;
            result = result.*hpf(uint8(4),result.NumberOfDigits);
            
          case 20
            % divide by 20 is most easily done as a multiply by 5 and a divide by 100
            % the sign is already set
            result.Exponent = result.Exponent - 2;
            result = result.*hpf(uint8(5),result.NumberOfDigits);
            
          case 40
            % divide by 40 is most easily done as a multiply by 0.025
            % the sign is already set
            result = result.*hpf('0.025',result.NumberOfDigits);
            
          case 50
            % divide by 50 is most easily done as a multiply by 0.02
            % the sign is already set
            result = result.*hpf('0.02',result.NumberOfDigits);
            
          case 3
            % form the inverse as an hpf directly
            deninv = hpf('0.3',result.NumberOfDigits);
            deninv.Migits = d2m(wrepvec(3, ...
              sum(result.NumberOfDigits)),result.DecimalBase);
            result = result.*deninv;
            
          case 6
            % form the inverse as an hpf directly
            deninv = hpf('0.1',result.NumberOfDigits);
            deninv.Migits = d2m([1,repmat(6,1, ...
              sum(result.NumberOfDigits) - 2),7],result.DecimalBase);
            result = result.*deninv;
            
          case 7
            % form the inverse as an hpf directly
            deninv = hpf('0.1',result.NumberOfDigits);
            deninv.Migits = d2m(wrepvec([1 4 2 8 5 7], ...
              sum(result.NumberOfDigits)),result.DecimalBase);
            result = result.*deninv;
            
          case 9
            % form the inverse as an hpf directly
            deninv = hpf('0.1',result.NumberOfDigits);
            deninv.Migits = d2m(wrepvec(1, ...
              sum(result.NumberOfDigits)),result.DecimalBase);
            result = result.*deninv;
            
          case 11
            % form the inverse as an hpf directly
            deninv = hpf('0.09',result.NumberOfDigits);
            M = wrepvec([9 0],sum(result.NumberOfDigits));
            if M(end) == 0
              M(end) = 1;
            end
            deninv.Migits(:) = d2m(M,result.DecimalBase);
            result = result.*deninv;
            
          case 12
            % form the inverse as an hpf directly
            deninv = hpf('0.08',result.NumberOfDigits);
            M = [8 wrepvec(3,sum(result.NumberOfDigits) - 2)];
            deninv.Migits = d2m(M,result.DecimalBase);
            result = result.*deninv;
            
          case 15
            % form the inverse as an hpf directly
            deninv = hpf('0.06',result.NumberOfDigits);
            M = wrepvec(6,sum(result.NumberOfDigits) - 1);
            M(end) = 7;
            deninv.Migits = d2m(M,result.DecimalBase);
            result = result.*deninv;
            
          otherwise
            % form the inverse as a multiply by the reciprocal
            result = result.*reciprocal(hpf(den,result.NumberOfDigits));
        end
      end
      
      
    end % result = rdivide(num,den)
    
    function R = reciprocal(D)
      % The reciprocal of the hpf number D
      % Is more efficient than 1./D, since less testing needs to be done.
      
      % Is this an array?
      if isempty(D)
        % empty begets empty
        R = D;
        return
      elseif numel(D) > 1
        R = D;
        for i = 1:numel(D)
          R(i) = reciprocal(D(i));
        end
        return
      end
      % if we drop through, then D was scalar
      
      % How many digits do we live in?
      NDig = D.NumberOfDigits;
      
      % trap for special cases
      if D.Sign == 0
        % 1/0 = inf
        R = hpf('inf',NDig);
        return
      elseif isnan(D.Exponent)
        % NaN begets NaN
        R = D;
        return
      elseif isinf(D.Exponent)
        % 1/inf = 0, 1/(-inf) = 0
        R = hpf('0',NDig);
        return
      elseif (D.Migits(1) == (D.Base/10)) && all(D.Migits(2:end) == 0)
        % +/-1, with some pure power of 10 applied to it
        R = D;
        
        % negate the exponent. the extra 2 in there is a reflection of the
        % scientific notation used by hpf.
        R.Exponent = 2-D.Exponent;
        
        % nothing else to do
        return
      end
      % we have dropped through, so we have eliminated [], inf, nan, 0, 1.
      % all that remains are generic scalars.
      
      % will need some simple hpf numbers for later
      two = hpf('2',NDig);
      one = hpf('1',NDig);
      
      % This will capture the sign of D, as well as creating an
      % hpf container for the reciprocal.
      R = D;
      
      % just negate the exponent
      R.Exponent = 2 - D.Exponent;
      % and stuff a unit set of migits
      R.Migits(:) = 0;
      R.Migits(1) = R.Base/10;
      
      % having done the above steps, we will reduce the problem
      % to one where D is in the interval (1,10). (Remember
      % that 1 is not a possibility at this point.)
      D.Exponent = 1;
      D.Sign = 1;
      
      % we need to reduce D to the interval (0.5,1), while getting
      % a very good initial estimate for 1/D. Do this by converting D
      % to a double.
      Ddouble = double(D);
      Rdouble = 1./Ddouble;
      % reduce it by a tiny amount
      Rdouble = Rdouble*(1 - 10*eps(Rdouble));
      
      % convert to hpf
      R0 = hpf(Rdouble,NDig);
      
      % by scaling both D and R by this estimate of the reciprocal,
      % we get a starting value for R that is within a tiny amount of
      % the final value. In fact, we should now have D approximately
      % equal to 1-10*eps
      R = R*R0;
      D = D*R0;
      
      % The inverse of D, for a number VERY near 1, but just slightly less
      % than 1. We will get a starting value from the Taylor series. This
      % will be accurate to about 30 digits of precision.
      u = one - D;
      if sum(NDig) < 13
        % for less than 30 digits required, this will be essentially exact
        Dinv = one + u;
      elseif sum(NDig) < 28
        % for less than 30 digits required, this will be essentially exact
        Dinv = one + (u + u*u);
      else
        % For more digits required than 30, this will give us roughly
        % 30 digits of precision to start with. The Newton's method
        % refinement that follows will be quadratically convergent.
        Dinv = one + (u + u.*u);
        
        % Newton's method, applied to f(x) = 1/x + D
        flag = true;
        Dold = Dinv;
        cycle = 0;
        while flag
          DinvD = Dinv.*D;
          Dinv = Dinv.*(two - DinvD);
          k = find(Dinv.Migits ~= Dold.Migits,1,'first');
          if isempty(k) || (cycle > 3)
            flag = false;
          else
            if k == numel(Dinv.Migits)
              cycle = cycle + 1;
            end
            Dold = Dinv;
          end
        end
      end
      
      % combine it all together
      R = R.*Dinv;
      
    end % function reciprocal
    
    function F = round(F)
      % Round towards the nearest integer
      
      for i = 1:numel(F)
        if isfinite(F.Numeric)
          % if F was not finite, then round is a no-op
          db = F(i).DecimalBase;
          
          % identify the decimal digit on which rounding will be based.
          % thus, the decimal place falls BEFORE this digit.
          decind = F(i).Exponent + 1;
          if decind <= 0
            % abs(x) <= 0.1
            F(i) = hpf('0',F(i).NumberOfDigits);
          elseif decind == 1
            % the number was in the half open interval [0.1,1)
            % extract that digit
            M = F(i).Migits(1);
            D = m2d(M,db);
            
            % was the first digit at least 5?
            F(i).Migits(:) = 0;
            if D(1) >= 5
              % we have either plus or minus 1 as a result
              F(i) = hpf(F(i).Sign,F(i).NumberOfDigits);
            else
              % we rounded down, and the result was zero.
              F(i) = hpf('0',F(i).NumberOfDigits);
            end
            
          elseif decind <= sum(F(i).NumberOfDigits)
            % there will be some rounding done, but it was NOT
            % on the first decimal digit of the number.
            % which Migit does it happen in?
            mind = ceil(decind./db);
            mindi = 1 + mod(decind - 1,db);
            
            % extract that migit
            M = F(i).Migits(mind);
            D = m2d(M,db);
            
            if D(mindi) >= 5
              % a round up. Add one to the previous digit, and set
              % the lower digits to zero.
              F(i).Migits((mind+1):end) = 0;
              D(mindi:end) = 0;
              M = d2m(D,db) + 10.^(db - mindi + 1);
              
              F(i).Migits(mind) = M;
              if M >= F(i).Base
                [F(i).Migits,exponentshift] = carryop(F(i).Migits,db,F(i).Base);
                if exponentshift ~= 0
                  F(i).Exponent = F(i).Exponent + exponentshift;
                end
              end
            else
              % a round down. set that digit, and all that follow to zero
              D(mindi:end) = 0;
              F(i).Migits((mind+1):end) = 0;
              F(i).Migits(mind) = d2m(D,db);
            end
          end
        end
      end
    end % F = round(F)
    
    function F = roundn(F,n)
      % Round towards the nearest multiple of 10.^n
      % i.e., roundn(X,-1) rounds to the nearest multiple of 0.1.
      % roundn(X,3) rounds to the nearest 1000.
      
      for i = 1:numel(F)
        if isfinite(F.Numeric)
          % if F was not finite, then roundn is a no-op
          dn = double(n);
          
          % Just shift the exponent. Admittedly, it would probably have
          % been slightly better to have written round to call roundn,
          % that to do it by roundn calling round. I'll do it that way the
          % next time around.
          F(i).Exponent = F(i).Exponent - dn;
          % round
          F(i) = round(F(i));
          % and then shift back
          
          if F(i).Sign ~= 0
            F(i).Exponent = F(i).Exponent + dn;
          end
        end
      end
    end % F = roundn(F)
    
    function F = sec(X)
      % evaluates sec(X) (sec(X) = 1/cos(X)) for an hpf number X in radians
      
      % a simple scheme suffices for now, costing only a spare divide
      F = reciprocal(cos(X));
      
    end % function F = sec(X)
    
    function F = secd(X)
      % secd(X) (secant function, where X is in degrees) for an hpf number X
      
      % a simple scheme suffices for now, costing only a spare divide
      F = reciprocal(cosd(X));
      
    end % function F = secd(X)
    
    function F = sech(X)
      % evaluates sech(X) (the hyperbolic secant) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology.
      F = X;
      for i = 1:numel(X)
        if isinf(X(i).Numeric)
          F(i) = hpf('0',X(i).NumberOfDigits);
        elseif X(i).Sign == 0
          F(i) = hpf('1',X(i).NumberOfDigits);
        elseif isnan(X(i).Numeric)
          F(i) = hpf('NaN',X(i).NumberOfDigits);
        else
          Fi = exp(X(i));
          F(i) = (2.*Fi)./(Fi.*Fi + 1);
        end
      end
      
    end % function F = sech(X)
    
    function S = sign(F)
      % Returns the sign of a hpf number F, as one of {-1, 0, +1}
      % This is often known as the signum function
      
      % Return the sign, UNLESS it was a nan.
      S = zeros(size(F));
      for i = 1:numel(F)
        if isnan(F.Numeric)
          S(i) = NaN;
        else
          S(i) = F(i).Sign;
        end
      end
    end % S = sign(F)
    
    function F = sin(X)
      % evaluates sin(X) for a hpf number X (in radians)
      
      % is X a scalar?
      if isempty(X)
        F = X;
        return
      elseif numel(X) > 1
        % vector or array
        F = X;
        for i = 1:numel(X)
          F(i) = sin(X(i));
        end
        return
      end
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if ~isfinite(X.Numeric)
        % sin(NaN) = sin(inf) = sin(-inf) = NaN
        F = hpf('NaN',NDig);
        return
      elseif X.Sign == 0
        % sin(0) == 0
        F = hpf('0',NDig);
        return
      end
      
      % get pi for that number of digits because sin(X) is a periodic
      % function.
      pie = hpf('pi',NDig);
      twopie = 2 .*pie;
      
      % We need to reduce X into the interval [-pi,pi].
      k = X./twopie;
      kfrac = fractionalpart(k);
      X = X - floor(k)*twopie;
      
      % convert kfrac to a double, to most easily determine the subinterval
      kfrac= double(kfrac);
      if kfrac < 0
        kfrac = 1 + kfrac;
      end
      
      % what sub-interval did k live in?
      % the break points for our bins are {1/8, 3/8, 5/8, 7/8}
      % so one trick here might have been to multiply by 8, then just take
      % the first digit pf the result. But it will probably be just as fast
      % to convert to double.
      bins = [1 3 5 7]/8;
      [~,kindex] = histc(kfrac,bins);
      switch kindex
        case 1
          % X was in the interval 2*pi*(n + [1/8 3/8]). This will be most
          % efficiently approximated using a shift, then a call to coscore
          F = coscore(X-pie./2);
        case 2
          % X was in the interval 2*pi*(n + [3/8 5/8]). This will be most
          % efficiently approximated using a shift, then a call to sincore
          F = uminus(sincore(X-pie));
        case 3
          % X was in the interval 2*pi*(n + [1/8 3/8]). This will be most
          % efficiently approximated using a shift, then a call to coscore
          % Note that 3/4 is exactly converted to a hpf number here.
          F = uminus(coscore(X-pie.*(3/2)));
        otherwise
          % X was in the interval 2*pi*(n + [-1/8 1/8]). This will be
          % simply done as a call to sincore.
          F = sincore(X);
      end
      
    end % function F = sin(X)
    
    function F = sind(X)
      % evaluates sin(X) for a hpf number X (in degrees)
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % simple conversion of degrees to radians
      F = sin(X.*hpf('pi',NDig)./180);
      
    end % function F = sind(X)
    
    function F = single(X)
      % Convert an hpf number into its single precision form.
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      if isempty(X)
        % empty begets empty
        F = single([]);
      else
        % F is an array
        F = zeros(size(X),'single');
        for i = 1:numel(X)
          % this is a simple solution, just using the conversion from
          % double to single. Since double is itself pretty efficient,
          % the cost is minimal.
          F(i) = single(double(X(i)));
        end
      end
    end % F = single(X)
    
    function F = sinh(X)
      % evaluates sinh(X) (the hyperbolic sine) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology.
      F = (exp(X) - exp(-X))./2;
      
    end % function F = sinh(X)
    
    function [Ns,tags] = sort(N,dim,sdir)
      % SORT: Sorts an hpf array in ascending or descending order.
      %
      %   For vectors, SORT(X) sorts the elements of X in ascending order.
      %   For matrices, SORT(X) sorts each column of X in ascending order.
      %   For N-D arrays, SORT(X) sorts the along the first non-singleton
      %   dimension of X.
      %
      %   Ns = SORT(N,DIM,SDIR)
      %
      %   SORT has two optional parameters:
      %   DIM selects a dimension along which to sort.
      %   SDIR selects the direction of the sort
      %      'ascend' results in ascending order
      %      'descend' results in descending order
      %
      %   The result is in Ns which has the same shape and type as X.
      %
      %   [Ns,tags] = SORT(X,DIM,MODE) also returns an index matrix tags.
      %
      %   If Ns is a vector, then Ns = N(tags).
      %   When N is an array, the elements in the specified dimension will
      %   be sorted, and tags will indicate the sort order as with a vector.
      %
      %   When more than one element has the same value, the order of the
      %   elements are preserved in the sorted result and the indexes of
      %   equal elements will be ascending in any index matrix.
      %
      %   The sorting methodology used is a merge sort:
      %   http://en.wikipedia.org/wiki/Merge_sort
      %
      % Example:
      %  See the help for the built-in sort
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release: 1.0
      %  Release date: 3/5/09
      
      
      if (nargin < 1) || (nargin > 3)
        error('HPF:SORT:arguments','Sort expects 1, 2, or 3 arguments')
      end
      
      % propagate empty
      if isempty(N)
        Ns = [];
        tags = [];
        return
      end
      
      % check for sort direction, or default
      if (nargin<3) || isempty(sdir)
        sdir = 'ascend';
      end
      if ~ischar(sdir)
        error('HPF:SORT:sortdirection', ...
          'Sort order must be either ''ascend'' or ''descend'' if supplied')
      end
      sdir = lower(sdir);
      if ~ismember(sdir,{'ascend','descend'})
        error('Sort order must be either ''ascend'' or ''descend'' if supplied')
      end
      
      % default dimension to sort on
      if (nargin<2) || isempty(dim)
        dim = find(size(N) > 1,1,'first');
        if isempty(dim)
          dim = 1;
        end
      end
      if ~isnumeric(dim) || (dim <= 0) || (dim > length(size(N))) || (dim~=round(dim))
        error('Dimension must be numeric, integer, positive, <= length(size(N))')
      end
      
      % sort over the dimension dim. Do so by
      % a permutation, then a vector sort in
      % a loop and an inverse permute.
      P = 1:length(size(N));
      P([1,dim]) = [dim,1];
      
      N = permute(N,P);
      sizeN = size(N);
      s1 = sizeN(1);
      Ns = reshape(N,s1,[]);
      tags = zeros(size(Ns));
      for j = 1:size(Ns,2)
        if strcmp(sdir,'ascend')
          [Ns(:,j),tags(:,j)] = sortvecA(Ns(:,j));
        else
          [Ns(:,j),tags(:,j)] = sortvecD(Ns(:,j));
        end
      end
      
      % do the inverse permutation
      Ns = ipermute(reshape(Ns,sizeN),P);
      tags = ipermute(reshape(tags,sizeN),P);
      
    end % function sort
    
    function xroot = sqrt(x)
      % Computes the square root of x, where x is an hpf number.
      % usage: xroot = sqrt(x)
      %
      % The method employed is more efficient than the Babylonian method,
      % since that solution uses divides. Instead, we apply Newton's method
      % to solving for 1/sqrt(x). This is a divide-free method, so is far
      % more efficient for hpf numbers.
      %
      % A derivation for this scheme can be found online, at:
      %  http://www.azillionmonkeys.com/qed/sqroot.html
      
      % is x a scalar?
      if isempty(x)
        xroot = x;
        return
      elseif numel(x) > 1
        % vector or array
        xroot = x;
        for i = 1:numel(x)
          xroot(i) = sqrt(x(i));
        end
        return
      elseif isnan(x) || (isinf(x) && (x.Sign >= 0))
        % catch the nans or +infs. -inf gets caught later
        xroot = x;
        return
      end
      
      % check for zero or negative
      if x.Sign == 0
        % a no-op
        xroot = x;
        return
      elseif x.Sign < 0
        % imaginary result
        xroot = hpf('NaN',x.NumberOfDigits);
        warning('HPF:sqrtofnegativenumber', ...
          'sqrt of negative value. Only real HPF numbers are currently supported')
        return
      end
      
      % how many digits do we need?
      NDig = x.NumberOfDigits;
      
      % reduce the problem to bring x into the interval [0.1,10)
      xroot = hpf('1',NDig);
      p = floor(x.Exponent/2);
      xroot.Exponent = p + 1;
      x.Exponent = x.Exponent - 2*p;
      % we have now reduced x to a reasonable number in case it would
      % have otherwise exceeded the dynamic range of a double. x must
      % now lie in the half open interval [1,10).
      
      % get a (very) good estimate for the inverse square root of x
      % using the double precision sqrt.
      xest = hpf(1./sqrt(double(x)),NDig);
      
      % Newton on the inverse sqrt, a divide-free method. This will
      % still double the number of correct digits for each iteration
      % of the process, and it is a divide-free iteration so very fast.
      % note that the division by 2 is actually accomplished in ./
      % by a divide by 10 and a multiply by 5.
      three = hpf('3',NDig);
      half = hpf('0.5',NDig);
      flag = true;
      cycle = 0;
      while flag
        xprior = xest.Migits;
        xest = half.*xest.*(three - x.*xest.*xest);
        
        k = find(xest.Migits ~= xprior,1,'first');
        if isempty(k) || (cycle > 2)
          flag = false;
        else
          if k == numel(xprior)
            cycle = cycle + 1;
          end
        end
      end
      % merge the pieces together. see that
      %   x.*1/sqrt(x) == sqrt(x)
      % so this recovers the desired root with nary a hpf divide required.
      xroot = xroot.*x.*xest;
      
    end % function xroot = sqrt(x)
    
    function Asum = sum(A,dim)
      % Asum: sum of a HPF array
      % usage: Asum = sum(A,dim);
      %
      % Arguments:
      %  A - an HPF object array
      %
      %  dim - (optional) dimension of A to sum over
      %      DEFAULT: dim is the first non-singleton
      %      dimension of A.
      %
      % Arguments: (output)
      %  Asum - the sum HPF object array
      %
      % Example:
      %  A = hpf(-3:2:25);
      %  sum(A)
      %
      %  ans =
      %      165
      %
      %  See also: prod, cumsum
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com

      if (nargin<1) || (nargin>2)
        error('sum takes one or two arguments')
      end
      
      if numel(A) == 1
        % the sum of a scalar is a no-op
        Asum = A;
      else
        % a vector or array
        
        % default for dim?
        if (nargin==1) || isempty(dim)
          dim = find(size(A)>1,1,'first');
          if isempty(dim)
            dim = 1;
          end
        end
        
        % sum over the dimension dim. Do so by
        % a permutation, then a sum and an
        % inverse permute.
        P = 1:length(size(A));
        P([1,dim]) = [dim,1];
        
        A = permute(A,P);
        NAsum = size(A);
        N1 = NAsum(1);
        NAsum(1) = 1;
        Asum = repmat(hpf(0,A(1).NumberOfDigits),NAsum);
        for j = 1:prod(NAsum(2:end))
          for i = 1:N1
            k = i + (j-1)*N1;
            Asum(j) = Asum(j) + A(k);
          end
        end
        
        % do the inverse permutation
        Asum = ipermute(Asum,P);
        
      end

    end % function Asum = sum(A,dim)
    
    function F = tan(X)
      % evaluates tan(X) for a hpf number X (in radians)
      % The scheme here will be to compute versin(X) internally,
      % after X has been brought into the correct period.
      % versin(X) = 1 - cos(X), thus we will not have subtractive
      % cancellation problems near 0. Then compute tan(X) from the
      % versin. All of this gains about 50% in speed over a simpler
      % scheme to compute the tangent.
      
      % is X a scalar?
      if isempty(X)
        F = X;
        return
      elseif numel(X) > 1
        % vector or array
        F = X;
        for i = 1:numel(X)
          F(i) = tan(X(i));
        end
        return
      end
      
      % how many digits do we need to carry?
      NDig = X.NumberOfDigits;
      
      % special cases
      if ~isfinite(X.Exponent)
        % tan(NaN) = tan(inf) = tan(-inf) = NaN
        F = hpf('NaN',NDig);
        return
      elseif X.Sign == 0
        % tan(0) == 1
        F = hpf('1',NDig);
        return
      end
      
      % get pi for that number of digits because tan(X) is a periodic
      % function, with period pi.
      pie = hpf('pi',NDig);
      
      % We need to reduce X into the interval [-pi/2,pi/2].
      k = X./pie;
      kint = round(k);
      kfrac = k - kint;
      X = X - kint*pie;
      % abs(kfrac) < 0.5, abs(X) <= pi/2
      
      % convert kfrac to a double, to most easily determine the subinterval
      kfrac= double(kfrac);
      
      % what sub-interval did X live in?
      % if X is near zero, then it is best to compute the versin
      % accurately using a series solution.
      % for X away from zero, then use versin(X) = 1-cos(X)
      if abs(kfrac) > 0.25
        % X was away from zero, so compute 1 - cos(X).
        % However, we know what cos will do here. It simply
        % converts the problem to a sincore problem anyway.
        if kfrac > 0
          % cos(X), for X in the interval [pi/4,pi/2]
          % is computed as a transformation into sin
          V = 1 + sincore(X - pie/2);
        else
          % cos(X), for X in the interval [-pi/2,-pi/4]
          % is computed as a transformation into a sin function
          V = 1 - sincore(X + pie/2);
        end
      else
        % X is near enough to zero that we will want to compute
        % the versin directly. This is simply the cos series, but
        % without the unit first term in that series.
        
        % How many terms in the Taylor series do we need to compute?
        % See coscore for details. There is a tweak, because the
        % first term in the series here is actually x^2/2
        %
        %  -1 + log(10^-NDig) + log(X^2/2) = log(X^m/m!) = m log(X) - log(m!)
        %
        %  -1 - NDig*log(10) - log(2) = -2*log(X) - m*log(X) - [1/2*log(2*pi) + 1/2*log(m) + m*log(m) - m]
        %
        %  1/2*(log(2*pi)-2) - NDig*log(10) + log(2) = 2*log(X) + m*(log(X) + 1) - log(m)*(m + 1/2)
        %
        % fun is a decreasing function. The zero crossing defines the number
        % of the last term we need to compute.
        dx = abs(double(X));
        fun = @(m) m.*(log(dx) + 1) - 2*log(dx) - log(m).*(m + 1/2) + ...
          sum(NDig)*log(10) - 1/2*(log(2*pi) - 2) + log(2);
        
        % bisection scheme.
        mterms = termsbisector(fun);
        
        % Each term in this series goes by two powers of X, so we will
        % really need only half as many terms for convergence.
        % Make sure that mterms is an even number, by rounding up to
        % an even number if necessary.
        mterms = mterms + double(mod(mterms,2));
        
        % run the versin series in reverse now.
        Xsq = X.*X;
        V = hpf('1',NDig);
        Fact = V;
        for m = mterms:-2:2
          Fact = uminus(Fact.*(m.*(m-1)));
          if m > 2
            V = Xsq*V + Fact;
          end
        end
        % because we ran the loop backwards, in the end we need to
        % divide by factorial(mterms). Don't forget that we are
        % computing 1-cos(X), so there is a sign change
        V = uminus(V).*Xsq./Fact;
      end
      
      % We now have V = versin(X). The tangent is now simply computed
      F = sqrt(V.*(2 - V))./(1-V);
      
      % make sure we get the correct sign on the result
      if X.Sign < 0
        F = uminus(F);
      end
      
    end % function F = tan(X)
    
    function F = tand(X)
      % evaluates tan(X) for a hpf number X (in degrees)
      
      % a simple scheme suffices for now
      F = sind(X)./cosd(X);
      
    end % function F = tand(X)
    
    function F = tanh(X)
      % evaluates tanh(X) (the hyperbolic tangent) for an hpf number X
      
      % a simple scheme suffices for now, since exp is efficient, with good
      % range reduction methodology. Just check for the special cases.
      F = X;
      for i = 1:numel(X)
        if isinf(X(i).Numeric)
          F(i) = hpf('1',X(i).NumberOfDigits);
          F(i).Sign = X(i).Sign;
        elseif X(i).Sign == 0
          F(i) = hpf('0',X(i).NumberOfDigits);
        else
          Fi = exp(2*X(i));
          F(i) = (Fi-1)./(Fi+1);
        end
      end
    end % function F = tanh(X)
    
    function result = times(A,B)
      % Multiplies two hpf numbers, or a hpf and another numeric class
      % Where the two numbers are represented with a different number of
      % digits, use the smaller of the two lengths.
      
      % at least one of these numbers must be a hpf, so if not A, then B.
      % if A is not a hpf, then swap the two numbers
      if ~isa(A,'hpf')
        % do the swap
        [B,A] = deal(A,B);
      end
      
      % Just copy A into result, to set up result as a hpf
      result = A;
      
      % are either A or B arrays of any size, or empty?
      if isempty(A)
        if numel(B) > 1
          error('HPF:TIMES:impropersize','Matrix dimensions must agree.')
        else
          result = hpf([]);
        end
        return
      elseif isempty(B)
        if numel(A) > 1
          error('HPF:TIMES:impropersize','Matrix dimensions must agree.')
        else
          result = hpf([]);
        end
        return
      elseif (numel(A) > 1) && (numel(B) == 1)
        % scalar expansion of B
        result = repmat(hpf(0),size(A));
        for i = 1:numel(A)
          result(i) = A(i).*B;
        end
        return
      elseif (numel(B) > 1) && (numel(A) == 1)
        % scalar expansion of A
        result = repmat(hpf(0),size(B));
        for i = 1:numel(B)
          result(i) = A.*B(i);
        end
        return
      elseif (numel(A) > 1) && (numel(B) > 1)
        % both A and B are arrays
        if ~isequal(size(A),size(B))
          error('HPF:TIMES:impropersize','Matrix dimensions must agree.')
        end
        
        result = repmat(hpf(0),size(B));
        for i = 1:numel(B)
          result(i) = A(i).*B(i);
        end
        return
      end
      % if we drop down, then both A and B are scalars
      
      % if A was zero, then we can quit now
      if isnan(A) || isnan(B)
        % nan times anything gives nan
        result = hpf('NaN',A.NumberOfDigits);
        
        if isa(B,'hpf')
          result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
        end
        return
      elseif isinf(A)
        % inf*anything is inf, except for inf*0
        if B == 0
          result = hpf('NaN',A.NumberOfDigits);
          if isa(B,'hpf')
            result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
          end
        else
          result = hpf('inf',A.NumberOfDigits);
          result.Sign = sign(A)*sign(B);
          if isa(B,'hpf')
            result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
          end
        end
        return
      elseif isinf(B)
        % inf*anything is inf, except for inf*0
        if A.Sign == 0
          result = hpf('NaN',A.NumberOfDigits);
          if isa(B,'hpf')
            result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
          end
        else
          result = hpf('inf',A.NumberOfDigits);
          result.Sign = sign(A)*sign(B);
          if isa(B,'hpf')
            result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
          end
        end
        return
      elseif (A.Sign == 0) || (B == 0)
        result = hpf('0',A.NumberOfDigits);
        if isa(B,'hpf')
          result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
        end
        return
      end
      % if we drop down here, neither A or B are 0, inf, nan, or empty.
      % they are just simple scalars. Is B an HPF, or is it something else?
      
      % check on B as a hpf, because if the the non-hpf number is a
      % simple multiply, then we can do things sometimes more efficiently.
      % otherwise, just make B an hpf too
      if isa(B,'hpf')
        % both are hpf numbers
        
        % do the numbers have the same decimal base? choose the smaller one
        % if they somehow differ. That would be unlikely though.
        if A.DecimalBase ~= B.DecimalBase
          if A.DecimalBase < B.DecimalBase
            B = adjustdecimalbase(B,A.DecimalBase);
          else
            A = adjustdecimalbase(A,B.DecimalBase);
          end
        end
        
        % set the known sign of the result
        result.Sign = A.Sign * B.Sign;
        
        result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
        
        % Be careful - the number 1 is stored in the form
        % 0.1000..., with an exponent of 1.
        result.Exponent = A.Exponent + B.Exponent - result.DecimalBase;
        
      elseif isa(B,'vpi')
        % just make B directly into an hpf
        B = hpf(B,A.NumberOfDigits);
        result.NumberOfDigits = combineNDig(A.NumberOfDigits,B.NumberOfDigits);
        
        result.Sign = A.Sign * B.Sign;
        result.Exponent = A.Exponent + B.Exponent - result.DecimalBase;
        
      else
        % B must be a standard numeric form
        
        % set the known sign of the result. the double resolves
        % the case of logical B. (sign fails on logicals for some
        % silly reason)
        result.Sign = A.Sign * sign(double(B));
        B = abs(B);
        
        % the number of digits in the result will be the same as A
        
        % check for the special cases
        switch B
          case {1 -1}
            % a no-op, with result already correct, even the sign.
            return
          case {10 -10}
            % just shift the exponent
            result.Exponent = result.Exponent + 1;
            return
          case {100 -100}
            % just shift the exponent
            result.Exponent = result.Exponent + 2;
            return
          case {1000 -1000}
            % just shift the exponent
            result.Exponent = result.Exponent + 3;
            return
          case {10000 -10000}
            % just shift the exponent
            result.Exponent = result.Exponent + 4;
            return
          case {100000 -100000}
            % just shift the exponent
            result.Exponent = result.Exponent + 5;
            return
          case {1000000 -1000000}
            % just shift the exponent
            result.Exponent = result.Exponent + 6;
            return
        end
        
        % if we drop through here, B is some other number
        % if it is a sufficiently small integer, then we still have a plan.
        % sicne any migit can never exceed 1e6, then as long as B is an
        % integer no larger than 1e9, the product will never exceed 1e15,
        % which is less than 2^53, so still a flint.
        if (B < 1e9) && (B == round(B))
          % B is small enough that we can do even better than conv
          % coupled with a carry. Just a scalar multiply, then carries.
          result.Migits = A.Migits.*double(B);
          % result.Exponent is already set to A.Exponent. The call to
          % carryop will shift the exponent as needed.
          [result.Migits,exponentshift] = carryop( ...
            result.Migits,result.DecimalBase,result.Base);
          result.Exponent = result.Exponent + exponentshift;
          
          return
        end
        
        % conv is blazingly fast anyway. So just make B into an
        % hpf with the proper number of digits.
        B = hpf(B,A.NumberOfDigits);
        result.Exponent = A.Exponent + B.Exponent - result.DecimalBase;

      end
      % Both numbers must now be hpf numbers, or we have already exited
      
      % the digit multiply is just a call to conv.
      multdigits = conv(A.Migits,B.Migits);
      
      % Do the carries...
      [multdigits,exponentshift] = carryop(multdigits, ...
        result.DecimalBase,result.Base);
      
      result.Exponent = result.Exponent + exponentshift;
      
      % stuff the digits into result
      result.Migits = multdigits(1:(sum(result.NumberOfDigits)/result.DecimalBase));
      
    end % result = times(A,B)
    
    function F = uint8(X)
      % Convert an hpf number into its uint8 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = uint8(uint64(X));
      
    end % F = uint8(X)
    
    function F = uint16(X)
      % Convert an hpf number into its uint16 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = uint16(uint64(X));
      
    end % F = uint32(X)
    
    function F = uint32(X)
      % Convert an hpf number into its uint32 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % simplest
      F = uint32(uint64(X));
      
    end % F = uint32(X)
    
    function F = uint64(X)
      % Convert an hpf number into its uint64 form (after floor)
      % Lower order digits will be lost.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      if isempty(X)
        % empty begets empty
        F = [];
      elseif numel(X) > 1
        % F is an array
        F = zeros(size(X));
        for i = 1:numel(X)
          F(i) = uint64(X(i));
        end
      else
        % X is scalar. any special cases?
        if ~isfinite(X.Numeric)
          % inf will map to 2^64-1. -inf and NaN will map to 0,
          % which is consistent with MATLAB, as uint64(NaN) = 0.
          F = uint64(X.Numeric.*X.Sign);
        elseif X.Sign <= 0
          F = uint64(0);
        elseif X.Exponent > 20
          % we can be sure this is an overflow.
          F = uint64(inf);
        else
          % something in the middle, although it may still
          % be an overflow for uint64
          X = floor(X);
          
          db = X.DecimalBase;
          E = X.Exponent;
          base = uint64(X.Base);
          migits = uint64(X.Migits);
          mlast = find(migits ~= 0,1,'last');
          F = uint64(0);
          for i = 1:(mlast-1)
            F = F.*base + migits(i);
            E = E - db;
          end
          
          % extract the high non-zero digits from the last non-zero migit
          Dlast = m2d(double(migits(mlast)),db);
          for i = 1:db
            F = F.*uint64(10) + Dlast(i);
            E = E - 1;
            if E <= 0
              break
            end
          end
          
          if E > 0
            % there were still a few powers of 10 in there
            F = F.*uint64(10).^E;
          end
        end
      end
    end % F = int64(X)
    
    function F = uminus(F)
      % unary minus for hpf numbers, thus: -F
      %
      % The syntax -F negates the hpf element F, calling uminus.
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      % this works for an array, vector, scalar, even empty, inf or nan
      for i = 1:numel(F)
        % leave the sign alone for NaNs
        if ~isnan(F(i).Numeric)
          F(i).Sign = -F(i).Sign;
        end
      end
      
    end % F = uminus(F)
    
    function F = uplus(F)
      % unary plus for a hpf number, thus: +F
      %
      % The syntax +F is an identity operation, returning F itself
      %
      %  See also:
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
    end % F = uplus(F)
    
    function N = vpi(F)
      % Convert a hpf number into its VPI form.
      % Fractional digits will be roounded.
      %
      %  See also: double
      %
      %  Author: John D'Errico
      %  e-mail: woodchips@rochester.rr.com
      %  Release date: 1/25/11
      
      if isempty(F)
        % empty begets empty
        N = vpi([]);
      elseif numel(F) > 1
        % F is an array
        N = vpi(zeros(size(F)));
        for i = 1:numel(X)
          N(I) = vpi(F(i));
        end
      else
        % a scalar hpf
        if ~isfinite(F.Exponent)
          error('HPF:VPI:nonfinite','VPI numbers do not allow infinite or NaN values')
        else
          % it is a real number
          
          % first, round it to the nearest integer
          F = round(F);
          
          % was it already zero, or became zero after rounding?
          if F.Sign == 0
            N = vpi(0);
            return
          end
          
          % extract the necessary digits
          digitlist = padz(mantissa(F),F.Exponent,'right');
          % as characters
          digitlist = char(digitlist + '0');
          
          % is it negative?
          if F.Sign < 0
            digitlist = ['-',digitlist];
          end
          
          % convert this digit string into a vpi
          N = vpi(digitlist);
        end
      end
    end % F = vpi(F)
    
  end % methods
  
end

% =============================================================
%      end mainline hpf, begin subfunctions
% =============================================================
function mterms = termsbisector(fun)
% mterms is the last term in the Taylor series that we will
% need to get a good approximation from that series. Always take
% at least 5 terms in the series.
%
% This function is used specifically by a few callers, mainly
% sincore, coscore, exp, etc.
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
end
end % function termsbisector
% =============================================================

function F = sincore(X)
% evaluates sin(X) for a scalar hpf number X, where -pi/4 <= X <= pi/4
% No test is made to assure that X is in the proper interval. X is
% assumed to be in radians.
%
% The methodology used is a simple Taylor series, but run in reverse.
% This minimizes all divisions except one divide at the very end.

% how many digits do we need to carry?
NDig = X.NumberOfDigits;

% sin(0) == 0
if X.Sign == 0
  F = hpf('0',NDig);
  return
end

% How many terms in the Taylor series do we need to compute?
% The sine series has as its general term z^m/factorial(m), with
% m = 2*n-1, n being the index of the term we choose to stop at.
%
% we can quit when the natural log of that expression is less
% than the log of our precision goal. (The -1 below gives us some
% protection for round-off).
%
%  -1 + log(10^-NDig) = log(X^m/m!) = m log(X) - log(m!)
%
% the logs here are natural logs of course. Using Stirling's
% approximation for m!, we get
%
%  -1 - NDig*log(10) = m*log(X) - [1/2*log(2*pi) + 1/2*log(m) + m*log(m) - m]
%
% Combining terms, this reduces to
%
%  1/2*(log(2*pi)-2) - NDig*log(10) = m*(log(X) + 1) - log(m)*(m + 1/2)
%
% fun is a decreasing function. The zero crossing defines the number
% of the last term we need to compute.
dx = double(X);
fun = @(m) m.*(log(dx) + 1) - log(m).*(m + 1/2) + sum(NDig)*log(10) - 1/2*(log(2*pi) - 2);

% this is a common thing to solve for in hpf, so use a simple
% bisection scheme.
mterms = termsbisector(fun);

% Each term in this series goes by two powers of X, so we will
% really need only half as many terms for convergence.
% Make sure that mterms is an odd number, by rounding up to
% an odd number if necessary.
mterms = mterms + double(mod(mterms,2) == 0);

% run the sine series in reverse now. This is why I wanted to know how
% many terms would be necessary. By running the loop backwards, I
% avoid divisions, and a divide is more expensive than a multiply for
% an hpf number.
Xsq = X.*X;
F = hpf('1',NDig);
Fact = F;
for m = mterms:-2:3
  Fact = uminus(Fact.*(m.*(m-1)));
  F = Xsq*F + Fact;
end
% because we ran the loop backwards, in the end we need to
% divide by factorial(mterms). If we started out with the wrong
% sign on F for the first term, that sign will be corrected now.
F = F.*X./Fact;

end % function F = sincore(X)

% =============================================================

function F = coscore(X)
% evaluates cos(X) for a scalar hpf number X, where -pi/4 <= X <= pi/4
% No test is made to assure that X is in the proper interval. X is
% assumed to be in radians.

% how many digits do we need to carry?
NDig = X.NumberOfDigits;

% sin(0) == 0
if X.Sign == 0
  F = hpf('1',NDig);
  return
end

% How many terms in the Taylor series do we need to compute?
% The cosine series has as its general term z^m/factorial(m), with
% m = 2*n, n being the index of the term we choose to stop at.
%
% we can quit when the natural log of that expression is less
% than the log of our precision goal. (The -1 below gives us some
% protection for round-off).
%
%  -1 + log(10^-NDig) = log(X^m/m!) = m log(X) - log(m!)
%
% the logs here are natural logs of course. Using Stirling's
% approximation for m!, we get
%
%  -1 - NDig*log(10) = m*log(X) - [1/2*log(2*pi) + 1/2*log(m) + m*log(m) - m]
%
% Combining terms, this reduces to
%
%  1/2*(log(2*pi)-2) - NDig*log(10) = m*(log(X) + 1) - log(m)*(m + 1/2)
%
% fun is a decreasing function. The zero crossing defines the number
% of the last term we need to compute.
dx = abs(double(X));
fun = @(m) m.*(log(dx) + 1) - log(m).*(m + 1/2) + sum(NDig)*log(10) - 1/2*(log(2*pi) - 2);

% this is a common thing to solve for in hpf, so use a simple
% bisection scheme.
mterms = termsbisector(fun);

% Each term in this series goes by two powers of X, so we will
% really need only half as many terms for convergence.
% Make sure that mterms is an even number, by rounding up to
% an even number if necessary.
mterms = mterms + double(mod(mterms,2));

% run the sine series in reverse now. This is why I wanted to know how
% many terms would be necessary. By running the loop backwards, I
% avoid divisions, and a divide is more expensive than a multiply for
% an hpf number.
Xsq = X.*X;
F = hpf('1',NDig);
Fact = F;
for m = mterms:-2:2
  Fact = uminus(Fact.*(m.*(m-1)));
  F = Xsq*F + Fact;
end
% because we ran the loop backwards, in the end we need to
% divide by factorial(mterms). If we started out with the wrong
% sign on F for the first term, that sign will be corrected now.
F = F./Fact;

end % function F = coscore(X)

% =============================================================
function [mant,exponentshift,signflag] = carryop(mant,decimalbase,base)
% carryop - does decimal carries, allowing for overflow
%
% arguments: (input)
% mant - a list (a row vector) of decimal migits to be checked for carries
%
% decimalbase - integer, from the set [1:6], defines the base of the
%        migit elements.
%
% base - 10.^decimalbase, precomputed to save time
%
% arguments: (output)
%  mant - vector of mantissa migit elements after carries are resolved
%        as well as any leading migits removed
%
%  exponentshift - any shift that was necessary to apply to the exponent,
%        due either to overflows in the leading migit, or because of
%        zero leading migits.
%
%  signflag - one of [-1, 0, +1], indicating if the value was actually
%        a negative number, zero, or a positive number. Normal operation
%        will have signflag always +1. 0 and -1 are the exceptions to watch
%        for.

% set these numbers as their expected values for normal cases
signflag = 1;
exponentshift = 0;
nmant = numel(mant);

% are there any zero leading migits?
k = find(mant ~= 0,1,'first');
if isempty(k)
  % the number was a true zero, with all zero migits
  signflag = 0;
  return
elseif (k > 1)
  % There were k-1 leading zero migits
  exponentshift = -(k-1)*decimalbase;
  
  % we need to trim off the leading zeros, padding with
  % the same number of trailing zeros. circshift does
  % this for us neatly. Note that we want to shift the
  % second dimension of mant, which is a row vector.
  mant = circshift(mant,[0,-(k-1)]);
end

carryindex = find((mant >= base) | (mant < 0));
while ~isempty(carryindex)
  % we don't want to do a carry past the highest
  % order migit right now
  if carryindex(1) == 1
    carryindex(1) = [];
  end
  remainder = mod(mant(carryindex),base);
  carry = (mant(carryindex) - remainder)/base;
  mant(carryindex) = remainder;
  carryindex = carryindex - 1;
  % do the carry itself
  mant(carryindex) = mant(carryindex) + carry;
  
  % Of the elements carried into, which ones are
  % still a problem?
  carryindex((mant(carryindex) < base) & (mant(carryindex) >= 0)) = [];
end

% just check again to make that no zero leading migits were created
% the case where all migits are zero should have been trapped already.
k = find(mant ~= 0,1,'first');
if isempty(k)
  % the number was a true zero, with all zero migits
  signflag = 0;
  return
elseif (k > 1)
  % There were k-1 leading zero migits
  exponentshift = exponentshift - (k-1)*decimalbase;
  
  % we need to trim off the leading zeros, padding with
  % the same number of trailing zeros. circshift does
  % this for us neatly. Note that we want to shift the
  % second dimension of mant, which is a row vector.
  mant = circshift(mant,[0,-(k-1)]);
end

% finally, look at the first migit itself. Is it negative?
% we should have caught that event before, but after the carries,
% check if the first migit is less than 0.
if mant(1) < 0
  % if so, the result was a negative number.
  signflag = -signflag;
  
  % so negate the migits, and redo the carries.
  mant = -mant;
  
  % do one more round of carries
  [mant,es] = carryop(mant,decimalbase,base);
  exponentshift = exponentshift + es;
end

% Is the first migit greater than base?
M1 = mant(1);
while (M1 >= base)
  % we do need to carry the first migit, appending one
  % or more new migits to the top end.
  remainder = mod(M1,base);
  carry = (M1 - remainder)./base;
  
  mant(1) = remainder;
  % since M1 was >= base to be in this loop, then carry
  % MUST be strictly greater than 0
  M1 = carry;
  % append a new migit to mant. later on we will strip off
  % any extraneous trailing migits.
  mant = [M1,mant]; %#ok
  exponentshift = exponentshift + decimalbase;
end

% does the first migit have leading zeros in the migit itself?
% this means that mant(1) is less than base/10. We need only
% bother with this if decimalbase is at least 2, so base is at
% least 100.
M1 = mant(1);
b10 = base/10;
if (decimalbase > 1) && (M1 < b10)
  % the first migit effectively has one or more leading zeros. remove
  % the zeros. a simple while loop is effective here, rather than
  % chasing log10 with a ceil and worrying about rounding issues
  % in log10.
  shift = 0;
  while M1 < b10
    M1 = M1*10;
    shift = shift + 1;
  end
  
  % multiply the migits to do the actual shift.
  mant = mant.*(10.^shift);
  
  exponentshift = exponentshift - shift;
  
  % we have probably now incurred a carry, but since we did a
  % carry before, the carries here won't persist past one step.
  % we also need not worry about the first migit exceeding base,
  % again because of where we are.
  carryindex = find((mant >= base) | (mant < 0));
  remainder = mod(mant(carryindex),base);
  carry = (mant(carryindex) - remainder)./base;
  
  mant(carryindex) = remainder;
  carryindex = carryindex - 1;
  mant(carryindex) = mant(carryindex) + carry;
  
end

% finally, have we appended new migits to the top end? If
% so, then they must be stripped off and rounding applied.
if numel(mant) > nmant
  % strip off trailing migits
  mtrail = mant(nmant+1);
  
  % I recall this is faster than selecting the first nmant
  % migits and overwriting mant
  mant((nmant+1):end) = [];
  
  % round on that last migit
  if (mtrail >= (base/2))
    mant(end) = mant(end) + 1;
    if mant(end) >= base
      % one final carry was just made necessary, but
      % this is a rare event. It also cannot change the sign
      % of the result.
      [mant,finalexpshift] = carryop(mant,decimalbase,base);
      
      exponentshift = exponentshift + finalexpshift;
    end
  end
end

end

% =============================================================

function [NDig,source] = combineNDig(NDig1,NDig2)
% Combined number of digits to be carried in an op betweeen two numbers
%
% arguments: (input)
% NDig1, NDig2 - pairs of vectors of length 2, containing the
%      number of reported digits, plus the number of shadow
%      digits for each operand
%
% NDig - the resulting number of digits as output

shadowdigits = min(NDig1(2),NDig2(2));

% The total number of digits to be carried is the smaller of the two
% totals.
[totaldigits,source] = min([sum(NDig1),sum(NDig2)]);

% combine the two
NDig = [totaldigits - shadowdigits , shadowdigits];

end

% =============================================================
function [D,nlz] = parsenumericdigits(dstr,opmode)
% parses out numeric digits from a character digit string
%
% arguments: (input)
%  dstr - character string of pure numeric digits. No test is
%        made to ensure that the string is composed of pure
%        digits.
%
%  opmode - controls whether leading zero digits are stripped
%        off of the string.
%
%        opmode = 0 --> (DEFAULT) Strip off any leading zeros,
%                     returning a count of how many zeros were
%                     taken off.
%
%        opmode = 1 --> Do not strip off leading zeros
%
% arguments: (output)
%  D   - Vector of integer digits, [0...9], as doubles
%        If the input was just a string of zeros, then
%        D will be a single zero.
%
%  nlz - the number of leading zero digits found

% if all are zeros, then just return 0
if isempty(dstr)
  D = [];
  nlz = 0;
elseif all(dstr == '0')
  D = 0;
  nlz = 0;
elseif opmode
  % just convert the digits to integers
  D = dstr - '0';
  nlz = 0;
else
  % we need to look for leading zeros. There must be
  % at least some
  k = find(dstr ~= '0',1,'first');
  if k == 1
    % no leading zeros, so just convert to numeric
    D = dstr - '0';
    nlz = 0;
  else
    % there were leading zeros, and we need to strip them out
    D = dstr(k:end) - '0';
    nlz = k - 1;
  end
end

end % D = parsenumericdigits(dstr,opmode)

% =============================================================

function vec = wrepvec(vec,finallength)
% replicates the elements in a vector to have the given length.
% finallength need not be an integer multiple of length(vec)

n0 = numel(vec);

% an index vector of the correct final length
ind = 1:finallength;

% do a clockmod on the index
ind = mod(ind,n0);
ind(ind == 0) = n0;

vec = vec(ind);

end

% =============================================================

function m = d2m(d,DBase)
% converts a character string of numbers from the set 0:9 into a vector of migits

% did d come in in character form?
if ischar(d)
  % check that the size of d is a multiple of DBase
  res = mod(numel(d),DBase);
  if res ~= 0
    d = [d,repmat('0',1,DBase - res)];
  end
  
  % create each row as one migit, in character form
  d = reshape(d,[],DBase);
  
  % convert to decimal integer, in base 10^DBase
  m = base2dec(d,10);
  
  % we want a row vec, and base2dec gives a column
  m = m.';
  
else
  % they were numeric digit values
  
  % check that the size of d is a multiple of DBase
  res = mod(numel(d),DBase);
  if res ~= 0
    d = [d,zeros(1,DBase - res)];
  end
  
  % create each row as one migit, in numeric form
  d = reshape(double(d),DBase,[]).';
  
  % convert to decimal integer, in base 10^DBase
  m = d(:,end);
  p = 1;
  for i = 2:DBase
    p = p*10;
    m = m + p*d(:,DBase - i  + 1);
  end
  
  % we want a row vec, and base2dec gives a column
  m = m.';
  
end

end % function d2m

% =============================================================

function D = m2d(M,DBase)
% Convert a single digit decimal string to migits of a given base
%
% M must be a list (row vector) of integers, from the set 0:(10^DBase - 1)
% DBase must be one of [1 2 3 4 5 6]

if DBase == 1
  % a no-op
  D = M;
else
  % D > 1
  D = dec2base(M,10,DBase).' - '0';
  D = D(:).';
end

end % function m2d

% =============================================================

function [Vs,tags] = sortvecA(V)
% sorts a single vector in ascending order
% The sorting scheme chosen is a basic merge sort.
n = length(V);
% special case on the length. otherwise, recurse
if n == 1
  % 1 is simple
  Vs = V;
  tags = 1;
elseif n == 2
  if V(1) <= V(2)
    Vs = V;
    tags = [1 2];
  else
    Vs = V([2 1]);
    tags = [2 1];
  end
else
  % at least 3
  n1 = floor(n/2);
  [Vs1,tags1] = sortvecA(V(1:n1));
  [Vs2,tags2] = sortvecA(V((n1+1):n));
  % shift the tags for the upper half
  tags2 = tags2 + n1;
  n2 = n - n1;
  
  % merge the two sorted subarrays
  i1 = 1;
  i2 = 1;
  k = 1;
  tags = zeros(n,1);
  while (i1<=n1) && (i2<=n2)
    if Vs1(i1) <= Vs2(i2)
      tags(k) = tags1(i1);
      i1 = i1 + 1;
      if (i1 > n1) && (k < n)
        tags((k+1):n) = tags2(i2:end);
        i2 = n2 + 1;
      end
    else
      tags(k) = tags2(i2);
      i2 = i2 + 1;
      if (i2 > n2) && (k < n)
        tags((k+1):n) = tags1(i1:end);
        i1 = n1 + 1;
      end
    end
    k = k + 1;
  end
  
  Vs = V(tags);
end

end

% =============================================================

function [Vs,tags] = sortvecD(V)
% sorts a single vector in descending order
% The sorting scheme chosen is a basic merge sort.
n = length(V);
% special case on the length. otherwise, recurse
if n == 1
  % 1 is simple
  Vs = V;
  tags = 1;
elseif n == 2
  if V(1) >= V(2)
    Vs = V;
    tags = [1 2];
  else
    Vs = V([2 1]);
    tags = [2 1];
  end
else
  % at least 3
  n1 = floor(n/2);
  [Vs1,tags1] = sortvecD(V(1:n1));
  [Vs2,tags2] = sortvecD(V((n1+1):n));
  % shift the tags for the upper half
  tags2 = tags2 + n1;
  n2 = n - n1;
  
  % merge the two sorted subarrays
  i1 = 1;
  i2 = 1;
  k = 1;
  tags = zeros(n,1);
  while (i1<=n1) && (i2<=n2)
    if Vs1(i1) >= Vs2(i2)
      tags(k) = tags1(i1);
      i1 = i1 + 1;
      if (i1 > n1) && (k < n)
        tags((k+1):n) = tags2(i2:end);
        i2 = n2 + 1;
      end
    else
      tags(k) = tags2(i2);
      i2 = i2 + 1;
      if (i2 > n2) && (k < n)
        tags((k+1):n) = tags1(i1:end);
        i1 = n1 + 1;
      end
    end
    k = k + 1;
  end
  
  Vs = V(tags);
end

end


