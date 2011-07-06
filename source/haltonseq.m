function H = haltonseq(NUMPTS,NDIMS);

%HALTONSEQ(NUMPTS,NDIMS,) Generate a Halton sequence in NDIMS dimensional space 
%   containing NUMPTS.  The output is between 0 and 1.  NUMPTS may be a vector
%   of integers in which they indicate which elements of the Halton sequence to
%   compute.  A vector NUMPTS is useful in sequential Halton sequence sampling.  
%
%Example:
%	H = haltonseq(1024,2);
%	subplot(2,2,1); plot(H(1:128,1),H(1:128,2),'bo');
%       title('points 1 to 128');
%	subplot(2,2,2); plot(H(129:512,1),H(129:512,2),'bo');
%       title('points 129 to 512');
%	subplot(2,2,3); plot(H(513:1024,1),H(513:1024,2),'bo');
%	title('points 513 to 1024');
%	subplot(2,2,4); plot(H(:,1),H(:,2),'bo');
%	title('points 1 to 1024');
%

if (NDIMS < 12)
	P = [2 3 5 7 11 13 17 19 23 29 31]; 
else
	P = primes(1.3*NDIMS*log(NDIMS));
	P = P(1:NDIMS);
end

%P = [2 3 5 7 11 13 17 19 23 29 31]; %The first N-primes up to 36.
				%36 is the largest base that 
				%dec2base can handle
% if (NDIMS > length(P))
% 	error(['Sorry, Halton sequences are only possible for up to ',num2str(length(P)),' dimensions.']);
% else  
% 	P = P(1:NDIMS); %Get the first NDIMS prime numbers.
% end

if isequal(size(NUMPTS),[1 1])
	int_pts = [1:NUMPTS];
else %User has put in the points to sample.
	int_pts = NUMPTS;
	NUMPTS = length(int_pts);
end

H = zeros(NUMPTS,NDIMS);



for i = 1:NDIMS %Generate the components for each dimension.
	%V = fliplr(dec2base(int_pts,P(i)));
	%V = V-'0'; %Converts string to a matrix of doubles with correct numeric 
				   %values. 
	V = fliplr(dec2bigbase(int_pts,P(i)));	
	pows = -repmat([1:size(V,2)],size(V,1),1);
	H(:,i) = sum(V.*(P(i).^pows),2);
end

function s = dec2bigbase(d,base,n)
%DEC2BIGBASE Convert decimal integer to base B vector.
%   DEC2BIGBASE(D,B) returns the representation of D as a vector of
%   digits in base B.  D must be a non-negative integer smaller than 2^52
%   and B must be an integer greater than 1.  
%
%   DEC2BIGBASE(D,B,N) produces a representation with at least N digits.
%
%   Examples
%       dec2bigbase(23,3) returns [2 1 2] 
%       dec2bigbase(23,3,5) returns [0 0 2 1 2]
%
%   See also DEC2BASE, BASE2DEC, DEC2HEX, DEC2BIN.

%   written by Douglas M. Schwarz
%   Eastman Kodak Company (on leave until 4 Jan 1999)
%   schwarz@kodak.com, schwarz@servtech.com
%   1 October 1998

error(nargchk(2,3,nargin));

if size(d,2) ~= 1, d = d(:); end

base = floor(base);
if base < 2, error('B must be greater than 1.'); end
if base == 2,
  [x,nreq] = log2(max(d));
else
  nreq = ceil(log2(max(d) + 1)/log2(base)); 
end

if nargin == 3
    nreq = max(nreq,1);
    n = max(n,nreq);
    last = n - nreq + 1;
else
    n = max(nreq,1);
    last = 1;
end

s(:,n) = rem(d,base);
while n ~= last
    n = n - 1;
    d = floor(d/base);
    s(:,n) = rem(d,base);
end
