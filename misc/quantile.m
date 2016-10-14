function y = quantile(x, p, dim)
%QUANTILE Empirical (sample) quantiles.
%
%   For vectors Q = QUANTILE(X, P) is the empirical quantiles of X for the
%   probabilities in P.  The smallest observation corresponds to P = 0 and
%   the largest to P = 1.  The length of Q is LENGTH(P).
%
%   For matrices, QUANTILE(X, P) is the empirical quantiles of each column.
%
%   In general, QUANTILE(X, P) is the empirical quantiles along the first
%   non-singleton dimension.
%
%   QUANTILE(X, P, DIM) returnes the quantiles along dimension DIM.
%
%   This is a MATLAB version of the R `quantile' function.

%   This is from the `R' documentation of the `quantile' function:
%
%   The generic function `quantile' produces sample quantiles corresponding
%   to the given probabilities. The smallest observation corresponds to a
%   probability of 0 and the largest to a probability of 1.
%
%   A vector of length `length(probs)' is returned; if `names = TRUE',
%   it has a `names' attribute.
%
%   `quantile(x,p)' as a function of `p' linearly interpolates the
%   points ( (i-1)/(n-1), ox[i] ), where `ox <- order(x)' (the ``order
%   statistics'') and `n <- length(x)'.
%
%   This gives `quantile(x, p) == (1-f)*ox[i] + f*ox[i+1]', where `r
%   <- 1 + (n-1)*p', `i <- floor(r)', `f <- r - i' and `ox[n+1] :=
%   ox[n]'.

%   See also MEAN, STD, MIN, MAX, COV.

%   Author:      Peter John Acklam
%   Time-stamp:  2004-09-22 19:14:03 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   narginchk(2, 3);

   sx = size(x);                % size of `x'
   dx = ndims(x);               % number of dimensions in `x'

   % make `p' a column vector
   p = p(:);
   np = length(p);              % number of elements in `p'

   % get first non-singleton dimension, or 1 if none found
   if nargsin < 3
      k = find(sx ~= 1);
         if isempty(k)
         dim = 1;
      else
         dim = k(1);
      end
   else
      if any(size(dim) ~= 1) || dim < 1 || dim ~= round(dim)
         error('Dimension must be a scalar positive integer.');
      end
   end

   n = size(x, dim);

   % special case when `x' is empty
   if isempty(x)
      y = zeros(sx);
      return;
   end

   % permute and reshape so DIM becomes the row dimension of a 2-D array
   perm = 1:dx;
   perm([1 dim]) = [dim 1];
   x = reshape(permute(x, perm), [n prod(sx)/n]);

   % compute the quantiles by linear interpolation
   y = interp1((0:n-1).'/(n-1), sort(x, 1), p);

   % reshape and permute back
   sy = sx;
   sy(dim) = np;
   y = permute(reshape(y, sy(perm)), perm);
