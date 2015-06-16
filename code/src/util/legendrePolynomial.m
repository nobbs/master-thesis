function pn = legendrePolynomial(x, n, span)
  % Evaluate a shifted Legendre polynomial.
  %
  % This is done by an implementation of the recursion formula of Bonnet, see
  % https://en.wikipedia.org/wiki/Legendre_polynomials#Recursive_definition,
  % and based on an algorithm from \cite Press:2007:NRE:1403886
  %
  % The used recursion formula for the interval `[-1, 1]` is of the form
  %   ``n P_n(x) = (2n - 1) x P_{n-1}(x) - (n - 1) P_{n-2}(x),``
  % where `P_0(x) = 1` und `P_1(x) = x`.
  %
  % Parameters:
  %   x: grid on which the polynomial is evaluated @type matrix
  %   n: degree of the Legendre polynomial to evaluate @type integer
  %   span: interval on which the polynomial should be shifted @type vector
  %     @default [-1 1]
  %
  % Return values:
  %   pn: function values of the Legendre polynomial @type vector
  %
  % See also:
  %   legendrePolynomialDerivative

  % set default values
  if nargin == 2
    span = [-1, 1];
  elseif nargin < 2
    error('Not enough arguments given!');
  end

  % Shift from the default interval [-1 1] to [a b]
  xs = 2 / (span(2) - span(1)) * x + (span(1) + span(2)) / (span(1) - span(2));

  % polynomials of degree 0, 1 and 2 are given, higher degrees are evaluated
  % through the recursion formula
  if n == 0
    pn = ones(size(xs, 1), size(xs, 2));
  elseif n == 1
    pn = xs;
  else
    pn1 = xs;
    pn  = (3 * xs.^2 - 1) / 2;
    if n > 2
      for m = 3:n
        tmp = pn;
        pn  = ((2 * m - 1) .* xs .* pn - (m - 1) * pn1) / m;
        pn1 = tmp;
      end
    end
  end

end
