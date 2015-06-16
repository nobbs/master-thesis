function pn = legendrePolynomialDerivative(x, n, span)
  % Evaluate the first derivative of the shifted Legendre polynomial.
  %
  % Analogue to legendrePolynomial this is done by a recursion, this time of the
  % form
  %   ``(n - 1) P_n(x) = (2n - 1) x P_{n-1}(x) - n P_{n-2}(x)``
  % where `P_0(x) = 0` und `P_1(x) = 1`, for the legendre polynomials on the
  % interval `[-1 1]`.
  %
  % Parameters:
  %   x: grid on which the polynomial is evaluated @type matrix
  %   n: degree of the Legendre polynomial to evaluate @type integer
  %   span: interval on which the polynomial should be shifted @type vector
  %     @default [-1 1]
  %
  % Return values:
  %   pn: function values of the first derivative of the Legendre polynomial
  %     @type vector
  %
  % See also:
  %   legendrePolynomial

  % set default values
  if nargin == 2
    span = [-1, 1];
  elseif nargin < 2
    error('Not enough arguments given!');
  end

  % Shift from the default interval [-1 1] to [a b]
  xs = 2 / (span(2) - span(1)) * x + (span(1) + span(2)) / (span(1) - span(2));
  dx = 2 / (span(2) - span(1));

  % polynomials of degree 0, 1 and 2 are given, higher degrees are evaluated
  % through the recursion formula
  if n == 0
    pn = dx * zeros(size(xs, 1), size(xs, 2));
  elseif n == 1
    pn = dx * ones(size(xs, 1), size(xs, 2));
  else
    pn1 = dx * ones(size(xs, 1), size(xs, 2));
    pn = dx * 3 * xs;
    if n > 2
      for m = 3:n
        tmp = pn;
        pn = ((2 * m - 1) .* xs .* pn - m * pn1) / (m - 1);
        pn1 = tmp;
      end
    end
  end

end
