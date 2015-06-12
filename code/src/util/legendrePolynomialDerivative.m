function pn = legendrePolynomialDerivative(x, n, span)
	% Wertet die erste Ableitung eines geshifteten Legendre-Polynoms aus.
	%
	% Ähnlich zu legendrePolynomial wird eine rekursive Darstellung zur Berechnung
	% verwendet. Die Rekursionsformel für das Intervall `[-1, 1]` lautet hierbei
	%   ``(n - 1) P_n(x) = (2n - 1) x P_{n-1}(x) - n P_{n-2}(x)``
	% mit `P_0(x) = 0` und `P_1(x) = 1`.
	%
	% Parameters:
	%   x: Gitter, auf dem das Polynom ausgewertet werden soll. @type matrix
	%   n: Index des auszuwertenden Polynoms. @type int
	%   span: Intervall, auf den das Polynom geshiftet werden soll. @type vector
	%
	% Return values:
	%   pn: Auswertung der Ableitung des Legendre-Polynoms. @type matrix
	%
	% See also:
	%   legendrePolynomial

	if nargin == 2
		span = [-1, 1];
	elseif nargin < 2
		error('Zu wenige Argumente übergeben!');
	end

	a = span(1);
	b = span(2);

	% Shift der Polynome auf den Intervall [0, 1] statt [-1, 1].
	xs = 2 / (b - a) * x + (a + b) / (a - b);
	dx = 2 / (b - a);

	% Die Grade 0, 1 und 2 werden ohne Rekursion bestimmt.
	if n == 0
		pn = dx * zeros(size(xs, 1), size(xs, 2));
	elseif n == 1
		pn = dx * ones(size(xs, 1), size(xs, 2));
	else
		pn1 = dx * ones(size(xs, 1), size(xs, 2));
		pn = dx * 3 * xs;
		if n > 2
			% Rekursionsformel von Bonet ausfÃ¼hren
			for m = 3:n
				tmp = pn;
				pn = ((2 * m - 1) .* xs .* pn - m * pn1) / (m - 1);
				pn1 = tmp;
			end
		end
	end

end
