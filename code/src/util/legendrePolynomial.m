function pn = legendrePolynomial(x, n, span)
	% Wertet ein geshiftetes Legendre-Polynom aus.
	%
	% Verwendet wird hierzu die Rekursionsformel von Bonnet, siehe beispielsweise
	% http://de.wikipedia.org/wiki/Legendre-Polynom#Rekursionsformeln, sowie eine an
	% \cite Press:2007:NRE:1403886 angelehnter rekursiver Algorithmus. Die
	% zugrundeliegende Rekursionsformel für die Legendre-Polynome auf dem Intevall
	% `[-1, 1]` lautet
	%   ``n P_n(x) = (2n - 1) x P_{n-1}(x) - (n - 1) P_{n-2}(x)``
	% mit `P_0(x) = 1` und `P_1(x) = x`.
	%
	% Parameters:
	%   x: Gitter, auf dem das Polynom ausgewertet werden soll. @type matrix
	%   n: Index des auszuwertenden Polynoms. @type int
	%   span: Intervall, auf den das Polynom geshiftet werden soll. @type vector
	%
	% Return values:
	%   pn: Auswertung des Legendre-Polynoms. @type matrix
	%
	% See also:
	%   legendrePolynomialDerivative

	if nargin == 2
		span = [-1, 1];
	elseif nargin < 2
		error('Zu wenige Argumente übergeben!');
	end

	a = span(1);
	b = span(2);

	% Shift der Polynome auf den Intervall [0, 1] statt [-1, 1].
	xs = 2 / (b - a) * x + (a + b) / (a - b);

	% Die Grade 0, 1 und 2 werden ohne Rekursion bestimmt.
	if n == 0
		pn = ones(size(xs, 1), size(xs, 2));
	elseif n == 1
		pn = xs;
	else
		pn1 = xs;
		pn  = (3 * xs.^2 - 1) / 2;
		if n > 2
			% Rekursionsformel von Bonnet ausführen
			for m = 3:n
				tmp = pn;
				pn  = ((2 * m - 1) .* xs .* pn - (m - 1) * pn1) / m;
				pn1 = tmp;
			end
		end
	end

end
