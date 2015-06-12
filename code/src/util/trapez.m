function Q = trapez(f, grid, is_periodic)
	% Zusammengesetzte Trapez-Quadratur.
	%
	% Berechnet das Integral `\int_{a}^{b} f(x) dx` mittels summierter Trapez-
	% Quadraturformel mit äquidistantem Gitter grid. Handelt es sich bei `f` um eine
	% periodische Funktion und wird das Integral über eine (oder mehrere) ganze
	% Perioden berechnet, dann sollte is_periodic auf true gesetzt werden und sowohl
	% f als auch grid sollten für das Intervall `[a, b)` bestimmt werden, das heißt,
	% der "doppelte" Wert am Periodenanfang und Periodenende darf *nicht* doppelt
	% vorkommen.
	%
	% Parameters:
	%   f: Funktionswerte an den Gitterpunkten @type vector
	%   grid: Gitter des Intervalls `[a, b]` @type vector
	%   is_periodic: Schalter, um zwischen periodischer / nicht-periodischer
	% 		Funktion umzuschalten @default true @type bool
	%
	% Return values:
	%   Q: Approximation des Integralwertes @type double

	% Default-Werte setzen
	if nargin == 2
		is_periodic = 1;
	end

	% Check, ob Anzahl Funktionsauswertungen gleich Anzahl der Gitterpunkte.
	assert(length(f) == length(grid));

	% Schrittweite bestimmen und numerische Quadratur durchführen
	h = grid(2) - grid(1);
	if is_periodic
		Q = h * (f(1) + sum(f(2:end)));
	else
		Q = (h / 2) * (f(1) + f(end) + 2 * sum(f(2:end-1)));
	end

end
