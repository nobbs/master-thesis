function Q = trapez(f, I, is_periodic)
% TRAPEZ  Trapez-Quadratur
%  Berechnet das Integral `\int_{I} f(x) dx` mittels summierter Trapez-Formel.
%  Parameter:
%    f  Funktionswerte in den Gitterpunkten
%    I  Gitter
%    is_periodic  Periodische Funktion, d.h. letzter Gitterpunkt "fehlt",
%                 da gleich erstem?

if nargin == 2
    is_periodic = 1;
end

% Check, ob Anzahl Funktionsauswertungen gleich Anzahl der Gitterpunkte.
assert(length(f) == length(I));

% Schrittweite bestimmen und numerische Quadratur durchführen
h = I(2) - I(1);
if is_periodic
    Q = h * (f(1) + sum(f(2:end)));
else
    Q = (h / 2) * (f(1) + f(end) + 2 * sum(f(2:end-1)));
end

end
