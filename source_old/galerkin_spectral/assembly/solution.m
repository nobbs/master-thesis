function [ y ] = solution(t, x, opt, coeffs)
% SOLUTION Baut aus dem berechneten Lösungsvektor eine Lösung. Dazu wird
% mittels der berechneten Koeffizienten eine Linearkombination der
% Ansatzfunktionen gebildet.

y = zeros(size(t, 1), size(t, 2));

for j = 1:opt.num.Xj
    for k = 1:opt.num.Xk
    	% Zur aktuellen Ansatzfunktion zugehörigen Koeffizienten bestimmen
    	cidx = (j - 1) * opt.num.Xk + k;
    	% Normierungsfaktor bestimmen
		cnrm = sqrt((1 + (pi * j)^2) / (2 * (2*(k - 1) + 1)) + ...
		            legendre_dP(1, k - 1));
		% Ansatzfunktion auswerten
        y = y + coeffs(cidx) * cnrm * sin(pi * j * x) .* legendre_P(t, k - 1);
    end
end

end

