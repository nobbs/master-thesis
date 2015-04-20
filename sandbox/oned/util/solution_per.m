function [ y ] = solution(t, x, opt, coeffs)
% SOLUTION Baut aus dem berechneten Lösungsvektor eine Lösung. Dazu wird
% mittels der berechneten Koeffizienten eine Linearkombination der
% Ansatzfunktionen gebildet.

y = zeros(size(t, 1), size(t, 2));
offset = (opt.num.M - 1) / 2;

function y = phi(x, k)
    if k == 0
        y = ones(size(x, 1), size(x, 2));
    elseif k < 0
        y = cos(2 * pi * k * x);
    else
        y = sin(2 * pi * k * x);
    end
end

for j = 1:opt.num.Xj
    jj = j - 1 - offset;
    for k = 1:opt.num.Xk
    	% Zur aktuellen Ansatzfunktion zugehörigen Koeffizienten bestimmen
    	cidx = (j - 1) * opt.num.Xk + k;
    	% Normierungsfaktor bestimmen
% 		cnrm = sqrt((1 + (pi * j)^2) / (2 * (2*(k - 1) + 1)) + ...
% 		            legendre_dP(1, k - 1));
        cnrm = 1;
		% Ansatzfunktion auswerten
        y = y + coeffs(cidx) * cnrm * phi(x, jj) .* legendre_P(t, k - 1);
    end
end

end

