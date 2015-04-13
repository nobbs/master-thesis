function pn = legendre_dP( x, n )
%LEGENDRE_DP Berechnet die geshiftete erste Ableitung der Legendre-Polynome.
% Berechnet die erste Ableitung des Legendre-Polynom n-ten Grades, geshiftet
% auf das Interval [0, 1], und wertet es in x aus. Die Implementierung basiert
% auf der verallgemeinerten Rekursionsformel von Bonet.

% Shift der Polynome auf den Intervall [0, 1] statt [-1, 1].
x = 2 * x - 1;
a = 2;

% Die Grade 0, 1 und 2 werden ohne Rekursion bestimmt.
if n == 0
    pn = a * zeros(size(x, 1), size(x, 2));
elseif n == 1
    pn = a * ones(size(x, 1), size(x, 2));
else
    pn1 = a * ones(size(x, 1), size(x, 2));
    pn2 = a * zeros(size(x, 1), size(x, 2));
    pn = a * 3 * x;
    if n > 2
        % Rekursionsformel von Bonet ausführen
        for m = 3:n
            tmp = pn;
            pn = ((2 * m - 1) .* x .* pn - m * pn1) / (m - 1);
            pn2 = pn1;
            pn1 = tmp;
        end
    end
end

end
