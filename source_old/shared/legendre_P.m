function pn = legendre_P( x, n )
%LEGENDRE_P Berechnet die auf [0, 1] geshifteten Legendre-Polynome.
% Berechnet das Legendre-Polynom n-ten Grades, geshiftet auf das Interval
% [0, 1], und wertet es in x aus. Die Implementierung basiert auf der
% Rekursionsformel von Bonet.

% Shift der Polynome auf den Intervall [0, 1] statt [-1, 1].
x = 2 * x - 1;

% Die Grade 0, 1 und 2 werden ohne Rekursion bestimmt.
if n == 0
    pn = ones(size(x, 1), size(x, 2));
elseif n == 1
    pn = x;
else
    pn1 = x;
    pn2 = ones(size(x, 1), size(x, 2));
    pn  = (3 * x.^2 - 1) / 2;
    if n > 2
        % Rekursionsformel von Bonet ausführen
        for m = 3:n
            tmp = pn;
            pn  = ((2 * m - 1) .* x .* pn - (m - 1) * pn1) / m;
            pn2 = pn1;
            pn1 = tmp;
        end
    end
end

end
