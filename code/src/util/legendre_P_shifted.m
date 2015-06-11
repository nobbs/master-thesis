function pn = legendre_P_shifted(x, n, span)
%LEGENDRE_P Berechnet die auf [0, 1] geshifteten Legendre-Polynome.
% Berechnet das Legendre-Polynom n-ten Grades, geshiftet auf das Interval
% [0, 1], und wertet es in x aus. Die Implementierung basiert auf der
% Rekursionsformel von Bonet.

if nargin == 2
  span = [-1, 1];
elseif nargin < 2
  error('Möp!');
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
    pn2 = ones(size(xs, 1), size(xs, 2));
    pn  = (3 * xs.^2 - 1) / 2;
    if n > 2
        % Rekursionsformel von Bonet ausführen
        for m = 3:n
            tmp = pn;
            pn  = ((2 * m - 1) .* xs .* pn - (m - 1) * pn1) / m;
            pn2 = pn1;
            pn1 = tmp;
        end
    end
end

end
