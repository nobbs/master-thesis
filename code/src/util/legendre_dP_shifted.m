function pn = legendre_dP_shifted(x, n, span)
%LEGENDRE_DP Berechnet die geshiftete erste Ableitung der Legendre-Polynome.
% Berechnet die erste Ableitung des Legendre-Polynom n-ten Grades, geshiftet
% auf das Interval [0, 1], und wertet es in x aus. Die Implementierung basiert
% auf der verallgemeinerten Rekursionsformel von Bonet.

if nargin == 2
  span = [-1, 1];
elseif nargin < 2
  error('Mˆp!');
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
    pn2 = dx * zeros(size(xs, 1), size(xs, 2));
    pn = dx * 3 * xs;
    if n > 2
        % Rekursionsformel von Bonet ausf√ºhren
        for m = 3:n
            tmp = pn;
            pn = ((2 * m - 1) .* xs .* pn - m * pn1) / (m - 1);
            pn2 = pn1;
            pn1 = tmp;
        end
    end
end

end
