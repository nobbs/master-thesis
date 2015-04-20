function [NX, NY] = normierungsmatrizen(opt)
% NORMIERUNGSMATRIZEN Berechnet Normierungsmatrizen.
% Diese werden auf die Steifigkeitsmatrix und den Lastvektor multipliziert um
% die Ansatz- und Testfunktionen zu normieren.

NXv = zeros(opt.dim.X, 1);
NYv = zeros(opt.dim.Y, 1);

% Normierung bezüglich X-Norm berechnen
for j = 1:opt.num.Xj
    for k = 1:opt.num.Xk
        pos = (j - 1) * opt.num.Xk + k;
        NXv(pos) = sqrt((1 + (pi * j)^2) / (2 * (2*(k - 1) + 1)) + ...
                        legendre_dP(1, k - 1));
    end
end

% Normierung bezüglich Y-Norm berechnen
for l = 1:opt.num.Yl
    for m = 1:opt.num.Ym
        pos = (l - 1) * opt.num.Ym + m;
        NYv(pos) = sqrt((1 + (pi * l)^2) / (2 * (2*(m - 1) + 1)));
    end
end
for n = 1:opt.num.Yn
    pos = opt.num.Yl * opt.num.Ym + n;
    NYv(pos) = sqrt(1/2);
end

% Normierungsmatrizen aufbauen
NX = spdiags(NXv, 0, opt.dim.X, opt.dim.X);
NY = spdiags(NYv, 0, opt.dim.Y, opt.dim.Y);

end
