function [F] = lastvektor_per(data, opt)
% LASTVEKTOR Berechnet den Lastvektor zu geg. Anfangs- und Randbedingungen.

offset = (opt.num.M - 1) / 2;
f = zeros(opt.dim.Y, 1);

% Wie zuvor betrachten wir die beiden Komponenten von $\mathcal Y_N$ getrennt.
% Zunächst also $v = (v_1, 0)$.
for l = 1:opt.num.Yl
    for m = 1:opt.num.Ym
        value = integral2(@(t, x) data.g(t, x) .* sin(pi*l*x) .* ...
                            legendre_P(t, m - 1), ...
                          data.tspan(1), data.tspan(2), ...
                          data.xspan(1), data.xspan(2));
        pos    = (l - 1) * opt.num.Ym + m;
        f(pos) = 0;
        f(pos) = value;
    end
end

% Und nun $v = (0, v_2)$.
for n = 1:opt.num.Yn
  nn = n - 1 - offset;
    % value = integral(@(x) data.u0(x) .* sin(pi*n*x), ...
    %                  data.xspan(1), data.xspan(2));
if nn == 0
  value = 1;
else
  value = 0;
end
    pos    = opt.num.Yl * opt.num.Ym + n;
    f(pos) = value;
end

F = sparse(f);

end
